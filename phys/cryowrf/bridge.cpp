/*
 * CRYOWRF: SNOWPACK + MeteoIO + WRF bridge
 *
 * This file provides the C interface between WRF (Fortran) and SNOWPACK (C++)
 * CRYOWRF serves as the coupling layer that integrates SNOWPACK snow physics
 * and MeteoIO meteorological data processing with the WRF atmospheric model
 */

#include "bridge.h"

#include <errno.h>
#include <pthread.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// Forward declarations for namespace functions

namespace SnowpackUtils {
void prepare_meteo_data(const MeteoInput& input, CurrentMeteo& Mdata,
                        const mio::Date& current_time,
                        double current_snow_depth,
                        const SnowpackConfig* config) {
  // Temperature sanity check
  const auto& constants = SnowpackConstants::get();
  const double safe_temp = std::max(
      constants.t_crazy_min, std::min(input.temp_air, constants.t_crazy_max));
  const double safe_rh = std::max(0.01, std::min(1.0, input.humidity));
  const double safe_wind = std::max(0.0, input.wind_speed);

  // Fill meteorological data structure using correct SNOWPACK API
  Mdata.date = current_time;
  Mdata.ta = safe_temp;  // Air temperature [K]
  Mdata.rh = safe_rh;    // Relative humidity [0-1]
  Mdata.rh_avg = safe_rh;
  Mdata.vw = safe_wind;  // Wind speed [m/s]
  Mdata.vw_avg = safe_wind;
  Mdata.vw_max = safe_wind;  // CRYOWRF src/coupler/main_coupler/Coupler.cpp:907
  Mdata.dw = input.wind_dir;  // Wind direction [degrees]
  Mdata.dw_drift = input.wind_dir;
  Mdata.vw_drift = safe_wind;
  Mdata.iswr = std::max(0.0, input.shortwave_in);  // Incoming shortwave [W/m²]
  Mdata.lw_net = input.longwave_in;  // Net longwave surrogate (no ilwr_v in
                                     // vendored Snowpack)
  Mdata.psum = std::max(0.0, input.precipitation);  // Precipitation [mm]
  Mdata.psum_ph =
      (safe_temp < constants.precip_phase_threshold)
          ? 0.0
          : 1.0;  // CRYOWRF src/coupler/main_coupler/Coupler.cpp:921
  Mdata.ea = std::max(
      0.6,
      std::min(1.0, 0.7 + 5.95e-5 * input.pressure *
                              std::exp(1500.0 * safe_rh /
                                       safe_temp)));  // Coupler atmospheric
                                                      // emissivity heuristic
  Mdata.tss =
      mio::IOUtils::nodata;  // Let SNOWPACK solve for surface temperature
  Mdata.ts0 =
      constants.surface_temp_guess;  // CRYOWRF
                                     // src/coupler/main_coupler/Coupler.cpp:909
                                     // (TSG=400K)
  Mdata.hs = std::max(0.0, current_snow_depth);  // Current snow depth [m]

  // Set roughness length and friction velocity for wind pumping calculations
  // Without these, wind pumping produces NaN thermal conductivity
  double roughness_length =
      constants.snow_roughness_length;  // Default fallback value

  if (config) {
    config->getValue("ROUGHNESS_LENGTH", "Snowpack", roughness_length,
                     mio::IOUtils::nothrow);
    if (roughness_length < 0.0) {
      roughness_length =
          constants.snow_roughness_length;  // Default if not in config
    }
  }

  // Adapt roughness based on snow presence (CRYOWRF pattern from Coupler.cpp
  // line 267)
  const double snow_depth_threshold =
      constants.snow_depth_roughness_threshold;  // 3cm snow depth threshold
  const double rough_len =
      (current_snow_depth > snow_depth_threshold)
          ? roughness_length
          : constants.bare_roughness_length;  // Bare soil z0 ~0.01
  Mdata.z0 = rough_len;

  // Compute friction velocity from wind speed using log profile (CRYOWRF
  // pattern from Coupler.cpp) u* = k * u / ln(z/z0), where k=0.4 (von Karman
  // constant)
  const double von_karman =
      0.4;  // Physical constant - matches CRYOWRF implementation
  const double z_wind = input.height;  // Measurement height [m]

  double log_argument =
      (z_wind > Mdata.z0) ? std::log(std::max(1.01, z_wind / Mdata.z0)) : 0.0;
  if (log_argument <= 0.0) {
    log_argument = 1.0;
  }
  Mdata.ustar = (z_wind > Mdata.z0)
                    ? (von_karman * std::max(0.0, Mdata.vw) / log_argument)
                    : 0.1 * std::max(0.0, Mdata.vw);  // Fallback if z <= z0

  const double rd = 287.05;  // Dry air gas constant [J/(kg K)]

  // Validate critical parameters that could cause crashes
  if (std::isnan(Mdata.ustar) || std::isinf(Mdata.ustar)) {
    printf(
        "SNOWPACK-FATAL [%d,%d]: ustar is NaN/Inf (%.6f)! This will crash "
        "SNOWPACK!\n",
        input.i_grid, input.j_grid, Mdata.ustar);
    printf("SNOWPACK-FATAL [%d,%d]: z_wind=%.6f, z0=%.6f, log(z/z0)=%.6f\n",
           input.i_grid, input.j_grid, z_wind, Mdata.z0, log_argument);
    std::abort();
  }
  if (std::isnan(Mdata.z0) || std::isinf(Mdata.z0) || Mdata.z0 <= 0.0) {
    printf(
        "SNOWPACK-FATAL [%d,%d]: z0 is invalid (%.6f)! This will crash "
        "SNOWPACK!\n",
        input.i_grid, input.j_grid, Mdata.z0);
    std::abort();
  }
}

void extract_surface_outputs(const SnowStation& station,
                             const SurfaceFluxes& fluxes, const BoundCond& bc,
                             const CurrentMeteo& meteo, SnowpackOutput& output,
                             double temp_air_fallback) {
  const auto& constants = SnowpackConstants::get();

  output.latent_flux_kg = 0.0;
  output.moisture_flux = 0.0;
  output.roughness_mom = std::max(1.0e-4, meteo.z0);
  output.roughness_heat = output.roughness_mom;
  output.friction_velocity = std::max(0.0, meteo.ustar);
  output.stability_param = meteo.psi_s;

  try {
    size_t ndata_size = station.Ndata.size();
    output.surface_temp = (ndata_size > 0)
                              ? station.Ndata.back().T
                              : temp_air_fallback;  // Surface temperature [K]
    output.snow_swe = station.swe;  // Snow water equivalent [mm]
    output.snow_depth =
        station.cH - station.Ground;  // Snow height [m] (CRYOWRF line 1129)
    output.heat_flux_sensible =
        -1.0 *
        fluxes.qs;  // Negative sign for WRF convention (CRYOWRF line 1123)
    const double latent_flux_w =
        -1.0 * bc.ql;  // Use boundary condition data (CRYOWRF line 1124)
    output.heat_flux_latent = latent_flux_w;
    output.albedo = station.Albedo;
    output.snow_coverage =
        1.0;  // Hardcoded to 1.0 following CRYOWRF pattern (line 1127)
    // Extract soil properties from SNOWPACK ground interface elements
    // SNOWPACK has full soil physics - this is why we use it instead of NoahMP
    // CRYOWRF: module_sf_snowpacklsm.F:364-366 should extract SNOWPACK soil
    // data

    // Initialize with safe defaults
    output.soil_temperature = constants.default_soil_temperature;
    output.soil_moisture_volumetric = constants.default_soil_moisture_vol;
    output.soil_moisture_liquid = constants.default_soil_moisture_liq;
    output.soil_moisture_avail = constants.default_soil_moisture_avail;
    output.soil_moisture_total = constants.default_soil_moisture_total;
    output.soil_density = constants.default_soil_density;
    output.soil_conductivity = constants.default_soil_conductivity;
    output.soil_heat_capacity = constants.default_soil_heat_capacity;

    // Extract real soil data from SNOWPACK elements
    // SNOWPACK soil elements are at the bottom of the Edata vector
    try {
      size_t n_elements = station.getNumberOfElements();
      output.soil_layer_count = 0;
      for (int idx = 0; idx < SnowpackOutput::MAX_SOIL_LAYERS; ++idx) {
        output.soil_temp_layers[idx] = output.soil_temperature;
        output.soil_moisture_vol_layers[idx] = 0.0;
        output.soil_moisture_liq_layers[idx] = 0.0;
      }

      double total_water_mm = 0.0;

      if (n_elements > 0) {
        size_t soil_node = static_cast<size_t>(station.SoilNode);

        if (soil_node > 0 && soil_node <= n_elements) {
          size_t layer_idx = 0;
          for (size_t elem_idx = soil_node - 1;
               elem_idx < n_elements &&
               layer_idx < SnowpackOutput::MAX_SOIL_LAYERS;
               ++elem_idx) {
            const ElementData& soil_elem = station.Edata[elem_idx];
            const double ice_fraction =
                (soil_elem.theta.size() > 0) ? soil_elem.theta[0] : 0.0;
            const double water_fraction =
                (soil_elem.theta.size() > 1) ? soil_elem.theta[1] : 0.0;
            const double air_fraction =
                (soil_elem.theta.size() > 2) ? soil_elem.theta[2] : 0.0;
            const double soil_fraction = std::max(
                0.0, 1.0 - ice_fraction - water_fraction - air_fraction);

            const double volumetric_liq = std::max(0.0, water_fraction);
            const double volumetric_total =
                std::max(0.0, std::min(1.0, water_fraction + ice_fraction));

            output.soil_temp_layers[layer_idx] = soil_elem.Te;
            output.soil_moisture_vol_layers[layer_idx] = volumetric_total;
            output.soil_moisture_liq_layers[layer_idx] = volumetric_liq;

            const double layer_depth_m = soil_elem.L;
            if (layer_depth_m > 0.0) {
              total_water_mm += volumetric_total * layer_depth_m * 1000.0;
            }

            if (layer_idx == 0) {
              const double availability_ref =
                  (constants.default_soil_moisture_liq > 1.0e-6)
                      ? constants.default_soil_moisture_liq
                      : 0.3;
              output.soil_temperature = soil_elem.Te;
              output.soil_moisture_volumetric = volumetric_total;
              output.soil_moisture_liquid = volumetric_liq;
              output.soil_moisture_avail =
                  std::min(1.0, (availability_ref > 0.0 && volumetric_liq > 0.0)
                                    ? volumetric_liq / availability_ref
                                    : 0.0);

              const double rho_ice = 917.0;
              const double rho_water = 1000.0;
              const double rho_air = 1.2;
              const double rho_soil = 2650.0;
              output.soil_density =
                  ice_fraction * rho_ice + water_fraction * rho_water +
                  air_fraction * rho_air + soil_fraction * rho_soil;

              const double k_ice = 2.2;
              const double k_water = 0.6;
              const double k_soil = 0.25;
              output.soil_conductivity = ice_fraction * k_ice +
                                         water_fraction * k_water +
                                         soil_fraction * k_soil;

              const double c_ice = 2100.0;
              const double c_water = 4186.0;
              const double c_soil = 800.0;
              output.soil_heat_capacity = ice_fraction * c_ice +
                                          water_fraction * c_water +
                                          soil_fraction * c_soil;
            }

            ++layer_idx;
            output.soil_layer_count = static_cast<int>(layer_idx);
          }
        }
      }

      output.soil_moisture_total = total_water_mm;

      if (output.soil_layer_count == 0) {
        output.soil_temperature = constants.default_soil_temperature;
        output.soil_moisture_volumetric = constants.default_soil_moisture_vol;
        output.soil_moisture_liquid = 0.0;
        output.soil_moisture_avail = 0.0;
        output.soil_moisture_total = 0.0;
        output.soil_density = constants.default_soil_density;
        output.soil_conductivity = constants.default_soil_conductivity;
        output.soil_heat_capacity = constants.default_soil_heat_capacity;
      }

      if (output.soil_layer_count > 0) {
        printf(
            "SNOWPACK-SOIL: Extracted %d layers; top layer T=%.1fK, "
            "theta=%.1f%%, density=%.0f kg/m³\n",
            output.soil_layer_count, output.soil_temperature,
            output.soil_moisture_volumetric * 100.0, output.soil_density);
      }
    } catch (const std::exception& e) {
      printf("SNOWPACK-WARNING: Could not extract soil properties: %s\n",
             e.what());
      // Keep default values if extraction fails
    }

    // Sanity checks for soil properties
    output.soil_moisture_volumetric =
        std::max(0.0, std::min(1.0, output.soil_moisture_volumetric));
    output.soil_moisture_liquid =
        std::max(0.0, std::min(1.0, output.soil_moisture_liquid));
    output.soil_moisture_avail =
        std::max(0.0, std::min(1.0, output.soil_moisture_avail));
    output.soil_moisture_total = std::max(0.0, output.soil_moisture_total);
    output.soil_temperature = std::max(
        150.0,
        std::min(350.0, output.soil_temperature));  // Reasonable temp range

    // Convert latent heat flux to mass flux (kg m-2 s-1) following CRYOWRF
    // convention
    const double latent_heat =
        (output.surface_temp > 273.15) ? 2.5104e6 : 2.8440e6;
    if (std::isfinite(latent_flux_w)) {
      output.latent_flux_kg = latent_flux_w / latent_heat;
    } else {
      output.latent_flux_kg = 0.0;
    }
    output.moisture_flux = output.latent_flux_kg;

    // Consistency check
    if (output.snow_depth > 0.001 && output.snow_swe <= 0.0) {
      // Use fallback density
      const auto& constants = SnowpackConstants::get();
      output.snow_swe = output.snow_depth * constants.snow_density_fallback;
    }

  } catch (const std::exception& e) {
    printf("SNOWPACK-FATAL: Exception in extract_surface_outputs: %s\n",
           e.what());
    std::abort();
  }
}

void extract_layer_data(const SnowStation& station,
                        SnowpackLayerData& layer_data) {
  try {
    size_t num_elements = station.getNumberOfElements();
    layer_data.n_layers = static_cast<int>(num_elements);

    size_t layers_to_extract = std::min(num_elements, size_t(100));

    for (size_t i = 0; i < layers_to_extract; i++) {
      const ElementData& elem = station.Edata[i];

      layer_data.layer_temp[i] = elem.Te;
      layer_data.layer_thick[i] = elem.L;
      layer_data.layer_vol_ice[i] = elem.theta[0] * 100.0;    // ICE index = 0
      layer_data.layer_vol_water[i] = elem.theta[1] * 100.0;  // WATER index = 1
      layer_data.layer_vol_air[i] = elem.theta[2] * 100.0;    // AIR index = 2

      layer_data.layer_grain_radius[i] = elem.rg;
      layer_data.layer_bond_radius[i] = elem.rb;
      layer_data.layer_dendricity[i] = elem.dd;
      layer_data.layer_sphericity[i] = elem.sp;
      layer_data.layer_cdot[i] = elem.CDot;
      layer_data.layer_meta[i] = elem.metamo;
      layer_data.layer_deposition_julian[i] = elem.depositionDate.getJulian();
      layer_data.layer_graintype[i] = static_cast<double>(elem.type);
      layer_data.layer_marker[i] = static_cast<double>(elem.mk);
    }

    // Clear unused layers
    for (size_t i = static_cast<size_t>(layer_data.n_layers); i < 100; i++) {
      layer_data.layer_temp[i] = 0.0;
      layer_data.layer_thick[i] = 0.0;
      layer_data.layer_vol_ice[i] = 0.0;
      layer_data.layer_vol_water[i] = 0.0;
      layer_data.layer_vol_air[i] = 0.0;
      layer_data.layer_grain_radius[i] = 0.0;
      layer_data.layer_bond_radius[i] = 0.0;
      layer_data.layer_dendricity[i] = 0.0;
      layer_data.layer_sphericity[i] = 0.0;
      layer_data.layer_cdot[i] = 0.0;
      layer_data.layer_meta[i] = 0.0;
      layer_data.layer_deposition_julian[i] = 0.0;
      layer_data.layer_graintype[i] = 0.0;
      layer_data.layer_marker[i] = 0.0;
    }

  } catch (const std::exception& e) {
    printf("SNOWPACK-FATAL: Exception in extract_layer_data: %s\n", e.what());
    layer_data.n_layers = 0;
    std::abort();
  }
}

void extract_budget_data(const SurfaceFluxes& fluxes, const BoundCond& bc,
                         double cumu_precip, const MeteoInput& input,
                         BudgetData& budget_data) {
  try {
    // Mass budgets (following CRYOWRF patterns) - SAFE ACCESS PATTERN
    budget_data.mass_precip = cumu_precip;  // Cumulative precipitation [kg/m²]
    budget_data.mass_sublim =
        fluxes.mass[SurfaceFluxes::MS_SUBLIMATION];  // Sublimation [kg/m²]
                                                     // (CRYOWRF line 1136)
    budget_data.mass_melt =
        fluxes
            .mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF];  // Melt runoff [kg/m²] -
                                                       // SNOWPACK v11.08 API
    budget_data.mass_swe =
        fluxes.mass[SurfaceFluxes::MS_SWE];  // Current SWE [kg/m²]

    // Calculate refreeze mass from ice base melting/freezing or negative runoff
    double ice_base_meltfreeze =
        fluxes.mass[SurfaceFluxes::MS_ICEBASE_MELTING_FREEZING];
    budget_data.mass_refreeze = (ice_base_meltfreeze > 0.0)
                                    ? ice_base_meltfreeze
                                    : 0.0;  // Positive values = freezing

    // Energy budgets
    budget_data.energy_lw_in = input.longwave_in;  // Incoming longwave [W/m²]
    budget_data.energy_lw_out =
        fluxes.lw_out;  // Outgoing longwave [W/m²] (CRYOWRF line 1087)
    budget_data.energy_sw_in = input.shortwave_in;  // Incoming shortwave [W/m²]
    budget_data.energy_sw_out = input.shortwave_in * (1.0) -
                                fluxes.sw_in;    // Reflected shortwave [W/m²]
    budget_data.energy_sensible = fluxes.qs;     // Sensible heat flux [W/m²]
    budget_data.energy_latent = bc.ql;           // Latent heat flux [W/m²]
    budget_data.energy_ground_flux = fluxes.qg;  // Ground heat flux
    budget_data.energy_rain =
        fluxes.qr;  // Rain heat flux [W/m²] (CRYOWRF pattern)
    budget_data.energy_total =
        fluxes.dIntEnergy / input.dt;  // Total energy [W/m²] (CRYOWRF pattern)

  } catch (const std::exception& e) {
    printf("SNOWPACK-FATAL: Exception in extract_budget_data: %s\n", e.what());
    std::abort();
  }
}

}  // namespace SnowpackUtils
void SnowpackBridge::initialize_time(int start_year, int start_month,
                                     int start_day, int start_hour,
                                     int start_minute) {
  if (time_initialized_.load(std::memory_order_acquire)) {
    return;
  }

  if (!config_initialized_.load(std::memory_order_acquire)) {
    throw std::runtime_error("Configuration must be initialized before time");
  }

  std::call_once(time_once_, [this, start_year, start_month, start_day,
                              start_hour, start_minute]() {
    try {
      current_simulation_date_ = mio::Date(start_year, start_month, start_day,
                                           start_hour, start_minute, 0.0, 0.0);
      simulation_start_date_ = current_simulation_date_;

      {
        std::lock_guard<std::mutex> guard(time_mutex_);
        station_time_state_.clear();
      }

      // Get calculation step length from configuration (in minutes)
      calculation_step_length_ =
          config_->get("CALCULATION_STEP_LENGTH", "Snowpack");
      time_initialized_.store(true, std::memory_order_release);

      printf("SNOWPACK-INFO: Time initialized to %04d-%02d-%02d %02d:%02d\n",
             start_year, start_month, start_day, start_hour, start_minute);

    } catch (const std::exception& e) {
      printf("SNOWPACK-FATAL: Time initialization failed: %s\n", e.what());
      throw;
    }
  });

  if (!time_initialized_.load(std::memory_order_acquire)) {
    throw std::runtime_error("SNOWPACK-ERROR: Time not initialized");
  }
}

Snowpack* SnowpackBridge::get_or_create_snowpack_instance(
    int i_grid, int j_grid, double wrf_lat, double wrf_lon, double wrf_alt) {
  const std::string station_key = SnowpackConstants::get().station_id_prefix +
                                  "_1_" + std::to_string(i_grid) + "_" +
                                  std::to_string(j_grid);

  std::lock_guard<std::mutex> lock(station_mutex_);
  auto it = grid_snowpack_instances_.find(station_key);
  if (it != grid_snowpack_instances_.end()) {
    return it->second.get();
  }

  auto snowpack = std::unique_ptr<Snowpack>(new Snowpack(*config_));
  Snowpack* snowpack_ptr = snowpack.get();
  grid_snowpack_instances_[station_key] = std::move(snowpack);
  return snowpack_ptr;
}

SnowStation* SnowpackBridge::get_or_create_snowstation(int i_grid, int j_grid,
                                                       int wrf_domain_id,
                                                       double wrf_lat,
                                                       double wrf_lon,
                                                       double wrf_alt) {
  static std::atomic<int> station_log_budget{10};
  const std::string station_key = SnowpackConstants::get().station_id_prefix +
                                  "_" + std::to_string(wrf_domain_id) + "_" +
                                  std::to_string(i_grid) + "_" +
                                  std::to_string(j_grid);
  bool use_canopy = false;
  bool use_soil = false;

  std::lock_guard<std::mutex> lock(station_mutex_);

  auto station_it = grid_snowstations_.find(station_key);
  if (station_it != grid_snowstations_.end()) {
    return station_it->second.get();
  }

  // Create new SnowStation for this grid point - read configuration for correct
  // parameters
  if (!config_) {
    printf(
        "SNOWPACK-FATAL: Configuration not initialized when creating "
        "SnowStation\n");
    std::abort();
  }

  config_->getValue("CANOPY", "Snowpack", use_canopy);
  config_->getValue("SNP_SOIL", "Snowpack", use_soil);
  const bool alpine3d = false;  // WRF integration, not Alpine3D
  const bool sea_ice = false;   // Standard snowpack mode
  auto new_station = std::unique_ptr<SnowStation>(
      new SnowStation(use_canopy, use_soil, alpine3d, sea_ice));

  // Set up basic station metadata using WRF coordinates
  mio::Coords position;

  try {
    position.setLatLon(wrf_lat, wrf_lon,
                       wrf_alt);  // Correct order: lat, lon, alt
    if (station_log_budget.load(std::memory_order_relaxed) > 0) {
      const int remaining = station_log_budget.fetch_sub(1, std::memory_order_relaxed);
      if (remaining > 0) {
        printf(
            "SNOWPACK-INFO: Successfully set geographic coordinates for station "
            "(%d,%d)\n",
            i_grid, j_grid);
        if (remaining == 1) {
          printf(
              "SNOWPACK-INFO: Further coordinate messages suppressed (quota reached).\n");
        }
      } else {
        station_log_budget.store(0, std::memory_order_relaxed);
      }
    }
  } catch (const mio::InvalidArgumentException& e) {
    printf("SNOWPACK-ERROR: Invalid coordinates for station (%d,%d): %s\n",
           i_grid, j_grid, e.what());
    printf(
        "SNOWPACK-ERROR: Received lat=%.6f, lon=%.6f - these may not be "
        "geographic coordinates\n",
        wrf_lat, wrf_lon);
    throw;
  }

  std::string stationID = "WRF_GRID_" + std::to_string(wrf_domain_id) + "_" +
                          std::to_string(i_grid) + "_" + std::to_string(j_grid);
  std::string stationName =
      "WRF Grid Point " + std::to_string(i_grid) + "," + std::to_string(j_grid);
  new_station->meta.setStationData(position, stationID, stationName);

  // Try to load existing .sno file state (CRYOWRF pattern)
  bool loaded_from_file = false;
  if (io_) {
    const std::string sno_filename = "snpack_" + std::to_string(wrf_domain_id) +
                                     "_" + std::to_string(i_grid) + "_" +
                                     std::to_string(j_grid) + ".sno";
    try {
      SN_SNOWSOIL_DATA ssdata;
      ZwischenData zdata;
      mio::Date profile_date;

      std::lock_guard<std::mutex> io_lock(io_mutex_);
      io_->readSnowCover(sno_filename, stationID, ssdata, zdata, false);

      ssdata.meta.position = position;
      ssdata.meta.stationID = stationID;
      ssdata.meta.stationName = stationName;
      new_station->initialize(ssdata, 0);

      for (size_t e = 0; e < new_station->getNumberOfElements(); e++) {
        if (new_station->Edata[e].k.size() < 3) {
          new_station->Edata[e].k.resize(3, 0.0);
        }
        if (new_station->Edata[e].c.size() < 3) {
          new_station->Edata[e].c.resize(3, 0.0);
        }

        for (size_t idx = 0; idx < 3; idx++) {
          if (std::isnan(new_station->Edata[e].k[idx]) ||
              new_station->Edata[e].k[idx] != new_station->Edata[e].k[idx]) {
            new_station->Edata[e].k[idx] = 0.0;
          }
          if (std::isnan(new_station->Edata[e].c[idx]) ||
              new_station->Edata[e].c[idx] != new_station->Edata[e].c[idx]) {
            new_station->Edata[e].c[idx] = 0.0;
          }
        }

        new_station->Edata[e].heatCapacity();
      }

      loaded_from_file = true;

    } catch (const std::exception& e) {
      static int init_failure_count = 0;
      if (init_failure_count < 5) {
        printf(
            "SNOWPACK-INFO: No existing state for grid (%d,%d), starting fresh "
            "(%s)\n",
            i_grid, j_grid, e.what());
        init_failure_count++;
      }
    }
  } else {
    printf(
        "SNOWPACK-INIT [%d,%d]: SnowpackIO is NULL - cannot load .sno files\n",
        i_grid, j_grid);
  }

  if (!loaded_from_file) {
    printf("SNOWPACK-INFO: Fresh SnowStation created for grid (%d,%d)\n",
           i_grid, j_grid);
  }

  SnowStation* station_ptr = new_station.get();
  grid_snowstations_[station_key] = std::move(new_station);
  return station_ptr;
}

SnowStation* SnowpackBridge::get_existing_snowstation(int i_grid, int j_grid) {
  std::string station_key = SnowpackConstants::get().station_id_prefix + "_1_" +
                            std::to_string(i_grid) + "_" +
                            std::to_string(j_grid);
  std::lock_guard<std::mutex> lock(station_mutex_);
  auto it = grid_snowstations_.find(station_key);
  if (it != grid_snowstations_.end()) {
    return it->second.get();
  }

  return nullptr;  // Station doesn't exist
}

void SnowpackBridge::save_snowstation_state(int i_grid, int j_grid) {
  if (!io_ || !config_initialized_.load(std::memory_order_acquire)) {
    return;
  }

  std::string station_key = SnowpackConstants::get().station_id_prefix + "_1_" +
                            std::to_string(i_grid) + "_" +
                            std::to_string(j_grid);
  SnowStation* station_ptr = nullptr;
  {
    std::lock_guard<std::mutex> lock(station_mutex_);
    auto station_it = grid_snowstations_.find(station_key);
    if (station_it != grid_snowstations_.end()) {
      station_ptr = station_it->second.get();
    }
  }

  if (station_ptr) {
    mio::Date snapshot_time = simulation_start_date_;
    {
      std::lock_guard<std::mutex> time_lock(time_mutex_);
      auto ts_it = station_time_state_.find(station_key);
      if (ts_it != station_time_state_.end() && ts_it->second.initialized) {
        snapshot_time = ts_it->second.current_time;
      }
    }
    try {
      std::lock_guard<std::mutex> io_lock(io_mutex_);
      ZwischenData zdata;  // Empty for basic usage
      io_->writeSnowCover(snapshot_time, *station_ptr, zdata, true);
    } catch (const std::exception& e) {
      printf("SNOWPACK-ERROR: Failed to save state for grid (%d,%d): %s\n",
             i_grid, j_grid, e.what());
    }
  }
}

void SnowpackBridge::save_all_snowpack_states() {
  if (!io_ || !config_initialized_.load(std::memory_order_acquire)) {
    printf("SNOWPACK-WARNING: Cannot save states - IO not initialized\n");
    return;
  }

  int saved_count = 0;
  std::vector<std::pair<std::string, SnowStation*>> stations_snapshot;
  {
    std::lock_guard<std::mutex> lock(station_mutex_);
    stations_snapshot.reserve(grid_snowstations_.size());
    for (const auto& entry : grid_snowstations_) {
      stations_snapshot.emplace_back(entry.first, entry.second.get());
    }
  }

  for (const auto& entry : stations_snapshot) {
    const std::string& key = entry.first;
    SnowStation* station = entry.second;
    if (!station) {
      continue;
    }
    try {
      mio::Date snapshot_time = simulation_start_date_;
      {
        std::lock_guard<std::mutex> time_lock(time_mutex_);
        auto ts_it = station_time_state_.find(key);
        if (ts_it != station_time_state_.end() && ts_it->second.initialized) {
          snapshot_time = ts_it->second.current_time;
        }
      }
      std::lock_guard<std::mutex> io_lock(io_mutex_);
      ZwischenData zdata;  // Empty for basic usage
      io_->writeSnowCover(snapshot_time, *station, zdata, true);
      saved_count++;
    } catch (const std::exception& e) {
      printf("SNOWPACK-ERROR: Failed to save state for %s: %s\n", key.c_str(),
             e.what());
    }
  }

  printf("SNOWPACK-INFO: Saved %d snowpack states to disk\n", saved_count);
}

void SnowpackBridge::execute_snowpack(const MeteoInput& input,
                                      SnowpackOutput& output,
                                      SnowpackLayerData* layer_data,
                                      BudgetData* budget_data) {
  static std::atomic<int> debug_trace_counter(500);
  // Validate grid coordinates first
  if (input.i_grid < 0 || input.i_grid > 10000 || input.j_grid < 0 ||
      input.j_grid > 10000) {
    printf(
        "SNOWPACK-ERROR: Invalid grid coordinates (%d,%d) - outside expected "
        "range [0-10000]\n",
        input.i_grid, input.j_grid);
    printf(
        "SNOWPACK-ERROR: This indicates a problem with WRF loop bounds "
        "initialization\n");
    return;
  }

  // Track physics calls for debugging
  int call_count = ++execute_call_count_;
  if (call_count <= 5 || (call_count % 1000 == 0)) {
    if (call_count <= 5) {
      printf("SNOWPACK-INFO: Executed snowpack call #%d - Grid (%d,%d)\n",
             call_count, input.i_grid, input.j_grid);
    } else {
      printf("SNOWPACK-INFO: Executed %d snowpack calls\n", call_count);
    }
  }

  // Validate initialization
  if (!config_initialized_.load(std::memory_order_acquire)) {
    printf(
        "SNOWPACK-FATAL: Configuration not initialized - call "
        "initialize_config() first\n");
    std::abort();
  }

  if (!time_initialized_.load(std::memory_order_acquire)) {
    printf(
        "SNOWPACK-FATAL: Time not initialized - call initialize_time() "
        "first\n");
    std::abort();
  }

  try {
    // 1. Get/create persistent SNOWPACK objects (singleton pattern only)
    auto objects = SnowpackObjects::create_station_objects(input, *this);

    if (execute_call_count_.load(std::memory_order_relaxed) < 5) {
      printf("SNOWPACK-TRACE exec-start [%d,%d]\n", input.i_grid, input.j_grid);
      fflush(stdout);
    }

    int trace_remaining = debug_trace_counter.load(std::memory_order_relaxed);
    if (trace_remaining > 0) {
      if (debug_trace_counter.compare_exchange_strong(
              trace_remaining, trace_remaining - 1,
              std::memory_order_relaxed)) {
        printf(
            "SNOWPACK-TRACE [%d,%d]: tid=%lu station=%p snowpack=%p call=%d\n",
            input.i_grid, input.j_grid,
            static_cast<unsigned long>(pthread_self()), (void*)objects.station,
            (void*)objects.instance, execute_call_count_.load());
      }
    }

    // 2. Ensure station is properly initialized (load from file or create
    // fresh)
    SnowpackObjects::ensure_station_initialized(objects.station, input, *this);

    // 3. Advance simulation time following CRYOWRF per-station counter pattern
    const std::string time_key = SnowpackConstants::get().station_id_prefix +
                                 "_1_" + std::to_string(input.i_grid) + "_" +
                                 std::to_string(input.j_grid);
    const double measured_height = std::max(1.0, input.height);
    bool height_changed = false;
    StationTimeState station_time_snapshot;
    {
      std::lock_guard<std::mutex> time_lock(time_mutex_);
      StationTimeState& state = station_time_state_[time_key];
      if (!state.initialized) {
        state.current_time = simulation_start_date_;
        state.call_counter = 0;
        state.initialized = true;
        state.last_height_m = measured_height;
        height_changed = true;
      }

      state.call_counter += 1;
      const double dt_seconds =
          (input.dt > 0.0) ? input.dt : calculation_step_length_ * 60.0;
      const double dt_days = dt_seconds / 86400.0;
      state.current_time += dt_days;
      if (!height_changed &&
          std::abs(state.last_height_m - measured_height) > 1.0e-6) {
        state.last_height_m = measured_height;
        height_changed = true;
      }
      station_time_snapshot = state;
    }

    current_simulation_date_ = station_time_snapshot.current_time;

    if (execute_call_count_.load(std::memory_order_relaxed) < 5) {
      printf("SNOWPACK-TRACE exec-time [%d,%d]: julian=%.6f height=%.3f\n",
             input.i_grid, input.j_grid,
             station_time_snapshot.current_time.getJulian(), measured_height);
      fflush(stdout);
    }

    if (height_changed && config_) {
      // CRYOWRF src/coupler/main_coupler/Coupler.cpp:832-858 updates
      // HEIGHT_OF_METEO_VALUES each call
      try {
        std::lock_guard<std::mutex> config_lock(config_mutex_);
        std::ostringstream height_stream;
        height_stream.setf(std::ios::fixed);
        height_stream.precision(6);
        height_stream << measured_height;
        config_->deleteKey("HEIGHT_OF_WIND_VALUE", "Snowpack");
        config_->deleteKey("HEIGHT_OF_METEO_VALUES", "Snowpack");
        const std::string height_value = height_stream.str();
        config_->addKey("HEIGHT_OF_WIND_VALUE", "Snowpack", height_value);
        config_->addKey("HEIGHT_OF_METEO_VALUES", "Snowpack", height_value);
      } catch (const std::exception& e) {
        printf(
            "SNOWPACK-WARNING [%d,%d]: Unable to update "
            "HEIGHT_OF_METEO_VALUES: %s\n",
            input.i_grid, input.j_grid, e.what());
      }
    }

    // 4. Prepare meteorological data using utility function
    auto Mdata = std::unique_ptr<CurrentMeteo>(new CurrentMeteo());
    double current_snow_depth =
        (execute_call_count_ > 1) ? objects.station->cH : 0.0;
    SnowpackUtils::prepare_meteo_data(input, *Mdata,
                                      station_time_snapshot.current_time,
                                      current_snow_depth, config_.get());
    int trace_after_prepare =
        debug_trace_counter.load(std::memory_order_relaxed);
    if (trace_after_prepare > 0 && execute_call_count_.load() <= 20) {
      printf(
          "SNOWPACK-TRACE prepare [%d,%d]: tid=%lu ta=%.2f wind=%.2f rh=%.3f "
          "dt=%.1f u*=%.3f z0=%.4f\n",
          input.i_grid, input.j_grid,
          static_cast<unsigned long>(pthread_self()), Mdata->ta, Mdata->vw,
          Mdata->rh, input.dt, Mdata->ustar, Mdata->z0);
      fflush(stdout);
    }

    // 5. Create surface fluxes and boundary condition objects
    auto surfFluxes = std::unique_ptr<SurfaceFluxes>(new SurfaceFluxes());
    auto sn_Bdata = std::unique_ptr<BoundCond>(new BoundCond());

    // 6. SINGLE SNOWPACK CALL (CRYOWRF pattern - line 1074 in original)
    double cumu_precip = 0.0;
    objects.instance->runSnowpackModel(*Mdata, *objects.station, cumu_precip,
                                       *sn_Bdata, *surfFluxes);

    if (execute_call_count_.load(std::memory_order_relaxed) < 5) {
      printf("SNOWPACK-TRACE exec-run [%d,%d]: completed Snowpack core\n",
             input.i_grid, input.j_grid);
      fflush(stdout);
    }

    // 7. Collect surface fluxes (CRYOWRF pattern - line 1078 in original)
    try {
      surfFluxes->collectSurfaceFluxes(*sn_Bdata, *objects.station, *Mdata);
    } catch (const std::exception& e) {
      printf("SNOWPACK-ERROR: Exception in collectSurfaceFluxes: %s\n",
             e.what());
      // Continue with partial data rather than aborting
    }

    // 8. Extract surface outputs using utility function (CRYOWRF pattern -
    // lines 1121-1130)
    SnowpackUtils::extract_surface_outputs(*objects.station, *surfFluxes,
                                           *sn_Bdata, *Mdata, output,
                                           input.temp_air);

    // 9. Extract detailed layer data if requested (CRYOWRF pattern - lines
    // 1153-1185)
    if (layer_data) {
      SnowpackUtils::extract_layer_data(*objects.station, *layer_data);
    }

    // 10. Extract budget data if requested (CRYOWRF pattern - lines 1134-1150)
    if (budget_data) {
      SnowpackUtils::extract_budget_data(*surfFluxes, *sn_Bdata, cumu_precip,
                                         input, *budget_data);
      budget_data->mass_swe = output.snow_swe;  // Ensure consistency
    }

    if (execute_call_count_.load() <= 20) {
      printf(
          "SNOWPACK-TRACE output [%d,%d]: tid=%lu Tsurf=%.2f SWE=%.2f "
          "depth=%.2f HFX=%.2f QFX=%.6f\n",
          input.i_grid, input.j_grid,
          static_cast<unsigned long>(pthread_self()), output.surface_temp,
          output.snow_swe, output.snow_depth, output.heat_flux_sensible,
          output.heat_flux_latent);
      fflush(stdout);
    }

    // 11. Validate outputs
    SnowpackObjects::validate_output_values(output, input.i_grid, input.j_grid);
    if (layer_data) {
      SnowpackObjects::validate_layer_data(*layer_data, input.i_grid,
                                           input.j_grid);
    }

  } catch (const std::exception& e) {
    printf("SNOWPACK-FATAL [%d,%d]: Exception in execute_snowpack: %s\n",
           input.i_grid, input.j_grid, e.what());
    printf(
        "SNOWPACK-ERROR: Grid point (%d,%d) - Ta=%.2fK, RH=%.2f, "
        "Precip=%.3fmm\n",
        input.i_grid, input.j_grid, input.temp_air, input.humidity,
        input.precipitation);

    // Provide a benign fallback so WRF can continue running.
    // Equivalent to the resilience path in the original coupler
    // (CRYOWRF/src/coupler/main_coupler/Coupler.cpp:1074 onwards).
    output.surface_temp = input.temp_air;
    output.snow_swe = 0.0;
    output.snow_depth = 0.0;
    output.heat_flux_sensible = 0.0;
    output.heat_flux_latent = 0.0;
    output.albedo = 0.5;
    output.snow_coverage = 0.0;
    output.friction_velocity = 0.0;
    output.stability_param = 0.0;
    output.soil_moisture_volumetric = 0.0;
    output.soil_temperature = input.temp_air;
    output.soil_density = 0.0;
    output.soil_conductivity = 0.0;
    output.soil_heat_capacity = 0.0;
    output.soil_moisture_liquid = 0.0;
    output.soil_moisture_avail = 0.0;
    output.soil_moisture_total = 0.0;
    output.roughness_mom = 0.002;
    output.roughness_heat = 0.002;
    output.latent_flux_kg = 0.0;
    output.moisture_flux = 0.0;

    if (layer_data) {
      layer_data->n_layers = 0;
      std::fill(std::begin(layer_data->layer_temp),
                std::end(layer_data->layer_temp), 0.0);
      std::fill(std::begin(layer_data->layer_thick),
                std::end(layer_data->layer_thick), 0.0);
      std::fill(std::begin(layer_data->layer_vol_ice),
                std::end(layer_data->layer_vol_ice), 0.0);
      std::fill(std::begin(layer_data->layer_vol_water),
                std::end(layer_data->layer_vol_water), 0.0);
      std::fill(std::begin(layer_data->layer_vol_air),
                std::end(layer_data->layer_vol_air), 0.0);
      std::fill(std::begin(layer_data->layer_grain_radius),
                std::end(layer_data->layer_grain_radius), 0.0);
      std::fill(std::begin(layer_data->layer_bond_radius),
                std::end(layer_data->layer_bond_radius), 0.0);
      std::fill(std::begin(layer_data->layer_dendricity),
                std::end(layer_data->layer_dendricity), 0.0);
      std::fill(std::begin(layer_data->layer_sphericity),
                std::end(layer_data->layer_sphericity), 0.0);
    }

    if (budget_data) {
      budget_data->mass_precip = 0.0;
      budget_data->mass_sublim = 0.0;
      budget_data->mass_melt = 0.0;
      budget_data->mass_swe = 0.0;
      budget_data->mass_refreeze = 0.0;
      budget_data->energy_lw_in = 0.0;
      budget_data->energy_lw_out = 0.0;
      budget_data->energy_sw_in = 0.0;
      budget_data->energy_sw_out = 0.0;
      budget_data->energy_sensible = 0.0;
      budget_data->energy_latent = 0.0;
      budget_data->energy_ground_flux = 0.0;
      budget_data->energy_rain = 0.0;
      budget_data->energy_total = 0.0;
    }

    return;
  }
}

// Global singleton instance for Fortran compatibility
SnowpackBridge& g_snowpack_bridge = SnowpackBridge::instance();

void SnowpackBridge::initialize_config(const std::string& ini_path) {
  if (config_initialized_.load(std::memory_order_acquire)) {
    return;
  }

  std::call_once(config_once_, [this, ini_path]() {
    // DEBUG: Track initialization attempts
    printf("SNOWPACK-DEBUG: initialize_config invoked by thread %lu\n",
           pthread_self());

    try {
      mio::Config file_config =
          SnowpackConfigManager::loadConfiguration(ini_path);
      SnowpackConfigManager::validateConfiguration(file_config);

      config_ =
          std::unique_ptr<SnowpackConfig>(new SnowpackConfig(file_config));

      // IMPORTANT: Ensure METEO_STEP_LENGTH is available in the final config
      // SnowpackConfig might re-read the ini file, so we need to inject
      // METEO_STEP_LENGTH again
      const double calculation_step_length =
          config_->get("CALCULATION_STEP_LENGTH", "Snowpack");
      const double meteo_step_length =
          calculation_step_length * 60.0;  // Convert minutes to seconds

      std::stringstream ss_meteo_length;
      ss_meteo_length << meteo_step_length;
      config_->addKey("METEO_STEP_LENGTH", "Snowpack", ss_meteo_length.str());

      // Verify it was added
      std::string verify_value;
      config_->getValue("METEO_STEP_LENGTH", "Snowpack", verify_value);
      printf(
          "SNOWPACK-CONFIG: Final METEO_STEP_LENGTH set to %s seconds (from "
          "CALCULATION_STEP_LENGTH=%.1f minutes)\n",
          verify_value.c_str(), calculation_step_length);

      // DEBUG: Check config pointer and content before creating SnowpackIO
      printf("SNOWPACK-DEBUG: About to create SnowpackIO with config_=0x%p\n",
             (void*)config_.get());
      if (config_) {
        printf("SNOWPACK-DEBUG: Config SNOW_WRITE exists: %s\n",
               config_->keyExists("SNOW_WRITE", "Output") ? "YES" : "NO");
        if (config_->keyExists("SNOW_WRITE", "Output")) {
          std::string snow_write_val = config_->get("SNOW_WRITE", "Output");
          printf("SNOWPACK-DEBUG: Config SNOW_WRITE value: '%s'\n",
                 snow_write_val.c_str());
        }
      } else {
        printf("SNOWPACK-DEBUG: ERROR: config_ pointer is NULL!\n");
      }

      SnowpackConstants::load_from_config(config_.get());

      io_ = std::unique_ptr<SnowpackIO>(new SnowpackIO(*config_));
      printf("SNOWPACK-DEBUG: SnowpackIO created successfully\n");

      config_file_path_ = ini_path;
      config_initialized_.store(true, std::memory_order_release);

      printf("SNOWPACK-INFO: Configuration initialized from %s\n",
             ini_path.c_str());
      printf("SNOWPACK-DEBUG: Thread %lu created singleton config_=0x%p\n",
             pthread_self(), (void*)config_.get());

    } catch (const std::exception& e) {
      printf("SNOWPACK-FATAL: Configuration failed for %s: %s\n",
             ini_path.c_str(), e.what());
      throw;
    }
  });

  if (!config_initialized_.load(std::memory_order_acquire)) {
    throw std::runtime_error("SNOWPACK-ERROR: Configuration not initialized");
  }
}
