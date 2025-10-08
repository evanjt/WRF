/*
 * Snowpack utilities - meteorological data preparation and time management
 * Extracted from snowpack_bridge.cpp to improve code organization
 */

#include <cmath>
#include <cstdio>
#include <algorithm>
#include "snowpack_bridge.h"
#include "snowpack_bridge_structs.h"
#include "meteoio/meteoio/MeteoIO.h"
#include "snowpack/snowpack/DataClasses.h"
#include "snowpack/snowpack/SnowpackConfig.h"

namespace {

// Initialize time step management following CRYOWRF pattern
void initialize_time_step_internal(const MeteoInput& input,
                                   mio::Date& current_time,
                                   double calculation_step_length,
                                   int& call_counter,
                                   int& compute_counter,
                                   bool& first_call) {

    if (first_call) {
        double wrf_dt = input.dt;  // WRF timestep in seconds
        double snowpack_dt = calculation_step_length * 60.0;  // SNOWPACK timestep in seconds

        // WRF calls physics every wrf_dt, SNOWPACK needs to run every snowpack_dt
        // So we run SNOWPACK every (wrf_dt/snowpack_dt) WRF calls
        compute_counter = (int)(wrf_dt / snowpack_dt);
        first_call = false;

        printf("SNOWPACK-INFO: Time management initialized - WRF dt=%.1fs, SNOWPACK dt=%.1fs, compute_counter=%d\n",
               wrf_dt, snowpack_dt, compute_counter);
    }

    call_counter++;

    // Only advance time when counter matches (CRYOWRF pattern)
    if (compute_counter > 0 && (call_counter % compute_counter) == 0) {
        current_time += (calculation_step_length / 1440.0);  // Convert minutes to days
        printf("SNOWPACK-INFO: Time advanced to %s (call %d)\n",
               current_time.toString().c_str(), call_counter);
    }
}

} // anonymous namespace

namespace SnowpackUtils {

void prepare_meteo_data(const MeteoInput& input,
                        CurrentMeteo& Mdata,
                        const mio::Date& current_time,
                        double current_snow_depth,
                        const SnowpackConfig* config) {

    // Temperature sanity check
    double safe_temp = std::max(SnowpackConstants::T_CRAZY_MIN_KELVIN,
                               std::min(input.temp_air, SnowpackConstants::T_CRAZY_MAX_KELVIN));

    // Fill meteorological data structure using correct SNOWPACK API
    Mdata.date = current_time;
    Mdata.ta = safe_temp;                                       // Air temperature [K]
    Mdata.rh = std::max(0.01, std::min(1.0, input.humidity));  // Relative humidity [0-1]
    Mdata.vw = std::max(0.1, input.wind_speed);               // Wind speed [m/s]
    Mdata.dw = input.wind_dir;                                 // Wind direction [degrees]
    Mdata.iswr = std::max(0.0, input.shortwave_in);           // Incoming shortwave [W/m²]
    Mdata.lw_net = std::max(0.0, input.longwave_in);         // Net longwave radiation [W/m²]
    Mdata.psum = std::max(0.0, input.precipitation);          // Precipitation [mm]

    // Calculate atmospheric emissivity using approximation (Stull 1988)
    // CRYOWRF SOURCE: Derived from atmospheric emissivity patterns in CRYOWRF coupler
    double atmospheric_emissivity = 0.7 + 5.95e-5 * input.pressure * exp(1500.0 * input.humidity / safe_temp);
    Mdata.ea = std::max(0.6, std::min(1.0, atmospheric_emissivity));  // Atmospheric emissivity [dimensionless]

    // Set ground temperature (CRYOWRF compatibility)
    // CRYOWRF SOURCE: ./CRYOWRF/src/coupler/main_coupler/Coupler.cpp - vecMyMeteo(MeteoData::TSG) = 400.00
    Mdata.ts0 = 400.0;  // Set ground temperature [K] (CRYOWRF pattern - fixed 400K ground temp)

    // Set maximum wind speed (CRYOWRF compatibility)
    // CRYOWRF SOURCE: ./CRYOWRF/src/coupler/main_coupler/Coupler.cpp - vecMyMeteo(MeteoData::VW_MAX) = l_VW_MAX
    Mdata.vw_max = input.wind_speed;  // Set max wind speed [m/s] (CRYOWRF pattern - track max wind)

    // Additional required meteorological parameters
    Mdata.psum_ph = (safe_temp < SnowpackConstants::PRECIP_PHASE_THRESHOLD_K) ? 0.0 : 1.0;  // Precipitation phase (0=snow, 1=rain)
    Mdata.tss = mio::IOUtils::nodata;                         // Surface temperature (let SNOWPACK compute)
    Mdata.ts0 = safe_temp - SnowpackConstants::BOTTOM_TEMP_OFFSET_K;  // Bottom temperature estimate
    Mdata.hs = current_snow_depth;                            // Current snow height [m]

    // Set roughness length and friction velocity for wind pumping calculations
    // Without these, wind pumping produces NaN thermal conductivity
    double roughness_length = 0.002;  // Default fallback value

    if (config) {
        config->getValue("ROUGHNESS_LENGTH", "Snowpack", roughness_length, mio::IOUtils::nothrow);
        if (roughness_length < 0.0) {
            roughness_length = 0.002;  // Default if not in config
        }
    }

    // Adapt roughness based on snow presence (CRYOWRF pattern from Coupler.cpp line 267)
    const double snow_depth_threshold = 0.03;  // 3cm snow depth threshold
    const double rough_len = (current_snow_depth > snow_depth_threshold) ? roughness_length : 0.01;  // Bare soil z0 ~0.01
    Mdata.z0 = rough_len;

    // Compute friction velocity from wind speed using log profile (CRYOWRF pattern from Coupler.cpp)
    // u* = k * u / ln(z/z0), where k=0.4 (von Karman constant)
    const double von_karman = 0.4;  // Physical constant - matches CRYOWRF implementation
    const double z_wind = input.height;   // Measurement height [m]

    Mdata.ustar = (z_wind > Mdata.z0) ?
                  (von_karman * Mdata.vw / std::log(z_wind / Mdata.z0)) :
                  (0.1 * Mdata.vw);  // Fallback if z <= z0

    // Validate critical parameters that could cause crashes
    if (std::isnan(Mdata.ustar) || std::isinf(Mdata.ustar)) {
        printf("SNOWPACK-FATAL [%d,%d]: ustar is NaN/Inf (%.6f)! This will crash SNOWPACK!\n",
               input.i_grid, input.j_grid, Mdata.ustar);
        printf("SNOWPACK-FATAL [%d,%d]: z_wind=%.6f, z0=%.6f, log(z/z0)=%.6f\n",
               input.i_grid, input.j_grid, z_wind, Mdata.z0, std::log(z_wind / Mdata.z0));
        std::abort();
    }
    if (std::isnan(Mdata.z0) || std::isinf(Mdata.z0) || Mdata.z0 <= 0.0) {
        printf("SNOWPACK-FATAL [%d,%d]: z0 is invalid (%.6f)! This will crash SNOWPACK!\n",
               input.i_grid, input.j_grid, Mdata.z0);
        std::abort();
    }
}

void extract_surface_outputs(const SnowStation& station,
                             const SurfaceFluxes& fluxes,
                             const BoundCond& bc,
                             SnowpackOutput& output,
                             double temp_air_fallback) {

    try {
        size_t ndata_size = station.Ndata.size();
        output.surface_temp = (ndata_size > 0) ? station.Ndata.back().T : temp_air_fallback;  // Surface temperature [K]
        output.snow_swe = station.swe;                                  // Snow water equivalent [mm]
        output.snow_depth = station.cH - station.Ground;                // Snow height [m] (CRYOWRF line 1129)
        output.heat_flux_sensible = -1.0 * fluxes.qs;                   // Negative sign for WRF convention (CRYOWRF line 1123)
        output.heat_flux_latent = -1.0 * bc.ql;                        // Use boundary condition data (CRYOWRF line 1124)
        output.albedo = station.Albedo;
        output.snow_coverage = 1.0;  // Hardcoded to 1.0 following CRYOWRF pattern (line 1127)
        output.friction_velocity = 0.0;  // Will be set from meteo data
        output.stability_param = 0.0;    // Will be calculated if needed

        // Consistency check
        if (output.snow_depth > 0.001 && output.snow_swe <= 0.0) {
            output.snow_swe = output.snow_depth * SnowpackConstants::SNOW_DENSITY_FALLBACK;  // Use fallback density
        }

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL: Exception in extract_surface_outputs: %s\n", e.what());
        std::abort();
    }
}

void extract_layer_data(const SnowStation& station, SnowpackLayerData& layer_data) {
    try {
        size_t num_elements = station.getNumberOfElements();
        layer_data.n_layers = static_cast<int>(num_elements);

        size_t layers_to_extract = std::min(num_elements, size_t(100));

        for (size_t i = 0; i < layers_to_extract; i++) {
            const ElementData& elem = station.Edata[i];

            layer_data.layer_temp[i] = elem.Te;
            layer_data.layer_thick[i] = elem.L;
            layer_data.layer_vol_ice[i] = elem.theta[0] * 100.0;    // ICE index = 0
            layer_data.layer_vol_water[i] = elem.theta[1] * 100.0; // WATER index = 1
            layer_data.layer_vol_air[i] = elem.theta[2] * 100.0;   // AIR index = 2

            layer_data.layer_grain_radius[i] = elem.rg;
            layer_data.layer_bond_radius[i] = elem.rb;
            layer_data.layer_dendricity[i] = elem.dd;
            layer_data.layer_sphericity[i] = elem.sp;
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
        }

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL: Exception in extract_layer_data: %s\n", e.what());
        layer_data.n_layers = 0;
        std::abort();
    }
}

void extract_budget_data(const SurfaceFluxes& fluxes,
                        const BoundCond& bc,
                        double cumu_precip,
                        const MeteoInput& input,
                        BudgetData& budget_data) {
    try {
        // Mass budgets (following CRYOWRF patterns) - SAFE ACCESS PATTERN
        budget_data.mass_precip = cumu_precip;                           // Cumulative precipitation [kg/m²]
        budget_data.mass_sublim = fluxes.mass[SurfaceFluxes::MS_SUBLIMATION];  // Sublimation [kg/m²] (CRYOWRF line 1136)
        budget_data.mass_melt = fluxes.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF]; // Melt runoff [kg/m²] - SNOWPACK v11.08 API
        budget_data.mass_swe = fluxes.mass[SurfaceFluxes::MS_SWE];        // Current SWE [kg/m²]

        // Calculate refreeze mass from ice base melting/freezing or negative runoff
        double ice_base_meltfreeze = fluxes.mass[SurfaceFluxes::MS_ICEBASE_MELTING_FREEZING];
        budget_data.mass_refreeze = (ice_base_meltfreeze > 0.0) ? ice_base_meltfreeze : 0.0;  // Positive values = freezing

        // Energy budgets
        budget_data.energy_lw_in = input.longwave_in;                           // Incoming longwave [W/m²]
        budget_data.energy_lw_out = fluxes.lw_out;                             // Outgoing longwave [W/m²] (CRYOWRF line 1087)
        budget_data.energy_sw_in = input.shortwave_in;                          // Incoming shortwave [W/m²]
        budget_data.energy_sw_out = input.shortwave_in * (1.0) - fluxes.sw_in; // Reflected shortwave [W/m²]
        budget_data.energy_sensible = fluxes.qs;                               // Sensible heat flux [W/m²]
        budget_data.energy_latent = bc.ql;                                     // Latent heat flux [W/m²]
        budget_data.energy_ground_flux = fluxes.qg;                            // Ground heat flux
        budget_data.energy_rain = fluxes.qr;                                   // Rain heat flux [W/m²] (CRYOWRF pattern)
        budget_data.energy_total = fluxes.dIntEnergy / input.dt;               // Total energy [W/m²] (CRYOWRF pattern)

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL: Exception in extract_budget_data: %s\n", e.what());
        std::abort();
    }
}

// Time management class for better organization
class TimeManager {
private:
    int call_counter_ = 0;
    int compute_counter_ = 0;
    bool first_call_ = true;

public:
    void advance_time(const MeteoInput& input, mio::Date& current_time, double calculation_step_length) {
        initialize_time_step_internal(input, current_time, calculation_step_length,
                                    call_counter_, compute_counter_, first_call_);
    }

    void reset() {
        call_counter_ = 0;
        compute_counter_ = 0;
        first_call_ = true;
    }
};

// Global time manager instance
static TimeManager g_time_manager;

void advance_simulation_time(const MeteoInput& input,
                           mio::Date& current_time,
                           double calculation_step_length) {
    g_time_manager.advance_time(input, current_time, calculation_step_length);
}

} // namespace SnowpackUtils