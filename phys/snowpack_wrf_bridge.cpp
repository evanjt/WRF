/*
 * SNOWPACK-WRF Bridge Implementation (Stateless Version)
 *
 * This file provides the C interface between WRF (Fortran) and SNOWPACK (C++)
 * Following CRYOWRF pattern: stateless calls with temporary objects
 * This eliminates persistent memory issues that cause WRF segfaults
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>

// SNOWPACK v11.08 headers - relative paths from phys/snowpack/
#include "meteoio/meteoio/MeteoIO.h"
#include "snowpack/snowpack/DataClasses.h"
#include "snowpack/snowpack/SnowpackConfig.h"
#include "snowpack/snowpack/snowpackCore/Snowpack.h"

// Global configuration (shared, read-only)
static std::unique_ptr<SnowpackConfig> global_config;
static bool config_initialized = false;

// SNOWPACK Configuration Constants
namespace SnowpackConstants {
  // Timestep configuration
  constexpr double CALCULATION_STEP_MINUTES = 15.0;  // WRF coupling timestep [minutes]
  
  // Meteorological measurement heights [m]
  constexpr int METEO_HEIGHT_METERS = 30;           // Height of meteorological measurements
  constexpr int WIND_HEIGHT_METERS = 30;            // Height of wind measurements
  
  // Surface properties
  constexpr double ROUGHNESS_LENGTH_METERS = 0.01;  // Surface roughness length [m]
  constexpr double GEO_HEAT_FLUX = 0.06;            // Geothermal heat flux [W/m²]
  constexpr int CANOPY_NONE = 0;                     // No canopy model
  
  // Temperature sanity checks [K] - prevents solver instabilities
  constexpr double T_CRAZY_MAX_KELVIN = 400.0;      // Maximum reasonable temperature (127°C)
  constexpr double T_CRAZY_MIN_KELVIN = 100.0;      // Minimum reasonable temperature (-173°C)
  
  // Default station metadata
  const std::string STATION_ID_PREFIX = "WRF_GRID";  // Station ID prefix for SNOWPACK
  constexpr double DEFAULT_LATITUDE = 46.0;          // Default latitude for physics [degrees N]
  constexpr double DEFAULT_LONGITUDE = 8.0;          // Default longitude for physics [degrees E]
}

// Initialize SNOWPACK configuration (once, globally)
void initialize_snowpack_config() {
  if (config_initialized) return;
  
  // Create empty config first, then add keys
  mio::Config base_config;

  // SNOWPACK configuration matching CRYOWRF exactly
  // [Snowpack] section - core settings
  base_config.addKey("CALCULATION_STEP_LENGTH", "Snowpack", std::to_string(SnowpackConstants::CALCULATION_STEP_MINUTES));
  base_config.addKey("METEO_STEP_LENGTH", "Snowpack", std::to_string(SnowpackConstants::CALCULATION_STEP_MINUTES));
  base_config.addKey("MEAS_TSS", "Snowpack", "FALSE");
  base_config.addKey("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack", "FALSE");
  base_config.addKey("SW_MODE", "Snowpack", "INCOMING");
  base_config.addKey("INCOMING_LONGWAVE", "Snowpack", "TRUE");
  base_config.addKey("HEIGHT_OF_WIND_VALUE", "Snowpack", std::to_string(SnowpackConstants::WIND_HEIGHT_METERS));
  base_config.addKey("HEIGHT_OF_METEO_VALUES", "Snowpack", std::to_string(SnowpackConstants::METEO_HEIGHT_METERS));
  base_config.addKey("ATMOSPHERIC_STABILITY", "Snowpack", "MO_HOLTSLAG");
  base_config.addKey("ROUGHNESS_LENGTH", "Snowpack", std::to_string(SnowpackConstants::ROUGHNESS_LENGTH_METERS));
  base_config.addKey("NUMBER_SLOPES", "Snowpack", "1");
  base_config.addKey("CHANGE_BC", "Snowpack", "FALSE");
  base_config.addKey("SNP_SOIL", "Snowpack", "FALSE");
  base_config.addKey("SOIL_FLUX", "Snowpack", "TRUE");
  base_config.addKey("GEO_HEAT", "Snowpack", std::to_string(SnowpackConstants::GEO_HEAT_FLUX));
  base_config.addKey("CANOPY", "Snowpack", std::to_string(SnowpackConstants::CANOPY_NONE));
  base_config.addKey("FORCING", "Snowpack", "ATMOS");
  
  // [SnowpackAdvanced] section - advanced settings
  base_config.addKey("VARIANT", "SnowpackAdvanced", "DEFAULT");
  base_config.addKey("RESEARCH_MODE", "SnowpackAdvanced", "TRUE");
  base_config.addKey("ALLOW_ADAPTIVE_TIMESTEPPING", "SnowpackAdvanced", "TRUE");
  base_config.addKey("SNOW_EROSION", "SnowpackAdvanced", "FALSE");
  base_config.addKey("DETECT_GRASS", "SnowpackAdvanced", "TRUE");
  base_config.addKey("HN_DENSITY", "SnowpackAdvanced", "PARAMETERIZED");
  base_config.addKey("HN_DENSITY_PARAMETERIZATION", "SnowpackAdvanced", "ZWART");
  base_config.addKey("AVG_METHOD_HYDRAULIC_CONDUCTIVITY", "SnowpackAdvanced", "ARITHMETICMEAN");
  base_config.addKey("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced", "RICHARDSEQUATION");
  base_config.addKey("WATERTRANSPORTMODEL_SOIL", "SnowpackAdvanced", "RICHARDSEQUATION");
  base_config.addKey("ENABLE_VAPOUR_TRANSPORT", "SnowpackAdvanced", "FALSE");

  // Create SnowpackConfig from mio::Config
  global_config = std::make_unique<SnowpackConfig>(base_config);
  config_initialized = true;

  printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: ✓ Timestep: %.1f min, SNP_SOIL: FALSE\n", 
         SnowpackConstants::CALCULATION_STEP_MINUTES);
}

// C interface for Fortran binding
extern "C" {

void snowpack_physics(double temp_air, double humidity, double wind_speed, double wind_dir,
                      double shortwave_in, double longwave_in, double precipitation, double pressure, double height, double dt,
                      int i_grid, int j_grid,
                      double* snow_swe, double* snow_depth, double* surface_temp,
                      double* heat_flux_sensible, double* heat_flux_latent, double* albedo, double* snow_coverage) {
  
  // Track physics calls for debugging
  static int call_count = 0;
  call_count++;
  
  // Periodic progress reporting - much less verbose
  if (call_count <= 5 || (call_count % 1000 == 0)) {  // First 5 calls, then every 1000th
    printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: Physics call #%d - Grid (%d,%d) - Stateless call\n", 
           call_count, i_grid, j_grid);
  }
  
  // Initialize configuration on first call (shared, read-only)
  try {
    initialize_snowpack_config();
  } catch (const std::exception& e) {
    printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: ❌ Configuration failed: %s\n", e.what());
    printf("SNOWPACK-FATAL: Aborting WRF run due to SNOWPACK configuration failure\n");
    std::abort();  // Abort instead of silent fallback
  }

  // CRYOWRF-style stateless approach: Create temporary objects for each call
  try {
    // Create temporary SNOWPACK instance (following CRYOWRF Coupler pattern)
    Snowpack snowpack_instance(*global_config);
    
    // Create temporary SnowStation (no canopy, no soil layers, no Alpine3D, no sea ice)
    SnowStation snow_station(false, false, false, false);
    
    // Create temporary meteorological data
    CurrentMeteo Mdata;
    SurfaceFluxes surfFluxes;
    BoundCond sn_Bdata;
    
    // Initialize snow station for this call
    SN_SNOWSOIL_DATA ssdata;
    ssdata.meta.stationID = SnowpackConstants::STATION_ID_PREFIX + "_" + 
                           std::to_string(i_grid) + "_" + std::to_string(j_grid);
    
    // Set position with default coordinates (SNOWPACK needs valid lat/lon)
    ssdata.meta.position.setLatLon(SnowpackConstants::DEFAULT_LATITUDE, 
                                   SnowpackConstants::DEFAULT_LONGITUDE,
                                   height);
    
    ssdata.Height = height;
    ssdata.nN = 1;       // Start with 1 node (ground only)  
    ssdata.nLayers = 0;  // No snow/soil layers initially (SNP_SOIL = FALSE)
    
    // Initialize surface properties
    ssdata.Albedo = 0.85;   // Default snow albedo
    ssdata.SoilAlb = 0.2;   // Default soil albedo
    ssdata.BareSoil_z0 = SnowpackConstants::ROUGHNESS_LENGTH_METERS;
    ssdata.HS_last = 0.0;   // No previous snow height
    ssdata.Ldata.clear();   // No layer data (matches nLayers = 0)
    
    // Initialize station
    snow_station.initialize(ssdata, 0);  // Initialize with sector 0
    
    // Set up current meteorology
    mio::Date current_time(2010, 7, 16, 0, 0, 0.0, 0.0);  // Dummy date with timezone
    
    // Temperature sanity check
    double safe_temp = std::max(SnowpackConstants::T_CRAZY_MIN_KELVIN, 
                               std::min(temp_air, SnowpackConstants::T_CRAZY_MAX_KELVIN));
    
    // Fill meteorological data structure using correct SNOWPACK API
    Mdata.date = current_time;
    Mdata.ta = safe_temp;                                       // Air temperature [K]
    Mdata.rh = std::max(0.01, std::min(1.0, humidity));       // Relative humidity [0-1]
    Mdata.vw = std::max(0.1, wind_speed);                     // Wind speed [m/s] 
    Mdata.dw = wind_dir;                                       // Wind direction [degrees]
    Mdata.iswr = std::max(0.0, shortwave_in);                 // Incoming shortwave [W/m²]
    Mdata.lw_net = std::max(0.0, longwave_in);               // Net longwave radiation [W/m²] 
    Mdata.psum = std::max(0.0, precipitation);                // Precipitation [mm]
    // Note: CurrentMeteo has no pressure field - pressure used internally by SNOWPACK
    
    // Additional required meteorological parameters
    Mdata.psum_ph = (safe_temp < 273.65) ? 0.0 : 1.0;        // Precipitation phase (0=snow, 1=rain)
    Mdata.tss = mio::IOUtils::nodata;                         // Surface temperature (let SNOWPACK compute)
    Mdata.ts0 = safe_temp - 5.0;                              // Bottom temperature estimate
    Mdata.hs = *snow_depth;                                   // Current snow height [m]
    
    // Run SNOWPACK model (temporary objects will auto-destruct)
    if (call_count <= 3) {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Calling SNOWPACK model for grid (%d,%d), precip=%.3f\n", 
             i_grid, j_grid, precipitation);
    }
    
    // Execute SNOWPACK physics (correct API with cumulative precipitation parameter)
    double cumu_precip = 0.0;  // Cumulative precipitation parameter  
    snowpack_instance.runSnowpackModel(Mdata, snow_station, cumu_precip, sn_Bdata, surfFluxes);
    
    if (call_count <= 3) {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: SNOWPACK model completed successfully for grid (%d,%d)\n", 
             i_grid, j_grid);
    }
    
    // Extract results from SNOWPACK (using correct member names)
    *surface_temp = (snow_station.Ndata.size() > 0) ? snow_station.Ndata.back().T : temp_air;  // Surface temperature [K]
    *snow_swe = snow_station.swe;                                  // Snow water equivalent [mm]
    *snow_depth = snow_station.cH;                                // Snow height [m]
    *heat_flux_sensible = surfFluxes.qs;                          // Sensible heat flux [W/m²]
    *heat_flux_latent = surfFluxes.ql;                            // Latent heat flux [W/m²]
    *albedo = snow_station.Albedo;                                // Surface albedo [0-1]
    *snow_coverage = (*snow_depth > 0.001) ? 1.0 : 0.0;          // Simple snow coverage [0-1]
    
    // Consistency checks and fallbacks
    if (*snow_depth > 0.001 && *snow_swe <= 0.0) {
        *snow_swe = *snow_depth * 100.0;  // Assume ~100kg/m³ density fallback
    }
    
    if (call_count <= 10) {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: SUCCESS #%d: Grid (%d,%d) T_sfc=%.1fK (%.1f°C), SWE=%.2fmm, depth=%.2fm\n",
             call_count, i_grid, j_grid, *surface_temp, *surface_temp - 273.15, *snow_swe, *snow_depth);
    }
    
    // Temporary objects automatically destroyed here - no memory leaks!
    
  } catch (const std::exception& e) {
    printf("SNOWPACK-ERROR [C++/snowpack_wrf_bridge.cpp]: ❌ Error in grid (%d,%d): %s\n", 
           i_grid, j_grid, e.what());
    printf("SNOWPACK-ERROR: Grid point (%d,%d) - Ta=%.2fK, RH=%.2f, Precip=%.3fmm\n", 
           i_grid, j_grid, temp_air, humidity, precipitation);
    printf("SNOWPACK-FATAL: Aborting WRF run due to SNOWPACK physics failure\n");
    std::abort();  // Abort instead of silent fallback
  }
}

// Enhanced interface with detailed layer extraction
void snowpack_physics_layers(double temp_air, double humidity, double wind_speed, double wind_dir,
                             double shortwave_in, double longwave_in, double precipitation, double pressure, double height, double dt,
                             int i_grid, int j_grid,
                             double* snow_swe, double* snow_depth, double* surface_temp,
                             double* heat_flux_sensible, double* heat_flux_latent, double* albedo, double* snow_coverage,
                             // Layer arrays (max 50 layers)
                             int* n_layers,
                             double* layer_temp, double* layer_thick, 
                             double* layer_vol_ice, double* layer_vol_water, double* layer_vol_air,
                             double* layer_grain_radius, double* layer_bond_radius,
                             double* layer_dendricity, double* layer_sphericity,
                             // Budget tracking
                             double* mass_precip, double* mass_sublim, double* mass_melt,
                             double* energy_lw_in, double* energy_sensible, double* energy_latent) {
  
  // Track physics calls for debugging
  static int call_count = 0;
  call_count++;
  
  // Initialize configuration on first call (shared, read-only)
  try {
    initialize_snowpack_config();
  } catch (const std::exception& e) {
    printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Configuration failed: %s\n", e.what());
    printf("SNOWPACK-FATAL: Aborting WRF run due to SNOWPACK configuration failure\n");
    std::abort();  // Abort instead of silent fallback
  }

  // Stateless approach: Create temporary objects for each call
  try {
    // Create temporary SNOWPACK instance
    Snowpack snowpack_instance(*global_config);
    
    // Create temporary SnowStation
    SnowStation snow_station(false, false, false, false);
    
    // Create temporary meteorological data
    CurrentMeteo Mdata;
    SurfaceFluxes surfFluxes;
    BoundCond sn_Bdata;
    
    // Initialize snow station for this call
    SN_SNOWSOIL_DATA ssdata;
    ssdata.meta.stationID = SnowpackConstants::STATION_ID_PREFIX + "_" + 
                           std::to_string(i_grid) + "_" + std::to_string(j_grid);
    
    // Set position with default coordinates
    ssdata.meta.position.setLatLon(SnowpackConstants::DEFAULT_LATITUDE, 
                                   SnowpackConstants::DEFAULT_LONGITUDE,
                                   height);
    
    ssdata.Height = height;
    ssdata.nN = 1;       // Start with 1 node (ground only)  
    ssdata.nLayers = 0;  // No snow/soil layers initially
    
    // Initialize surface properties
    ssdata.Albedo = 0.85;   // Default snow albedo
    ssdata.SoilAlb = 0.2;   // Default soil albedo
    ssdata.BareSoil_z0 = SnowpackConstants::ROUGHNESS_LENGTH_METERS;
    ssdata.HS_last = 0.0;   // No previous snow height
    ssdata.Ldata.clear();   // No layer data
    
    // Initialize station
    snow_station.initialize(ssdata, 0);
    
    // Set up current meteorology
    mio::Date current_time(2010, 7, 16, 0, 0, 0.0, 0.0);
    
    // Temperature sanity check
    double safe_temp = std::max(SnowpackConstants::T_CRAZY_MIN_KELVIN, 
                               std::min(temp_air, SnowpackConstants::T_CRAZY_MAX_KELVIN));
    
    // Fill meteorological data structure
    Mdata.date = current_time;
    Mdata.ta = safe_temp;
    Mdata.rh = std::max(0.01, std::min(1.0, humidity));
    Mdata.vw = std::max(0.1, wind_speed);
    Mdata.dw = wind_dir;
    Mdata.iswr = std::max(0.0, shortwave_in);
    Mdata.lw_net = std::max(0.0, longwave_in);
    Mdata.psum = std::max(0.0, precipitation);
    Mdata.psum_ph = (safe_temp < 273.65) ? 0.0 : 1.0;
    Mdata.tss = mio::IOUtils::nodata;
    Mdata.ts0 = safe_temp - 5.0;
    Mdata.hs = *snow_depth;
    
    // Execute SNOWPACK physics
    double cumu_precip = 0.0;
    snowpack_instance.runSnowpackModel(Mdata, snow_station, cumu_precip, sn_Bdata, surfFluxes);
    
    // Extract basic results
    *surface_temp = (snow_station.Ndata.size() > 0) ? snow_station.Ndata.back().T : temp_air;
    *snow_swe = snow_station.swe;
    *snow_depth = snow_station.cH;
    *heat_flux_sensible = surfFluxes.qs;
    *heat_flux_latent = surfFluxes.ql;
    *albedo = snow_station.Albedo;
    *snow_coverage = (*snow_depth > 0.001) ? 1.0 : 0.0;
    
    // Extract detailed layer information from SnowStation
    size_t num_elements = snow_station.getNumberOfElements();
    *n_layers = static_cast<int>(num_elements);
    
    // Limit to max 50 layers for WRF arrays
    size_t layers_to_extract = std::min(num_elements, size_t(50));
    
    for (size_t i = 0; i < layers_to_extract; i++) {
      const ElementData& elem = snow_station.Edata[i];
      
      // Layer properties
      layer_temp[i] = elem.Te;                              // Element temperature [K]
      layer_thick[i] = elem.L;                              // Element thickness [m]
      layer_vol_ice[i] = elem.theta[ICE] * 100.0;          // Ice volume fraction [%]
      layer_vol_water[i] = elem.theta[WATER] * 100.0;      // Water volume fraction [%]  
      layer_vol_air[i] = elem.theta[AIR] * 100.0;          // Air volume fraction [%]
      
      // Grain properties
      layer_grain_radius[i] = elem.rg;                      // Grain radius [mm]
      layer_bond_radius[i] = elem.rb;                       // Bond radius [mm]
      layer_dendricity[i] = elem.dd;                        // Dendricity [-]
      layer_sphericity[i] = elem.sp;                        // Sphericity [-]
    }
    
    // Clear unused layers
    for (size_t i = layers_to_extract; i < 50; i++) {
      layer_temp[i] = 0.0;
      layer_thick[i] = 0.0;
      layer_vol_ice[i] = 0.0;
      layer_vol_water[i] = 0.0;
      layer_vol_air[i] = 0.0;
      layer_grain_radius[i] = 0.0;
      layer_bond_radius[i] = 0.0;
      layer_dendricity[i] = 0.0;
      layer_sphericity[i] = 0.0;
    }
    
    // Extract mass budget information from SurfaceFluxes
    *mass_precip = surfFluxes.mass[SurfaceFluxes::MS_RAIN] + 
                   surfFluxes.mass[SurfaceFluxes::MS_SNOW];     // Total precipitation [kg/m²]
    *mass_sublim = surfFluxes.mass[SurfaceFluxes::MS_SUBLIMATION]; // Sublimation [kg/m²]
    *mass_melt = surfFluxes.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF]; // Melt runoff [kg/m²]
    
    // Extract energy budget information
    *energy_lw_in = Mdata.lw_net;                           // Net LW radiation [W/m²]
    *energy_sensible = surfFluxes.qs;                       // Sensible heat flux [W/m²]
    *energy_latent = surfFluxes.ql;                         // Latent heat flux [W/m²]
    
    if (call_count <= 5) {
      printf("SNOWPACK-LAYERS [C++]: Grid (%d,%d) - %d layers, T_sfc=%.1fK, SWE=%.2fmm, depth=%.2fm\n",
             i_grid, j_grid, *n_layers, *surface_temp, *snow_swe, *snow_depth);
    }
    
  } catch (const std::exception& e) {
    printf("SNOWPACK-ERROR [C++]: Error in grid (%d,%d): %s\n", i_grid, j_grid, e.what());
    printf("SNOWPACK-ERROR: Grid point (%d,%d) - Ta=%.2fK, RH=%.2f, Precip=%.3fmm\n", 
           i_grid, j_grid, temp_air, humidity, precipitation);
    printf("SNOWPACK-FATAL: Aborting WRF run due to SNOWPACK physics failure\n");
    std::abort();  // Abort instead of silent fallback
  }
}

} // extern "C"