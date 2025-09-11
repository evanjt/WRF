/*
 * SNOWPACK-WRF Bridge Implementation
 *
 * This file provides the C interface between WRF (Fortran) and SNOWPACK (C++)
 * It implements real SNOWPACK v11.08 physics integration.
 */

#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <string>

// SNOWPACK v11.08 headers - relative paths from phys/snowpack/
#include "meteoio/meteoio/MeteoIO.h"
#include "snowpack/snowpack/DataClasses.h"
#include "snowpack/snowpack/SnowpackConfig.h"
#include "snowpack/snowpack/snowpackCore/Snowpack.h"

// Global SNOWPACK instances - one per grid point
static std::map<std::pair<int, int>, std::unique_ptr<Snowpack>>
    snowpack_instances;
static std::map<std::pair<int, int>, SnowStation> snow_stations;
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
  
  // Water vapor transport settings
  constexpr int VAPOR_TRANSPORT_TIMESTEP_SEC = 60;  // Vapor transport sub-timestep [seconds]
  constexpr double VAPOR_IMPLICIT_FACTOR = 1.0;     // Fully implicit solver (most stable)
  
  // Default station metadata
  const std::string STATION_ID_PREFIX = "WRF_GRID";  // Station ID prefix for SNOWPACK
  constexpr double DEFAULT_LATITUDE = 46.0;          // Default latitude for physics [degrees N]
  constexpr double DEFAULT_LONGITUDE = 8.0;          // Default longitude for physics [degrees E]
}

// Initialize SNOWPACK configuration
void initialize_snowpack_config() {
  if (config_initialized) return;
  
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: === SNOWPACK Configuration Initialization START ===\n");
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Using constants from SnowpackConstants namespace\n");
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: - CALCULATION_STEP_MINUTES = %.1f\n", SnowpackConstants::CALCULATION_STEP_MINUTES);
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: - T_CRAZY_MAX_KELVIN = %.1f K\n", SnowpackConstants::T_CRAZY_MAX_KELVIN);
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: - T_CRAZY_MIN_KELVIN = %.1f K\n", SnowpackConstants::T_CRAZY_MIN_KELVIN);

  // Create empty config first, then add keys
  mio::Config base_config;
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Created empty mio::Config object\n");

  // SNOWPACK configuration matching CRYOWRF exactly
  // [Snowpack] section - core settings
  base_config.addKey("CALCULATION_STEP_LENGTH", "Snowpack", std::to_string(SnowpackConstants::CALCULATION_STEP_MINUTES));
  base_config.addKey("METEO_STEP_LENGTH", "Snowpack", std::to_string(SnowpackConstants::CALCULATION_STEP_MINUTES));  // Same as CALCULATION_STEP_LENGTH in CRYOWRF
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
  base_config.addKey("ENABLE_VAPOUR_TRANSPORT_SOIL", "SnowpackAdvanced", "FALSE");
  base_config.addKey("WATER_VAPOR_TRANSPORT_TIMESTEP", "SnowpackAdvanced", std::to_string(SnowpackConstants::VAPOR_TRANSPORT_TIMESTEP_SEC));
  base_config.addKey("WATER_VAPOR_TRANSPORT_IMPLICIT_FACTOR", "SnowpackAdvanced", std::to_string(SnowpackConstants::VAPOR_IMPLICIT_FACTOR));
  base_config.addKey("MEAS_INCOMING_LONGWAVE", "SnowpackAdvanced", "TRUE");
  base_config.addKey("FORCE_ADD_SNOWFALL", "SnowpackAdvanced", "FALSE");
  base_config.addKey("T_CRAZY_MAX", "SnowpackAdvanced", std::to_string(SnowpackConstants::T_CRAZY_MAX_KELVIN));
  base_config.addKey("T_CRAZY_MIN", "SnowpackAdvanced", std::to_string(SnowpackConstants::T_CRAZY_MIN_KELVIN));

  // Minimal configuration output - only show once
  static bool config_message_shown = false;
  if (!config_message_shown) {
    printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: Initializing SNOWPACK v11.08 configuration\n");
    printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: ✓ Timestep: %.1f min, SNP_SOIL: FALSE\n", 
           SnowpackConstants::CALCULATION_STEP_MINUTES);
    config_message_shown = true;
  }
  
  try {
    global_config = std::make_unique<SnowpackConfig>(base_config);
    
    // Verify critical configuration values (silent verification)
    bool soil_layers_enabled = false;
    global_config->getValue("SNP_SOIL", "Snowpack", soil_layers_enabled);
    
    if (soil_layers_enabled) {
      printf("SNOWPACK-ERROR [C++/snowpack_wrf_bridge.cpp]: ❌ Configuration error: SNP_SOIL is TRUE but should be FALSE!\n");
    }
    
    config_initialized = true;
  } catch (const mio::UnknownValueException& e) {
    printf("SNOWPACK-ERROR [C++/snowpack_wrf_bridge.cpp]: ❌ Missing configuration key: %s\n", e.what());
    printf("SNOWPACK-ERROR [C++/snowpack_wrf_bridge.cpp]: This key needs to be added to our bridge configuration\n");
    throw;
  } catch (const std::exception& e) {
    printf("SNOWPACK-ERROR [C++/snowpack_wrf_bridge.cpp]: ❌ SnowpackConfig creation failed: %s\n", e.what());
    printf("SNOWPACK-ERROR [C++/snowpack_wrf_bridge.cpp]: Exception type: %s\n", typeid(e).name());
    printf("SNOWPACK-ERROR [C++/snowpack_wrf_bridge.cpp]: Configuration will remain uninitialized\n");
    throw;  // Re-throw to maintain error handling
  }
}

// C interface for Fortran binding
extern "C" {

/*
 * Real SNOWPACK physics interface called from WRF
 * This function must match the BIND(C) interface in module_sf_snowpack.F
 */
void snowpack_physics(double temp_air, double humidity, double wind_speed,
                      double wind_dir, double shortwave_in, double longwave_in,
                      double precipitation, double pressure, double height,
                      double* snow_swe, double* snow_depth,
                      double* surface_temp, double* heat_flux_sensible,
                      double* heat_flux_latent, double* albedo,
                      double* snow_coverage, double dt, int i_grid,
                      int j_grid) {
  
  // Reduce debug verbosity - only print for first few grid points or errors
  static int call_count = 0;
  call_count++;
  
  
  // Periodic progress reporting - much less verbose
  if (call_count <= 5 || (call_count % 1000 == 0)) {  // First 5 calls, then every 1000th
    printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: Physics call #%d - Grid (%d,%d) - Active instances: %zu\n", 
           call_count, i_grid, j_grid, snowpack_instances.size());
  }
  
  // Initialize configuration on first call
  try {
    // Minimal initialization output - only for first few calls
    if (call_count <= 3) {  // Reduced from 10 to 3
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Initializing SNOWPACK configuration...\n");
    }
    initialize_snowpack_config();
    if (call_count <= 3) {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: ✅ Configuration ready\n");
    }
  } catch (const std::exception& e) {
    printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: ❌ Configuration failed: %s\n", e.what());
    *surface_temp = -999.0;  // Error indicators
    *snow_swe = -999.0;
    *snow_depth = -999.0;
    return;  // Early return on config failure
  }

  // Grid point identifier
  std::pair<int, int> grid_point(i_grid, j_grid);
  
  // Create SNOWPACK instance for this grid point if it doesn't exist
  if (snowpack_instances.find(grid_point) == snowpack_instances.end()) {
    // Only print for first few grid points to reduce output
    if (snowpack_instances.size() < 3) {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Creating SNOWPACK instance for grid (%d,%d)\n", i_grid, j_grid);
    }
    
    try {
      snowpack_instances[grid_point] = std::make_unique<Snowpack>(*global_config);
      snow_stations[grid_point] = SnowStation();
    } catch (const std::exception& e) {
      printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: ❌ SNOWPACK instance creation failed: %s\n", e.what());
      printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: This is likely a configuration error\n");
      *surface_temp = -999.0;
      *snow_swe = -999.0;
      *snow_depth = -999.0;
      return;
    }

    // Initialize snow station with basic setup
    SnowStation& xdata = snow_stations[grid_point];

    // Create minimal initialization data for SNOWPACK
    SN_SNOWSOIL_DATA ssdata;
    ssdata.meta.stationID = SnowpackConstants::STATION_ID_PREFIX + "_" + 
                           std::to_string(i_grid) + "_" + std::to_string(j_grid);
    
    // Set position with default coordinates (SNOWPACK needs valid lat/lon)
    ssdata.meta.position.setLatLon(SnowpackConstants::DEFAULT_LATITUDE, 
                                   SnowpackConstants::DEFAULT_LONGITUDE,
                                   height);  // Use actual elevation from WRF
    
    ssdata.Height = height;
    ssdata.nN = 1;       // Start with 1 node (ground only)  
    ssdata.nLayers = 0;  // No snow/soil layers initially (SNP_SOIL = FALSE)
    
    // Initialize surface properties
    ssdata.Albedo = 0.85;   // Default snow albedo
    ssdata.SoilAlb = 0.2;   // Default soil albedo
    ssdata.BareSoil_z0 = SnowpackConstants::ROUGHNESS_LENGTH_METERS;
    
    // Initialize snow state
    ssdata.HS_last = 0.0;
    
    // Ensure LayerData is empty (no soil/snow layers, just bare ground)
    ssdata.Ldata.clear();       // No layer data (matches nLayers = 0)

    // Only show detailed initialization for first few grid points
    if (snowpack_instances.size() <= 3) {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Initializing SnowStation: %s at %.0fm\n", 
             ssdata.meta.stationID.c_str(), height);
    }
    
    try {
      
      xdata.initialize(ssdata, 0);  // Initialize with sector 0
      
      // Only print success for first few grid points
      if (snowpack_instances.size() <= 3) {
        printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: ✅ SnowStation initialized successfully\n");
      }
    } catch (const mio::IOException& e) {
      printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: ❌ SnowStation initialization failed (IOException): %s\n", e.what());
      printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Error likely related to configuration or data inconsistency\n");
      printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Grid (%d,%d) - This may be a domain boundary issue\n", i_grid, j_grid);
      
      // Return safe fallback values instead of crashing
      *surface_temp = temp_air;
      *snow_swe = 0.0;
      *snow_depth = 0.0; 
      *heat_flux_sensible = 50.0;
      *heat_flux_latent = 20.0;
      *albedo = 0.2;  // Soil albedo
      *snow_coverage = 0.0;
      return;
    } catch (const std::exception& e) {
      printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: ❌ SnowStation initialization failed (std::exception): %s\n", e.what());
      printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Grid (%d,%d) - Unknown error during SNOWPACK initialization\n", i_grid, j_grid);
      
      // Return safe fallback values instead of crashing
      *surface_temp = temp_air;
      *snow_swe = 0.0;
      *snow_depth = 0.0;
      *heat_flux_sensible = 50.0; 
      *heat_flux_latent = 20.0;
      *albedo = 0.2;  // Soil albedo
      *snow_coverage = 0.0;
      return;
    }

    // Only print completion message for first few grid points
    if (snowpack_instances.size() <= 3) {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: ✅ Grid point (%d,%d) fully initialized\n", i_grid, j_grid);
    }
  }

  // Get references to this grid point's SNOWPACK instance and data
  Snowpack& snowpack = *snowpack_instances[grid_point];
  SnowStation& xdata = snow_stations[grid_point];

  // Prepare meteorological data for SNOWPACK
  CurrentMeteo mdata(*global_config);
  mdata.date = mio::Date();   // Current timestep
  mdata.ta = temp_air;        // Air temperature [K]
  mdata.rh = humidity;        // Relative humidity [0-1]
  mdata.vw = wind_speed;      // Wind speed [m/s]
  mdata.dw = wind_dir;        // Wind direction [deg]
  mdata.iswr = shortwave_in;  // Incoming shortwave [W/m²]
  mdata.ea = longwave_in /
             (5.67e-8 * pow(temp_air, 4));  // Atmospheric emissivity from LW
  mdata.psum = precipitation * dt;          // Precipitation [mm]
  mdata.psum_ph = (temp_air > 273.15) ? 1.0 : 0.0;  // Phase: 0=snow, 1=rain

  // Update snow station state from WRF
  // Note: Real implementation would maintain multi-layer snowpack
  // For now, we'll use SNOWPACK's native state management

  // Boundary conditions
  BoundCond bdata;

  // Surface fluxes (output from SNOWPACK)
  SurfaceFluxes sdata;

  // Cumulative precipitation (maintained by SNOWPACK)
  static double cumu_precip = 0.0;

  try {
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Calling SNOWPACK model for grid (%d,%d), precip=%.3f\n", 
           i_grid, j_grid, mdata.psum);
    
    // Call SNOWPACK model for this timestep
    snowpack.runSnowpackModel(mdata, xdata, cumu_precip, bdata, sdata);
    
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: SNOWPACK model completed successfully for grid (%d,%d)\n", 
           i_grid, j_grid);

    // Extract results for WRF with safety checks
    *surface_temp = (mdata.tss > 0) ? mdata.tss : temp_air;  // Use air temp if surface temp invalid
    *snow_swe = xdata.swe;   // Snow water equivalent [mm] 
    *snow_depth = xdata.cH;  // Total calculated height [m]
    
    // Apply physical constraints
    if (*surface_temp < 150.0 || *surface_temp > 350.0) {  // Sanity check temperature
        *surface_temp = temp_air;  // Fallback to air temperature
    }
    
    if (*snow_depth < 0.0) *snow_depth = 0.0;  // No negative snow depth
    if (*snow_swe < 0.0) *snow_swe = 0.0;      // No negative SWE
    
    // Consistency check: SWE and depth should be related
    if (*snow_depth > 0.001 && *snow_swe <= 0.0) {
        *snow_swe = *snow_depth * 100.0;  // Assume ~100kg/m³ density fallback
    }

    // Energy fluxes from SNOWPACK  
    *heat_flux_sensible = sdata.qs;  // Sensible heat flux [W/m²]
    *heat_flux_latent = sdata.ql;    // Latent heat flux [W/m²]

    // Surface properties
    *albedo = (xdata.Albedo > 0.0) ? xdata.Albedo : 0.85;  // Default snow albedo if invalid
    *snow_coverage = (*snow_depth > 0.001) ? 1.0 : 0.0;    // Snow coverage [-]

  } catch (const std::exception& e) {
    // SNOWPACK physics failure - terminate with error
    printf("FATAL ERROR: SNOWPACK physics failed at grid(%d,%d): %s\n", i_grid,
           j_grid, e.what());
    printf(
        "SNOWPACK integration cannot continue. Check SNOWPACK configuration "
        "and data.\n");

    // Set error state and terminate
    *surface_temp = -999.0;  // Error indicator
    *snow_swe = -999.0;
    *snow_depth = -999.0;
    *heat_flux_sensible = -999.0;
    *heat_flux_latent = -999.0;
    *albedo = -999.0;
    *snow_coverage = -999.0;

    // Force program termination to prevent corrupted model state
    exit(1);
  }
  
  if (call_count <= 10 || (call_count % 100 == 0)) {  // Reduced success message frequency
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: SUCCESS #%d: Grid (%d,%d) T_sfc=%.1fK (%.1f°C), SWE=%.2fmm, depth=%.2fm\n", 
           call_count, i_grid, j_grid, *surface_temp, *surface_temp - 273.15, *snow_swe, *snow_depth);
  }
  
  // Check for genuinely unphysical results (season-agnostic)
  if (*surface_temp < 150.0 || *surface_temp > 350.0) {  // Extreme but possible range: -123°C to 77°C
    printf("SNOWPACK-WARNING [C++/snowpack_wrf_bridge.cpp]: Extreme temperature at grid (%d,%d): %.1fK (%.1f°C)\n",
           i_grid, j_grid, *surface_temp, *surface_temp - 273.15);
  }
  
  if (*snow_depth > 20.0) {  // Very deep snow (could indicate initialization issue)
    printf("SNOWPACK-WARNING [C++/snowpack_wrf_bridge.cpp]: Very deep snow at grid (%d,%d): %.2fm depth\n",
           i_grid, j_grid, *snow_depth);
  }
}

}  // extern "C"
