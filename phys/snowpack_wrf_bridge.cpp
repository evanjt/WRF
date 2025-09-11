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

  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: === Creating SnowpackConfig ===\n");
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Key sections added: [Snowpack], [SnowpackAdvanced]\n");
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Critical keys configured:\n");
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   ✓ CALCULATION_STEP_LENGTH = %.1f min\n", SnowpackConstants::CALCULATION_STEP_MINUTES);
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   ✓ METEO_STEP_LENGTH = %.1f min\n", SnowpackConstants::CALCULATION_STEP_MINUTES);
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   ✓ T_CRAZY_MAX = %.1f K\n", SnowpackConstants::T_CRAZY_MAX_KELVIN);
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   ✓ T_CRAZY_MIN = %.1f K\n", SnowpackConstants::T_CRAZY_MIN_KELVIN);
  
  try {
    global_config = std::make_unique<SnowpackConfig>(base_config);
    config_initialized = true;
    
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: ✅ SnowpackConfig created successfully!\n");
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: === SNOWPACK Configuration Initialization COMPLETE ===\n");
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
  
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: === SNOWPACK PHYSICS CALL START ===\n");
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Grid point: (%d,%d)\n", i_grid, j_grid);
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Input meteorology:\n");
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   - Temperature: %.2f K (%.2f°C)\n", temp_air, temp_air - 273.15);
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   - Humidity: %.3f\n", humidity);
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   - Wind speed: %.2f m/s\n", wind_speed);
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   - Precipitation: %.3f mm\n", precipitation);
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   - SW radiation: %.2f W/m²\n", shortwave_in);
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   - LW radiation: %.2f W/m²\n", longwave_in);
  
  // Initialize configuration on first call
  try {
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Calling initialize_snowpack_config()...\n");
    initialize_snowpack_config();
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Configuration initialization completed\n");
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
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Creating new SNOWPACK instance for grid (%d,%d)\n", i_grid, j_grid);
    
    try {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Constructing Snowpack object with global_config...\n");
      snowpack_instances[grid_point] = std::make_unique<Snowpack>(*global_config);
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: ✅ Snowpack object created successfully\n");
      
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Creating SnowStation object...\n");
      snow_stations[grid_point] = SnowStation();
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: ✅ SnowStation object created successfully\n");
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
    
    // Initialize basic snow/soil properties
    ssdata.Albedo = 0.85;  // Fresh snow albedo
    ssdata.SoilAlb = 0.2;  // Soil albedo 
    ssdata.BareSoil_z0 = SnowpackConstants::ROUGHNESS_LENGTH_METERS;
    
    // Ensure LayerData is empty (no soil/snow layers, just bare ground)
    ssdata.Ldata.clear();       // No layer data (matches nLayers = 0)

    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Initializing SnowStation with:\n");
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   - Station ID: %s\n", ssdata.meta.stationID.c_str());
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   - Position: %.2f°N, %.2f°E, %.0fm\n", 
           SnowpackConstants::DEFAULT_LATITUDE, SnowpackConstants::DEFAULT_LONGITUDE, height);
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   - Structure: nN=%zu, nLayers=%zu, Ldata.size=%zu\n", 
           ssdata.nN, ssdata.nLayers, ssdata.Ldata.size());
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]:   - Configuration: SNP_SOIL=FALSE (bare ground only)\n");
    
    try {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Calling xdata.initialize()...\n");
      xdata.initialize(ssdata, 0);  // Initialize with sector 0
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: ✅ SnowStation initialized successfully\n");
    } catch (const mio::IOException& e) {
      printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: ❌ SnowStation initialization failed: %s\n", e.what());
      printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: This is likely a soil layer configuration issue\n");
      *surface_temp = -999.0;
      *snow_swe = -999.0;
      *snow_depth = -999.0;
      return;
    } catch (const std::exception& e) {
      printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: ❌ SnowStation initialization failed: %s\n", e.what());
      *surface_temp = -999.0;
      *snow_swe = -999.0;
      *snow_depth = -999.0;
      return;
    }

    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: ✅ Grid point (%d,%d) fully initialized\n", i_grid, j_grid);
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

    // Extract results for WRF
    *surface_temp =
        mdata.tss;  // Surface temperature from SNOWPACK (in CurrentMeteo)
    *snow_swe = xdata.swe;   // Snow water equivalent [mm] (direct member)
    *snow_depth = xdata.cH;  // Total calculated height [m] (direct member)

    // Energy fluxes from SNOWPACK
    *heat_flux_sensible = sdata.qs;  // Sensible heat flux [W/m²]
    *heat_flux_latent = sdata.ql;    // Latent heat flux [W/m²]

    // Surface properties
    *albedo = xdata.Albedo;                              // Snow albedo [-]
    *snow_coverage = (*snow_depth > 0.001) ? 1.0 : 0.0;  // Snow coverage [-]

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
  
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: === SNOWPACK PHYSICS CALL SUCCESS ===\n");
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Grid (%d,%d) completed successfully\n", i_grid, j_grid);
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Output: T_sfc=%.2fK, SWE=%.3fmm, depth=%.3fm\n", 
         *surface_temp, *snow_swe, *snow_depth);
}

}  // extern "C"
