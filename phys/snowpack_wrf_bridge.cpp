/*
 * SNOWPACK-WRF Bridge Implementation (Stateless Version)
 *
 * This file provides the C interface between WRF (Fortran) and SNOWPACK (C++)
 * Following CRYOWRF pattern: stateless calls with temporary objects
 * This eliminates persistent memory issues that cause WRF segfaults
 */

#include <cmath>
#include <cstdio>
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
    *surface_temp = temp_air;  // Fallback to air temperature
    *snow_swe = 0.0;
    *snow_depth = 0.0;
    return;
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
    mio::Date current_time;  // Default time is fine for physics
    current_time.setFromComponents(2010, 7, 16, 0, 0, 0.0);  // Dummy date
    
    // Temperature sanity check
    double safe_temp = std::max(SnowpackConstants::T_CRAZY_MIN_KELVIN, 
                               std::min(temp_air, SnowpackConstants::T_CRAZY_MAX_KELVIN));
    
    // Fill meteorological data structure
    Mdata.date = current_time;
    Mdata(MeteoData::TA) = safe_temp;          // Air temperature [K]
    Mdata(MeteoData::RH) = std::max(0.01, std::min(1.0, humidity));  // Relative humidity [0-1]
    Mdata(MeteoData::VW) = std::max(0.1, wind_speed);                // Wind speed [m/s] 
    Mdata(MeteoData::DW) = wind_dir;                                 // Wind direction [degrees]
    Mdata(MeteoData::ISWR) = std::max(0.0, shortwave_in);           // Incoming shortwave [W/m²]
    Mdata(MeteoData::ILWR) = std::max(0.0, longwave_in);            // Incoming longwave [W/m²] 
    Mdata(MeteoData::PSUM) = std::max(0.0, precipitation);          // Precipitation [mm]
    Mdata(MeteoData::P) = pressure;                                  // Pressure [Pa]
    
    // Additional required meteorological parameters
    Mdata(MeteoData::PSUM_PH) = (safe_temp < 273.65) ? 0.0 : 1.0;  // Precipitation phase (0=snow, 1=rain)
    Mdata(MeteoData::TSS) = mio::IOUtils::nodata;                   // Surface temperature (let SNOWPACK compute)
    Mdata(MeteoData::TSG) = safe_temp - 5.0;                        // Ground temperature estimate
    Mdata(MeteoData::HS) = *snow_depth;                             // Current snow height [m]
    
    // Run SNOWPACK model (temporary objects will auto-destruct)
    if (call_count <= 3) {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Calling SNOWPACK model for grid (%d,%d), precip=%.3f\n", 
             i_grid, j_grid, precipitation);
    }
    
    // Execute SNOWPACK physics
    snowpack_instance.runSnowpackModel(Mdata, snow_station, sn_Bdata, surfFluxes);
    
    if (call_count <= 3) {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: SNOWPACK model completed successfully for grid (%d,%d)\n", 
             i_grid, j_grid);
    }
    
    // Extract results from SNOWPACK
    *surface_temp = snow_station.Sdata.front().T;                    // Surface temperature [K]
    *snow_swe = snow_station.SWE;                                    // Snow water equivalent [mm]
    *snow_depth = snow_station.getHS();                             // Snow height [m]
    *heat_flux_sensible = surfFluxes.qs;                           // Sensible heat flux [W/m²]
    *heat_flux_latent = surfFluxes.ql;                             // Latent heat flux [W/m²]
    *albedo = snow_station.Albedo;                                 // Surface albedo [0-1]
    *snow_coverage = (*snow_depth > 0.001) ? 1.0 : 0.0;           // Simple snow coverage [0-1]
    
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
    
    // Return safe fallback values
    *surface_temp = temp_air;  
    *snow_swe = 0.0;
    *snow_depth = 0.0;
    *heat_flux_sensible = 0.0;
    *heat_flux_latent = 0.0;
    *albedo = 0.2;  // Bare soil albedo
    *snow_coverage = 0.0;
  }
}

} // extern "C"