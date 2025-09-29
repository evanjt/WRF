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
#include <cstring>
#include <vector>

// SNOWPACK v11.08 headers - relative paths from phys/snowpack/
#include "meteoio/meteoio/MeteoIO.h"
#include "snowpack/snowpack/DataClasses.h"
#include "snowpack/snowpack/SnowpackConfig.h"
#include "snowpack/snowpack/snowpackCore/Snowpack.h"
#include "snowpack/snowpack/plugins/SnowpackIO.h"

// Configuration management
class SnowpackConfigManager {
public:
    static mio::Config loadConfiguration(const std::string& ini_file_path);
    static void validateConfiguration(const mio::Config& cfg);
    static std::string getDefaultConfigPath();
};

// Global configuration and I/O (shared, read-only)
static std::unique_ptr<SnowpackConfig> global_config;
static std::unique_ptr<SnowpackIO> global_snowpack_io;
static bool config_initialized = false;
static std::string config_file_path;
static bool use_state_persistence = false;  // Enable .sno file persistence

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

// SnowpackConfigManager implementation
mio::Config SnowpackConfigManager::loadConfiguration(const std::string& ini_file_path) {
    try {
        // Load configuration from file
        mio::Config config(ini_file_path);
        printf("SNOWPACK-INFO [C++/SnowpackConfigManager]: Loaded configuration from %s\n", ini_file_path.c_str());
        return config;
    } catch (const std::exception& e) {
        printf("SNOWPACK-ERROR [C++/SnowpackConfigManager]: Failed to load %s: %s\n", ini_file_path.c_str(), e.what());
        throw;
    }
}

void SnowpackConfigManager::validateConfiguration(const mio::Config& cfg) {
    // Check for essential SNOWPACK parameters
    std::vector<std::pair<std::string, std::string>> required_params = {
        {"CALCULATION_STEP_LENGTH", "Snowpack"},
        {"FORCING", "Snowpack"},
        {"SNP_SOIL", "Snowpack"},
        {"SOIL_FLUX", "Snowpack"},
        {"VARIANT", "SnowpackAdvanced"}
    };
    
    for (const auto& param : required_params) {
        try {
            std::string value;
            cfg.getValue(param.first, param.second, value);
            printf("SNOWPACK-VALIDATE [C++]: %s::%s = %s\n", param.second.c_str(), param.first.c_str(), value.c_str());
        } catch (const std::exception& e) {
            printf("SNOWPACK-ERROR [C++/SnowpackConfigManager]: Missing required parameter %s::%s\n", 
                   param.second.c_str(), param.first.c_str());
            throw std::runtime_error("Configuration validation failed: missing " + param.first);
        }
    }
    
    printf("SNOWPACK-INFO [C++/SnowpackConfigManager]: Configuration validation passed\n");
}

std::string SnowpackConfigManager::getDefaultConfigPath() {
    return "./io.ini";  // Default path in WRF run directory
}

// Helper functions for SNOWPACK state persistence
std::string generateSnoFilename(int i, int j, int domain = 1) {
    // Follow CRYOWRF naming convention: snpack_domain_j_i.sno
    return "./snowpack_states/snpack_" + std::to_string(domain) + "_" + 
           std::to_string(j) + "_" + std::to_string(i) + ".sno";
}

std::string generateStationID(int i, int j, int domain = 1) {
    // Follow CRYOWRF station ID convention
    return "snpack_" + std::to_string(domain) + "_" + std::to_string(j) + "_" + std::to_string(i);
}

// Function declarations
void initialize_snowpack_config_with_path(const std::string& ini_path);

// Initialize SNOWPACK configuration (once, globally)
void initialize_snowpack_config() {
    initialize_snowpack_config_with_path(SnowpackConfigManager::getDefaultConfigPath());
}

// Initialize SNOWPACK configuration with specific file path
void initialize_snowpack_config_with_path(const std::string& ini_path) {
    if (config_initialized) return;
    
    try {
        // Load configuration from file
        mio::Config file_config = SnowpackConfigManager::loadConfiguration(ini_path);
        
        // Validate configuration
        SnowpackConfigManager::validateConfiguration(file_config);
        
        // Create SnowpackConfig from file
        global_config = std::make_unique<SnowpackConfig>(file_config);
        
        // Create SnowpackIO instance for state persistence
        global_snowpack_io = std::make_unique<SnowpackIO>(*global_config);
        
        config_file_path = ini_path;
        config_initialized = true;
        
        // Extract and report key settings
        std::string calc_step, snp_soil;
        file_config.getValue("CALCULATION_STEP_LENGTH", "Snowpack", calc_step);
        file_config.getValue("SNP_SOIL", "Snowpack", snp_soil);
        
        printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: Configured from %s - Timestep: %s min, SNP_SOIL: %s\n", 
               ini_path.c_str(), calc_step.c_str(), snp_soil.c_str());
               
    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Configuration failed for %s: %s\n", 
               ini_path.c_str(), e.what());
        printf("SNOWPACK-FATAL: Unable to load SNOWPACK configuration - WRF run will abort\n");
        throw;
    }
}

// C interface for Fortran binding
extern "C" {

// Initialize configuration with specific path (called from Fortran)
void initialize_snowpack_config_c(const char* ini_file_path) {
    std::string path_str(ini_file_path);
    initialize_snowpack_config_with_path(path_str);
}

// Get current configuration file path
void get_snowpack_config_path_c(char* path_buffer, int buffer_size) {
    if (config_initialized) {
        strncpy(path_buffer, config_file_path.c_str(), buffer_size - 1);
        path_buffer[buffer_size - 1] = '\0';
    } else {
        strncpy(path_buffer, "NOT_INITIALIZED", buffer_size - 1);
        path_buffer[buffer_size - 1] = '\0';
    }
}

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
    printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Configuration failed: %s\n", e.what());
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
    printf("SNOWPACK-ERROR [C++/snowpack_wrf_bridge.cpp]: Error in grid (%d,%d): %s\n", 
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
                   surfFluxes.mass[SurfaceFluxes::MS_HNW];      // Total precipitation [kg/m²] (rain + solid)
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

// Extract detailed layer data from SNOWPACK for CRYOWRF compatibility
void extract_snowpack_layers_c(SnowpackInterface* interface_data,
                               float layer_temps[50], float layer_thick[50],
                               float layer_voli[50], float layer_volw[50], float layer_volv[50],
                               float layer_rg[50], float layer_rb[50], 
                               float layer_dd[50], float layer_sp[50],
                               int* n_layers) {
  try {
    // Initialize all arrays to zero
    for (int i = 0; i < 50; i++) {
      layer_temps[i] = 0.0f;
      layer_thick[i] = 0.0f;
      layer_voli[i] = 0.0f;
      layer_volw[i] = 0.0f;
      layer_volv[i] = 0.0f;
      layer_rg[i] = 0.0f;
      layer_rb[i] = 0.0f;
      layer_dd[i] = 0.0f;
      layer_sp[i] = 1.0f;  // Default sphericity
    }
    
    // Access SnowStation from the interface
    // Note: This is simplified - in reality we'd need to access the actual
    // SNOWPACK SnowStation object from the interface_data structure
    // For now, populate with placeholder values based on surface conditions
    
    if (interface_data->surface.snow_depth > 0.001) {
      // Simple snow layers based on snow depth
      *n_layers = std::min(10, (int)(interface_data->surface.snow_depth * 100));  // 1 layer per cm
      
      for (int i = 0; i < *n_layers; i++) {
        // Simple layer properties (to be replaced with actual SNOWPACK extraction)
        layer_temps[i] = std::min(273.15f, interface_data->surface.surface_temp);  // Snow temp <= 0°C
        layer_thick[i] = interface_data->surface.snow_depth / *n_layers;           // Equal thickness layers
        layer_voli[i] = 0.3f;     // 30% ice volume fraction
        layer_volw[i] = 0.0f;     // No liquid water initially
        layer_volv[i] = 0.7f;     // 70% air volume fraction
        layer_rg[i] = 0.4f;       // Grain radius [mm]
        layer_rb[i] = 0.1f;       // Bond radius [mm]
        layer_dd[i] = 0.0f;       // Dendricity
        layer_sp[i] = 1.0f;       // Sphericity
      }
    } else {
      *n_layers = 0;
    }
    
  } catch (const std::exception& e) {
    printf("SNOWPACK-ERROR: Layer extraction failed: %s\n", e.what());
    *n_layers = 0;
  }
}

} // extern "C"