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
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>

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

// Global time management (CRYOWRF pattern)
static mio::Date current_simulation_date;  // Current WRF simulation time
static bool time_initialized = false;      // Track initialization
static double calculation_step_length = 0.0;  // Read from SNOWPACK config (minutes)
static bool use_state_persistence = true;   // Enable .sno file persistence (CRYOWRF pattern)

// Persistent SnowStation storage per grid point (CRYOWRF pattern)
// Key format: "i_j" (e.g., "125_67" for grid point i=125, j=67)
static std::map<std::string, std::unique_ptr<SnowStation>> grid_snowstations;
static std::map<std::string, std::unique_ptr<Snowpack>> grid_snowpack_instances;
static bool persistent_objects_initialized = false;

namespace SnowpackConstants {
  // Temperature sanity checks [K] - prevents solver instabilities
  constexpr double T_CRAZY_MAX_KELVIN = 400.0;  // 127°C
  constexpr double T_CRAZY_MIN_KELVIN = 100.0;  // -173°C
  
  // Default station metadata
  const std::string STATION_ID_PREFIX = "WRF_GRID";  // Station ID prefix for SNOWPACK
}

// SnowpackConfigManager implementation
mio::Config SnowpackConfigManager::loadConfiguration(const std::string& ini_file_path) {
    try {
        // Check if file exists first
        std::ifstream test_file(ini_file_path);
        if (!test_file.good()) {
            char cwd[1024];
            getcwd(cwd, sizeof(cwd));
            printf("SNOWPACK-ERROR: io.ini file not found at: %s\n", ini_file_path.c_str());
            printf("SNOWPACK-ERROR: Current working directory: %s\n", cwd);
            throw std::runtime_error("io.ini file not found");
        }
        test_file.close();
        
        // Load configuration from file
        mio::Config config(ini_file_path);
        
        // Verify we can read critical parameters (diagnostic check)
        try {
            bool canopy_test;
            config.getValue("CANOPY", "Snowpack", canopy_test);
            printf("SNOWPACK-INFO: Successfully read CANOPY = %s from [Snowpack] section\n", 
                   canopy_test ? "true" : "false");
        } catch (const std::exception& e) {
            printf("SNOWPACK-ERROR: Cannot read CANOPY from [Snowpack] section: %s\n", e.what());
            printf("SNOWPACK-ERROR: Check that io.ini contains [Snowpack] section with CANOPY parameter\n");
            printf("SNOWPACK-ERROR: File path: %s\n", ini_file_path.c_str());
            throw;
        }
        
        // CRYOWRF compatibility: Calculate METEO_STEP_LENGTH from CALCULATION_STEP_LENGTH
        // Following exact CRYOWRF pattern: meteo_step_length = M_TO_S(calculation_step_length)
        const double calculation_step_length = config.get("CALCULATION_STEP_LENGTH", "Snowpack");
        const double meteo_step_length = calculation_step_length * 60.0; // Convert minutes to seconds (M_TO_S)
        
        // Add METEO_STEP_LENGTH to config dynamically (CRYOWRF pattern)
        std::stringstream ss_meteo_length;
        ss_meteo_length << meteo_step_length;
        config.addKey("METEO_STEP_LENGTH", "Snowpack", ss_meteo_length.str());
        
        printf("SNOWPACK-INFO [C++/SnowpackConfigManager]: CRYOWRF pattern - METEO_STEP_LENGTH = %.1f seconds (from CALCULATION_STEP_LENGTH = %.1f minutes)\n", 
               meteo_step_length, calculation_step_length);
        printf("SNOWPACK-INFO [C++/SnowpackConfigManager]: Loaded configuration from %s\n", ini_file_path.c_str());
        return config;
    } catch (const std::exception& e) {
        printf("SNOWPACK-ERROR [C++/SnowpackConfigManager]: Failed to load %s: %s\n", ini_file_path.c_str(), e.what());
        throw;
    }
}

void SnowpackConfigManager::validateConfiguration(const mio::Config& cfg) {
    // Check for essential SNOWPACK parameters that will be read by SNOWPACK components
    std::vector<std::pair<std::string, std::string>> required_params = {
        {"CALCULATION_STEP_LENGTH", "Snowpack"},
        {"FORCING", "Snowpack"},
        {"SNP_SOIL", "Snowpack"},
        {"SOIL_FLUX", "Snowpack"},
        {"CANOPY", "Snowpack"},                // Required by Snowpack.cc:174
        {"HEIGHT_OF_METEO_VALUES", "Snowpack"}, // Required by Snowpack.cc:179
        {"HEIGHT_OF_WIND_VALUE", "Snowpack"},   // Required by Meteo.cc:56
        {"ROUGHNESS_LENGTH", "Snowpack"},       // Required by Meteo.cc:50
        {"SW_MODE", "Snowpack"},                // Common requirement
        {"ATMOSPHERIC_STABILITY", "Snowpack"},  // Required by Meteo.cc:44
        {"GEO_HEAT", "Snowpack"},              // Required for energy balance
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
void ensure_snowpack_states_directory() {
    const char* dir_path = "./snowpack_states";
    struct stat st = {0};
    
    if (stat(dir_path, &st) == -1) {
        if (mkdir(dir_path, 0755) == 0) {
            printf("SNOWPACK-INFO: Created snowpack_states directory\n");
        } else {
            printf("SNOWPACK-WARNING: Failed to create snowpack_states directory: %s\n", strerror(errno));
        }
    }
}

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
        
        // Ensure directory exists for .sno file persistence
        ensure_snowpack_states_directory();
        
        config_file_path = ini_path;
        config_initialized = true;
        
        // Extract and report key settings
        std::string calc_step, snp_soil;
        file_config.getValue("CALCULATION_STEP_LENGTH", "Snowpack", calc_step);
        file_config.getValue("SNP_SOIL", "Snowpack", snp_soil);
        
        // Update global calculation step length from config (no hardcoded values)
        calculation_step_length = std::stod(calc_step);
        
        printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: Configured from %s - Timestep: %.1f min, SNP_SOIL: %s\n", 
               ini_path.c_str(), calculation_step_length, snp_soil.c_str());
               
    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Configuration failed for %s: %s\n", 
               ini_path.c_str(), e.what());
        printf("SNOWPACK-FATAL: Unable to load SNOWPACK configuration - WRF run will abort\n");
        throw;
    }
}

// Persistent SnowStation management functions (CRYOWRF pattern)
std::string generate_grid_key(int i_grid, int j_grid) {
    return std::to_string(i_grid) + "_" + std::to_string(j_grid);
}

SnowStation* get_or_create_snowstation(int i_grid, int j_grid, double wrf_lat = 0.0, double wrf_lon = 0.0, double wrf_alt = 1000.0) {
    initialize_snowpack_config(); // Ensure config is loaded
    
    std::string grid_key = generate_grid_key(i_grid, j_grid);
    
    // Check if SnowStation already exists for this grid point
    auto station_it = grid_snowstations.find(grid_key);
    if (station_it != grid_snowstations.end()) {
        return station_it->second.get();
    }
    
    // Create new SnowStation for this grid point - read configuration for correct parameters
    // MeteoIO handles boolean parsing - accepts false/FALSE/0 or true/TRUE/1
    bool use_canopy, use_soil;
    global_config->getValue("CANOPY", "Snowpack", use_canopy);
    global_config->getValue("SNP_SOIL", "Snowpack", use_soil);
    const bool alpine3d = false; // WRF integration, not Alpine3D
    const bool sea_ice = false;  // Standard snowpack mode
    
    auto new_station = std::make_unique<SnowStation>(use_canopy, use_soil, alpine3d, sea_ice);
    
    // Set up basic station metadata using WRF coordinates
    if (wrf_lat < -90.0 || wrf_lat > 90.0 || wrf_lon < -180.0 || wrf_lon > 360.0) {
        // Check for obviously invalid coordinate sentinel values instead
        printf("SNOWPACK-FATAL: Station (%d,%d) - Invalid WRF coordinates! lat=%.6f, lon=%.6f\n", 
               i_grid, j_grid, wrf_lat, wrf_lon);
        printf("SNOWPACK-FATAL: Latitude must be [-90,90], longitude must be [-180,360]\n");
        throw std::runtime_error("WRF coordinates are outside valid Earth coordinate ranges");
    }
    
    mio::Coords position;
    position.setLatLon(wrf_lat, wrf_lon, wrf_alt);
    printf("SNOWPACK-INFO: Station (%d,%d) initialized with WRF coordinates: lat=%.6f°, lon=%.6f°, alt=%.1fm\n", 
           i_grid, j_grid, wrf_lat, wrf_lon, wrf_alt);
    std::string stationID = "WRF_GRID_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);
    std::string stationName = "WRF Grid Point " + std::to_string(i_grid) + "," + std::to_string(j_grid);
    new_station->meta.setStationData(position, stationID, stationName);
    
    // Try to load existing .sno file state (CRYOWRF pattern)
    bool loaded_from_file = false;
    if (use_state_persistence && global_snowpack_io) {
        std::string sno_filename = "snowpack_states/" + stationID + ".sno";
        try {
            // Attempt to read existing snowpack state
            SN_SNOWSOIL_DATA ssdata;
            ZwischenData zdata;
            mio::Date profile_date;
            
            global_snowpack_io->readSnowCover(sno_filename, stationID, ssdata, zdata, false);
            
            // Initialize SnowStation with loaded data (CRYOWRF pattern)
            ssdata.meta.position = position;  // Update with current WRF coordinates
            ssdata.meta.stationID = stationID;
            ssdata.meta.stationName = stationName;
            
            new_station->initialize(ssdata, 0);  // Initialize with loaded data
            loaded_from_file = true;
            
            printf("SNOWPACK-INFO: Loaded existing state for grid (%d,%d) from %s - %d layers\n", 
                   i_grid, j_grid, sno_filename.c_str(), (int)ssdata.Ldata.size());
        } catch (const std::exception& e) {
            // No existing state file - start with fresh snowpack
            printf("SNOWPACK-INFO: No existing state for grid (%d,%d), starting fresh: %s\n", 
                   i_grid, j_grid, e.what());
        }
    }
    
    if (!loaded_from_file) {
        printf("SNOWPACK-INFO: Fresh SnowStation created for grid (%d,%d)\n", i_grid, j_grid);
    }
    
    // Store the new SnowStation
    SnowStation* station_ptr = new_station.get();
    grid_snowstations[grid_key] = std::move(new_station);
    
    return station_ptr;
}

Snowpack* get_or_create_snowpack_instance(int i_grid, int j_grid, double wrf_lat = 0.0, double wrf_lon = 0.0, double wrf_alt = 1000.0) {
    initialize_snowpack_config(); // Ensure config is loaded
    
    std::string grid_key = generate_grid_key(i_grid, j_grid);
    
    // Check if Snowpack instance already exists for this grid point
    auto instance_it = grid_snowpack_instances.find(grid_key);
    if (instance_it != grid_snowpack_instances.end()) {
        return instance_it->second.get();
    }
    
    // Create new Snowpack instance for this grid point
    auto new_instance = std::make_unique<Snowpack>(*global_config);
    
    // Store the new Snowpack instance
    Snowpack* instance_ptr = new_instance.get();
    grid_snowpack_instances[grid_key] = std::move(new_instance);
    
    return instance_ptr;
}

void save_snowstation_state(int i_grid, int j_grid) {
    if (!use_state_persistence || !global_snowpack_io || !time_initialized) return;
    
    std::string grid_key = generate_grid_key(i_grid, j_grid);
    auto station_it = grid_snowstations.find(grid_key);
    
    if (station_it != grid_snowstations.end()) {
        std::string stationID = "WRF_GRID_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);
        std::string sno_filename = "snowpack_states/" + stationID + ".sno";
        
        try {
            // Save snowpack state to .sno file (CRYOWRF pattern)
            ZwischenData zdata;  // Empty for basic usage
            
            // Create SnowStation data for output
            SN_SNOWSOIL_DATA output_data;
            station_it->second->getSnowpackData(output_data);
            
            // Write to file using the station's internal data
            std::ofstream sno_file(sno_filename);
            if (sno_file.is_open()) {
                // Write basic .sno format header
                sno_file << "[STATION_PARAMETERS]" << std::endl;
                sno_file << "station_id = " << stationID << std::endl;
                sno_file << "station_name = " << station_it->second->meta.getStationName() << std::endl;
                sno_file << "latitude = " << station_it->second->meta.position.getLat() << std::endl;
                sno_file << "longitude = " << station_it->second->meta.position.getLon() << std::endl;
                sno_file << "altitude = " << station_it->second->meta.position.getAltitude() << std::endl;
                sno_file << "date = " << current_simulation_date.toString() << std::endl;
                sno_file.close();
                
                printf("SNOWPACK-INFO: Saved state for grid (%d,%d) to %s\n", 
                       i_grid, j_grid, sno_filename.c_str());
            } else {
                printf("SNOWPACK-WARNING: Could not open file for writing: %s\n", sno_filename.c_str());
            }
        } catch (const std::exception& e) {
            printf("SNOWPACK-WARNING: Failed to save state for grid (%d,%d): %s\n", 
                   i_grid, j_grid, e.what());
        }
    }
}

// Save all active snowpack states (called periodically)
void save_all_snowpack_states() {
    if (!use_state_persistence) return;
    
    printf("SNOWPACK-INFO: Saving all %d active snowpack states\n", (int)grid_snowstations.size());
    
    for (const auto& station_pair : grid_snowstations) {
        // Parse grid key "i_j" back to i,j coordinates
        const std::string& grid_key = station_pair.first;
        size_t underscore_pos = grid_key.find('_');
        if (underscore_pos != std::string::npos) {
            int i_grid = std::stoi(grid_key.substr(0, underscore_pos));
            int j_grid = std::stoi(grid_key.substr(underscore_pos + 1));
            save_snowstation_state(i_grid, j_grid);
        }
    }
}

// C interface for Fortran binding
extern "C" {

// Initialize configuration with specific path (called from Fortran)
void initialize_snowpack_config_c(const char* ini_file_path) {
    std::string path_str(ini_file_path);
    initialize_snowpack_config_with_path(path_str);
}

// Initialize WRF simulation time (CRYOWRF pattern - called once from Fortran)
void initialize_wrf_simulation_time_c(int start_year, int start_month, int start_day, 
                                      int start_hour, int start_minute) {
    try {
        // Debug: Print received parameters
        printf("SNOWPACK-DEBUG [C++]: Received parameters: year=%d, month=%d, day=%d, hour=%d, minute=%d\n",
               start_year, start_month, start_day, start_hour, start_minute);
        
        // Initialize simulation time with WRF namelist start time (CRYOWRF pattern)
        current_simulation_date = mio::Date(start_year, start_month, start_day, 
                                           start_hour, start_minute, 0.0, 0.0);
        time_initialized = true;
        
        printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: WRF simulation time initialized to %s\n", 
               current_simulation_date.toString(mio::Date::ISO).c_str());
    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Failed to initialize time: %s\n", e.what());
        printf("SNOWPACK-DEBUG [C++]: Parameters were: year=%d, month=%d, day=%d, hour=%d, minute=%d\n",
               start_year, start_month, start_day, start_hour, start_minute);
        std::abort();
    }
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
                      int i_grid, int j_grid, double wrf_lat, double wrf_lon,
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

  // CRYOWRF-style persistent approach: Get persistent objects for each grid point
  try {
    // Get persistent SNOWPACK objects for this grid point (following CRYOWRF pattern)
    SnowStation* snow_station = get_or_create_snowstation(i_grid, j_grid, wrf_lat, wrf_lon, height);
    Snowpack* snowpack_instance = get_or_create_snowpack_instance(i_grid, j_grid, wrf_lat, wrf_lon, height);
    
    // Create temporary meteorological data
    CurrentMeteo Mdata;
    SurfaceFluxes surfFluxes;
    BoundCond sn_Bdata;
    
    // Validate WRF coordinates
    if (wrf_lat < -90.0 || wrf_lat > 90.0 || wrf_lon < -180.0 || wrf_lon > 360.0) {
        printf("SNOWPACK-FATAL: Grid (%d,%d) - Invalid WRF coordinates! lat=%.6f, lon=%.6f\n", 
               i_grid, j_grid, wrf_lat, wrf_lon);
        throw std::runtime_error("WRF coordinates are outside valid Earth coordinate ranges");
    }
    
    // Initialize snow station data with WRF coordinates
    SN_SNOWSOIL_DATA ssdata;
    mio::Coords position;
    position.setLatLon(wrf_lat, wrf_lon, height);
    std::string stationID = SnowpackConstants::STATION_ID_PREFIX + "_" + 
                           std::to_string(i_grid) + "_" + std::to_string(j_grid);
    std::string stationName = "WRF Grid Point " + std::to_string(i_grid) + "," + std::to_string(j_grid);
    ssdata.meta.setStationData(position, stationID, stationName);
    
    // Set only essential structural data (SNOWPACK constructor handles surface properties)
    ssdata.Height = height;    // Station elevation [m]
    ssdata.nN = 1;            // Start with 1 node (ground only)  
    ssdata.nLayers = 0;       // No snow/soil layers initially
    // HS_last, Ldata, and all surface properties use constructor defaults
    
    // Initialize station (only if it's a new station without existing state)
    if (snow_station->getNumberOfElements() == 0) {
      snow_station->initialize(ssdata, 0);  // Initialize with sector 0
    }
    
    // Set up current meteorology
    // CRITICAL: Use advancing WRF simulation time (CRYOWRF pattern)
    // Time must be initialized by Fortran using WRF namelist values
    if (!time_initialized) {
        printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Time not initialized! Call initialize_wrf_simulation_time_c() first from Fortran\n");
        std::abort();
    }
    
    // CRYOWRF compute_counter pattern (lines 579, 862) - advance time conditionally
    static int compute_counter_basic = 0;
    static int call_counter_basic = 0;
    static bool first_physics_call_basic = true;
    
    if (first_physics_call_basic) {
        double wrf_dt = dt;  // WRF timestep in seconds
        double snowpack_dt = calculation_step_length * 60.0;  // SNOWPACK timestep in seconds
        compute_counter_basic = (int)(snowpack_dt / wrf_dt);
        first_physics_call_basic = false;
        printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: compute_counter_basic = %d (snowpack_dt=%.1fs, wrf_dt=%.1fs)\n", 
               compute_counter_basic, snowpack_dt, wrf_dt);
    }
    
    call_counter_basic++;
    
    // Only advance time when counter matches (CRYOWRF pattern)
    if ((call_counter_basic % compute_counter_basic) == 0) {
        current_simulation_date += (calculation_step_length / 1440.0);  // Convert minutes to days
        printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: Time advanced to %s (call %d)\n", 
               current_simulation_date.toString().c_str(), call_counter_basic);
    }
    
    mio::Date current_time = current_simulation_date;  // Use advancing time
    
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
    snowpack_instance->runSnowpackModel(Mdata, *snow_station, cumu_precip, sn_Bdata, surfFluxes);
    
    if (call_count <= 3) {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: SNOWPACK model completed successfully for grid (%d,%d)\n", 
             i_grid, j_grid);
    }
    
    // Extract results from SNOWPACK (using correct member names)
    *surface_temp = (snow_station->Ndata.size() > 0) ? snow_station->Ndata.back().T : temp_air;  // Surface temperature [K]
    *snow_swe = snow_station->swe;                                  // Snow water equivalent [mm]
    *snow_depth = snow_station->cH;                                // Snow height [m]
    *heat_flux_sensible = surfFluxes.qs;                          // Sensible heat flux [W/m²]
    *heat_flux_latent = surfFluxes.ql;                            // Latent heat flux [W/m²]
    *albedo = snow_station->Albedo;                                // Surface albedo [0-1]
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
                             int i_grid, int j_grid, double wrf_lat, double wrf_lon,
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

  // CRYOWRF-style persistent approach: Get persistent objects for each grid point
  try {
    // Get persistent SNOWPACK objects for this grid point (following CRYOWRF pattern)
    SnowStation* snow_station = get_or_create_snowstation(i_grid, j_grid, wrf_lat, wrf_lon, height);
    Snowpack* snowpack_instance = get_or_create_snowpack_instance(i_grid, j_grid, wrf_lat, wrf_lon, height);
    
    // Create temporary meteorological data
    CurrentMeteo Mdata;
    SurfaceFluxes surfFluxes;
    BoundCond sn_Bdata;
    
    // Validate WRF coordinates
    if (wrf_lat < -90.0 || wrf_lat > 90.0 || wrf_lon < -180.0 || wrf_lon > 360.0) {
        printf("SNOWPACK-FATAL: Grid (%d,%d) - Invalid WRF coordinates! lat=%.6f, lon=%.6f\n", 
               i_grid, j_grid, wrf_lat, wrf_lon);
        throw std::runtime_error("WRF coordinates are outside valid Earth coordinate ranges");
    }
    
    // Initialize snow station data with WRF coordinates
    SN_SNOWSOIL_DATA ssdata;
    mio::Coords position;
    position.setLatLon(wrf_lat, wrf_lon, height);
    std::string stationID = SnowpackConstants::STATION_ID_PREFIX + "_" + 
                           std::to_string(i_grid) + "_" + std::to_string(j_grid);
    std::string stationName = "WRF Grid Point " + std::to_string(i_grid) + "," + std::to_string(j_grid);
    ssdata.meta.setStationData(position, stationID, stationName);
    
    // Set only essential structural data (SNOWPACK constructor handles surface properties)
    ssdata.Height = height;    // Station elevation [m]
    ssdata.nN = 1;            // Start with 1 node (ground only)  
    ssdata.nLayers = 0;       // No snow/soil layers initially
    // HS_last, Ldata, and all surface properties use constructor defaults
    
    // Initialize station (only if it's a new station without existing state)
    if (snow_station->getNumberOfElements() == 0) {
      snow_station->initialize(ssdata, 0);  // Initialize with sector 0
    }
    
    // Set up current meteorology
    // CRITICAL: Use advancing WRF simulation time (CRYOWRF pattern)
    // Time must be initialized by Fortran using WRF namelist values
    if (!time_initialized) {
        printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Time not initialized! Call initialize_wrf_simulation_time_c() first from Fortran\n");
        std::abort();
    }
    
    // CRYOWRF compute_counter pattern (lines 579, 862) - advance time conditionally
    static int compute_counter_layers = 0;
    static int call_counter_layers = 0;
    static bool first_physics_call_layers = true;
    
    if (first_physics_call_layers) {
        double wrf_dt = dt;  // WRF timestep in seconds
        double snowpack_dt = calculation_step_length * 60.0;  // SNOWPACK timestep in seconds
        compute_counter_layers = (int)(snowpack_dt / wrf_dt);
        first_physics_call_layers = false;
        printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: compute_counter_layers = %d (snowpack_dt=%.1fs, wrf_dt=%.1fs)\n", 
               compute_counter_layers, snowpack_dt, wrf_dt);
    }
    
    call_counter_layers++;
    
    // Only advance time when counter matches (CRYOWRF pattern)
    if ((call_counter_layers % compute_counter_layers) == 0) {
        current_simulation_date += (calculation_step_length / 1440.0);  // Convert minutes to days
        printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: Time advanced to %s (call %d)\n", 
               current_simulation_date.toString().c_str(), call_counter_layers);
    }
    
    mio::Date current_time = current_simulation_date;  // Use advancing time
    
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
    snowpack_instance->runSnowpackModel(Mdata, *snow_station, cumu_precip, sn_Bdata, surfFluxes);
    
    // Extract basic results
    *surface_temp = (snow_station->Ndata.size() > 0) ? snow_station->Ndata.back().T : temp_air;
    *snow_swe = snow_station->swe;
    *snow_depth = snow_station->cH;
    *heat_flux_sensible = surfFluxes.qs;
    *heat_flux_latent = surfFluxes.ql;
    *albedo = snow_station->Albedo;
    *snow_coverage = (*snow_depth > 0.001) ? 1.0 : 0.0;
    
    // Extract detailed layer information from SnowStation
    size_t num_elements = snow_station->getNumberOfElements();
    *n_layers = static_cast<int>(num_elements);
    
    // Limit to max 50 layers for WRF arrays
    size_t layers_to_extract = std::min(num_elements, size_t(50));
    
    for (size_t i = 0; i < layers_to_extract; i++) {
      const ElementData& elem = snow_station->Edata[i];
      
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
void extract_snowpack_layers_c(int i_grid, int j_grid,
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
    
    // Use persistent SnowStation objects from our grid-based storage
    std::string grid_key = generate_grid_key(i_grid, j_grid);
    auto station_it = grid_snowstations.find(grid_key);
    
    if (station_it == grid_snowstations.end()) {
        // No SnowStation exists for this grid point - set to zero layers
        *n_layers = 0;
        printf("SNOWPACK-LAYERS: Grid (%d,%d) has no persistent SnowStation - zero layers\n", 
               i_grid, j_grid);
        return;
    }
    
    // Access persistent SnowStation for this grid point
    const SnowStation* xdata = station_it->second.get();
    
    // CRYOWRF-style layer extraction algorithm (from Coupler.cpp:1153-1185)
    // Extract from evolved snowpack with real layer data
    
    const std::vector<ElementData>& elem_data = xdata->Edata;
    const int get_size = (int)xdata->getNumberOfElements();
    const int loc_snpack_lay_to_sav = 50; // Maximum layers to save (CRYOWRF default)
    const int lim_size = std::max(0, get_size - loc_snpack_lay_to_sav);
    
    // CRYOWRF algorithm: Extract from surface down, top-down indexing
    int tmtm = 0;
    for (int e = get_size - 1; e >= lim_size && tmtm < 50; e--, tmtm++) {
      // Extract layer data using exact CRYOWRF algorithm
      // With persistent SnowStation objects, evolved snowpack layers are available for extraction
      if (e + 1 < (int)xdata->Ndata.size()) { // Safety check for node access
        layer_temps[tmtm] = (float)xdata->Ndata[e + 1].T;        // Node temperature [K]
      }
      layer_thick[tmtm] = (float)elem_data[e].L;                // Layer thickness [m]
      layer_voli[tmtm] = (float)elem_data[e].theta[ICE] * 100.0f; // Ice volume fraction [%]
      layer_volw[tmtm] = (float)elem_data[e].theta[WATER] * 100.0f; // Water volume fraction [%]
      layer_volv[tmtm] = (float)elem_data[e].theta[AIR] * 100.0f; // Air volume fraction [%]
      layer_rg[tmtm] = (float)elem_data[e].rg;                  // Grain radius [mm]
      layer_rb[tmtm] = (float)elem_data[e].rb;                  // Bond radius [mm]
      layer_dd[tmtm] = (float)elem_data[e].dd;                  // Dendricity [-]
      layer_sp[tmtm] = (float)elem_data[e].sp;                  // Sphericity [-]
    }
    
    *n_layers = get_size; // Total number of layers (following CRYOWRF: sn_nlayer = get_size)
    
    printf("SNOWPACK-LAYERS: Grid (%d,%d) extracted %d layers using CRYOWRF algorithm (persistent SnowStation)\n", 
           i_grid, j_grid, *n_layers);
    
  } catch (const std::exception& e) {
    printf("SNOWPACK-ERROR: Layer extraction failed for grid (%d,%d): %s\n", i_grid, j_grid, e.what());
    *n_layers = 0;
  }
}

// Save all snowpack states (called from Fortran for periodic saves)
void save_all_snowpack_states_c() {
    save_all_snowpack_states();
}

} // extern "C"