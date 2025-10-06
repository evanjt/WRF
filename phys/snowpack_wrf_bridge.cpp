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
#include <algorithm>
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
// State persistence is always enabled for WRF-SNOWPACK integration

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

// Global call counters for debugging - track ALL function calls from Fortran
static int config_init_calls = 0;
static int time_init_calls = 0;
static int physics_calls = 0;
static int physics_layers_calls = 0;
static int structured_calls = 0;
static int station_creation_calls = 0;

// Initialize SNOWPACK configuration (once, globally)
void initialize_snowpack_config() {
    printf("SNOWPACK-DEBUG: ======== SNOWPACK INITIALIZATION STARTING ========\n");
    printf("SNOWPACK-DEBUG: initialize_snowpack_config() called - loading default config\n");
    initialize_snowpack_config_with_path(SnowpackConfigManager::getDefaultConfigPath());
    printf("SNOWPACK-DEBUG: ======== SNOWPACK INITIALIZATION COMPLETED ========\n");
}

// Initialize SNOWPACK configuration with specific file path
void initialize_snowpack_config_with_path(const std::string& ini_path) {
    printf("SNOWPACK-DEBUG: initialize_snowpack_config_with_path() called with '%s'\n", ini_path.c_str());
    printf("SNOWPACK-DEBUG: config_initialized=%s\n", config_initialized ? "true" : "false");

    if (config_initialized) {
        printf("SNOWPACK-DEBUG: Configuration already initialized, skipping\n");
        return;
    }

    printf("SNOWPACK-DEBUG: Starting configuration loading process...\n");
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

        config_initialized = true;
        printf("SNOWPACK-DEBUG: Configuration successfully loaded and initialized\n");

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

SnowStation* get_or_create_snowstation(int i_grid, int j_grid, int wrf_domain_id = 1, double wrf_lat = 0.0, double wrf_lon = 0.0, double wrf_alt = 1000.0) {
    // wrf_domain_id: 1=first domain (most common), 2=nested domain #2, etc.

    station_creation_calls++;
    // DEBUG: Track who is calling this function and with what coordinates
    printf("SNOWPACK-DEBUG [CALL #%d]: get_or_create_snowstation(%d,%d,domain=%d,lat=%.6f,lon=%.6f,alt=%.1f)\n",
           station_creation_calls, i_grid, j_grid, wrf_domain_id, wrf_lat, wrf_lon, wrf_alt);

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
    // WRF's XLAT and XLONG contain geographic coordinates in degrees (regardless of map projection)
    
    mio::Coords position;
    
    // Debug: Log what we receive from WRF
    printf("SNOWPACK-DEBUG [C++]: Station (%d,%d) - Received from WRF: lat=%.6f°, lon=%.6f°, alt=%.1fm\n", 
           i_grid, j_grid, wrf_lat, wrf_lon, wrf_alt);
    
    // WRF passes geographic coordinates (latitude/longitude in degrees)
    // CRYOWRF uses setLatLon with these values directly
    // setLatLon expects (latitude, longitude, altitude) - standard geographic order
    try {
        position.setLatLon(wrf_lat, wrf_lon, wrf_alt);  // Correct order: lat, lon, alt
        printf("SNOWPACK-INFO: Successfully set geographic coordinates for station (%d,%d)\n", i_grid, j_grid);
    } catch (const mio::InvalidArgumentException& e) {
        printf("SNOWPACK-ERROR: Invalid coordinates for station (%d,%d): %s\n", i_grid, j_grid, e.what());
        printf("SNOWPACK-ERROR: Received lat=%.6f, lon=%.6f - these may not be geographic coordinates\n", wrf_lat, wrf_lon);
        // If coordinates are invalid, there's likely a data passing issue
        // Don't try to use setXY as a fallback - that's for projected coordinates in meters
        throw;
    }
    std::string stationID = "WRF_GRID_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);
    std::string stationName = "WRF Grid Point " + std::to_string(i_grid) + "," + std::to_string(j_grid);
    new_station->meta.setStationData(position, stationID, stationName);
    
    // Try to load existing .sno file state (CRYOWRF pattern)
    bool loaded_from_file = false;
    if (global_snowpack_io) {
        // CRYOWRF C++ naming pattern: snpack_{grid_id}_{I}_{J}.sno (from Coupler.cpp line 613)
        // Let SNOWPACK handle the path through SNOWPATH configuration
        std::string sno_filename = "snpack_" + std::to_string(wrf_domain_id) + "_" + std::to_string(i_grid) + "_" + std::to_string(j_grid) + ".sno";

        printf("SNOWPACK-INIT [%d,%d]: Attempting to load .sno file: %s\n", i_grid, j_grid, sno_filename.c_str());

        try {
            // Attempt to read existing snowpack state
            SN_SNOWSOIL_DATA ssdata;
            ZwischenData zdata;
            mio::Date profile_date;

            printf("SNOWPACK-INIT [%d,%d]: Calling readSnowCover()...\n", i_grid, j_grid);
            global_snowpack_io->readSnowCover(sno_filename, stationID, ssdata, zdata, false);
            printf("SNOWPACK-INIT [%d,%d]: readSnowCover() returned %zu layers\n", i_grid, j_grid, ssdata.nLayers);

            // DEBUG: Check what coordinates were loaded from .sno file
            printf("SNOWPACK-INIT [%d,%d]: Coordinates loaded from .sno file: lat=%.6f, lon=%.6f, alt=%.1f\n",
                   i_grid, j_grid, ssdata.meta.position.getLat(), ssdata.meta.position.getLon(), ssdata.meta.position.getAltitude());

            // IMPORTANT FIX: Don't overwrite valid coordinates from .sno file with default (0.0,0.0) values
            // Only use WRF coordinates if .sno file has invalid coordinates (both lat and lon are 0)
            double sno_lat = ssdata.meta.position.getLat();
            double sno_lon = ssdata.meta.position.getLon();
            double sno_alt = ssdata.meta.position.getAltitude();

            if (std::abs(sno_lat) > 0.001 && std::abs(sno_lon) > 0.001) {
                // .sno file has valid coordinates - use them!
                printf("SNOWPACK-INIT [%d,%d]: Using valid coordinates from .sno file\n", i_grid, j_grid);
                // position = ssdata.meta.position;  // Keep the loaded coordinates
            } else {
                // .sno file has invalid coordinates - use WRF coordinates
                printf("SNOWPACK-INIT [%d,%d]: .sno file has invalid coordinates, using WRF coordinates\n", i_grid, j_grid);
                ssdata.meta.position = position;  // Update with WRF coordinates
            }

            ssdata.meta.stationID = stationID;
            ssdata.meta.stationName = stationName;

            printf("SNOWPACK-INIT [%d,%d]: Final coordinates: lat=%.6f, lon=%.6f, alt=%.1f\n",
                   i_grid, j_grid, ssdata.meta.position.getLat(), ssdata.meta.position.getLon(), ssdata.meta.position.getAltitude());

            printf("SNOWPACK-INIT [%d,%d]: Calling SnowStation::initialize()...\n", i_grid, j_grid);
            new_station->initialize(ssdata, 0);  // Initialize with loaded data
            printf("SNOWPACK-INIT [%d,%d]: SnowStation::initialize() completed - %zu elements created\n",
                   i_grid, j_grid, new_station->getNumberOfElements());

            // CRITICAL FIX: ASCII SMET .sno files don't contain k/c thermal property vectors
            // initialize() should set them to 0, but explicitly verify and fix if needed
            printf("SNOWPACK-INIT [%d,%d]: Verifying k/c thermal property vectors...\n", i_grid, j_grid);

            size_t nan_k_count = 0, nan_c_count = 0, resize_k_count = 0, resize_c_count = 0;

            for (size_t e = 0; e < new_station->getNumberOfElements(); e++) {
                // Check initial state before fixes
                bool had_k_issue = (new_station->Edata[e].k.size() < 3);
                bool had_c_issue = (new_station->Edata[e].c.size() < 3);

                if (had_k_issue) {
                    printf("SNOWPACK-INIT [%d,%d] Layer %zu: k vector size=%zu (WRONG, should be 3) - RESIZING\n",
                           i_grid, j_grid, e, new_station->Edata[e].k.size());
                    resize_k_count++;
                }
                if (had_c_issue) {
                    printf("SNOWPACK-INIT [%d,%d] Layer %zu: c vector size=%zu (WRONG, should be 3) - RESIZING\n",
                           i_grid, j_grid, e, new_station->Edata[e].c.size());
                    resize_c_count++;
                }

                // Ensure k and c vectors have proper size (3 elements: TEMPERATURE, SEEPAGE, SETTLEMENT)
                if (new_station->Edata[e].k.size() < 3) {
                    new_station->Edata[e].k.resize(3, 0.0);  // Resize and initialize to 0
                }
                if (new_station->Edata[e].c.size() < 3) {
                    new_station->Edata[e].c.resize(3, 0.0);  // Resize and initialize to 0
                }

                // Force initialize to 0.0 if NaN or uninitialized
                for (size_t i = 0; i < 3; i++) {
                    if (std::isnan(new_station->Edata[e].k[i]) || new_station->Edata[e].k[i] != new_station->Edata[e].k[i]) {
                        printf("SNOWPACK-INIT [%d,%d] Layer %zu: k[%zu]=NaN detected - setting to 0.0\n",
                               i_grid, j_grid, e, i);
                        new_station->Edata[e].k[i] = 0.0;
                        nan_k_count++;
                    }
                    if (std::isnan(new_station->Edata[e].c[i]) || new_station->Edata[e].c[i] != new_station->Edata[e].c[i]) {
                        printf("SNOWPACK-INIT [%d,%d] Layer %zu: c[%zu]=NaN detected - setting to 0.0\n",
                               i_grid, j_grid, e, i);
                        new_station->Edata[e].c[i] = 0.0;
                        nan_c_count++;
                    }
                }

                // Recompute heat capacity (c[TEMPERATURE]) from layer properties
                new_station->Edata[e].heatCapacity();

                // Verify after fixes (only print first 3 layers to avoid spam)
                if (e < 3) {
                    printf("SNOWPACK-INIT [%d,%d] Layer %zu after fixes: k=[%.6f,%.6f,%.6f], c=[%.3f,%.3f,%.3f]\n",
                           i_grid, j_grid, e,
                           new_station->Edata[e].k[0], new_station->Edata[e].k[1], new_station->Edata[e].k[2],
                           new_station->Edata[e].c[0], new_station->Edata[e].c[1], new_station->Edata[e].c[2]);
                }
            }

            if (resize_k_count > 0 || resize_c_count > 0 || nan_k_count > 0 || nan_c_count > 0) {
                printf("SNOWPACK-INIT [%d,%d]: ISSUES FIXED - k_resize=%zu, c_resize=%zu, k_nan=%zu, c_nan=%zu\n",
                       i_grid, j_grid, resize_k_count, resize_c_count, nan_k_count, nan_c_count);
            } else {
                printf("SNOWPACK-INIT [%d,%d]: All k/c vectors verified OK (no fixes needed)\n", i_grid, j_grid);
            }

            loaded_from_file = true;
            printf("SNOWPACK-INFO: Loaded existing state for grid (%d,%d) from %s - %zu layers\n",
                   i_grid, j_grid, sno_filename.c_str(), new_station->getNumberOfElements());
            if (new_station->getNumberOfElements() > 0) {
                printf("SNOWPACK-DEBUG: First layer after load - T=%.2fK, L=%.3fm, theta_i=%.3f, k[0]=%.6f, c[0]=%.3f\n",
                       new_station->Edata[0].Te, new_station->Edata[0].L,
                       new_station->Edata[0].theta[ICE],
                       new_station->Edata[0].k[0], new_station->Edata[0].c[0]);
            }
        } catch (const std::exception& e) {
            // No existing state file - start with fresh snowpack
            printf("SNOWPACK-INIT [%d,%d]: Failed to load .sno file - %s\n", i_grid, j_grid, e.what());
            printf("SNOWPACK-INFO: No existing state for grid (%d,%d), starting fresh\n", i_grid, j_grid);
        }
    } else {
        printf("SNOWPACK-INIT [%d,%d]: global_snowpack_io is NULL - cannot load .sno files\n", i_grid, j_grid);
    }
    
    if (!loaded_from_file) {
        printf("SNOWPACK-INFO: Fresh SnowStation created for grid (%d,%d)\n", i_grid, j_grid);
    }
    
    // Store the new SnowStation
    SnowStation* station_ptr = new_station.get();
    grid_snowstations[grid_key] = std::move(new_station);

    printf("SNOWPACK-DEBUG [CALL #%d]: get_or_create_snowstation() completed successfully for grid (%d,%d)\n",
           station_creation_calls, i_grid, j_grid);
    printf("SNOWPACK-DEBUG: Total stations created so far: %zu\n", grid_snowstations.size());

    return station_ptr;
}

Snowpack* get_or_create_snowpack_instance(int i_grid, int j_grid, double wrf_lat = 0.0, double wrf_lon = 0.0, double wrf_alt = 1000.0) {
    printf("SNOWPACK-DEBUG: Creating Snowpack instance for grid (%d,%d), lat=%.6f, lon=%.6f\n",
           i_grid, j_grid, wrf_lat, wrf_lon);
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
    if (!global_snowpack_io || !time_initialized) return;
    
    std::string grid_key = generate_grid_key(i_grid, j_grid);
    auto station_it = grid_snowstations.find(grid_key);
    
    if (station_it != grid_snowstations.end()) {
        std::string stationID = "WRF_GRID_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);
        // CRYOWRF output naming: SNOWPATH directory with stationID.sno format
        // This follows SNOWPACK configuration SNOWPATH = ./snowpack_states
        std::string sno_filename = "./snowpack_states/" + stationID + ".sno";
        
        try {
            // Save snowpack state to .sno file (CRYOWRF pattern)
            ZwischenData zdata;  // Empty for basic usage
            
            // Use SNOWPACK's official writeSnowCover method (CRYOWRF pattern)
            
            global_snowpack_io->writeSnowCover(current_simulation_date, *(station_it->second), zdata, true);
            printf("SNOWPACK-INFO: Saved .sno state for grid (%d,%d) using SnowpackIO\n", i_grid, j_grid);
        } catch (const std::exception& e) {
            printf("SNOWPACK-WARNING: Failed to save .sno state for grid (%d,%d): %s\n", 
                   i_grid, j_grid, e.what());
        }
    }
}

// Save all active snowpack states (called periodically)
void save_all_snowpack_states() {
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
    config_init_calls++;
    printf("SNOWPACK-DEBUG [CALL #%d]: initialize_snowpack_config_c() called with ini_file_path='%s'\n",
           config_init_calls, ini_file_path);

    std::string path_str(ini_file_path);
    initialize_snowpack_config_with_path(path_str);

    printf("SNOWPACK-DEBUG [CALL #%d]: initialize_snowpack_config_c() completed successfully\n", config_init_calls);
}

// Initialize WRF simulation time (CRYOWRF pattern - called once from Fortran)
void initialize_wrf_simulation_time_c(int start_year, int start_month, int start_day,
                                      int start_hour, int start_minute) {
    time_init_calls++;
    printf("SNOWPACK-DEBUG [CALL #%d]: initialize_wrf_simulation_time_c() called with %d-%02d-%02d %02d:%02d\n",
           time_init_calls, start_year, start_month, start_day, start_hour, start_minute);

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
        printf("SNOWPACK-DEBUG [CALL #%d]: initialize_wrf_simulation_time_c() completed successfully\n", time_init_calls);
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
    // Pass domain_id=1 explicitly to avoid parameter mismatch
    printf("SNOWPACK-DEBUG [%d,%d]: About to create/get SNOWPACK objects (lat=%.6f, lon=%.6f, height=%.1f)...\n",
           i_grid, j_grid, wrf_lat, wrf_lon, height);
    SnowStation* snow_station = get_or_create_snowstation(i_grid, j_grid, 1, wrf_lat, wrf_lon, height);
    printf("SNOWPACK-DEBUG [%d,%d]: snow_station created/retrieved: %p\n", i_grid, j_grid, (void*)snow_station);

    Snowpack* snowpack_instance = get_or_create_snowpack_instance(i_grid, j_grid, wrf_lat, wrf_lon, height);
    printf("SNOWPACK-DEBUG [%d,%d]: snowpack_instance created/retrieved: %p\n", i_grid, j_grid, (void*)snowpack_instance);

    // Validate pointers before use
    if (!snow_station) {
        printf("SNOWPACK-FATAL [%d,%d]: snow_station pointer is null!\n", i_grid, j_grid);
        std::abort();
    }
    if (!snowpack_instance) {
        printf("SNOWPACK-FATAL [%d,%d]: snowpack_instance pointer is null!\n", i_grid, j_grid);
        std::abort();
    }

    // Create temporary meteorological data
    printf("SNOWPACK-DEBUG [%d,%d]: Creating temporary data objects...\n", i_grid, j_grid);
    auto Mdata = std::make_unique<CurrentMeteo>();
    auto surfFluxes = std::make_unique<SurfaceFluxes>();
    auto sn_Bdata = std::make_unique<BoundCond>();
    printf("SNOWPACK-DEBUG [%d,%d]: Temporary objects created successfully\n", i_grid, j_grid);
    
    // Note: WRF coordinates may be projected (grid coordinates in meters) rather than lat/lon
    
    // Initialize snow station data with WRF coordinates
    SN_SNOWSOIL_DATA ssdata;
    mio::Coords position;
    // WRF always passes geographic coordinates (XLAT/XLONG in degrees)
    // Use setLatLon with proper order: latitude, longitude, altitude
    position.setLatLon(wrf_lat, wrf_lon, height);  // Geographic coordinates from WRF
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
        // WRF calls physics every wrf_dt, SNOWPACK needs to run every snowpack_dt
        // So we run SNOWPACK every (wrf_dt/snowpack_dt) WRF calls
        // Following CRYOWRF pattern exactly
        compute_counter_basic = (int)(wrf_dt / snowpack_dt);
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
    printf("SNOWPACK-DEBUG [%d,%d]: Setting up Mdata meteorological structure...\n", i_grid, j_grid);
    Mdata->date = current_time;
    Mdata->ta = safe_temp;                                       // Air temperature [K]
    Mdata->rh = std::max(0.01, std::min(1.0, humidity));       // Relative humidity [0-1]
    Mdata->vw = std::max(0.1, wind_speed);                     // Wind speed [m/s]
    Mdata->dw = wind_dir;                                       // Wind direction [degrees]
    Mdata->iswr = std::max(0.0, shortwave_in);                 // Incoming shortwave [W/m²]
    Mdata->lw_net = std::max(0.0, longwave_in);               // Net longwave radiation [W/m²]
    Mdata->psum = std::max(0.0, precipitation);                // Precipitation [mm]

    printf("SNOWPACK-DEBUG [%d,%d]: Basic Mdata set - ta=%.2f, rh=%.3f, vw=%.2f, dw=%.1f, iswr=%.1f, lw_net=%.1f, psum=%.3f\n",
           i_grid, j_grid, Mdata->ta, Mdata->rh, Mdata->vw, Mdata->dw, Mdata->iswr, Mdata->lw_net, Mdata->psum);
    printf("SNOWPACK-DEBUG [%d,%d]: >>> WARNING: Received pressure=%.1f hPa from WRF but CurrentMeteo has no pressure field! <<<\n",
           i_grid, j_grid, pressure);
    printf("SNOWPACK-DEBUG [%d,%d]: This could cause incorrect air density, humidity, and turbulence calculations\n",
           i_grid, j_grid);
    // Note: CurrentMeteo has no pressure field - pressure used internally by SNOWPACK

    // Additional required meteorological parameters
    Mdata->psum_ph = (safe_temp < 273.65) ? 0.0 : 1.0;        // Precipitation phase (0=snow, 1=rain)
    Mdata->tss = mio::IOUtils::nodata;                         // Surface temperature (let SNOWPACK compute)
    Mdata->ts0 = safe_temp - 5.0;                              // Bottom temperature estimate
    Mdata->hs = *snow_depth;                                   // Current snow height [m]

    // CRITICAL: Set roughness length and friction velocity for wind pumping calculations
    // Without these, wind pumping produces NaN thermal conductivity
    // Read ROUGHNESS_LENGTH from config (SNOWPACK's Meteo class also does this)
    static double roughness_length = -1.0;  // Cache the value
    if (roughness_length < 0.0) {
        if (global_config) {
            global_config->getValue("ROUGHNESS_LENGTH", "Snowpack", roughness_length, mio::IOUtils::nothrow);
            if (roughness_length < 0.0) {
                roughness_length = 0.002;  // Default if not in config
                printf("SNOWPACK-INFO: ROUGHNESS_LENGTH not found in config, using default %.4f m\n", roughness_length);
            } else {
                printf("SNOWPACK-INFO: Read ROUGHNESS_LENGTH=%.4f m from config\n", roughness_length);
            }
        } else {
            roughness_length = 0.002;  // Fallback default
        }
    }

    // Adapt roughness based on snow presence (same logic as Meteo::MicroMet line 267)
    const double rough_len = (*snow_depth > 0.03) ? roughness_length : 0.01;  // BareSoil_z0 typically ~0.01
    Mdata->z0 = rough_len;

    // Compute friction velocity from wind speed using log profile (same as Meteo::MicroMet line 316)
    // u* = k * u / ln(z/z0), where k=0.4 (von Karman constant)
    const double von_karman = 0.4;
    const double z_wind = height;  // Measurement height [m]

    printf("SNOWPACK-DEBUG [%d,%d]: Computing turbulence - z_wind=%.2f, z0=%.6f, snow_depth=%.3f\n",
           i_grid, j_grid, z_wind, Mdata->z0, *snow_depth);
    printf("SNOWPACK-DEBUG [%d,%d]: ustar calculation: (z_wind > z0)? %d, log(z/z0)=%.6f\n",
           i_grid, j_grid, (z_wind > Mdata->z0), (z_wind > Mdata->z0) ? std::log(z_wind / Mdata->z0) : 0.0);

    Mdata->ustar = (z_wind > Mdata->z0) ?
                  (von_karman * Mdata->vw / std::log(z_wind / Mdata->z0)) :
                  (0.1 * Mdata->vw);  // Fallback if z <= z0

    printf("SNOWPACK-DEBUG [%d,%d]: Computed ustar=%.6f m/s using wind speed %.2f m/s\n",
           i_grid, j_grid, Mdata->ustar, Mdata->vw);

    // Check for invalid values that could cause crashes
    if (std::isnan(Mdata->ustar) || std::isinf(Mdata->ustar)) {
        printf("SNOWPACK-FATAL [%d,%d]: ustar is NaN/Inf (%.6f)! This will crash SNOWPACK!\n", i_grid, j_grid, Mdata->ustar);
        printf("SNOWPACK-FATAL [%d,%d]: z_wind=%.6f, z0=%.6f, log(z/z0)=%.6f\n", i_grid, j_grid, z_wind, Mdata->z0, std::log(z_wind / Mdata->z0));
        std::abort();
    }
    if (std::isnan(Mdata->z0) || std::isinf(Mdata->z0) || Mdata->z0 <= 0.0) {
        printf("SNOWPACK-FATAL [%d,%d]: z0 is invalid (%.6f)! This will crash SNOWPACK!\n", i_grid, j_grid, Mdata->z0);
        std::abort();
    }
    
    // DEBUG: Print z0 and ustar values
    if (call_count <= 5) {
        printf("SNOWPACK-DEBUG [%d,%d]: z0=%.6f m, ustar=%.6f m/s, VW=%.2f m/s, height=%.2f m\n",
               i_grid, j_grid, Mdata->z0, Mdata->ustar, Mdata->vw, height);
    }
    
    // Run SNOWPACK model (temporary objects will auto-destruct)
    if (call_count <= 3) {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Calling SNOWPACK model for grid (%d,%d), precip=%.3f\n",
             i_grid, j_grid, precipitation);
    }

    // Debug: Print Mdata values before SNOWPACK call
    printf("SNOWPACK-PHYSICS [%d,%d]: BEFORE runSnowpackModel - Ta=%.2fK, RH=%.4f, VW=%.2fm/s, ISWR=%.1fW/m2, psum=%.3fmm\n",
           i_grid, j_grid, Mdata->ta, Mdata->rh, Mdata->vw, Mdata->iswr, Mdata->psum);

    // Debug: Check snow station state before calling SNOWPACK
    printf("SNOWPACK-PHYSICS [%d,%d]: SnowStation state - nElements=%zu, nNodes=%zu, cH=%.3fm, swe=%.2fmm\n",
           i_grid, j_grid, snow_station->getNumberOfElements(), snow_station->getNumberOfNodes(),
           snow_station->cH, snow_station->swe);

    if (snow_station->getNumberOfElements() > 0) {
        printf("SNOWPACK-PHYSICS [%d,%d]: First element BEFORE runSnowpackModel - Te=%.2fK, L=%.3fm, k=[%.6f,%.6f,%.6f], c=[%.3f,%.3f,%.3f]\n",
               i_grid, j_grid, snow_station->Edata[0].Te, snow_station->Edata[0].L,
               snow_station->Edata[0].k[0], snow_station->Edata[0].k[1], snow_station->Edata[0].k[2],
               snow_station->Edata[0].c[0], snow_station->Edata[0].c[1], snow_station->Edata[0].c[2]);
    }

    // Execute SNOWPACK physics (correct API with cumulative precipitation parameter)
    double cumu_precip = 0.0;  // Cumulative precipitation parameter

    // Detailed input parameter logging before SNOWPACK call
    printf("SNOWPACK-PHYSICS [%d,%d]: >>> INPUT TO SNOWPACK <<<\n", i_grid, j_grid);
    printf("SNOWPACK-PHYSICS [%d,%d]:   Mdata->ta=%.2f°C (%.1fK), rh=%.4f, vw=%.2fm/s\n",
           i_grid, j_grid, Mdata->ta - 273.15, Mdata->ta, Mdata->rh, Mdata->vw);
    printf("SNOWPACK-PHYSICS [%d,%d]:   Mdata->iswr=%.1fW/m2, Mdata->lw_net=%.1fW/m2, Mdata->psum=%.6fmm\n",
           i_grid, j_grid, Mdata->iswr, Mdata->lw_net, Mdata->psum);
    printf("SNOWPACK-PHYSICS [%d,%d]:   Mdata->tss=%.2f°C, hs=%.3fm, cumu_precip=%.6fmm\n",
           i_grid, j_grid, (Mdata->tss != mio::IOUtils::nodata) ? Mdata->tss - 273.15 : -999.0, Mdata->hs, cumu_precip);
    printf("SNOWPACK-PHYSICS [%d,%d]:   Station: z=%.1fm, azi=%.1f°, slope=%.1f°\n",
           i_grid, j_grid, snow_station->meta.position.getAltitude(),
           snow_station->meta.getAzimuth(), snow_station->meta.getSlopeAngle());
    printf("SNOWPACK-PHYSICS [%d,%d]: >>> CALLING snowpack_instance->runSnowpackModel() <<<\n", i_grid, j_grid);
    snowpack_instance->runSnowpackModel(*Mdata, *snow_station, cumu_precip, *sn_Bdata, *surfFluxes);
    printf("SNOWPACK-PHYSICS [%d,%d]: >>> runSnowpackModel() returned successfully <<<\n", i_grid, j_grid);

    // Detailed output logging after SNOWPACK call
    printf("SNOWPACK-PHYSICS [%d,%d]: <<< OUTPUT FROM SNOWPACK <<<\n", i_grid, j_grid);
    printf("SNOWPACK-PHYSICS [%d,%d]:   Station data: SWE=%.3fmm, depth=%.3fm, elements=%zu\n",
           i_grid, j_grid, snow_station->swe, snow_station->cH, snow_station->getNumberOfElements());
    printf("SNOWPACK-PHYSICS [%d,%d]:   Surface data: Ts=%.2f°C\n",
           i_grid, j_grid, (snow_station->Cdata.temp != mio::IOUtils::nodata) ? snow_station->Cdata.temp - 273.15 : -999.0);

    if (call_count <= 3) {
      printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: SNOWPACK model completed successfully for grid (%d,%d)\n",
             i_grid, j_grid);
    }

    // Debug: Check what changed after SNOWPACK ran
    if (snow_station->getNumberOfElements() > 0) {
        printf("SNOWPACK-PHYSICS [%d,%d]: First element AFTER runSnowpackModel - Te=%.2fK, L=%.3fm, k=[%.6f,%.6f,%.6f], c=[%.3f,%.3f,%.3f]\n",
               i_grid, j_grid, snow_station->Edata[0].Te, snow_station->Edata[0].L,
               snow_station->Edata[0].k[0], snow_station->Edata[0].k[1], snow_station->Edata[0].k[2],
               snow_station->Edata[0].c[0], snow_station->Edata[0].c[1], snow_station->Edata[0].c[2]);
    }
    
    // Extract results from SNOWPACK (using correct member names)
    *surface_temp = (snow_station->Ndata.size() > 0) ? snow_station->Ndata.back().T : temp_air;  // Surface temperature [K]
    *snow_swe = snow_station->swe;                                  // Snow water equivalent [mm]
    *snow_depth = snow_station->cH;                                // Snow height [m]
    *heat_flux_sensible = surfFluxes->qs;                          // Sensible heat flux [W/m²]
    *heat_flux_latent = surfFluxes->ql;                            // Latent heat flux [W/m²]
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

// Enhanced interface with detailed layer extraction (internal C++ function)
void snowpack_physics_layers_internal(double temp_air, double humidity, double wind_speed, double wind_dir,
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
                             // Budget tracking - mass budgets
                             double* mass_precip, double* mass_sublim, double* mass_melt, double* mass_swe, double* mass_refreeze,
                             // Budget tracking - energy budgets
                             double* energy_lw_in, double* energy_lw_out, double* energy_sw_in, double* energy_sw_out,
                             double* energy_sensible, double* energy_latent, double* energy_ground_flux, double* energy_rain, double* energy_total) {
  
  // Track physics calls for debugging
  static int call_count = 0;
  call_count++;

  printf("SNOWPACK-DEBUG: ======== STARTING snowpack_physics_layers_internal CALL #%d ========\n", call_count);
    printf("SNOWPACK-DEBUG: Previous call completed successfully, starting new call...\n");
  printf("SNOWPACK-DEBUG [%d,%d]: Input parameters - temp=%.2f, humidity=%.3f, wind=%.2f, dir=%.1f\n",
         i_grid, j_grid, temp_air, humidity, wind_speed, wind_dir);
  printf("SNOWPACK-DEBUG [%d,%d]: Grid coords - lat=%.6f, lon=%.6f, height=%.1f, dt=%.1f\n",
         i_grid, j_grid, wrf_lat, wrf_lon, height, dt);
  printf("SNOWPACK-DEBUG [%d,%d]: Radiation - sw=%.1f, lw=%.1f, precip=%.3f, pressure=%.1f\n",
         i_grid, j_grid, shortwave_in, longwave_in, precipitation, pressure);
  printf("SNOWPACK-DEBUG [%d,%d]: Current snow state - swe=%.3f, depth=%.3f, temp=%.2f\n",
         i_grid, j_grid, *snow_swe, *snow_depth, *surface_temp);

  // Initialize configuration on first call (shared, read-only)
  try {
    printf("SNOWPACK-DEBUG [%d,%d]: Initializing SNOWPACK configuration...\n", i_grid, j_grid);
    initialize_snowpack_config();
    printf("SNOWPACK-DEBUG [%d,%d]: SNOWPACK configuration initialized successfully\n", i_grid, j_grid);
  } catch (const std::exception& e) {
    printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Configuration failed: %s\n", e.what());
    printf("SNOWPACK-FATAL: Aborting WRF run due to SNOWPACK configuration failure\n");
    std::abort();  // Abort instead of silent fallback
  }

  // CRYOWRF-style persistent approach: Get persistent objects for each grid point
  try {
    // Get persistent SNOWPACK objects for this grid point (following CRYOWRF pattern)
    // Pass domain_id=1 explicitly to avoid parameter mismatch
    printf("SNOWPACK-DEBUG [%d,%d]: About to create/get SNOWPACK objects (lat=%.6f, lon=%.6f, height=%.1f)...\n",
           i_grid, j_grid, wrf_lat, wrf_lon, height);
    SnowStation* snow_station = get_or_create_snowstation(i_grid, j_grid, 1, wrf_lat, wrf_lon, height);
    printf("SNOWPACK-DEBUG [%d,%d]: snow_station created/retrieved: %p\n", i_grid, j_grid, (void*)snow_station);

    Snowpack* snowpack_instance = get_or_create_snowpack_instance(i_grid, j_grid, wrf_lat, wrf_lon, height);
    printf("SNOWPACK-DEBUG [%d,%d]: snowpack_instance created/retrieved: %p\n", i_grid, j_grid, (void*)snowpack_instance);

    // Validate pointers before use
    if (!snow_station) {
        printf("SNOWPACK-FATAL [%d,%d]: snow_station pointer is null!\n", i_grid, j_grid);
        std::abort();
    }
    if (!snowpack_instance) {
        printf("SNOWPACK-FATAL [%d,%d]: snowpack_instance pointer is null!\n", i_grid, j_grid);
        std::abort();
    }

    // Create temporary meteorological data
    printf("SNOWPACK-DEBUG [%d,%d]: Creating temporary data objects...\n", i_grid, j_grid);
    auto Mdata = std::make_unique<CurrentMeteo>();
    auto surfFluxes = std::make_unique<SurfaceFluxes>();
    auto sn_Bdata = std::make_unique<BoundCond>();
    printf("SNOWPACK-DEBUG [%d,%d]: Temporary objects created successfully\n", i_grid, j_grid);
    
    // Note: WRF coordinates may be projected (grid coordinates in meters) rather than lat/lon
    
    // Initialize snow station data with WRF coordinates
    SN_SNOWSOIL_DATA ssdata;
    mio::Coords position;
    // WRF always passes geographic coordinates (XLAT/XLONG in degrees)
    // Use setLatLon with proper order: latitude, longitude, altitude
    position.setLatLon(wrf_lat, wrf_lon, height);  // Geographic coordinates from WRF
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
        // WRF calls physics every wrf_dt, SNOWPACK needs to run every snowpack_dt
        // So we run SNOWPACK every (wrf_dt/snowpack_dt) WRF calls
        // Following CRYOWRF pattern exactly
        compute_counter_layers = (int)(wrf_dt / snowpack_dt);
        first_physics_call_layers = false;
        printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: compute_counter_layers = %d (snowpack_dt=%.1fs, wrf_dt=%.1fs)\n", 
               compute_counter_layers, snowpack_dt, wrf_dt);
    }
    
    call_counter_layers++;
    
    // Only advance time when counter matches (CRYOWRF pattern)
    if (compute_counter_layers > 0 && (call_counter_layers % compute_counter_layers) == 0) {
        current_simulation_date += (calculation_step_length / 1440.0);  // Convert minutes to days
        printf("SNOWPACK-INFO [C++/snowpack_wrf_bridge.cpp]: Time advanced to %s (call %d)\n", 
               current_simulation_date.toString().c_str(), call_counter_layers);
    }
    
    mio::Date current_time = current_simulation_date;  // Use advancing time
    
    // Temperature sanity check
    double safe_temp = std::max(SnowpackConstants::T_CRAZY_MIN_KELVIN, 
                               std::min(temp_air, SnowpackConstants::T_CRAZY_MAX_KELVIN));
    
    // Fill meteorological data structure
    Mdata->date = current_time;
    Mdata->ta = safe_temp;
    Mdata->rh = std::max(0.01, std::min(1.0, humidity));
    Mdata->vw = std::max(0.1, wind_speed);
    Mdata->dw = wind_dir;
    Mdata->iswr = std::max(0.0, shortwave_in);
    Mdata->lw_net = std::max(0.0, longwave_in);
    Mdata->psum = std::max(0.0, precipitation);
    Mdata->psum_ph = (safe_temp < 273.65) ? 0.0 : 1.0;
    Mdata->tss = mio::IOUtils::nodata;
    Mdata->ts0 = safe_temp - 5.0;

    // Initialize z0 and ustar (required for wind pumping in thermal conductivity)
    // Read ROUGHNESS_LENGTH from config (same as SNOWPACK's Meteo::MicroMet)
    static double roughness_length = -1.0;
    if (roughness_length < 0.0) {
        if (global_config) {
            global_config->getValue("ROUGHNESS_LENGTH", "Snowpack", roughness_length, mio::IOUtils::nothrow);
            if (roughness_length < 0.0) {
                roughness_length = 0.002;  // Default if not in config
                printf("SNOWPACK-INFO: ROUGHNESS_LENGTH not found in config, using default %.4f m\n", roughness_length);
            } else {
                printf("SNOWPACK-INFO: Read ROUGHNESS_LENGTH=%.4f m from config\n", roughness_length);
            }
        }
    }
    
    // Adapt roughness based on snow presence (same logic as Meteo::MicroMet line 267)
    const double rough_len = (*snow_depth > 0.03) ? roughness_length : 0.01;  // BareSoil_z0 typically ~0.01
    Mdata->z0 = rough_len;
    
    // Compute friction velocity from wind speed using log profile (same as Meteo::MicroMet line 316)
    // u* = k * u / ln(z/z0), where k=0.4 (von Karman constant)
    const double von_karman = 0.4;
    const double z_wind = height;  // Measurement height [m]
    Mdata->ustar = (z_wind > Mdata->z0) ?
                  (von_karman * Mdata->vw / std::log(z_wind / Mdata->z0)) :
                  (0.1 * Mdata->vw);  // Fallback if z <= z0
    
    // DEBUG: Print z0 and ustar values for first few calls
    if (call_count <= 5) {
        printf("SNOWPACK-DEBUG [%d,%d]: z0=%.6f m, ustar=%.6f m/s, VW=%.2f m/s, height=%.2f m, snow_depth=%.3f m\n",
               i_grid, j_grid, Mdata->z0, Mdata->ustar, Mdata->vw, height, *snow_depth);
    }
    Mdata->hs = *snow_depth;
    
    // Execute SNOWPACK physics
    double cumu_precip = 0.0;
    printf("SNOWPACK-DEBUG [%d,%d]: >>> About to call runSnowpackModel <<<\n", i_grid, j_grid);
    printf("SNOWPACK-DEBUG [%d,%d]:   snowpack_instance=%p, snow_station=%p\n", i_grid, j_grid, (void*)snowpack_instance, (void*)snow_station);
    printf("SNOWPACK-DEBUG [%d,%d]:   Mdata address: %p, sn_Bdata address: %p, surfFluxes address: %p\n",
           i_grid, j_grid, (void*)Mdata.get(), (void*)sn_Bdata.get(), (void*)surfFluxes.get());
    printf("SNOWPACK-DEBUG [%d,%d]:   cumu_precip=%.6f, Mdata->vw=%.2f, Mdata->ta=%.2f, Mdata->rh=%.2f\n",
           i_grid, j_grid, cumu_precip, Mdata->vw, Mdata->ta, Mdata->rh);
    printf("SNOWPACK-DEBUG [%d,%d]:   Mdata->z0=%.6f, Mdata->ustar=%.6f, Mdata->hs=%.3f\n",
           i_grid, j_grid, Mdata->z0, Mdata->ustar, Mdata->hs);

    printf("SNOWPACK-DEBUG [%d,%d]: >>> CALLING snowpack_instance->runSnowpackModel() <<<\n", i_grid, j_grid);
    snowpack_instance->runSnowpackModel(*Mdata, *snow_station, cumu_precip, *sn_Bdata, *surfFluxes);
    printf("SNOWPACK-DEBUG [%d,%d]: >>> runSnowpackModel completed successfully <<<\n", i_grid, j_grid);

    // CRITICAL: Populate SurfaceFluxes with actual SNOWPACK results
    // In SNOWPACK v11.08, collectSurfaceFluxes must be called to populate surfFluxes
    printf("SNOWPACK-DEBUG [%d,%d]: Collecting surface fluxes from SNOWPACK results...\n", i_grid, j_grid);
    printf("SNOWPACK-DEBUG [%d,%d]:   BEFORE collectSurfaceFluxes call:\n", i_grid, j_grid);
    printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes->qs = %.2f (should be 0.0 before collection)\n", i_grid, j_grid, surfFluxes->qs);
    printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes->ql = %.2f (should be 0.0 before collection)\n", i_grid, j_grid, surfFluxes->ql);
    printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes->lw_out = %.2f (should be 0.0 before collection)\n", i_grid, j_grid, surfFluxes->lw_out);
    printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes object address: %p\n", i_grid, j_grid, surfFluxes.get());

    printf("SNOWPACK-DEBUG [%d,%d]:   About to call collectSurfaceFluxes(sn_Bdata, *snow_station, Mdata)...\n", i_grid, j_grid);
    printf("SNOWPACK-DEBUG [%d,%d]:   sn_Bdata address: %p\n", i_grid, j_grid, sn_Bdata.get());
    printf("SNOWPACK-DEBUG [%d,%d]:   snow_station address: %p\n", i_grid, j_grid, snow_station);
    printf("SNOWPACK-DEBUG [%d,%d]:   Mdata address: %p\n", i_grid, j_grid, Mdata.get());

    try {
        surfFluxes->collectSurfaceFluxes(*sn_Bdata, *snow_station, *Mdata);
        printf("SNOWPACK-DEBUG [%d,%d]:   collectSurfaceFluxes call completed without exception\n", i_grid, j_grid);
    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [%d,%d]: Exception in collectSurfaceFluxes: %s\n", i_grid, j_grid, e.what());
        std::abort();
    } catch (...) {
        printf("SNOWPACK-FATAL [%d,%d]: Unknown exception in collectSurfaceFluxes\n", i_grid, j_grid);
        std::abort();
    }

    printf("SNOWPACK-DEBUG [%d,%d]:   AFTER collectSurfaceFluxes call:\n", i_grid, j_grid);
    printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes->qs = %.2f (should be populated now)\n", i_grid, j_grid, surfFluxes->qs);
    printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes->ql = %.2f (should be populated now)\n", i_grid, j_grid, surfFluxes->ql);
    printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes->lw_out = %.2f (THIS IS THE CRITICAL VALUE)\n", i_grid, j_grid, surfFluxes->lw_out);
    printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes->qg = %.2f (ground heat flux)\n", i_grid, j_grid, surfFluxes->qg);
    printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes->sw_in = %.2f (shortwave in)\n", i_grid, j_grid, surfFluxes->sw_in);
    printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes->dIntEnergy = %.2f (internal energy change)\n", i_grid, j_grid, surfFluxes->dIntEnergy);

    printf("SNOWPACK-DEBUG [%d,%d]: Surface fluxes collection analysis complete\n", i_grid, j_grid);

    // CRITICAL: Debug post-SNOWPACK processing step by step
    printf("SNOWPACK-DEBUG [%d,%d]: ======== STARTING POST-PROCESSING VALIDATION ========\n", i_grid, j_grid);

    // STEP 1: Validate object integrity before data extraction
    printf("SNOWPACK-DEBUG [%d,%d]: STEP 1: Validating object integrity...\n", i_grid, j_grid);
    if (!snow_station) {
        printf("SNOWPACK-FATAL [%d,%d]: STEP 1 FAILED - snow_station is NULL!\n", i_grid, j_grid);
        std::abort();
    }
    if (!snowpack_instance) {
        printf("SNOWPACK-FATAL [%d,%d]: STEP 1 FAILED - snowpack_instance is NULL!\n", i_grid, j_grid);
        std::abort();
    }
    printf("SNOWPACK-DEBUG [%d,%d]: STEP 1 PASSED - Objects are valid\n", i_grid, j_grid);

    // STEP 2: Basic result extraction with validation
    printf("SNOWPACK-DEBUG [%d,%d]: STEP 2: Extracting basic results...\n", i_grid, j_grid);
    try {
        printf("SNOWPACK-DEBUG [%d,%d]:   Accessing snow_station->Ndata.size()...\n", i_grid, j_grid);
        size_t ndata_size = snow_station->Ndata.size();
        printf("SNOWPACK-DEBUG [%d,%d]:   snow_station->Ndata.size() = %zu\n", i_grid, j_grid, ndata_size);

        printf("SNOWPACK-DEBUG [%d,%d]:   Extracting surface temperature...\n", i_grid, j_grid);
        *surface_temp = (ndata_size > 0) ? snow_station->Ndata.back().T : temp_air;
        printf("SNOWPACK-DEBUG [%d,%d]:   surface_temp extracted: %.2f\n", i_grid, j_grid, *surface_temp);

        printf("SNOWPACK-DEBUG [%d,%d]:   Extracting SWE and depth...\n", i_grid, j_grid);
        *snow_swe = snow_station->swe;
        *snow_depth = snow_station->cH - snow_station->Ground;  // Subtract ground height (CRYOWRF line 1129)
        printf("SNOWPACK-DEBUG [%d,%d]:   SWE=%.3f, depth=%.3f\n", i_grid, j_grid, *snow_swe, *snow_depth);

        printf("SNOWPACK-DEBUG [%d,%d]:   Extracting surface fluxes...\n", i_grid, j_grid);
        *heat_flux_sensible = -1.0 * surfFluxes->qs;      // Negative sign for WRF convention (CRYOWRF line 1123)
        *heat_flux_latent = -1.0 * sn_Bdata->ql;         // Use boundary condition data (CRYOWRF line 1124)
        printf("SNOWPACK-DEBUG [%d,%d]:   Sensible=%.2f, Latent=%.2f\n", i_grid, j_grid, *heat_flux_sensible, *heat_flux_latent);

        printf("SNOWPACK-DEBUG [%d,%d]:   Extracting albedo and snow coverage...\n", i_grid, j_grid);
        *albedo = snow_station->Albedo;
        *snow_coverage = 1.0;  // Hardcoded to 1.0 following CRYOWRF pattern (line 1127)
        printf("SNOWPACK-DEBUG [%d,%d]:   Albedo=%.3f, Coverage=%.3f\n", i_grid, j_grid, *albedo, *snow_coverage);

        printf("SNOWPACK-DEBUG [%d,%d]: STEP 2 COMPLETED - Basic results extracted\n", i_grid, j_grid);

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [%d,%d]: STEP 2 CRASHED - Exception: %s\n", i_grid, j_grid, e.what());
        std::abort();
    }

    printf("SNOWPACK-DEBUG [%d,%d]: ======== POST-PROCESSING STEP 2 COMPLETE ========\n", i_grid, j_grid);

    // STEP 3: Layer data extraction with validation
    printf("SNOWPACK-DEBUG [%d,%d]: STEP 3: Extracting detailed layer information...\n", i_grid, j_grid);
    try {
        printf("SNOWPACK-DEBUG [%d,%d]:   Calling snow_station->getNumberOfElements()...\n", i_grid, j_grid);
        size_t num_elements = snow_station->getNumberOfElements();
        printf("SNOWPACK-DEBUG [%d,%d]:   Number of elements: %zu\n", i_grid, j_grid, num_elements);
        *n_layers = static_cast<int>(num_elements);
        printf("SNOWPACK-DEBUG [%d,%d]:   n_layers set to: %d\n", i_grid, j_grid, *n_layers);

        // Limit to max 50 layers for WRF arrays
        printf("SNOWPACK-DEBUG [%d,%d]:   Calculating layers to extract (max 50)...\n", i_grid, j_grid);
        size_t layers_to_extract = std::min(num_elements, size_t(50));
        printf("SNOWPACK-DEBUG [%d,%d]:   Will extract %zu layers\n", i_grid, j_grid, layers_to_extract);

        printf("SNOWPACK-DEBUG [%d,%d]:   Starting layer extraction loop...\n", i_grid, j_grid);
        for (size_t i = 0; i < layers_to_extract; i++) {
            printf("SNOWPACK-DEBUG [%d,%d]:   Processing layer %zu/%zu...\n", i_grid, j_grid, i, layers_to_extract);

            printf("SNOWPACK-DEBUG [%d,%d]:     Accessing Edata[%zu]...\n", i_grid, j_grid, i);
            const ElementData& elem = snow_station->Edata[i];

            printf("SNOWPACK-DEBUG [%d,%d]:     Extracting basic properties...\n", i_grid, j_grid);
            layer_temp[i] = elem.Te;
            layer_thick[i] = elem.L;
            layer_vol_ice[i] = elem.theta[ICE] * 100.0;
            layer_vol_water[i] = elem.theta[WATER] * 100.0;
            layer_vol_air[i] = elem.theta[AIR] * 100.0;

            printf("SNOWPACK-DEBUG [%d,%d]:     Extracting grain properties...\n", i_grid, j_grid);
            layer_grain_radius[i] = elem.rg;
            layer_bond_radius[i] = elem.rb;
            layer_dendricity[i] = elem.dd;
            layer_sphericity[i] = elem.sp;

            if (i < 3) {  // Only print first few layers to avoid spam
                printf("SNOWPACK-DEBUG [%d,%d]:     Layer %zu: T=%.1f, L=%.3f, ice=%.1f%%, rg=%.3f\n",
                       i_grid, j_grid, i, layer_temp[i], layer_thick[i], layer_vol_ice[i], layer_grain_radius[i]);
            }
        }
        printf("SNOWPACK-DEBUG [%d,%d]:   Layer extraction completed successfully\n", i_grid, j_grid);
        printf("SNOWPACK-DEBUG [%d,%d]: STEP 3 COMPLETED - Layer data extracted\n", i_grid, j_grid);

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [%d,%d]: STEP 3 CRASHED - Exception: %s\n", i_grid, j_grid, e.what());
        printf("SNOWPACK-FATAL [%d,%d]: This indicates layer extraction is the crash point!\n", i_grid, j_grid);
        std::abort();
    }

    printf("SNOWPACK-DEBUG [%d,%d]: ======== POST-PROCESSING STEP 3 COMPLETE ========\n", i_grid, j_grid);

    // STEP 4: Budget calculations and final output validation
    printf("SNOWPACK-DEBUG [%d,%d]: STEP 4: Computing budget calculations...\n", i_grid, j_grid);
    try {
        printf("SNOWPACK-DEBUG [%d,%d]:   Calculating mass budgets...\n", i_grid, j_grid);

        // Mass budgets (following CRYOWRF patterns) - SAFE ACCESS PATTERN
        *mass_precip = cumu_precip;                           // Cumulative precipitation [kg/m²]

        printf("SNOWPACK-DEBUG [%d,%d]:   About to access surfFluxes->mass[MS_SUBLIMATION]...\n", i_grid, j_grid);
        printf("SNOWPACK-DEBUG [%d,%d]:   surfFluxes->mass array address: %p\n", i_grid, j_grid, surfFluxes->mass.data());
        *mass_sublim = surfFluxes->mass[SurfaceFluxes::MS_SUBLIMATION];  // Sublimation [kg/m²] (CRYOWRF line 1136)
        printf("SNOWPACK-DEBUG [%d,%d]:   MS_SUBLIMATION extracted: %.6f\n", i_grid, j_grid, *mass_sublim);

        // Calculate melt mass from SNOWPACK v11.08 SurfaceFluxes mass balance
        printf("SNOWPACK-DEBUG [%d,%d]:   Calculating melt mass from SurfaceFluxes...\n", i_grid, j_grid);
        *mass_melt = surfFluxes->mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF]; // Melt runoff [kg/m²] - SNOWPACK v11.08 API
        printf("SNOWPACK-DEBUG [%d,%d]:   MS_SNOWPACK_RUNOFF extracted: %.6f\n", i_grid, j_grid, *mass_melt);

        // Check some basic flux values that should be populated by collectSurfaceFluxes
        printf("SNOWPACK-DEBUG [%d,%d]:   Checking basic SurfaceFluxes values:\n", i_grid, j_grid);
        printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes->qs = %.2f (should be non-zero if populated)\n", i_grid, j_grid, surfFluxes->qs);
        printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes->ql = %.2f (should be non-zero if populated)\n", i_grid, j_grid, surfFluxes->ql);
        printf("SNOWPACK-DEBUG [%d,%d]:     surfFluxes->sw_in = %.2f (should match input if populated)\n", i_grid, j_grid, surfFluxes->sw_in);

        *mass_swe = snow_station->swe;                        // Current SWE [kg/m²]

        // Calculate refreeze mass from ice base melting/freezing or negative runoff
        // MS_ICEBASE_MELTING_FREEZING: mass gain/loss of ice base due to melting-freezing
        // Also consider negative MS_SNOWPACK_RUNOFF which indicates refreezing
        double ice_base_meltfreeze = surfFluxes->mass[SurfaceFluxes::MS_ICEBASE_MELTING_FREEZING];
        *mass_refreeze = (ice_base_meltfreeze > 0.0) ? ice_base_meltfreeze : 0.0;  // Positive values = freezing
        printf("SNOWPACK-DEBUG [%d,%d]:   Refreeze mass calculated: %.6f (ice_base_meltfreeze=%.6f)\n",
               i_grid, j_grid, *mass_refreeze, ice_base_meltfreeze);

        printf("SNOWPACK-DEBUG [%d,%d]:   Mass: precip=%.6f, sublim=%.6f, melt=%.6f, swe=%.6f\n",
               i_grid, j_grid, *mass_precip, *mass_sublim, *mass_melt, *mass_swe);

        printf("SNOWPACK-DEBUG [%d,%d]:   Calculating energy budgets...\n", i_grid, j_grid);

        // Energy budgets - SAFE ACCESS PATTERN
        *energy_lw_in = longwave_in;                           // Incoming longwave [W/m²]

        printf("SNOWPACK-DEBUG [%d,%d]:   About to access surfFluxes->lw_out...\n", i_grid, j_grid);
        printf("SNOWPACK-DEBUG [%d,%d]:   surfFluxes object address: %p\n", i_grid, j_grid, surfFluxes.get());
        printf("SNOWPACK-DEBUG [%d,%d]:   surfFluxes->lw_out address: %p\n", i_grid, j_grid, &(surfFluxes->lw_out));
        printf("SNOWPACK-DEBUG [%d,%d]:   surfFluxes object size: %zu bytes\n", i_grid, j_grid, sizeof(surfFluxes));

        // Memory validation - check if surfFluxes object is properly aligned
        uintptr_t addr = reinterpret_cast<uintptr_t>(surfFluxes.get());
        printf("SNOWPACK-DEBUG [%d,%d]:   Memory alignment check: address=%p, alignment=%zu\n", i_grid, j_grid, (void*)addr, alignof(SurfaceFluxes));
        bool alignment_ok = (addr % alignof(SurfaceFluxes) == 0);
        printf("SNOWPACK-DEBUG [%d,%d]:   Alignment status: %s\n", i_grid, j_grid, alignment_ok ? "OK" : "BAD");

        // Validate stack frame integrity
        void* stack_ptr = __builtin_frame_address(0);
        printf("SNOWPACK-DEBUG [%d,%d]:   Stack frame address: %p\n", i_grid, j_grid, stack_ptr);
        printf("SNOWPACK-DEBUG [%d,%d]:   Distance from stack to surfFluxes: %ld bytes\n", i_grid, j_grid, (long)stack_ptr - (long)surfFluxes.get());

        // Check other surfFluxes members to see if object is partially corrupted
        printf("SNOWPACK-DEBUG [%d,%d]:   Checking surfFluxes member integrity:\n", i_grid, j_grid);
        printf("SNOWPACK-DEBUG [%d,%d]:     qs (sensible): addr=%p, value=%.2f\n", i_grid, j_grid, &(surfFluxes->qs), surfFluxes->qs);
        printf("SNOWPACK-DEBUG [%d,%d]:     ql (latent): addr=%p, value=%.2f\n", i_grid, j_grid, &(surfFluxes->ql), surfFluxes->ql);
        printf("SNOWPACK-DEBUG [%d,%d]:     qg (ground): addr=%p, value=%.2f\n", i_grid, j_grid, &(surfFluxes->qg), surfFluxes->qg);
        printf("SNOWPACK-DEBUG [%d,%d]:     sw_in: addr=%p, value=%.2f\n", i_grid, j_grid, &(surfFluxes->sw_in), surfFluxes->sw_in);
        printf("SNOWPACK-DEBUG [%d,%d]:     qr (rain): addr=%p, value=%.2f\n", i_grid, j_grid, &(surfFluxes->qr), surfFluxes->qr);
        printf("SNOWPACK-DEBUG [%d,%d]:     dIntEnergy: addr=%p, value=%.2f\n", i_grid, j_grid, &(surfFluxes->dIntEnergy), surfFluxes->dIntEnergy);

        // Check memory around lw_out to detect corruption
        printf("SNOWPACK-DEBUG [%d,%d]:   Memory analysis around lw_out:\n", i_grid, j_grid);
        double* lw_out_ptr = &(surfFluxes->lw_out);
        printf("SNOWPACK-DEBUG [%d,%d]:     lw_out pointer value: %p\n", i_grid, j_grid, lw_out_ptr);

        // Try to read memory as raw bytes to check if it's readable
        volatile bool memory_readable = true;
        try {
            volatile double test_val = *lw_out_ptr;
            printf("SNOWPACK-DEBUG [%d,%d]:     Raw memory read successful: test_val=%.2f\n", i_grid, j_grid, test_val);
        } catch (...) {
            memory_readable = false;
            printf("SNOWPACK-FATAL [%d,%d]:     Raw memory read FAILED - memory corruption detected!\n", i_grid, j_grid);
        }

        if (!memory_readable) {
            printf("SNOWPACK-FATAL [%d,%d]: CRITICAL: Memory at lw_out location is not readable!\n", i_grid, j_grid);
            printf("SNOWPACK-FATAL [%d,%d]: This indicates severe memory corruption or stack overflow!\n", i_grid, j_grid);
            std::abort();
        }

    // CRITICAL: Fix memory access pattern for SurfaceFluxes
    // Copy SurfaceFluxes to heap-allocated memory to avoid stack corruption
    printf("SNOWPACK-DEBUG [%d,%d]: Creating heap-allocated SurfaceFluxes backup...\n", i_grid, j_grid);
    std::unique_ptr<SurfaceFluxes> safe_fluxes_backup = std::make_unique<SurfaceFluxes>();

    // Immediately copy the populated SurfaceFluxes to safe storage
    *safe_fluxes_backup = *surfFluxes;
    printf("SNOWPACK-DEBUG [%d,%d]: SurfaceFluxes safely copied to heap memory\n", i_grid, j_grid);

    // Validate the copy worked before proceeding
    printf("SNOWPACK-DEBUG [%d,%d]: Validating copied SurfaceFluxes...\n", i_grid, j_grid);
    printf("SNOWPACK-DEBUG [%d,%d]:   Backup qs=%.2f, ql=%.2f, lw_out=%.2f\n",
           i_grid, j_grid, safe_fluxes_backup->qs, safe_fluxes_backup->ql, safe_fluxes_backup->lw_out);
    printf("SNOWPACK-DEBUG [%d,%d]:   Original qs=%.2f, ql=%.2f, lw_out=%.2f\n",
           i_grid, j_grid, surfFluxes->qs, surfFluxes->ql, surfFluxes->lw_out);

    printf("SNOWPACK-DEBUG [%d,%d]:   Attempting to read safe_fluxes_backup->lw_out...\n", i_grid, j_grid);
    printf("SNOWPACK-DEBUG [%d,%d]:   safe_fluxes_backup address: %p\n", i_grid, j_grid, safe_fluxes_backup.get());
    printf("SNOWPACK-DEBUG [%d,%d]:   safe_fluxes_backup->lw_out address: %p\n", i_grid, j_grid, &safe_fluxes_backup->lw_out);
    printf("SNOWPACK-DEBUG [%d,%d]:   safe_fluxes_backup object size: %zu bytes\n", i_grid, j_grid, sizeof(*safe_fluxes_backup));

    // Memory validation - check if safe_fluxes_backup object is properly aligned
    uintptr_t safe_addr = reinterpret_cast<uintptr_t>(safe_fluxes_backup.get());
    printf("SNOWPACK-DEBUG [%d,%d]:   Memory alignment check: address=%p, alignment=%zu\n", i_grid, j_grid, (void*)safe_addr, alignof(SurfaceFluxes));
    bool safe_alignment_ok = (safe_addr % alignof(SurfaceFluxes) == 0);
    printf("SNOWPACK-DEBUG [%d,%d]:   Alignment status: %s\n", i_grid, j_grid, safe_alignment_ok ? "OK" : "BAD");

    // Check memory around lw_out to detect corruption
    printf("SNOWPACK-DEBUG [%d,%d]:   Memory analysis around safe lw_out:\n", i_grid, j_grid);
    double* safe_lw_out_ptr = &safe_fluxes_backup->lw_out;
    printf("SNOWPACK-DEBUG [%d,%d]:     safe lw_out pointer value: %p\n", i_grid, j_grid, safe_lw_out_ptr);

    // Try to read memory as raw bytes to check if it's readable
    volatile bool safe_memory_readable = true;
    try {
        volatile double safe_test_val = *safe_lw_out_ptr;
        printf("SNOWPACK-DEBUG [%d,%d]:     Raw memory read successful: test_val=%.2f\n", i_grid, j_grid, safe_test_val);
    } catch (...) {
        safe_memory_readable = false;
        printf("SNOWPACK-FATAL [%d,%d]:     Raw memory read FAILED - memory corruption detected!\n", i_grid, j_grid);
    }

    if (!safe_memory_readable) {
        printf("SNOWPACK-FATAL [%d,%d]: CRITICAL: Heap memory at lw_out location is not readable!\n", i_grid, j_grid);
        std::abort();
    }

    // Now safely extract from the heap-allocated copy instead of stack memory
    try {
        printf("SNOWPACK-DEBUG [%d,%d]:   About to execute: *energy_lw_out = safe_fluxes_backup->lw_out\n", i_grid, j_grid);
        *energy_lw_out = safe_fluxes_backup->lw_out;                   // Outgoing longwave [W/m²] (CRYOWRF line 1087)
        printf("SNOWPACK-DEBUG [%d,%d]:   Assignment completed: energy_lw_out = %.2f\n", i_grid, j_grid, *energy_lw_out);
        printf("SNOWPACK-DEBUG [%d,%d]:   SUCCESS: lw_out extracted from safe heap backup: %.2f\n", i_grid, j_grid, *energy_lw_out);
    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [%d,%d]: C++ Exception accessing safe_fluxes_backup->lw_out: %s\n", i_grid, j_grid, e.what());
        printf("SNOWPACK-FATAL [%d,%d]: Exception type: %s\n", i_grid, j_grid, typeid(e).name());
        std::abort();
    } catch (...) {
        printf("SNOWPACK-FATAL [%d,%d]: Unknown non-C++ exception accessing safe_fluxes_backup->lw_out\n", i_grid, j_grid);
        std::abort();
    }

        *energy_sw_in = shortwave_in;                          // Incoming shortwave [W/m²]
        *energy_sw_out = shortwave_in * (1.0 - *albedo);      // Reflected shortwave [W/m²]

        // Turbulent fluxes - use safe heap copy to avoid stack corruption
        *energy_sensible = safe_fluxes_backup->qs;                      // Sensible heat flux [W/m²]
        *energy_latent = safe_fluxes_backup->ql;                        // Latent heat flux [W/m²]

        // Ground heat flux from surface fluxes (following CRYOWRF pattern)
        *energy_ground_flux = safe_fluxes_backup->qg;

        // Rain energy (following CRYOWRF line 1149)
        printf("SNOWPACK-DEBUG [%d,%d]:   About to access safe_fluxes_backup->qr...\n", i_grid, j_grid);
        *energy_rain = safe_fluxes_backup->qr;                          // Rain heat flux [W/m²] (CRYOWRF pattern)
        printf("SNOWPACK-DEBUG [%d,%d]:   qr extracted: %.2f\n", i_grid, j_grid, *energy_rain);

        // Total energy balance (following CRYOWRF line 1150)
        printf("SNOWPACK-DEBUG [%d,%d]:   About to access safe_fluxes_backup->dIntEnergy...\n", i_grid, j_grid);
        *energy_total = safe_fluxes_backup->dIntEnergy / dt;           // Total energy [W/m²] (CRYOWRF pattern)
        printf("SNOWPACK-DEBUG [%d,%d]:   dIntEnergy extracted: %.2f, total=%.2f\n", i_grid, j_grid, safe_fluxes_backup->dIntEnergy, *energy_total);

        printf("SNOWPACK-DEBUG [%d,%d]:   Energy: SW_in=%.1f, LW_in=%.1f, Sensible=%.1f, Latent=%.1f, Ground=%.1f\n",
               i_grid, j_grid, *energy_sw_in, *energy_lw_in, *energy_sensible, *energy_latent, *energy_ground_flux);
        printf("SNOWPACK-DEBUG [%d,%d]: STEP 4 COMPLETED - Budget calculations done\n", i_grid, j_grid);

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [%d,%d]: STEP 4 CRASHED - Exception: %s\n", i_grid, j_grid, e.what());
        printf("SNOWPACK-FATAL [%d,%d]: This indicates budget calculation is the crash point!\n", i_grid, j_grid);
        std::abort();
    }

    printf("SNOWPACK-DEBUG [%d,%d]: ======== POST-PROCESSING STEP 4 COMPLETE ========\n", i_grid, j_grid);

    // STEP 5: Final output validation before return to Fortran
    printf("SNOWPACK-DEBUG [%d,%d]: STEP 5: Final output validation...\n", i_grid, j_grid);
    try {
        printf("SNOWPACK-DEBUG [%d,%d]:   All calculations completed, preparing to return to Fortran\n", i_grid, j_grid);
        printf("SNOWPACK-DEBUG [%d,%d]: STEP 5 COMPLETED - Ready to return to Fortran\n", i_grid, j_grid);

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [%d,%d]: STEP 5 CRASHED - Exception: %s\n", i_grid, j_grid, e.what());
        printf("SNOWPACK-FATAL [%d,%d]: This indicates final validation is the crash point!\n", i_grid, j_grid);
        std::abort();
    }

    printf("SNOWPACK-DEBUG [%d,%d]: ======== ALL POST-PROCESSING STEPS COMPLETED ========\n", i_grid, j_grid);

    // Clear unused layers (initialize remaining array elements)
    printf("SNOWPACK-DEBUG [%d,%d]: Clearing unused layer array elements...\n", i_grid, j_grid);
    for (size_t i = static_cast<size_t>(*n_layers); i < 50; i++) {
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
    printf("SNOWPACK-DEBUG [%d,%d]: Array cleanup completed\n", i_grid, j_grid);

    // STEP 6: Final summary before returning to Fortran
    printf("SNOWPACK-DEBUG [%d,%d]: STEP 6: Final summary before return to Fortran...\n", i_grid, j_grid);
    if (call_count <= 5) {
      printf("SNOWPACK-LAYERS [C++]: Grid (%d,%d) - %d layers, T_sfc=%.1fK, SWE=%.2fmm, depth=%.2fm\n",
             i_grid, j_grid, *n_layers, *surface_temp, *snow_swe, *snow_depth);
    }

    printf("SNOWPACK-DEBUG [%d,%d]: FINAL OUTPUT VALIDATION BEFORE RETURN TO FORTRAN:\n", i_grid, j_grid);
    printf("SNOWPACK-DEBUG [%d,%d]:   surface_temp=%.2f (was %.2f), snow_swe=%.3f (was %.3f), snow_depth=%.3f (was %.3f)\n",
           i_grid, j_grid, *surface_temp, temp_air, *snow_swe, *snow_swe, *snow_depth, *snow_depth);
    printf("SNOWPACK-DEBUG [%d,%d]:   heat_flux_sensible=%.2f, heat_flux_latent=%.2f, albedo=%.3f, snow_coverage=%.3f\n",
           i_grid, j_grid, *heat_flux_sensible, *heat_flux_latent, *albedo, *snow_coverage);
    printf("SNOWPACK-DEBUG [%d,%d]:   n_layers=%d, mass_precip=%.6f, energy_sensible=%.2f\n",
           i_grid, j_grid, *n_layers, *mass_precip, *energy_sensible);

    // Check for invalid output values that could cause Fortran crashes
    if (std::isnan(*surface_temp) || std::isinf(*surface_temp)) {
        printf("SNOWPACK-FATAL [%d,%d]: surface_temp is NaN/Inf (%.6f)!\n", i_grid, j_grid, *surface_temp);
        std::abort();
    }
    if (std::isnan(*snow_swe) || std::isinf(*snow_swe)) {
        printf("SNOWPACK-FATAL [%d,%d]: snow_swe is NaN/Inf (%.6f)!\n", i_grid, j_grid, *snow_swe);
        std::abort();
    }
    if (*n_layers < 0 || *n_layers > 50) {
        printf("SNOWPACK-FATAL [%d,%d]: n_layers is invalid (%d)!\n", i_grid, j_grid, *n_layers);
        std::abort();
    }

    printf("SNOWPACK-DEBUG [%d,%d]: ======== snowpack_physics_layers_internal CALL #%d COMPLETED SUCCESSFULLY ========\n",
           i_grid, j_grid, call_count);

    // DEBUG: Memory status check after processing (helps identify corruption between calls)
    printf("SNOWPACK-DEBUG [%d,%d]: Memory status check:\n", i_grid, j_grid);
    printf("SNOWPACK-DEBUG [%d,%d]:   snowpack_instance pointer: %p (still valid)\n", i_grid, j_grid, snowpack_instance);
    printf("SNOWPACK-DEBUG [%d,%d]:   snow_station pointer: %p (still valid)\n", i_grid, j_grid, snow_station);

    // DEBUG: Validate object integrity after use
    printf("SNOWPACK-DEBUG [%d,%d]: Object integrity validation:\n", i_grid, j_grid);
    if (snowpack_instance) {
      printf("SNOWPACK-DEBUG [%d,%d]:   snowpack_instance appears valid\n", i_grid, j_grid);
    } else {
      printf("SNOWPACK-FATAL [%d,%d]:   snowpack_instance is NULL after processing!\n", i_grid, j_grid);
    }

    if (snow_station) {
      printf("SNOWPACK-DEBUG [%d,%d]:   snow_station appears valid, swe=%.3f, depth=%.3f\n",
             i_grid, j_grid, snow_station->swe, snow_station->cH);
    } else {
      printf("SNOWPACK-FATAL [%d,%d]:   snow_station is NULL after processing!\n", i_grid, j_grid);
    }

    printf("SNOWPACK-DEBUG [%d,%d]: Next grid point will be: ", i_grid, j_grid);
    if (j_grid < 100) {
      printf("(%d,%d)\n", i_grid, j_grid+1);
    } else if (i_grid < 200) {
      printf("(%d,%d)\n", i_grid+1, 1);
    } else {
      printf("END OF GRID\n");
    }
    printf("SNOWPACK-DEBUG: ==============================================================================\n");
    printf("SNOWPACK-DEBUG: \n");

  } catch (const std::exception& e) {
    // DEBUG: Enhanced exception handling
    printf("SNOWPACK-FATAL [%d,%d]: >>> C++ EXCEPTION CAUGHT <<<\n", i_grid, j_grid);
    printf("SNOWPACK-FATAL [%d,%d]: Exception type: %s\n", i_grid, j_grid, typeid(e).name());
    printf("SNOWPACK-FATAL [%d,%d]: Exception message: %s\n", i_grid, j_grid, e.what());
    printf("SNOWPACK-FATAL [%d,%d]: This indicates a critical error in SNOWPACK processing!\n", i_grid, j_grid);

    // DEBUG: Try to identify what caused the exception
    printf("SNOWPACK-FATAL [%d,%d]: Context at time of crash:\n", i_grid, j_grid);
    printf("SNOWPACK-FATAL [%d,%d]:   temp_air=%.2f, wind_speed=%.2f, n_layers_requested=%d\n",
           i_grid, j_grid, temp_air, wind_speed, *n_layers);
    printf("SNOWPACK-FATAL [%d,%d]:   Note: Exception occurred during physics processing\n",
           i_grid, j_grid);

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

// Extern C wrapper for Fortran BIND(C) interface
extern "C" void snowpack_physics_layers(double temp_air, double humidity, double wind_speed, double wind_dir,
                             double shortwave_in, double longwave_in, double precipitation, double pressure, double height, double dt,
                             int i_grid, int j_grid, double wrf_lat, double wrf_lon,
                             double* snow_swe, double* snow_depth, double* surface_temp,
                             double* heat_flux_sensible, double* heat_flux_latent, double* albedo, double* snow_coverage,
                             int* n_layers,
                             double* layer_temp, double* layer_thick,
                             double* layer_vol_ice, double* layer_vol_water, double* layer_vol_air,
                             double* layer_grain_radius, double* layer_bond_radius,
                             double* layer_dendricity, double* layer_sphericity,
                             double* mass_precip, double* mass_sublim, double* mass_melt, double* mass_swe, double* mass_refreeze,
                             double* energy_lw_in, double* energy_lw_out, double* energy_sw_in, double* energy_sw_out,
                             double* energy_sensible, double* energy_latent, double* energy_ground_flux, double* energy_rain, double* energy_total) {
    physics_layers_calls++;
    printf("SNOWPACK-DEBUG [CALL #%d]: snowpack_physics_layers() called for grid (%d,%d), lat=%.6f, lon=%.6f\n",
           physics_layers_calls, i_grid, j_grid, wrf_lat, wrf_lon);
    printf("SNOWPACK-DEBUG [CALL #%d]:   temp=%.2f, wind=%.2f, precip=%.3f, dt=%.1f\n",
           physics_layers_calls, temp_air, wind_speed, precipitation, dt);

    // CRITICAL: Add detailed debug checkpoints to isolate segfault location
    printf("SNOWPACK-DEBUG [CALL #%d]: CHECKPOINT 1 - Starting pointer validation\n", physics_layers_calls);

    // CRITICAL: Stack protection for Fortran-C++ bridge
    // Validate only essential output pointers that must not be NULL
    printf("SNOWPACK-DEBUG [CALL #%d]: CHECKPOINT 2 - Pointer validation complete\n", physics_layers_calls);

    if (!snow_swe || !snow_depth || !surface_temp || !energy_lw_out) {
        printf("SNOWPACK-FATAL [CALL #%d]: Essential NULL pointer detected - snow_swe=%p, snow_depth=%p, surface_temp=%p, energy_lw_out=%p\n",
               physics_layers_calls, (void*)snow_swe, (void*)snow_depth, (void*)surface_temp, (void*)energy_lw_out);
        std::abort();
    }

    printf("SNOWPACK-DEBUG [CALL #%d]: CHECKPOINT 3 - Starting defensive initialization\n", physics_layers_calls);

    // Initialize all output arrays to safe values before processing (defensive initialization)
    try {
        printf("SNOWPACK-DEBUG [CALL #%d]: CHECKPOINT 4 - Starting essential output initialization\n", physics_layers_calls);

        // Essential outputs (these should never be NULL)
        *snow_swe = 0.0;
        printf("SNOWPACK-DEBUG [CALL #%d]: CHECKPOINT 5 - After snow_swe assignment\n", physics_layers_calls);
        *snow_depth = 0.0;
        printf("SNOWPACK-DEBUG [CALL #%d]: CHECKPOINT 6 - After snow_depth assignment\n", physics_layers_calls);
        *surface_temp = 273.15;  // Initialize to freezing point
        printf("SNOWPACK-DEBUG [CALL #%d]: CHECKPOINT 7 - After surface_temp assignment\n", physics_layers_calls);
        *energy_lw_out = 0.0;
        printf("SNOWPACK-DEBUG [CALL #%d]: CHECKPOINT 8 - After energy_lw_out assignment\n", physics_layers_calls);

        printf("SNOWPACK-DEBUG [CALL #%d]: CHECKPOINT 9 - Starting optional output initialization\n", physics_layers_calls);

        // Optional outputs (initialize only if pointers are not NULL)
        if (heat_flux_sensible) *heat_flux_sensible = 0.0;
        if (heat_flux_latent) *heat_flux_latent = 0.0;
        if (albedo) *albedo = 0.3;
        if (snow_coverage) *snow_coverage = 0.0;
        if (n_layers) *n_layers = 0;
        if (mass_precip) *mass_precip = 0.0;
        if (mass_sublim) *mass_sublim = 0.0;
        if (mass_melt) *mass_melt = 0.0;
        if (mass_swe) *mass_swe = 0.0;
        if (mass_refreeze) *mass_refreeze = 0.0;
        if (energy_lw_in) *energy_lw_in = longwave_in;
        if (energy_sw_in) *energy_sw_in = shortwave_in;
        if (energy_sw_out) *energy_sw_out = 0.0;
        if (energy_sensible) *energy_sensible = 0.0;
        if (energy_latent) *energy_latent = 0.0;
        if (energy_ground_flux) *energy_ground_flux = 0.0;
        if (energy_rain) *energy_rain = 0.0;
        if (energy_total) *energy_total = 0.0;

        // Initialize layer arrays to safe values (only if pointers are not NULL)
        if (layer_temp && layer_thick && layer_vol_ice && layer_vol_water && layer_vol_air &&
            layer_grain_radius && layer_bond_radius && layer_dendricity && layer_sphericity) {
            for (int i = 0; i < 50; i++) {
                layer_temp[i] = temp_air;
                layer_thick[i] = 0.0;
                layer_vol_ice[i] = 0.0;
                layer_vol_water[i] = 0.0;
                layer_vol_air[i] = 0.0;
                layer_grain_radius[i] = 0.0;
                layer_bond_radius[i] = 0.0;
                layer_dendricity[i] = 0.0;
                layer_sphericity[i] = 1.0;
            }
        }
    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [CALL #%d]: Exception initializing output arrays: %s\n", physics_layers_calls, e.what());
        std::abort();
    }

    // Call internal C++ function with proper coordinates from WRF
    snowpack_physics_layers_internal(temp_air, humidity, wind_speed, wind_dir,
                            shortwave_in, longwave_in, precipitation, pressure, height, dt,
                            i_grid, j_grid, wrf_lat, wrf_lon,  // Use actual WRF coordinates
                            snow_swe, snow_depth, surface_temp,
                            heat_flux_sensible, heat_flux_latent, albedo, snow_coverage,
                            n_layers,
                            layer_temp, layer_thick,
                            layer_vol_ice, layer_vol_water, layer_vol_air,
                            layer_grain_radius, layer_bond_radius,
                            layer_dendricity, layer_sphericity,
                            mass_precip, mass_sublim, mass_melt, mass_swe, mass_refreeze,
                            energy_lw_in, energy_lw_out, energy_sw_in, energy_sw_out,
                            energy_sensible, energy_latent, energy_ground_flux, energy_rain, energy_total);

    // Final output summary before returning to Fortran
    printf("SNOWPACK-DEBUG [CALL #%d]: ====== FINAL OUTPUT SUMMARY ======\n", physics_layers_calls);
    printf("SNOWPACK-DEBUG [CALL #%d]: PRIMARY OUTPUTS (essential):\n", physics_layers_calls);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *snow_swe = %.6f mm\n", physics_layers_calls, snow_swe ? *snow_swe : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *snow_depth = %.6f m\n", physics_layers_calls, snow_depth ? *snow_depth : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *surface_temp = %.2f°C\n", physics_layers_calls, surface_temp ? *surface_temp - 273.15 : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *energy_lw_out = %.2f W/m² (CRITICAL VALUE)\n", physics_layers_calls, energy_lw_out ? *energy_lw_out : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]: SECONDARY OUTPUTS:\n", physics_layers_calls);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *heat_flux_sensible = %.2f W/m²\n", physics_layers_calls, heat_flux_sensible ? *heat_flux_sensible : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *heat_flux_latent = %.2f W/m²\n", physics_layers_calls, heat_flux_latent ? *heat_flux_latent : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *albedo = %.3f\n", physics_layers_calls, albedo ? *albedo : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *snow_coverage = %.3f\n", physics_layers_calls, snow_coverage ? *snow_coverage : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *n_layers = %d\n", physics_layers_calls, n_layers ? *n_layers : -999);
    printf("SNOWPACK-DEBUG [CALL #%d]: MASS BUDGET OUTPUTS:\n", physics_layers_calls);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *mass_precip = %.6f kg/m²\n", physics_layers_calls, mass_precip ? *mass_precip : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *mass_sublim = %.6f kg/m²\n", physics_layers_calls, mass_sublim ? *mass_sublim : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *mass_melt = %.6f kg/m²\n", physics_layers_calls, mass_melt ? *mass_melt : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *mass_swe = %.6f kg/m²\n", physics_layers_calls, mass_swe ? *mass_swe : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]:   *mass_refreeze = %.6f kg/m²\n", physics_layers_calls, mass_refreeze ? *mass_refreeze : -999.0);
    printf("SNOWPACK-DEBUG [CALL #%d]: ====== END OUTPUT SUMMARY ======\n", physics_layers_calls);

    printf("SNOWPACK-DEBUG [CALL #%d]: snowpack_physics_layers() completed successfully for grid (%d,%d)\n",
           physics_layers_calls, i_grid, j_grid);

    // Print summary of calls so far
    printf("SNOWPACK-DEBUG: ======== CALL SUMMARY ========\n");
    printf("SNOWPACK-DEBUG: Config init calls: %d\n", config_init_calls);
    printf("SNOWPACK-DEBUG: Time init calls: %d\n", time_init_calls);
    printf("SNOWPACK-DEBUG: Station creation calls: %d\n", station_creation_calls);
    printf("SNOWPACK-DEBUG: Physics layers calls: %d\n", physics_layers_calls);
    printf("SNOWPACK-DEBUG: Total stations in memory: %zu\n", grid_snowstations.size());
    printf("SNOWPACK-DEBUG: ================================\n");
}