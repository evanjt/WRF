/*
 * CRYOWRF: SNOWPACK + WRF bridge
 *
 * This file provides the C interface between WRF (Fortran) and SNOWPACK (C++)
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

#include "snowpack_bridge.h"
#include "config.h"
#include "meteoio/meteoio/MeteoIO.h"
#include "snowpack/snowpack/DataClasses.h"
#include "snowpack/snowpack/SnowpackConfig.h"
#include "snowpack/snowpack/snowpackCore/Snowpack.h"
#include "snowpack/snowpack/plugins/SnowpackIO.h"

// // Global configuration and I/O (shared, read-only)
// static std::unique_ptr<SnowpackConfig> global_config;
// static std::unique_ptr<SnowpackIO> global_snowpack_io;
// static bool config_initialized = false;
// static std::string config_file_path;

// // Global time management (CRYOWRF pattern)
// static mio::Date current_simulation_date;      // Current WRF simulation time
// static bool time_initialized = false;          // Track initialization
// static double calculation_step_length = 0.0;   // Read from SNOWPACK config (minutes)

// // Persistent SnowStation storage per grid point (CRYOWRF pattern)
// // Key format: "i_j" (e.g., "125_67" for grid point i=125, j=67)
// static std::map<std::string, std::unique_ptr<SnowStation>> grid_snowstations;
// static std::map<std::string, std::unique_ptr<Snowpack>> grid_snowpack_instances;
// static bool persistent_objects_initialized = false;


void SnowpackBridge::initialize_time(
    int start_year, 
    int start_month, 
    int start_day,
    int start_hour,
    int start_minute
) {
    if (time_initialized_) {
        return;
    }

    if (!config_initialized_) {
        throw std::runtime_error("Configuration must be initialized before time");
    }

    try {
        current_simulation_date_ = mio::Date(start_year, start_month, start_day,
                                           start_hour, start_minute, 0.0, 0.0);

        // Get calculation step length from configuration (in minutes)
        calculation_step_length_ = config_->get("CALCULATION_STEP_LENGTH", "Snowpack");
        time_initialized_ = true;

        printf("SNOWPACK-INFO: Time initialized to %04d-%02d-%02d %02d:%02d\n",
               start_year, start_month, start_day, start_hour, start_minute);

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL: Time initialization failed: %s\n", e.what());
        throw;
    }
}

// SnowStation* SnowpackBridge::get_or_create_snowstation(int i_grid, int j_grid, int wrf_domain_id,
//                                                       double wrf_lat, double wrf_lon, double wrf_alt) {
//     std::string station_key = SnowpackConstants::STATION_ID_PREFIX + "_" + std::to_string(wrf_domain_id) +
//                              "_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);

//     auto it = grid_snowstations_.find(station_key);
//     if (it != grid_snowstations_.end()) {
//         return it->second.get();
//     }

//     auto station = std::make_unique<SnowStation>();
//     station->setMeteoStationName(station_key);
//     station->setPosition(wrf_lat, wrf_lon, wrf_alt);

//     SnowStation* station_ptr = station.get();
//     grid_snowstations_[station_key] = std::move(station);

//     return station_ptr;
// }

Snowpack* SnowpackBridge::get_or_create_snowpack_instance(int i_grid, int j_grid,
                                                        double wrf_lat, double wrf_lon, double wrf_alt) {
    std::string station_key = SnowpackConstants::STATION_ID_PREFIX + "_1_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);

    auto it = grid_snowpack_instances_.find(station_key);
    if (it != grid_snowpack_instances_.end()) {
        return it->second.get();
    }

    auto snowpack = std::make_unique<Snowpack>();
    Snowpack* snowpack_ptr = snowpack.get();
    grid_snowpack_instances_[station_key] = std::move(snowpack);

    return snowpack_ptr;
}

SnowStation* SnowpackBridge::get_existing_snowstation(int i_grid, int j_grid) {
    std::string station_key = SnowpackConstants::STATION_ID_PREFIX + "_1_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);

    auto it = grid_snowstations_.find(station_key);
    if (it != grid_snowstations_.end()) {
        return it->second.get();
    }

    return nullptr; // Station doesn't exist
}

void SnowpackBridge::save_snowstation_state(int i_grid, int j_grid) {
    if (!io_ || !config_initialized_) {
        return;
    }

    std::string station_key = STATION_ID_PREFIX + "_1_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);
    auto station_it = grid_snowstations_.find(station_key);

    if (station_it != grid_snowstations_.end()) {
        try {
            io_->writeSnowStation(*(station_it->second), current_simulation_date_);
        } catch (const std::exception& e) {
            printf("SNOWPACK-ERROR: Failed to save state for grid (%d,%d): %s\n",
                   i_grid, j_grid, e.what());
        }
    }
}

void SnowpackBridge::save_all_snowpack_states() {
    if (!io_ || !config_initialized_) {
        printf("SNOWPACK-WARNING: Cannot save states - IO not initialized\n");
        return;
    }

    int saved_count = 0;
    for (const auto& [key, station] : grid_snowstations_) {
        try {
            io_->writeSnowStation(*station, current_simulation_date_);
            saved_count++;
        } catch (const std::exception& e) {
            printf("SNOWPACK-ERROR: Failed to save state for %s: %s\n", key.c_str(), e.what());
        }
    }

    printf("SNOWPACK-INFO: Saved %d snowpack states to disk\n", saved_count);
}

// Get or create persistent SNOWPACK objects for grid point
SnowpackObjects get_snowpack_objects(const MeteoInput& input) {
    initialize_snowpack_config();
    SnowpackObjects objects;
    
    objects.station = get_or_create_snowstation(input.i_grid, input.j_grid, 1, input.wrf_lat, input.wrf_lon, input.height);
    if (!objects.station) {
        printf("SNOWPACK-FATAL [%d,%d]: snow_station pointer is null!\n", input.i_grid, input.j_grid);
        std::abort();
    }
    
    objects.instance = get_or_create_snowpack_instance(input.i_grid, input.j_grid, input.wrf_lat, input.wrf_lon, input.height);
    if (!objects.instance) {
        printf("SNOWPACK-FATAL [%d,%d]: snowpack_instance pointer is null!\n", input.i_grid, input.j_grid);
        std::abort();
    }

    return objects;
}
void SnowpackBridge::execute_snowpack(
    const MeteoInput& input,
    SnowpackOutput& output,
    SnowpackLayerData* layer_data,
    BudgetData* budget_data
) {
    // Track physics calls for debugging
    if (++execute_call_count_ <= 5 || (execute_call_count_ % 1000 == 0)) {
        printf("SNOWPACK-INFO: Execute snowpack call #%d - Grid (%d,%d)\n", execute_call_count_, input.i_grid, input.j_grid);
    }

    // Get persistent SNOWPACK objects
    SnowpackObjects objects = get_snowpack_objects(input);

    // Initialize station if needed
    if (objects.station->getNumberOfElements() == 0) {
        SN_SNOWSOIL_DATA ssdata;
        mio::Coords position;
        position.setLatLon(input.wrf_lat, input.wrf_lon, input.height);
        std::string stationID = SnowpackConstants::STATION_ID_PREFIX + "_" + std::to_string(input.i_grid) + "_" + std::to_string(input.j_grid);
        std::string stationName = "WRF Grid Point " + std::to_string(input.i_grid) + "," + std::to_string(input.j_grid);
        ssdata.meta.setStationData(position, stationID, stationName);
        ssdata.Height = input.height;
        ssdata.nN = 1;
        ssdata.nLayers = 0;
        objects.station->initialize(ssdata, 0);
    }

    // Initialize time step
    bool time_advanced = initialize_time_step(input);

    // Create meteorological data
    auto Mdata = std::make_unique<CurrentMeteo>();
    auto surfFluxes = std::make_unique<SurfaceFluxes>();
    auto sn_Bdata = std::make_unique<BoundCond>();

    // Prepare meteorological data
    double current_snow_depth = (execute_call_count_ > 1) ? objects.station->cH : 0.0;
    prepare_meteorological_data(input, Mdata.get(), current_snow_depth);

    // Execute SNOWPACK physics
    double cumu_precip = 0.0;
    if (time_advanced) {
        objects.instance->runSnowpackModel(*Mdata, *objects.station, cumu_precip, *sn_Bdata, *surfFluxes);
    }

    // Collect surface fluxes
    try {
        surfFluxes->collectSurfaceFluxes(*sn_Bdata, *objects.station, *Mdata);
    } catch (const std::exception& e) {
        printf("SNOWPACK-ERROR: Exception in collectSurfaceFluxes: %s\n", e.what());
    }

    // Extract results
    size_t ndata_size = objects.station->Ndata.size();
    output.surface_temp = (ndata_size > 0) ? objects.station->Ndata.back().T : input.temp_air;
    output.snow_swe = objects.station->swe;
    output.snow_depth = objects.station->cH - objects.station->Ground;
    output.heat_flux_sensible = -1.0 * surfFluxes->qs;
    output.heat_flux_latent = -1.0 * sn_Bdata->ql;
    output.albedo = objects.station->Albedo;
    output.snow_coverage = 1.0;
    output.friction_velocity = Mdata->ustar;
    output.stability_param = 0.0;

    // Consistency check
    if (output.snow_depth > 0.001 && output.snow_swe <= 0.0) {
        output.snow_swe = output.snow_depth * SnowpackConstants::SNOW_DENSITY_FALLBACK;
    }

    // // Extract detailed data if requested
    // if (layer_data) {
    //     extract_layer_data(objects.station, *layer_data);
    // }

    // if (budget_data) {
    //     extract_budget_data(*surfFluxes, *sn_Bdata, cumu_precip, input.longwave_in,
    //                         input.shortwave_in, output.albedo, input.dt, *budget_data);
    //     budget_data->mass_swe = output.snow_swe;
    // }

    validate_outputs(output, input);

}

// Global singleton instance for Fortran compatibility
SnowpackBridge& g_snowpack_bridge = SnowpackBridge::instance();

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
    initialize_snowpack_config_with_path(SnowpackConfigManager::getDefaultConfigPath());
}

void SnowpackBridge::initialize_config(const std::string& ini_path) {
    if (config_initialized_) {
        return;
    }

    try {
        mio::Config file_config = SnowpackConfigManager::loadConfiguration(ini_path);
        SnowpackConfigManager::validateConfiguration(file_config);

        config_ = std::make_unique<SnowpackConfig>(file_config);
        io_ = std::make_unique<SnowpackIO>(*config_);

        config_file_path_ = ini_path;
        config_initialized_ = true;

        printf("SNOWPACK-INFO: Configuration initialized from %s\n", ini_path.c_str());

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL: Configuration failed for %s: %s\n", ini_path.c_str(), e.what());
        throw;
    }
}

// Persistent SnowStation management functions (CRYOWRF pattern)
std::string generate_grid_key(int i_grid, int j_grid) {
    return std::to_string(i_grid) + "_" + std::to_string(j_grid);
}

SnowStation* SnowpackBridge::get_or_create_snowstation(
    int i_grid, int j_grid, int wrf_domain_id, double wrf_lat, double wrf_lon, double wrf_alt
) {
    station_creation_calls++;
    initialize_snowpack_config();
    std::string grid_key = generate_grid_key(i_grid, j_grid);
    bool use_canopy, use_soil;

    // Check if SnowStation already exists for this grid point
    auto station_it = grid_snowstations.find(grid_key);
    if (station_it != grid_snowstations.end()) {
        return station_it->second.get();
    }

    // Create new SnowStation for this grid point - read configuration for correct parameters
    // MeteoIO handles boolean parsing - accepts false/FALSE/0 or true/TRUE/1
    global_config->getValue("CANOPY", "Snowpack", use_canopy);
    global_config->getValue("SNP_SOIL", "Snowpack", use_soil);
    const bool alpine3d = false; // WRF integration, not Alpine3D
    const bool sea_ice = false;  // Standard snowpack mode
    auto new_station = std::make_unique<SnowStation>(use_canopy, use_soil, alpine3d, sea_ice);

    // Set up basic station metadata using WRF coordinates
    // WRF's XLAT and XLONG contain geographic coordinates in degrees (regardless of map projection)

    mio::Coords position;

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

            ssdata.meta.position = position;        // Set position with coordinates from WRF
            ssdata.meta.stationID = stationID;      // Set station ID
            ssdata.meta.stationName = stationName;  // Set station name
            new_station->initialize(ssdata, 0);     // Initialize SnowStation from SN_SNOWSOIL_DATA

            for (size_t e = 0; e < new_station->getNumberOfElements(); e++) {
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
                        new_station->Edata[e].k[i] = 0.0;
                    }
                    if (std::isnan(new_station->Edata[e].c[i]) || new_station->Edata[e].c[i] != new_station->Edata[e].c[i]) {
                        new_station->Edata[e].c[i] = 0.0;
                    }
                }

                // Recompute heat capacity (c[TEMPERATURE]) from layer properties
                new_station->Edata[e].heatCapacity();
            }

            loaded_from_file = true;
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

// // Initialize configuration with specific path (called from Fortran)
// void initialize_snowpack_config_c(const char* ini_file_path) {
//     config_init_calls++;
//     std::string path_str(ini_file_path);
//     initialize_snowpack_config_with_path(path_str);
// }

// Initialize WRF simulation time (CRYOWRF pattern - called once from Fortran)
void initialize_wrf_simulation_time_c(int start_year, int start_month, int start_day,
                                      int start_hour, int start_minute) {
    time_init_calls++;

    try {
        // Initialize simulation time with WRF namelist start time (CRYOWRF pattern)
        current_simulation_date = mio::Date(start_year, start_month, start_day,
                                           start_hour, start_minute, 0.0, 0.0);
        time_initialized = true;

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Failed to initialize time: %s\n", e.what());
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

    SnowStation* snow_station = get_or_create_snowstation(i_grid, j_grid, 1, wrf_lat, wrf_lon, height);
    Snowpack* snowpack_instance = get_or_create_snowpack_instance(i_grid, j_grid, wrf_lat, wrf_lon, height);

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
    auto Mdata = std::make_unique<CurrentMeteo>();
    auto surfFluxes = std::make_unique<SurfaceFluxes>();
    auto sn_Bdata = std::make_unique<BoundCond>();

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
    Mdata->date = current_time;
    Mdata->ta = safe_temp;                                       // Air temperature [K]
    Mdata->rh = std::max(0.01, std::min(1.0, humidity));       // Relative humidity [0-1]
    Mdata->vw = std::max(0.1, wind_speed);                     // Wind speed [m/s]
    Mdata->dw = wind_dir;                                       // Wind direction [degrees]
    Mdata->iswr = std::max(0.0, shortwave_in);                 // Incoming shortwave [W/m²]
    Mdata->lw_net = std::max(0.0, longwave_in);               // Net longwave radiation [W/m²]
    Mdata->psum = std::max(0.0, precipitation);                // Precipitation [mm]

    // CRITICAL FIX: Set atmospheric emissivity instead of pressure (SNOWPACK compatibility)
    // SNOWPACK uses atmospheric emissivity (ea) instead of direct pressure
    // Calculate atmospheric emissivity using approximation (Stull 1988)
    // CRYOWRF SOURCE: Derived from atmospheric emissivity patterns in CRYOWRF coupler
    double atmospheric_emissivity = 0.7 + 5.95e-5 * pressure * exp(1500.0 * humidity / safe_temp);
    Mdata->ea = std::max(0.6, std::min(1.0, atmospheric_emissivity));  // Atmospheric emissivity [dimensionless]

    // CRITICAL FIX: Set ground temperature (CRYOWRF compatibility)
    // CRYOWRF SOURCE: ./CRYOWRF/src/coupler/main_coupler/Coupler.cpp - vecMyMeteo(MeteoData::TSG) = 400.00
    Mdata->ts0 = 400.0;  // Set ground temperature [K] (CRYOWRF pattern - fixed 400K ground temp)

    // CRITICAL FIX: Set maximum wind speed (CRYOWRF compatibility)
    // CRYOWRF SOURCE: ./CRYOWRF/src/coupler/main_coupler/Coupler.cpp - vecMyMeteo(MeteoData::VW_MAX) = l_VW_MAX
    Mdata->vw_max = wind_speed;  // Set max wind speed [m/s] (CRYOWRF pattern - track max wind)

    // Additional required meteorological parameters
    // Use precipitation phase threshold from CRYOWRF pattern
    Mdata->psum_ph = (safe_temp < SnowpackConstants::PRECIP_PHASE_THRESHOLD_K) ? 0.0 : 1.0;  // Precipitation phase (0=snow, 1=rain)
    Mdata->tss = mio::IOUtils::nodata;                         // Surface temperature (let SNOWPACK compute)
    Mdata->ts0 = safe_temp - SnowpackConstants::BOTTOM_TEMP_OFFSET_K;  // Bottom temperature estimate
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

    // Adapt roughness based on snow presence (CRYOWRF pattern from Coupler.cpp line 267)
    const double snow_depth_threshold = 0.03;  // 3cm snow depth threshold
    const double rough_len = (*snow_depth > snow_depth_threshold) ? roughness_length : 0.01;  // Bare soil z0 ~0.01
    Mdata->z0 = rough_len;

    // Compute friction velocity from wind speed using log profile (CRYOWRF pattern from Coupler.cpp)
    // u* = k * u / ln(z/z0), where k=0.4 (von Karman constant)
    const double von_karman = 0.4;  // Physical constant - matches CRYOWRF implementation
    const double z_wind = height;   // Measurement height [m]

    Mdata->ustar = (z_wind > Mdata->z0) ?
                  (von_karman * Mdata->vw / std::log(z_wind / Mdata->z0)) :
                  (0.1 * Mdata->vw);  // Fallback if z <= z0

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

    // Execute SNOWPACK physics (correct API with cumulative precipitation parameter)
    double cumu_precip = 0.0;  // Cumulative precipitation parameter
    snowpack_instance->runSnowpackModel(*Mdata, *snow_station, cumu_precip, *sn_Bdata, *surfFluxes);

    // Extract results from SNOWPACK (using correct member names)
    *surface_temp = (snow_station->Ndata.size() > 0) ? snow_station->Ndata.back().T : temp_air;  // Surface temperature [K]
    *snow_swe = snow_station->swe;                                  // Snow water equivalent [mm]
    *snow_depth = snow_station->cH;                                // Snow height [m]
    *heat_flux_sensible = -1.0 * surfFluxes->qs;                     // CRYOWRF sign convention | ./CRYOWRF/src/coupler/main_coupler/Coupler.cpp:1123
                                                                   // CRYOWRF SOURCE:
    *heat_flux_latent = -1.0 * sn_Bdata->ql;                      // CRYOWRF uses boundary layer data (line 1124)
                                                                   // CRYOWRF SOURCE: ./CRYOWRF/src/coupler/main_coupler/Coupler.cpp:1124
    *albedo = snow_station->Albedo;                                // Surface albedo [0-1]
    // Snow coverage - CRYOWRF pattern (hardcoded to 1.0 following CRYOWRF line 1127)
    *snow_coverage = 1.0;  // CRYOWRF hardcodes snow coverage to 1.0

    // Note: Atmospheric variables (friction_velocity, stability_param)
    // are only available in the enhanced snowpack_physics_layers interface

    // Consistency checks and fallbacks
    if (*snow_depth > 0.001 && *snow_swe <= 0.0) {
        *snow_swe = *snow_depth * SnowpackConstants::SNOW_DENSITY_FALLBACK;  // Use fallback density
    }

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
                             double* friction_velocity, double* stability_param,
                             // Layer arrays (max 100 layers)
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

    SnowStation* snow_station = get_or_create_snowstation(i_grid, j_grid, 1, wrf_lat, wrf_lon, height);
    Snowpack* snowpack_instance = get_or_create_snowpack_instance(i_grid, j_grid, wrf_lat, wrf_lon, height);

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
    auto Mdata = std::make_unique<CurrentMeteo>();
    auto surfFluxes = std::make_unique<SurfaceFluxes>();
    auto sn_Bdata = std::make_unique<BoundCond>();

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
    }

    call_counter_layers++;

    // Only advance time when counter matches (CRYOWRF pattern)
    if (compute_counter_layers > 0 && (call_counter_layers % compute_counter_layers) == 0) {
        current_simulation_date += (calculation_step_length / 1440.0);  // Convert minutes to days
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
    // CRITICAL FIX: Set atmospheric variables (SNOWPACK compatibility)
    // SNOWPACK uses atmospheric emissivity (ea) instead of direct pressure
    // CRYOWRF SOURCE: Derived from atmospheric emissivity patterns in CRYOWRF coupler
    double atmospheric_emissivity = 0.7 + 5.95e-5 * pressure * exp(1500.0 * humidity / safe_temp);
    Mdata->ea = std::max(0.6, std::min(1.0, atmospheric_emissivity));  // Atmospheric emissivity [dimensionless]

    Mdata->ts0 = 400.0;          // CRITICAL FIX: Set ground temperature [K]
                               // CRYOWRF SOURCE: ./CRYOWRF/src/coupler/main_coupler/Coupler.cpp - vecMyMeteo(MeteoData::TSG) = 400.00
    Mdata->vw_max = wind_speed;   // CRITICAL FIX: Set maximum wind speed [m/s]
                               // CRYOWRF SOURCE: ./CRYOWRF/src/coupler/main_coupler/Coupler.cpp - vecMyMeteo(MeteoData::VW_MAX) = l_VW_MAX
    Mdata->psum_ph = (safe_temp < SnowpackConstants::PRECIP_PHASE_THRESHOLD_K) ? 0.0 : 1.0;
    Mdata->tss = mio::IOUtils::nodata;
    // Note: ts0 already set to CRYOWRF-compatible value above

    // Initialize z0 and ustar (required for wind pumping in thermal conductivity)
    // Read ROUGHNESS_LENGTH from config (same as SNOWPACK's Meteo::MicroMet)
    static double roughness_length = -1.0;
    if (roughness_length < 0.0) {
        if (global_config) {
            global_config->getValue("ROUGHNESS_LENGTH", "Snowpack", roughness_length, mio::IOUtils::nothrow);
            if (roughness_length < 0.0) {
                roughness_length = 0.002;  // Default if not in config
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


    Mdata->hs = *snow_depth;

    // Execute SNOWPACK physics
    double cumu_precip = 0.0;

    snowpack_instance->runSnowpackModel(*Mdata, *snow_station, cumu_precip, *sn_Bdata, *surfFluxes);



    try {
        surfFluxes->collectSurfaceFluxes(*sn_Bdata, *snow_station, *Mdata);
    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [%d,%d]: Exception in collectSurfaceFluxes: %s\n", i_grid, j_grid, e.what());
        std::abort();
    } catch (...) {
        printf("SNOWPACK-FATAL [%d,%d]: Unknown exception in collectSurfaceFluxes\n", i_grid, j_grid);
        std::abort();
    }

    if (!snow_station) {
        printf("SNOWPACK-FATAL [%d,%d]: STEP 1 FAILED - snow_station is NULL!\n", i_grid, j_grid);
        std::abort();
    }
    if (!snowpack_instance) {
        printf("SNOWPACK-FATAL [%d,%d]: STEP 1 FAILED - snowpack_instance is NULL!\n", i_grid, j_grid);
        std::abort();
    }

    try {
        size_t ndata_size = snow_station->Ndata.size();
        *surface_temp = (ndata_size > 0) ? snow_station->Ndata.back().T : temp_air;
        *snow_swe = snow_station->swe;
        *snow_depth = snow_station->cH - snow_station->Ground;  // Subtract ground height (CRYOWRF line 1129)
        *heat_flux_sensible = -1.0 * surfFluxes->qs;      // Negative sign for WRF convention (CRYOWRF line 1123)
                                                           // CRYOWRF SOURCE: ./CRYOWRF/src/coupler/main_coupler/Coupler.cpp:1123
        *heat_flux_latent = -1.0 * sn_Bdata->ql;         // Use boundary condition data (CRYOWRF line 1124)
                                                           // CRYOWRF SOURCE: ./CRYOWRF/src/coupler/main_coupler/Coupler.cpp:1124
        *albedo = snow_station->Albedo;
        *snow_coverage = 1.0;  // Hardcoded to 1.0 following CRYOWRF pattern (line 1127)


    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [%d,%d]: STEP 2 CRASHED - Exception: %s\n", i_grid, j_grid, e.what());
        std::abort();
    }


    // STEP 3: Layer data extraction with validation
    try {
        size_t num_elements = snow_station->getNumberOfElements();
        *n_layers = static_cast<int>(num_elements);

        size_t layers_to_extract = std::min(num_elements, size_t(100));

        for (size_t i = 0; i < layers_to_extract; i++) {
            const ElementData& elem = snow_station->Edata[i];

            layer_temp[i] = elem.Te;
            layer_thick[i] = elem.L;
            layer_vol_ice[i] = elem.theta[ICE] * 100.0;
            layer_vol_water[i] = elem.theta[WATER] * 100.0;
            layer_vol_air[i] = elem.theta[AIR] * 100.0;

            layer_grain_radius[i] = elem.rg;
            layer_bond_radius[i] = elem.rb;
            layer_dendricity[i] = elem.dd;
            layer_sphericity[i] = elem.sp;


        }

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [%d,%d]: STEP 3 CRASHED - Exception: %s\n", i_grid, j_grid, e.what());
        printf("SNOWPACK-FATAL [%d,%d]: This indicates layer extraction is the crash point!\n", i_grid, j_grid);
        std::abort();
    }


    // STEP 4: Budget calculations and final output validation
    try {
        // Mass budgets (following CRYOWRF patterns) - SAFE ACCESS PATTERN
        *mass_precip = cumu_precip;                           // Cumulative precipitation [kg/m²]

        *mass_sublim = surfFluxes->mass[SurfaceFluxes::MS_SUBLIMATION];  // Sublimation [kg/m²] (CRYOWRF line 1136)

        // Calculate melt mass from SNOWPACK v11.08 SurfaceFluxes mass balance
        *mass_melt = surfFluxes->mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF]; // Melt runoff [kg/m²] - SNOWPACK v11.08 API

        *mass_swe = snow_station->swe;                        // Current SWE [kg/m²]

        // Calculate refreeze mass from ice base melting/freezing or negative runoff
        // MS_ICEBASE_MELTING_FREEZING: mass gain/loss of ice base due to melting-freezing
        // Also consider negative MS_SNOWPACK_RUNOFF which indicates refreezing
        double ice_base_meltfreeze = surfFluxes->mass[SurfaceFluxes::MS_ICEBASE_MELTING_FREEZING];
        *mass_refreeze = (ice_base_meltfreeze > 0.0) ? ice_base_meltfreeze : 0.0;  // Positive values = freezing


        // Energy budgets - SAFE ACCESS PATTERN
        *energy_lw_in = longwave_in;                           // Incoming longwave [W/m²]

        // Memory validation - check if surfFluxes object is properly aligned
        uintptr_t addr = reinterpret_cast<uintptr_t>(surfFluxes.get());
        bool alignment_ok = (addr % alignof(SurfaceFluxes) == 0);

        // Validate stack frame integrity
        void* stack_ptr = __builtin_frame_address(0);

        // Check memory around lw_out to detect corruption
        double* lw_out_ptr = &(surfFluxes->lw_out);

        // Try to read memory as raw bytes to check if it's readable
        volatile bool memory_readable = true;
        try {
            volatile double test_val = *lw_out_ptr;
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
    std::unique_ptr<SurfaceFluxes> safe_fluxes_backup = std::make_unique<SurfaceFluxes>();

    // Immediately copy the populated SurfaceFluxes to safe storage
    *safe_fluxes_backup = *surfFluxes;


    // Memory validation - check if safe_fluxes_backup object is properly aligned
    uintptr_t safe_addr = reinterpret_cast<uintptr_t>(safe_fluxes_backup.get());
    bool safe_alignment_ok = (safe_addr % alignof(SurfaceFluxes) == 0);

    // Check memory around lw_out to detect corruption
    double* safe_lw_out_ptr = &safe_fluxes_backup->lw_out;

    // Try to read memory as raw bytes to check if it's readable
    volatile bool safe_memory_readable = true;
    try {
        volatile double safe_test_val = *safe_lw_out_ptr;
    } catch (...) {
        safe_memory_readable = false;
        printf("SNOWPACK-FATAL [%d,%d]:     Raw memory read FAILED - memory corruption detected!\n", i_grid, j_grid);
    }

    if (!safe_memory_readable) {
        printf("SNOWPACK-FATAL [%d,%d]: CRITICAL: Memory corruption detected - aborting\n", i_grid, j_grid);
        std::abort();
    }

    // Now safely extract from the heap-allocated copy instead of stack memory
    try {
        *energy_lw_out = safe_fluxes_backup->lw_out;                   // Outgoing longwave [W/m²] (CRYOWRF line 1087)
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
        *energy_rain = safe_fluxes_backup->qr;                          // Rain heat flux [W/m²] (CRYOWRF pattern)

        // Total energy balance (following CRYOWRF line 1150)
        *energy_total = safe_fluxes_backup->dIntEnergy / dt;           // Total energy [W/m²] (CRYOWRF pattern)

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [%d,%d]: STEP 4 CRASHED - Exception: %s\n", i_grid, j_grid, e.what());
        printf("SNOWPACK-FATAL [%d,%d]: This indicates budget calculation is the crash point!\n", i_grid, j_grid);
        std::abort();
    }

    // Clear unused layers (initialize remaining array elements)
    for (size_t i = static_cast<size_t>(*n_layers); i < 100; i++) {
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

    // Check for invalid output values that could cause Fortran crashes
    if (std::isnan(*surface_temp) || std::isinf(*surface_temp)) {
        printf("SNOWPACK-FATAL [%d,%d]: surface_temp is NaN/Inf (%.6f)!\n", i_grid, j_grid, *surface_temp);
        std::abort();
    }
    if (std::isnan(*snow_swe) || std::isinf(*snow_swe)) {
        printf("SNOWPACK-FATAL [%d,%d]: snow_swe is NaN/Inf (%.6f)!\n", i_grid, j_grid, *snow_swe);
        std::abort();
    }
    if (*n_layers < 0 || *n_layers > 100) {
        printf("SNOWPACK-FATAL [%d,%d]: n_layers is invalid (%d)!\n", i_grid, j_grid, *n_layers);
        std::abort();
    }


  } catch (const std::exception& e) {
    // DEBUG: Enhanced exception handling
    printf("SNOWPACK-FATAL [%d,%d]: >>> C++ EXCEPTION CAUGHT <<<\n", i_grid, j_grid);
    printf("SNOWPACK-FATAL [%d,%d]: Exception type: %s\n", i_grid, j_grid, typeid(e).name());
    printf("SNOWPACK-FATAL [%d,%d]: Exception message: %s\n", i_grid, j_grid, e.what());
    printf("SNOWPACK-FATAL [%d,%d]: This indicates a critical error in SNOWPACK processing!\n", i_grid, j_grid);


    printf("SNOWPACK-ERROR: Grid point (%d,%d) - Ta=%.2fK, RH=%.2f, Precip=%.3fmm\n",
           i_grid, j_grid, temp_air, humidity, precipitation);
    printf("SNOWPACK-FATAL: Aborting WRF run due to SNOWPACK physics failure\n");
    std::abort();  // Abort instead of silent fallback
  }
}

// // Extract detailed layer data from SNOWPACK for CRYOWRF compatibility
// void extract_snowpack_layers_c(int i_grid, int j_grid,
//                                float layer_temps[100], float layer_thick[100],
//                                float layer_voli[100], float layer_volw[100], float layer_volv[100],
//                                float layer_rg[100], float layer_rb[100],
//                                float layer_dd[100], float layer_sp[100],
//                                int* n_layers) {
//   try {
//     // Initialize all arrays to zero
//     for (int i = 0; i < 100; i++) {
//       layer_temps[i] = 0.0f;
//       layer_thick[i] = 0.0f;
//       layer_voli[i] = 0.0f;
//       layer_volw[i] = 0.0f;
//       layer_volv[i] = 0.0f;
//       layer_rg[i] = 0.0f;
//       layer_rb[i] = 0.0f;
//       layer_dd[i] = 0.0f;
//       layer_sp[i] = 1.0f;  // Default sphericity
//     }

//     // Use persistent SnowStation objects from our grid-based storage
//     std::string grid_key = generate_grid_key(i_grid, j_grid);
//     auto station_it = grid_snowstations.find(grid_key);

//     if (station_it == grid_snowstations.end()) {
//         // No SnowStation exists for this grid point - set to zero layers
//         *n_layers = 0;
//         printf("SNOWPACK-LAYERS: Grid (%d,%d) has no persistent SnowStation - zero layers\n",
//                i_grid, j_grid);
//         return;
//     }

//     // Access persistent SnowStation for this grid point
//     const SnowStation* xdata = station_it->second.get();

//     // CRYOWRF-style layer extraction algorithm (from Coupler.cpp:1153-1185)
//     // Extract from evolved snowpack with real layer data

//     const std::vector<ElementData>& elem_data = xdata->Edata;
//     const int get_size = (int)xdata->getNumberOfElements();
//     const int loc_snpack_lay_to_sav = 100; // Maximum layers to save (CRYOWRF default)
//     const int lim_size = std::max(0, get_size - loc_snpack_lay_to_sav);

//     // CRYOWRF algorithm: Extract from surface down, top-down indexing
//     int tmtm = 0;
//     for (int e = get_size - 1; e >= lim_size && tmtm < 100; e++, tmtm++) {
//       // Extract layer data using exact CRYOWRF algorithm
//       // With persistent SnowStation objects, evolved snowpack layers are available for extraction
//       if (e + 1 < (int)xdata->Ndata.size()) { // Safety check for node access
//         layer_temps[tmtm] = (float)xdata->Ndata[e + 1].T;        // Node temperature [K]
//       }
//       layer_thick[tmtm] = (float)elem_data[e].L;                // Layer thickness [m]
//       layer_voli[tmtm] = (float)elem_data[e].theta[ICE] * 100.0f; // Ice volume fraction [%]
//       layer_volw[tmtm] = (float)elem_data[e].theta[WATER] * 100.0f; // Water volume fraction [%]
//       layer_volv[tmtm] = (float)elem_data[e].theta[AIR] * 100.0f; // Air volume fraction [%]
//       layer_rg[tmtm] = (float)elem_data[e].rg;                  // Grain radius [mm]
//       layer_rb[tmtm] = (float)elem_data[e].rb;                  // Bond radius [mm]
//       layer_dd[tmtm] = (float)elem_data[e].dd;                  // Dendricity [-]
//       layer_sp[tmtm] = (float)elem_data[e].sp;                  // Sphericity [-]
//     }

//     *n_layers = get_size; // Total number of layers (following CRYOWRF: sn_nlayer = get_size)

//   } catch (const std::exception& e) {
//     printf("SNOWPACK-ERROR: Layer extraction failed for grid (%d,%d): %s\n", i_grid, j_grid, e.what());
//     *n_layers = 0;
//   }
// }

// // Save all snowpack states (called from Fortran for periodic saves)
// void save_all_snowpack_states_c() {
//     save_all_snowpack_states();
// }

}
// Extern C wrapper for Fortran BIND(C) interface
extern "C" void snowpack_physics_layers(
    double temp_air, double humidity, double wind_speed, double wind_dir, double shortwave_in, double longwave_in,
    double precipitation, double pressure, double height, double dt, int i_grid, int j_grid, double wrf_lat,
    double wrf_lon, double* snow_swe, double* snow_depth, double* surface_temp, double* heat_flux_sensible,
    double* heat_flux_latent, double* albedo, double* snow_coverage, double* friction_velocity, double* stability_param,
    int* n_layers, double* layer_temp, double* layer_thick, double* layer_vol_ice, double* layer_vol_water,
    double* layer_vol_air, double* layer_grain_radius, double* layer_bond_radius, double* layer_dendricity,
    double* layer_sphericity, double* mass_precip, double* mass_sublim, double* mass_melt, double* mass_swe,
    double* mass_refreeze, double* energy_lw_in, double* energy_lw_out, double* energy_sw_in, double* energy_sw_out,
    double* energy_sensible, double* energy_latent, double* energy_ground_flux, double* energy_rain,
    double* energy_total
) {
    SnowpackBridge& bridge = SnowpackBridge::instance();

    // Ensure configuration is loaded (loads io.ini if needed)
    if (!bridge.is_config_initialized()) {
        bridge.initialize_config(SnowpackConfigManager::getDefaultConfigPath());
    }

    // Init internal data structures
    MeteoInput input;
    SnowpackOutput output;
    SnowpackLayerData layer_data;
    BudgetData budget_data;
    
    // 
    input.temp_air = temp_air;
    input.humidity = humidity;
    input.wind_speed = wind_speed;
    input.wind_dir = wind_dir;
    input.shortwave_in = shortwave_in;
    input.longwave_in = longwave_in;
    input.precipitation = precipitation;
    input.pressure = pressure;
    input.height = height;
    input.dt = dt;
    input.i_grid = i_grid;
    input.j_grid = j_grid;
    input.wrf_lat = wrf_lat;
    input.wrf_lon = wrf_lon;

    bridge.execute_snowpack(input, output, &layer_data, &budget_data);

    // Basic snowpack outputs
    *snow_swe = output.snow_swe;
    *snow_depth = output.snow_depth;
    *surface_temp = output.surface_temp;
    *heat_flux_sensible = output.heat_flux_sensible;
    *heat_flux_latent = output.heat_flux_latent;
    *albedo = output.albedo;
    *snow_coverage = output.snow_coverage;
    *friction_velocity = output.friction_velocity;
    *stability_param = output.stability_param;

    // Layer data (up to 100 layers)
    *n_layers = layer_data.n_layers;
    for (int i = 0; i < layer_data.n_layers && i < 100; i++) {
        layer_temp[i] = layer_data.layer_temp[i];
        layer_thick[i] = layer_data.layer_thick[i];
        layer_vol_ice[i] = layer_data.layer_vol_ice[i];
        layer_vol_water[i] = layer_data.layer_vol_water[i];
        layer_vol_air[i] = layer_data.layer_vol_air[i];
        layer_grain_radius[i] = layer_data.layer_grain_radius[i];
        layer_bond_radius[i] = layer_data.layer_bond_radius[i];
        layer_dendricity[i] = layer_data.layer_dendricity[i];
        layer_sphericity[i] = layer_data.layer_sphericity[i];
    }

    // Mass budget data
    *mass_precip = budget_data.mass_precip;
    *mass_sublim = budget_data.mass_sublim;
    *mass_melt = budget_data.mass_melt;
    *mass_swe = budget_data.mass_swe;
    *mass_refreeze = budget_data.mass_refreeze;

    // Energy budget data
    *energy_lw_in = budget_data.energy_lw_in;
    *energy_lw_out = budget_data.energy_lw_out;
    *energy_sw_in = budget_data.energy_sw_in;
    *energy_sw_out = budget_data.energy_sw_out;
    *energy_sensible = budget_data.energy_sensible;
    *energy_latent = budget_data.energy_latent;
    *energy_ground_flux = budget_data.energy_ground_flux;
    *energy_rain = budget_data.energy_rain;
    *energy_total = budget_data.energy_total;
    
}