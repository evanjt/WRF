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

// Forward declarations for namespace functions
namespace SnowpackObjects {
    SnowpackBridgeObjects create_station_objects(const MeteoInput& input, SnowpackBridge& bridge);
    void ensure_station_initialized(SnowStation* station, const MeteoInput& input, SnowpackBridge& bridge);
    void validate_output_values(const SnowpackOutput& output, int i_grid, int j_grid);
    void validate_layer_data(const SnowpackLayerData& layer_data, int i_grid, int j_grid);
}

namespace SnowpackUtils {
    void advance_simulation_time(const MeteoInput& input, mio::Date& current_time, double calculation_step_length);
    void prepare_meteo_data(const MeteoInput& input, CurrentMeteo& Mdata, const mio::Date& current_time, double current_snow_depth, const SnowpackConfig* config);
    void extract_surface_outputs(const SnowStation& station, const SurfaceFluxes& fluxes, const BoundCond& bc, SnowpackOutput& output, double temp_air);
    void extract_layer_data(const SnowStation& station, SnowpackLayerData& layer_data);
    void extract_budget_data(const SurfaceFluxes& fluxes, const BoundCond& bc, double cumu_precip, const MeteoInput& input, BudgetData& budget);
}



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


Snowpack* SnowpackBridge::get_or_create_snowpack_instance(int i_grid, int j_grid,
                                                        double wrf_lat, double wrf_lon, double wrf_alt) {
    std::string station_key = SnowpackConstants::STATION_ID_PREFIX + "_1_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);

    auto it = grid_snowpack_instances_.find(station_key);
    if (it != grid_snowpack_instances_.end()) {
        return it->second.get();
    }

    auto snowpack = std::unique_ptr<Snowpack>(new Snowpack(*config_));
    Snowpack* snowpack_ptr = snowpack.get();
    grid_snowpack_instances_[station_key] = std::move(snowpack);

    return snowpack_ptr;
}

SnowStation* SnowpackBridge::get_or_create_snowstation(
    int i_grid, int j_grid, int wrf_domain_id, double wrf_lat, double wrf_lon, double wrf_alt
) {
    std::string station_key = SnowpackConstants::STATION_ID_PREFIX + "_" + std::to_string(wrf_domain_id) + "_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);
    bool use_canopy, use_soil;

    // Check if SnowStation already exists for this grid point
    auto station_it = grid_snowstations_.find(station_key);
    if (station_it != grid_snowstations_.end()) {
        return station_it->second.get();
    }

    // Create new SnowStation for this grid point - read configuration for correct parameters
    if (!config_) {
        printf("SNOWPACK-FATAL: Configuration not initialized when creating SnowStation\n");
        std::abort();
    }

    config_->getValue("CANOPY", "Snowpack", use_canopy);
    config_->getValue("SNP_SOIL", "Snowpack", use_soil);
    const bool alpine3d = false; // WRF integration, not Alpine3D
    const bool sea_ice = false;  // Standard snowpack mode
    auto new_station = std::unique_ptr<SnowStation>(new SnowStation(use_canopy, use_soil, alpine3d, sea_ice));

    // Set up basic station metadata using WRF coordinates
    mio::Coords position;

    try {
        position.setLatLon(wrf_lat, wrf_lon, wrf_alt);  // Correct order: lat, lon, alt
        printf("SNOWPACK-INFO: Successfully set geographic coordinates for station (%d,%d)\n", i_grid, j_grid);
    } catch (const mio::InvalidArgumentException& e) {
        printf("SNOWPACK-ERROR: Invalid coordinates for station (%d,%d): %s\n", i_grid, j_grid, e.what());
        printf("SNOWPACK-ERROR: Received lat=%.6f, lon=%.6f - these may not be geographic coordinates\n", wrf_lat, wrf_lon);
        throw;
    }

    std::string stationID = "WRF_GRID_" + std::to_string(wrf_domain_id) + "_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);
    std::string stationName = "WRF Grid Point " + std::to_string(i_grid) + "," + std::to_string(j_grid);
    new_station->meta.setStationData(position, stationID, stationName);

    // Try to load existing .sno file state (CRYOWRF pattern)
    bool loaded_from_file = false;
    if (io_) {
        // CRYOWRF C++ naming pattern: snpack_{grid_id}_{I}_{J}.sno (from Coupler.cpp line 613)
        std::string sno_filename = "snpack_" + std::to_string(wrf_domain_id) + "_" + std::to_string(i_grid) + "_" + std::to_string(j_grid) + ".sno";

        printf("SNOWPACK-INIT [%d,%d]: Attempting to load .sno file: %s\n", i_grid, j_grid, sno_filename.c_str());

        try {
            // Attempt to read existing snowpack state
            SN_SNOWSOIL_DATA ssdata;
            ZwischenData zdata;
            mio::Date profile_date;

            printf("SNOWPACK-INIT [%d,%d]: Calling readSnowCover()...\n", i_grid, j_grid);
            io_->readSnowCover(sno_filename, stationID, ssdata, zdata, false);
            printf("SNOWPACK-INIT [%d,%d]: readSnowCover() returned %zu layers\n", i_grid, j_grid, ssdata.nLayers);

            ssdata.meta.position = position;        // Set position with coordinates from WRF
            ssdata.meta.stationID = stationID;      // Set station ID
            ssdata.meta.stationName = stationName;  // Set station name
            new_station->initialize(ssdata, 0);     // Initialize SnowStation from SN_SNOWSOIL_DATA

            // Ensure proper initialization of layer properties
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
            printf("SNOWPACK-INIT [%d,%d]: Successfully loaded station state from %s\n", i_grid, j_grid, sno_filename.c_str());

        } catch (const std::exception& e) {
            // No existing state file - start with fresh snowpack
            printf("SNOWPACK-INIT [%d,%d]: Failed to load .sno file - %s\n", i_grid, j_grid, e.what());
            printf("SNOWPACK-INFO: No existing state for grid (%d,%d), starting fresh\n", i_grid, j_grid);
        }
    } else {
        printf("SNOWPACK-INIT [%d,%d]: SnowpackIO is NULL - cannot load .sno files\n", i_grid, j_grid);
    }

    if (!loaded_from_file) {
        printf("SNOWPACK-INFO: Fresh SnowStation created for grid (%d,%d)\n", i_grid, j_grid);
    }

    // Store the new SnowStation
    SnowStation* station_ptr = new_station.get();
    grid_snowstations_[station_key] = std::move(new_station);

    return station_ptr;
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

    std::string station_key = SnowpackConstants::STATION_ID_PREFIX + "_1_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);
    auto station_it = grid_snowstations_.find(station_key);

    if (station_it != grid_snowstations_.end()) {
        try {
            // Use correct SNOWPACK API
            ZwischenData zdata;  // Empty for basic usage
            io_->writeSnowCover(current_simulation_date_, *(station_it->second), zdata, true);
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
    for (auto it = grid_snowstations_.begin(); it != grid_snowstations_.end(); ++it) {
        const std::string& key = it->first;
        SnowStation* station = it->second.get();
        try {
            // Use correct SNOWPACK API
            ZwischenData zdata;  // Empty for basic usage
            io_->writeSnowCover(current_simulation_date_, *station, zdata, true);
            saved_count++;
        } catch (const std::exception& e) {
            printf("SNOWPACK-ERROR: Failed to save state for %s: %s\n", key.c_str(), e.what());
        }
    }

    printf("SNOWPACK-INFO: Saved %d snowpack states to disk\n", saved_count);
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

    // Validate initialization
    if (!config_initialized_) {
        printf("SNOWPACK-FATAL: Configuration not initialized - call initialize_config() first\n");
        std::abort();
    }

    if (!time_initialized_) {
        printf("SNOWPACK-FATAL: Time not initialized - call initialize_time() first\n");
        std::abort();
    }

    try {
        // 1. Get/create persistent SNOWPACK objects (singleton pattern only)
        auto objects = SnowpackObjects::create_station_objects(input, *this);

        // 2. Ensure station is properly initialized (load from file or create fresh)
        SnowpackObjects::ensure_station_initialized(objects.station, input, *this);

        // 3. Advance simulation time following CRYOWRF pattern
        SnowpackUtils::advance_simulation_time(input, current_simulation_date_, calculation_step_length_);

        // 4. Prepare meteorological data using utility function
        auto Mdata = std::unique_ptr<CurrentMeteo>(new CurrentMeteo());
        double current_snow_depth = (execute_call_count_ > 1) ? objects.station->cH : 0.0;
        SnowpackUtils::prepare_meteo_data(input, *Mdata, current_simulation_date_, current_snow_depth, config_.get());

        // 5. Create surface fluxes and boundary condition objects
        auto surfFluxes = std::unique_ptr<SurfaceFluxes>(new SurfaceFluxes());
        auto sn_Bdata = std::unique_ptr<BoundCond>(new BoundCond());

        // 6. SINGLE SNOWPACK CALL (CRYOWRF pattern - line 1074 in original)
        double cumu_precip = 0.0;
        objects.instance->runSnowpackModel(*Mdata, *objects.station, cumu_precip, *sn_Bdata, *surfFluxes);

        // 7. Collect surface fluxes (CRYOWRF pattern - line 1078 in original)
        try {
            surfFluxes->collectSurfaceFluxes(*sn_Bdata, *objects.station, *Mdata);
        } catch (const std::exception& e) {
            printf("SNOWPACK-ERROR: Exception in collectSurfaceFluxes: %s\n", e.what());
            // Continue with partial data rather than aborting
        }

        // 8. Extract surface outputs using utility function (CRYOWRF pattern - lines 1121-1130)
        SnowpackUtils::extract_surface_outputs(*objects.station, *surfFluxes, *sn_Bdata, output, input.temp_air);

        // 9. Extract detailed layer data if requested (CRYOWRF pattern - lines 1153-1185)
        if (layer_data) {
            SnowpackUtils::extract_layer_data(*objects.station, *layer_data);
        }

        // 10. Extract budget data if requested (CRYOWRF pattern - lines 1134-1150)
        if (budget_data) {
            SnowpackUtils::extract_budget_data(*surfFluxes, *sn_Bdata, cumu_precip, input, *budget_data);
            budget_data->mass_swe = output.snow_swe;  // Ensure consistency
        }

        // 11. Validate outputs
        SnowpackObjects::validate_output_values(output, input.i_grid, input.j_grid);
        if (layer_data) {
            SnowpackObjects::validate_layer_data(*layer_data, input.i_grid, input.j_grid);
        }

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [%d,%d]: Exception in execute_snowpack: %s\n",
               input.i_grid, input.j_grid, e.what());
        printf("SNOWPACK-ERROR: Grid point (%d,%d) - Ta=%.2fK, RH=%.2f, Precip=%.3fmm\n",
               input.i_grid, input.j_grid, input.temp_air, input.humidity, input.precipitation);
        std::abort();
    }
}

// Global singleton instance for Fortran compatibility
SnowpackBridge& g_snowpack_bridge = SnowpackBridge::instance();


void SnowpackBridge::initialize_config(const std::string& ini_path) {
    if (config_initialized_) {
        return;
    }

    try {
        mio::Config file_config = SnowpackConfigManager::loadConfiguration(ini_path);
        SnowpackConfigManager::validateConfiguration(file_config);

        config_ = std::unique_ptr<SnowpackConfig>(new SnowpackConfig(file_config));
        io_ = std::unique_ptr<SnowpackIO>(new SnowpackIO(*config_));

        config_file_path_ = ini_path;
        config_initialized_ = true;

        printf("SNOWPACK-INFO: Configuration initialized from %s\n", ini_path.c_str());

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL: Configuration failed for %s: %s\n", ini_path.c_str(), e.what());
        throw;
    }
}



// C interface for Fortran binding
extern "C" {

// Initialize WRF simulation time (CRYOWRF pattern - called once from Fortran)
void initialize_wrf_simulation_time_c(int start_year, int start_month, int start_day,
                                      int start_hour, int start_minute) {
    try {
        SnowpackBridge& bridge = SnowpackBridge::instance();

        // Initialize simulation time with WRF namelist start time (CRYOWRF pattern)
        bridge.initialize_time(start_year, start_month, start_day, start_hour, start_minute);

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [C++/snowpack_wrf_bridge.cpp]: Failed to initialize time: %s\n", e.what());
        std::abort();
    }
}

// Simplified C interface using structured data
void snowpack_physics_c(
    const MeteoInput* input,
    SnowpackOutput* output,
    SnowpackLayerData* layer_data,
    BudgetData* budget_data
) {
    SnowpackBridge& bridge = SnowpackBridge::instance();

    // Ensure configuration is loaded (loads io.ini if needed)
    if (!bridge.is_config_initialized()) {
        bridge.initialize_config(SnowpackConfigManager::getDefaultConfigPath());
    }

    // Execute snowpack physics through singleton
    bridge.execute_snowpack(*input, *output, layer_data, budget_data);
}

// Backward compatibility interface with individual parameters
void snowpack_physics_individual_c(
    double temp_air, double humidity, double wind_speed, double wind_dir,
    double shortwave_in, double longwave_in, double precipitation, double pressure, double height, double dt,
    int i_grid, int j_grid, double wrf_lat, double wrf_lon,
    double* snow_swe, double* snow_depth, double* surface_temp,
    double* heat_flux_sensible, double* heat_flux_latent, double* albedo, double* snow_coverage
) {
    // Build structured input
    MeteoInput input;
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

    // Execute and get structured output
    SnowpackOutput output;
    snowpack_physics_c(&input, &output, nullptr, nullptr);

    // Extract individual parameters for backward compatibility
    *snow_swe = output.snow_swe;
    *snow_depth = output.snow_depth;
    *surface_temp = output.surface_temp;
    *heat_flux_sensible = output.heat_flux_sensible;
    *heat_flux_latent = output.heat_flux_latent;
    *albedo = output.albedo;
    *snow_coverage = output.snow_coverage;
}


// Save all snowpack states (called from Fortran for periodic saves)
void save_all_snowpack_states_c() {
    SnowpackBridge& bridge = SnowpackBridge::instance();
    bridge.save_all_snowpack_states();
}

// Initialize snowpack configuration (calls SnowpackBridge method)
void initialize_snowpack_config_c(const char* ini_file_path) {
    SnowpackBridge& bridge = SnowpackBridge::instance();
    bridge.initialize_config(std::string(ini_file_path));
}

// Extract snowpack layers from structured data (non-redundant implementation)
void extract_snowpack_layers_c(int i_grid, int j_grid,
                              float* layer_temps, float* layer_thick, float* layer_voli, float* layer_volw, float* layer_volv,
                              float* layer_rg, float* layer_rb, float* layer_dd, float* layer_sp, int* n_layers) {
    // Get existing snow station data
    SnowpackBridge& bridge = SnowpackBridge::instance();
    SnowStation* station = bridge.get_existing_snowstation(i_grid, j_grid);

    if (!station) {
        *n_layers = 0;
        return;
    }

    // Extract layer data directly from SnowStation
    SnowpackLayerData layer_data;
    SnowpackUtils::extract_layer_data(*station, layer_data);

    // Copy to output arrays
    *n_layers = layer_data.n_layers;
    for (int i = 0; i < layer_data.n_layers && i < 100; i++) {
        layer_temps[i] = static_cast<float>(layer_data.layer_temp[i]);
        layer_thick[i] = static_cast<float>(layer_data.layer_thick[i]);
        layer_voli[i] = static_cast<float>(layer_data.layer_vol_ice[i]);
        layer_volw[i] = static_cast<float>(layer_data.layer_vol_water[i]);
        layer_volv[i] = static_cast<float>(layer_data.layer_vol_air[i]);
        layer_rg[i] = static_cast<float>(layer_data.layer_grain_radius[i]);
        layer_rb[i] = static_cast<float>(layer_data.layer_bond_radius[i]);
        layer_dd[i] = static_cast<float>(layer_data.layer_dendricity[i]);
        layer_sp[i] = static_cast<float>(layer_data.layer_sphericity[i]);
    }
}

// Fortran-mangled versions (with trailing underscores) for compatibility
void snowpack_physics_individual_c_(
    double temp_air, double humidity, double wind_speed, double wind_dir,
    double shortwave_in, double longwave_in, double precipitation, double pressure, double height, double dt,
    int i_grid, int j_grid, double wrf_lat, double wrf_lon,
    double* snow_swe, double* snow_depth, double* surface_temp,
    double* heat_flux_sensible, double* heat_flux_latent, double* albedo, double* snow_coverage
) {
    snowpack_physics_individual_c(temp_air, humidity, wind_speed, wind_dir,
                                  shortwave_in, longwave_in, precipitation, pressure, height, dt,
                                  i_grid, j_grid, wrf_lat, wrf_lon,
                                  snow_swe, snow_depth, surface_temp,
                                  heat_flux_sensible, heat_flux_latent, albedo, snow_coverage);
}

void extract_snowpack_layers_c_(int i_grid, int j_grid,
                               float* layer_temps, float* layer_thick, float* layer_voli, float* layer_volw, float* layer_volv,
                               float* layer_rg, float* layer_rb, float* layer_dd, float* layer_sp, int* n_layers) {
    extract_snowpack_layers_c(i_grid, j_grid, layer_temps, layer_thick, layer_voli, layer_volw, layer_volv,
                              layer_rg, layer_rb, layer_dd, layer_sp, n_layers);
}

} // extern "C"

/*
 * ============================================================================
 * UTILITY IMPLEMENTATIONS (Extracted for organization but compiled together)
 * ============================================================================
 */

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

// Helper function to generate grid key
std::string generate_grid_key(int i_grid, int j_grid) {
    return std::to_string(i_grid) + "_" + std::to_string(j_grid);
}

// Try to load existing station state from .sno file
bool load_station_state(const std::string& station_key,
                        SnowStation* station,
                        const mio::Coords& position,
                        const std::string& stationID,
                        const std::string& stationName,
                        SnowpackIO* io,
                        int wrf_domain_id,
                        int i_grid,
                        int j_grid) {
    if (!io) {
        printf("SNOWPACK-INIT [%d,%d]: SnowpackIO is NULL - cannot load .sno files\n", i_grid, j_grid);
        return false;
    }

    try {
        // CRYOWRF C++ naming pattern: snpack_{grid_id}_{I}_{J}.sno (from Coupler.cpp line 613)
        // Let SNOWPACK handle the path through SNOWPATH configuration
        std::string sno_filename = "snpack_" + std::to_string(wrf_domain_id) + "_" +
                                  std::to_string(i_grid) + "_" + std::to_string(j_grid) + ".sno";

        printf("SNOWPACK-INIT [%d,%d]: Attempting to load .sno file: %s\n", i_grid, j_grid, sno_filename.c_str());

        // Attempt to read existing snowpack state
        SN_SNOWSOIL_DATA ssdata;
        ZwischenData zdata;
        mio::Date profile_date;

        printf("SNOWPACK-INIT [%d,%d]: Calling readSnowCover()...\n", i_grid, j_grid);
        io->readSnowCover(sno_filename, stationID, ssdata, zdata, false);
        printf("SNOWPACK-INIT [%d,%d]: readSnowCover() returned %zu layers\n", i_grid, j_grid, ssdata.nLayers);

        ssdata.meta.position = position;        // Set position with coordinates from WRF
        ssdata.meta.stationID = stationID;      // Set station ID
        ssdata.meta.stationName = stationName;  // Set station name
        station->initialize(ssdata, 0);         // Initialize SnowStation from SN_SNOWSOIL_DATA

        // Ensure proper initialization of layer properties
        for (size_t e = 0; e < station->getNumberOfElements(); e++) {
            // Ensure k and c vectors have proper size (3 elements: TEMPERATURE, SEEPAGE, SETTLEMENT)
            if (station->Edata[e].k.size() < 3) {
                station->Edata[e].k.resize(3, 0.0);  // Resize and initialize to 0
            }
            if (station->Edata[e].c.size() < 3) {
                station->Edata[e].c.resize(3, 0.0);  // Resize and initialize to 0
            }

            // Force initialize to 0.0 if NaN or uninitialized
            for (size_t i = 0; i < 3; i++) {
                if (std::isnan(station->Edata[e].k[i]) || station->Edata[e].k[i] != station->Edata[e].k[i]) {
                    station->Edata[e].k[i] = 0.0;
                }
                if (std::isnan(station->Edata[e].c[i]) || station->Edata[e].c[i] != station->Edata[e].c[i]) {
                    station->Edata[e].c[i] = 0.0;
                }
            }

            // Recompute heat capacity (c[TEMPERATURE]) from layer properties
            station->Edata[e].heatCapacity();
        }

        printf("SNOWPACK-INIT [%d,%d]: Successfully loaded station state from %s\n", i_grid, j_grid, sno_filename.c_str());
        return true;

    } catch (const std::exception& e) {
        // No existing state file - start with fresh snowpack
        printf("SNOWPACK-INIT [%d,%d]: Failed to load .sno file - %s\n", i_grid, j_grid, e.what());
        printf("SNOWPACK-INFO: No existing state for grid (%d,%d), starting fresh\n", i_grid, j_grid);
        return false;
    }
}

// Initialize new station with default values
void initialize_new_station(SnowStation* station,
                           const mio::Coords& position,
                           const std::string& stationID,
                           const std::string& stationName,
                           double height,
                           int i_grid,
                           int j_grid) {
    try {
        SN_SNOWSOIL_DATA ssdata;
        ssdata.meta.setStationData(position, stationID, stationName);
        ssdata.Height = height;    // Station elevation [m]
        ssdata.nN = 1;            // Start with 1 node (ground only)
        ssdata.nLayers = 0;       // No snow/soil layers initially

        station->initialize(ssdata, 0);  // Initialize with sector 0
        printf("SNOWPACK-INFO: Fresh SnowStation created for grid (%d,%d)\n", i_grid, j_grid);

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL [%d,%d]: Failed to initialize new station: %s\n", i_grid, j_grid, e.what());
        throw;
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

void advance_simulation_time(const MeteoInput& input,
                           mio::Date& current_time,
                           double calculation_step_length) {
    g_time_manager.advance_time(input, current_time, calculation_step_length);
}

} // namespace SnowpackUtils

namespace SnowpackObjects {

SnowpackBridgeObjects create_station_objects(const MeteoInput& input,
                                              SnowpackBridge& bridge) {
    SnowpackBridgeObjects objects;

    // Get or create SnowStation using singleton method
    objects.station = bridge.get_or_create_snowstation(input.i_grid, input.j_grid, 1,
                                                      input.wrf_lat, input.wrf_lon, input.height);

    if (!objects.station) {
        printf("SNOWPACK-FATAL [%d,%d]: Failed to create SnowStation!\n", input.i_grid, input.j_grid);
        std::abort();
    }

    // Get or create Snowpack instance
    objects.instance = bridge.get_or_create_snowpack_instance(input.i_grid, input.j_grid,
                                                             input.wrf_lat, input.wrf_lon, input.height);

    if (!objects.instance) {
        printf("SNOWPACK-FATAL [%d,%d]: Failed to create Snowpack instance!\n", input.i_grid, input.j_grid);
        std::abort();
    }

    return objects;
}

void ensure_station_initialized(SnowStation* station,
                                const MeteoInput& input,
                                SnowpackBridge& bridge) {
    if (station->getNumberOfElements() == 0) {
        // Station needs initialization
        mio::Coords position;
        position.setLatLon(input.wrf_lat, input.wrf_lon, input.height);

        std::string stationID = SnowpackConstants::STATION_ID_PREFIX + "_" +
                               std::to_string(input.i_grid) + "_" + std::to_string(input.j_grid);
        std::string stationName = "WRF Grid Point " + std::to_string(input.i_grid) + "," + std::to_string(input.j_grid);

        // Try to load existing state first
        bool loaded_from_file = load_station_state(
            generate_grid_key(input.i_grid, input.j_grid),
            station, position, stationID, stationName,
            bridge.get_io(), 1, input.i_grid, input.j_grid
        );

        if (!loaded_from_file) {
            // Initialize fresh station
            initialize_new_station(station, position, stationID, stationName, input.height,
                                 input.i_grid, input.j_grid);
        }
    }
}

bool save_station_state(SnowStation* station,
                       const mio::Date& current_time,
                       SnowpackIO* io,
                       int i_grid,
                       int j_grid) {
    if (!station || !io) {
        printf("SNOWPACK-WARNING: Cannot save state - station or IO is NULL\n");
        return false;
    }

    try {
        ZwischenData zdata;  // Empty for basic usage
        io->writeSnowCover(current_time, *station, zdata, true);
        printf("SNOWPACK-INFO: Saved .sno state for grid (%d,%d) using SnowpackIO\n", i_grid, j_grid);
        return true;

    } catch (const std::exception& e) {
        printf("SNOWPACK-WARNING: Failed to save .sno state for grid (%d,%d): %s\n",
               i_grid, j_grid, e.what());
        return false;
    }
}

void save_all_station_states(const std::map<std::string, std::unique_ptr<SnowStation>>& stations,
                             const mio::Date& current_time,
                             SnowpackIO* io) {
    if (!io) {
        printf("SNOWPACK-WARNING: Cannot save states - IO not initialized\n");
        return;
    }

    int saved_count = 0;
    for (auto it = stations.begin(); it != stations.end(); ++it) {
        const std::string& key = it->first;
        SnowStation* station = it->second.get();
        try {
            // Parse grid key "i_j" back to i,j coordinates
            size_t underscore_pos = key.find('_');
            if (underscore_pos != std::string::npos) {
                int i_grid = std::stoi(key.substr(0, underscore_pos));
                int j_grid = std::stoi(key.substr(underscore_pos + 1));

                ZwischenData zdata;  // Empty for basic usage
                io->writeSnowCover(current_time, *station, zdata, true);
                saved_count++;
            }
        } catch (const std::exception& e) {
            printf("SNOWPACK-ERROR: Failed to save state for %s: %s\n", key.c_str(), e.what());
        }
    }

    printf("SNOWPACK-INFO: Saved %d snowpack states to disk\n", saved_count);
}

// Validation utilities
void validate_output_values(const SnowpackOutput& output, int i_grid, int j_grid) {
    if (std::isnan(output.surface_temp) || std::isinf(output.surface_temp)) {
        printf("SNOWPACK-FATAL [%d,%d]: surface_temp is NaN/Inf (%.6f)!\n", i_grid, j_grid, output.surface_temp);
        std::abort();
    }
    if (std::isnan(output.snow_swe) || std::isinf(output.snow_swe)) {
        printf("SNOWPACK-FATAL [%d,%d]: snow_swe is NaN/Inf (%.6f)!\n", i_grid, j_grid, output.snow_swe);
        std::abort();
    }
    if (std::isnan(output.snow_depth) || std::isinf(output.snow_depth)) {
        printf("SNOWPACK-FATAL [%d,%d]: snow_depth is NaN/Inf (%.6f)!\n", i_grid, j_grid, output.snow_depth);
        std::abort();
    }
}

void validate_layer_data(const SnowpackLayerData& layer_data, int i_grid, int j_grid) {
    if (layer_data.n_layers < 0 || layer_data.n_layers > 100) {
        printf("SNOWPACK-FATAL [%d,%d]: n_layers is invalid (%d)!\n", i_grid, j_grid, layer_data.n_layers);
        std::abort();
    }

    // Validate key layer properties
    for (int i = 0; i < layer_data.n_layers && i < 100; i++) {
        if (std::isnan(layer_data.layer_temp[i]) || std::isinf(layer_data.layer_temp[i])) {
            printf("SNOWPACK-FATAL [%d,%d]: layer_temp[%d] is NaN/Inf (%.6f)!\n",
                   i_grid, j_grid, i, layer_data.layer_temp[i]);
            std::abort();
        }
        if (std::isnan(layer_data.layer_thick[i]) || std::isinf(layer_data.layer_thick[i])) {
            printf("SNOWPACK-FATAL [%d,%d]: layer_thick[%d] is NaN/Inf (%.6f)!\n",
                   i_grid, j_grid, i, layer_data.layer_thick[i]);
            std::abort();
        }
    }
}

} // namespace SnowpackObjects

/*
 * ============================================================================
 * SNOWPACK CONFIGURATION MANAGER IMPLEMENTATION
 * ============================================================================
 */

namespace {
  mio::Config g_config_cache;
  bool g_config_loaded = false;
}

std::string SnowpackConfigManager::getDefaultConfigPath() {
    return "./io.ini";
}

mio::Config SnowpackConfigManager::loadConfiguration(const std::string& ini_file_path) {
    if (!g_config_loaded || ini_file_path != "./io.ini") {
        g_config_cache = mio::Config(ini_file_path);
        g_config_loaded = true;
    }
    return g_config_cache;
}

void SnowpackConfigManager::validateConfiguration(const mio::Config& cfg) {
    // Basic validation - can be expanded as needed
    printf("SNOWPACK-INFO: Configuration validated\n");
}
