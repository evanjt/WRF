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
    bool use_canopy = false;
    bool use_soil = false;

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

        // printf("SNOWPACK-INIT [%d,%d]: Attempting to load .sno file: %s\n", i_grid, j_grid, sno_filename.c_str());

        try {
            // Attempt to read existing snowpack state
            SN_SNOWSOIL_DATA ssdata;
            ZwischenData zdata;
            mio::Date profile_date;

            // printf("SNOWPACK-INIT [%d,%d]: Calling readSnowCover()...\n", i_grid, j_grid);
            io_->readSnowCover(sno_filename, stationID, ssdata, zdata, false);
            // printf("SNOWPACK-INIT [%d,%d]: readSnowCover() returned %zu layers\n", i_grid, j_grid, ssdata.nLayers);

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
            // printf("SNOWPACK-INIT [%d,%d]: Successfully loaded station state from %s\n", i_grid, j_grid, sno_filename.c_str());

        } catch (const std::exception& e) {
            // No existing state file - start with fresh snowpack
            // Only show first few initialization failures to avoid spam
            static int init_failure_count = 0;
            if (init_failure_count < 5) {
                printf("SNOWPACK-INFO: No existing state for grid (%d,%d), starting fresh\n", i_grid, j_grid);
                init_failure_count++;
            }
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
    // Validate grid coordinates first
    if (input.i_grid < 0 || input.i_grid > 10000 || input.j_grid < 0 || input.j_grid > 10000) {
        printf("SNOWPACK-ERROR: Invalid grid coordinates (%d,%d) - outside expected range [0-10000]\n",
               input.i_grid, input.j_grid);
        printf("SNOWPACK-ERROR: This indicates a problem with WRF loop bounds initialization\n");
        return;
    }

    // Track physics calls for debugging
    execute_call_count_++;
    if (execute_call_count_ <= 5 || (execute_call_count_ % 1000 == 0)) {
        if (execute_call_count_ <= 5) {
            printf("SNOWPACK-INFO: Executed snowpack call #%d - Grid (%d,%d)\n", execute_call_count_, input.i_grid, input.j_grid);
        } else {
            printf("SNOWPACK-INFO: Executed %d snowpack calls\n", execute_call_count_);
        }
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

        // Time advancement logging removed to avoid verbosity

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

        // IMPORTANT: Ensure METEO_STEP_LENGTH is available in the final config
        // SnowpackConfig might re-read the ini file, so we need to inject METEO_STEP_LENGTH again
        const double calculation_step_length = config_->get("CALCULATION_STEP_LENGTH", "Snowpack");
        const double meteo_step_length = calculation_step_length * 60.0; // Convert minutes to seconds

        std::stringstream ss_meteo_length;
        ss_meteo_length << meteo_step_length;
        config_->addKey("METEO_STEP_LENGTH", "Snowpack", ss_meteo_length.str());

        // Verify it was added
        std::string verify_value;
        config_->getValue("METEO_STEP_LENGTH", "Snowpack", verify_value);
        printf("SNOWPACK-CONFIG: Final METEO_STEP_LENGTH set to %s seconds (from CALCULATION_STEP_LENGTH=%.1f minutes)\n",
               verify_value.c_str(), calculation_step_length);

        io_ = std::unique_ptr<SnowpackIO>(new SnowpackIO(*config_));

        config_file_path_ = ini_path;
        config_initialized_ = true;

        printf("SNOWPACK-INFO: Configuration initialized from %s\n", ini_path.c_str());

    } catch (const std::exception& e) {
        printf("SNOWPACK-FATAL: Configuration failed for %s: %s\n", ini_path.c_str(), e.what());
        throw;
    }
}
