/*
 * CRYOWRF object factory and management for SNOWPACK integration
 * CRYOWRF provides the bridge between WRF and SNOWPACK/MeteoIO
 * Extracted from bridge.cpp to improve code organization
 */

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include "bridge.h"
#include "structs.h"
#include "meteoio/meteoio/MeteoIO.h"
#include "snowpack/snowpack/DataClasses.h"
#include "snowpack/snowpack/SnowpackConfig.h"
#include "snowpack/snowpack/snowpackCore/Snowpack.h"
#include "snowpack/snowpack/plugins/SnowpackIO.h"

namespace {

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
                        std::mutex& io_mutex,
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

        // Attempt to read existing snowpack state
        SN_SNOWSOIL_DATA ssdata;
        ZwischenData zdata;
        mio::Date profile_date;

        std::lock_guard<std::mutex> io_lock(io_mutex);
        io->readSnowCover(sno_filename, stationID, ssdata, zdata, false);

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

} // anonymous namespace

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

        const auto& constants = SnowpackConstants::get();
        std::string stationID = constants.station_id_prefix + "_" +
                               std::to_string(input.i_grid) + "_" + std::to_string(input.j_grid);
        std::string stationName = "WRF Grid Point " + std::to_string(input.i_grid) + "," + std::to_string(input.j_grid);

        // Try to load existing state first
        bool loaded_from_file = load_station_state(
            generate_grid_key(input.i_grid, input.j_grid),
            station, position, stationID, stationName,
            bridge.get_io(), bridge.io_mutex(), 1, input.i_grid, input.j_grid
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
                       std::mutex& io_mutex,
                       int i_grid,
                       int j_grid) {
    if (!station || !io) {
        printf("SNOWPACK-WARNING: Cannot save state - station or IO is NULL\n");
        return false;
    }

    try {
        ZwischenData zdata;  // Empty for basic usage
        std::lock_guard<std::mutex> io_lock(io_mutex);
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
                             SnowpackIO* io,
                             std::mutex& io_mutex) {
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
                std::lock_guard<std::mutex> io_lock(io_mutex);
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
    if (std::isnan(output.friction_velocity) || std::isinf(output.friction_velocity)) {
        printf("SNOWPACK-FATAL [%d,%d]: friction_velocity is NaN/Inf (%.6f)!\n", i_grid, j_grid, output.friction_velocity);
        std::abort();
    }
    if (output.roughness_mom <= 0.0 || std::isnan(output.roughness_mom) || std::isinf(output.roughness_mom)) {
        printf("SNOWPACK-FATAL [%d,%d]: roughness_mom invalid (%.6f)!\n", i_grid, j_grid, output.roughness_mom);
        std::abort();
    }
    if (output.roughness_heat <= 0.0 || std::isnan(output.roughness_heat) || std::isinf(output.roughness_heat)) {
        printf("SNOWPACK-FATAL [%d,%d]: roughness_heat invalid (%.6f)!\n", i_grid, j_grid, output.roughness_heat);
        std::abort();
    }
    if (std::isnan(output.latent_flux_kg) || std::isinf(output.latent_flux_kg)) {
        printf("SNOWPACK-FATAL [%d,%d]: latent_flux_kg is NaN/Inf (%.6f)!\n", i_grid, j_grid, output.latent_flux_kg);
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
