#include "snowpack_bridge.h"

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


// // Save all snowpack states (called from Fortran for periodic saves)
// void save_all_snowpack_states_c() {
//     SnowpackBridge& bridge = SnowpackBridge::instance();
//     bridge.save_all_snowpack_states();
// }

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
void extract_snowpack_layers_c_(int i_grid, int j_grid,
                               float* layer_temps, float* layer_thick, float* layer_voli, float* layer_volw, float* layer_volv,
                               float* layer_rg, float* layer_rb, float* layer_dd, float* layer_sp, int* n_layers) {
    extract_snowpack_layers_c(i_grid, j_grid, layer_temps, layer_thick, layer_voli, layer_volw, layer_volv,
                              layer_rg, layer_rb, layer_dd, layer_sp, n_layers);
}

} // extern "C"
