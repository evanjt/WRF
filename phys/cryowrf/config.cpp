
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <limits.h>

#include "meteoio/meteoio/MeteoIO.h"
#include "snowpack/SnowpackConfig.h"
#include "config.h"

// SnowpackConfigManager implementation
mio::Config SnowpackConfigManager::loadConfiguration(const std::string& ini_file_path) {
    try {
        // Use SnowpackConfig instead of mio::Config to get SNOWPACK defaults

        // Get absolute path
        char abs_path[PATH_MAX];
        if (realpath(ini_file_path.c_str(), abs_path) != NULL) {
            printf("SNOWPACK-DEBUG: Loading SnowpackConfig from %s (absolute: %s)\n", ini_file_path.c_str(), abs_path);
        } else {
            printf("SNOWPACK-DEBUG: Loading SnowpackConfig from %s (could not resolve absolute path)\n", ini_file_path.c_str());
        }

        SnowpackConfig config(ini_file_path);

        // Check if SNOW_WRITE exists in the configuration
        if (config.keyExists("SNOW_WRITE", "Output")) {
            printf("SNOWPACK-DEBUG: SNOW_WRITE found in Output section\n");
            std::string snow_write_value = config.get("SNOW_WRITE", "Output");
            printf("SNOWPACK-DEBUG: SNOW_WRITE value = '%s'\n", snow_write_value.c_str());
        } else {
            printf("SNOWPACK-DEBUG: SNOW_WRITE NOT found in Output section\n");
        }

        // Check for other missing parameters
        printf("SNOWPACK-DEBUG: Checking for missing parameters...\n");
        printf("SNOWPACK-DEBUG: RIME_INDEX exists: %s\n", config.keyExists("RIME_INDEX", "SnowpackAdvanced") ? "YES" : "NO");
        printf("SNOWPACK-DEBUG: WATER_LAYER exists: %s\n", config.keyExists("WATER_LAYER", "SnowpackAdvanced") ? "YES" : "NO");
        printf("SNOWPACK-DEBUG: TIME_ZONE exists: %s\n", config.keyExists("TIME_ZONE", "Input") ? "YES" : "NO");

        const double calculation_step_length = config.get("CALCULATION_STEP_LENGTH", "Snowpack");
        const double meteo_step_length = calculation_step_length * 60.0; // Convert minutes to seconds

        // Add METEO_STEP_LENGTH to config dynamically (CRYOWRF pattern)
        std::stringstream ss_meteo_length;
        ss_meteo_length << meteo_step_length;
        config.addKey("METEO_STEP_LENGTH", "Snowpack", ss_meteo_length.str());

        // Verify the key was added correctly
        std::string verify_value;
        config.getValue("METEO_STEP_LENGTH", "Snowpack", verify_value);
        printf("SNOWPACK-CONFIG: METEO_STEP_LENGTH set to %s seconds (from CALCULATION_STEP_LENGTH=%.1f minutes)\n",
               verify_value.c_str(), calculation_step_length);

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
        {"CANOPY", "Snowpack"},                 // Required by Snowpack.cc:174
        {"HEIGHT_OF_METEO_VALUES", "Snowpack"}, // Required by Snowpack.cc:179
        {"HEIGHT_OF_WIND_VALUE", "Snowpack"},   // Required by Meteo.cc:56
        {"ROUGHNESS_LENGTH", "Snowpack"},       // Required by Meteo.cc:50
        {"SW_MODE", "Snowpack"},                // Common requirement
        {"ATMOSPHERIC_STABILITY", "Snowpack"},  // Required by Meteo.cc:44
        {"GEO_HEAT", "Snowpack"},               // Required for energy balance
        {"VARIANT", "SnowpackAdvanced"}
    };

    for (const auto& param : required_params) {
        try {
            std::string value;
            cfg.getValue(param.first, param.second, value);
        } catch (const std::exception& e) {
            printf("SNOWPACK-ERROR [C++/SnowpackConfigManager]: Missing required parameter %s::%s\n",
                   param.second.c_str(), param.first.c_str());
            throw std::runtime_error("Configuration validation failed: missing " + param.first);
        }
    }
    }

std::string SnowpackConfigManager::getDefaultConfigPath() {
    // Allow environment variable override for config path
    const char* config_env = std::getenv("SNOWPACK_CONFIG_PATH");
    if (config_env) {
        printf("SNOWPACK-INFO: Using config path from environment: %s\n", config_env);
        return std::string(config_env);
    }
    return "./io.ini";  // Default path in WRF run directory
}