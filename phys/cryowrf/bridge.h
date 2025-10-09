
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "structs.h"
#include "utils.h"
#include "config.h"
#include "meteoio/meteoio/MeteoIO.h"
#include "snowpack/snowpack/DataClasses.h"
#include "snowpack/snowpack/SnowpackConfig.h"
#include "snowpack/snowpack/snowpackCore/Snowpack.h"
#include "snowpack/snowpack/plugins/SnowpackIO.h"

namespace SnowpackConstants {
  constexpr double T_CRAZY_MAX_KELVIN = 400.0;  //  127°C - configurable via T_CRAZY_MAX
  constexpr double T_CRAZY_MIN_KELVIN = 100.0;  // -173°C - configurable via T_CRAZY_MIN

  // Physical thresholds - can be overridden by configuration
  constexpr double PRECIP_PHASE_THRESHOLD_K = 273.65;  // 0.5°C threshold for snow/rain
  constexpr double BOTTOM_TEMP_OFFSET_K = 5.0;         // Offset below air temp for ground temp
  constexpr double SNOW_DENSITY_FALLBACK = 100.0;      // kg/m³ fallback density

  // Default station metadata
  const std::string STATION_ID_PREFIX = "WRF_GRID";  // Station ID prefix for SNOWPACK
}

struct SnowpackBridgeObjects {
    SnowStation* station;
    Snowpack* instance;
};


// SnowpackBridge singleton class for state management
class SnowpackBridge {
private:
    // Private constructor for singleton pattern
    SnowpackBridge() = default;

    // Delete copy constructor and assignment operator
    SnowpackBridge(const SnowpackBridge&) = delete;
    SnowpackBridge& operator=(const SnowpackBridge&) = delete;

    // Configuration and I/O
    std::unique_ptr<SnowpackConfig> config_;
    std::unique_ptr<SnowpackIO> io_;
    bool config_initialized_ = false;
    std::string config_file_path_;

    // Time management
    mio::Date current_simulation_date_;
    bool time_initialized_ = false;
    double calculation_step_length_ = 0.0;

    // Persistent object storage per grid point
    std::map<std::string, std::unique_ptr<SnowStation>> grid_snowstations_;
    std::map<std::string, std::unique_ptr<Snowpack>> grid_snowpack_instances_;

    // Call tracking for debugging
    int execute_call_count_ = 0;

public:
    // Get singleton instance
    static SnowpackBridge& instance() {
        static SnowpackBridge bridge;
        return bridge;
    }

    // Configuration management
    void initialize_config(const std::string& ini_path);
    bool is_config_initialized() const { return config_initialized_; }
    SnowpackConfig* get_config() const { return config_.get(); }
    SnowpackIO* get_io() const { return io_.get(); }

    // Time management
    void initialize_time(int start_year, int start_month, int start_day,
                        int start_hour, int start_minute);
    bool is_time_initialized() const { return time_initialized_; }
    const mio::Date& get_current_time() const {
        return current_simulation_date_;
    }

    // Object management
    SnowStation* get_or_create_snowstation(
        int i_grid, int j_grid, int wrf_domain_id,
        double wrf_lat, double wrf_lon, double wrf_alt
    );
    Snowpack* get_or_create_snowpack_instance(
        int i_grid, int j_grid, double wrf_lat, double wrf_lon, double wrf_alt
    );
    SnowStation* get_existing_snowstation(int i_grid, int j_grid);

    // State persistence
    void save_snowstation_state(int i_grid, int j_grid);
    void save_all_snowpack_states();

    // Core physics execution
    void execute_snowpack(
        const MeteoInput& input,
        SnowpackOutput& output,
        SnowpackLayerData* layer_data = nullptr,
        BudgetData* budget_data = nullptr
    );
};

// Forward declarations for utility functions
void prepare_meteo_data(const MeteoInput& input,
                        CurrentMeteo& Mdata,
                        const mio::Date& current_time,
                        double current_snow_depth,
                        const SnowpackConfig* config);

void extract_surface_outputs(const SnowStation& station,
                             const SurfaceFluxes& fluxes,
                             const BoundCond& bc,
                             SnowpackOutput& output,
                             double temp_air_fallback);

void extract_layer_data(const SnowStation& station, SnowpackLayerData& layer_data);

void extract_budget_data(const SurfaceFluxes& fluxes,
                        const BoundCond& bc,
                        double cumu_precip,
                        const MeteoInput& input,
                        BudgetData& budget_data);

void advance_simulation_time(const MeteoInput& input,
                           mio::Date& current_time,
                           double calculation_step_length);

SnowpackBridgeObjects create_station_objects(const MeteoInput& input,
                                           SnowpackBridge& bridge);

void ensure_station_initialized(SnowStation* station,
                                const MeteoInput& input,
                                SnowpackBridge& bridge);

bool save_station_state(SnowStation* station,
                       const mio::Date& current_time,
                       SnowpackIO* io,
                       int i_grid,
                       int j_grid);

void save_all_station_states(const std::map<std::string, std::unique_ptr<SnowStation>>& stations,
                             const mio::Date& current_time,
                             SnowpackIO* io);

void validate_output_values(const SnowpackOutput& output, int i_grid, int j_grid);

void validate_layer_data(const SnowpackLayerData& layer_data, int i_grid, int j_grid);
