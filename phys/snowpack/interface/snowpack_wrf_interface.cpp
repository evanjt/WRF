/**
 * Modern SNOWPACK-WRF Interface 
 * Direct coupling to SNOWPACK v11.08 using clean APIs
 * Eliminates complex CRYOWRF coupling system
 */

#include <snowpack/snowpackCore/Snowpack.h>
#include <snowpack/DataClasses.h>
#include <snowpack/SnowpackConfig.h>
#include <meteoio/MeteoIO.h>

#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <cmath>

using namespace std;
using namespace mio;

// Grid point storage for SNOWPACK instances
struct GridPointData {
    std::unique_ptr<Snowpack> snowpack;
    std::unique_ptr<SnowpackConfig> config;
    SnowStation station;
    CurrentMeteo meteo;
    SurfaceFluxes fluxes;
    BoundCond boundaries;
    double cumulative_precip;
    bool initialized;
    
    GridPointData() : cumulative_precip(0.0), initialized(false) {}
};

// Global storage for grid points (could be improved with better memory management)
static std::map<std::pair<int,int>, GridPointData> grid_storage;

/**
 * Initialize SNOWPACK configuration for a grid point
 */
bool initialize_snowpack_point(int i, int j, double latitude, double longitude, double altitude) {
    auto grid_key = std::make_pair(i, j);
    auto& data = grid_storage[grid_key];
    
    if (data.initialized) {
        return true; // Already initialized
    }
    
    try {
        // Create minimal SNOWPACK configuration
        std::ostringstream config_stream;
        config_stream << "[General]" << std::endl;
        config_stream << "BUFFER_SIZE = 370" << std::endl;
        config_stream << "BUFF_BEFORE = 1.5" << std::endl; 
        config_stream << std::endl;
        
        config_stream << "[Input]" << std::endl;
        config_stream << "COORDSYS = proj4" << std::endl;
        config_stream << "TIME_ZONE = 0" << std::endl;
        config_stream << std::endl;
        
        config_stream << "[Interpolations1D]" << std::endl;
        config_stream << "WINDOW_SIZE = 86400" << std::endl;
        config_stream << std::endl;
        
        config_stream << "[SnowpackAdvanced]" << std::endl;
        config_stream << "VARIANT = DEFAULT" << std::endl;
        config_stream << "CALCULATION_STEP_LENGTH = 15" << std::endl;
        config_stream << "CHANGE_BC = false" << std::endl;
        config_stream << "THRESH_CHANGE_BC = -1" << std::endl;
        config_stream << "SNP_SOIL = false" << std::endl;
        config_stream << "SOIL_FLUX = false" << std::endl;
        config_stream << "GEO_HEAT = 0.06" << std::endl;
        config_stream << "CANOPY = false" << std::endl;
        config_stream << std::endl;
        
        // Create config from string stream
        data.config = std::make_unique<SnowpackConfig>(config_stream.str(), "WRF_COUPLING");
        
        // Create SNOWPACK instance
        data.snowpack = std::make_unique<Snowpack>(*data.config);
        
        // Initialize station data
        data.station.meta.position.setLatLon(latitude, longitude, altitude);
        data.station.meta.stationID = "WRF_" + std::to_string(i) + "_" + std::to_string(j);
        data.station.meta.stationName = data.station.meta.stationID;
        
        // Initialize with minimal snowpack
        data.station.resize(1); // Start with soil layer only
        data.station.Edata[0].Te = 273.15; // Soil temperature
        data.station.Edata[0].L = 0.2;     // Soil layer thickness  
        data.station.Edata[0].theta[SOLID] = 0.7;  // Soil solid fraction
        data.station.Edata[0].theta[LIQUID] = 0.2; // Soil water content
        data.station.Edata[0].theta[GAS] = 0.1;    // Soil air content
        
        // Initialize cumulative precipitation
        data.cumulative_precip = 0.0;
        
        data.initialized = true;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error initializing SNOWPACK for grid point (" << i << "," << j << "): " 
                  << e.what() << std::endl;
        return false;
    }
}

/**
 * Convert relative humidity to specific humidity (approximation)
 */
double rh_to_specific_humidity(double rh, double temp, double pressure) {
    // Saturation vapor pressure (Magnus formula)
    const double es = 610.78 * exp(17.27 * (temp - 273.15) / (temp - 35.86));
    
    // Actual vapor pressure  
    const double e = rh * es;
    
    // Specific humidity
    const double epsilon = 0.622;
    return epsilon * e / (pressure - (1.0 - epsilon) * e);
}

extern "C" {

/**
 * Main SNOWPACK physics interface called from Fortran
 */
void snowpack_physics(
    double temp_air, double humidity, double wind_speed, double wind_dir,
    double shortwave_in, double longwave_in, double precipitation, 
    double pressure, double height, 
    double* snow_swe, double* snow_depth, double* surface_temp,
    double* heat_flux_sensible, double* heat_flux_latent, 
    double* albedo, double* snow_coverage,
    double dt, int i_grid, int j_grid) {

    try {
        // Initialize if needed (use dummy coordinates for now)
        if (!initialize_snowpack_point(i_grid, j_grid, 46.0, 8.0, 1500.0)) {
            std::cerr << "Failed to initialize SNOWPACK for grid (" << i_grid << "," << j_grid << ")" << std::endl;
            return;
        }
        
        auto grid_key = std::make_pair(i_grid, j_grid);
        auto& data = grid_storage[grid_key];
        
        // Set current date/time (simplified - would need proper calendar integration)
        Date current_date;
        current_date.setFromComponents(2024, 1, 1, 12, 0, 0); // Placeholder
        
        // Fill meteorological data
        data.meteo.reset(*data.config);
        data.meteo.date = current_date;
        data.meteo.ta = temp_air;                    // Air temperature [K] 
        data.meteo.rh = std::max(0.1, std::min(1.0, humidity)); // Relative humidity [0-1]
        data.meteo.vw = std::max(0.1, wind_speed);   // Wind speed [m/s] (min for stability)
        data.meteo.dw = wind_dir;                    // Wind direction [deg]
        data.meteo.vw_max = wind_speed;              // Max wind speed
        data.meteo.iswr = std::max(0.0, shortwave_in);  // Incoming SW [W/m2]
        data.meteo.ilwr = std::max(50.0, longwave_in);  // Incoming LW [W/m2] (min for stability)  
        data.meteo.psum = std::max(0.0, precipitation * dt / 3600.0); // Precipitation [mm/h]
        data.meteo.psum_tech = 0.0;                  // Technical snow
        data.meteo.hs = *snow_depth;                 // Snow height [m]
        
        // Set measurement height
        data.meteo.z_wind = height;
        
        // Update cumulative precipitation
        data.cumulative_precip += data.meteo.psum;
        
        // Set boundary conditions  
        data.boundaries.lw_out = 0.0;      // Will be calculated
        data.boundaries.qs = 0.0;          // Will be calculated  
        data.boundaries.ql = 0.0;          // Will be calculated
        data.boundaries.hoar = 0.0;        // No surface hoar
        data.boundaries.sw_hor = shortwave_in;
        data.boundaries.sw_in = shortwave_in;
        data.boundaries.sw_out = 0.0;      // Will be calculated
        
        // Call SNOWPACK physics
        data.snowpack->runSnowpackModel(data.meteo, data.station, data.cumulative_precip, 
                                       data.boundaries, data.fluxes);
        
        // Extract results
        *surface_temp = data.station.Edata.back().Te;  // Surface temperature [K]
        *snow_depth = data.station.cH - data.station.Ground; // Snow depth [m] 
        *snow_swe = data.station.SWE;                   // Snow water equivalent [mm]
        
        // Surface fluxes
        *heat_flux_sensible = data.fluxes.qs;          // Sensible heat [W/m2]
        *heat_flux_latent = data.fluxes.ql;            // Latent heat [W/m2] 
        *albedo = data.fluxes.albedo;                   // Albedo [0-1]
        
        // Snow coverage (simple approximation)
        *snow_coverage = (*snow_depth > 0.001) ? 1.0 : 0.0;
        
    } catch (const std::exception& e) {
        std::cerr << "SNOWPACK physics error at grid (" << i_grid << "," << j_grid << "): " 
                  << e.what() << std::endl;
        
        // Return safe default values on error
        *heat_flux_sensible = 0.0;
        *heat_flux_latent = 0.0; 
        *albedo = 0.3;  // Default snow/soil albedo
        *snow_coverage = (*snow_depth > 0.001) ? 1.0 : 0.0;
        // Don't modify snow states on error
    }
}

/**
 * Cleanup function (optional - called at end of simulation)
 */
void snowpack_cleanup() {
    grid_storage.clear();
}

} // extern "C"