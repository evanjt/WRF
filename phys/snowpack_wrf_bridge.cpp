/*
 * SNOWPACK-WRF Bridge Implementation
 *
 * This file provides the C interface between WRF (Fortran) and SNOWPACK (C++)
 * It implements real SNOWPACK v11.08 physics integration.
 */

#include <cmath>
#include <cstdio>
#include <map>
#include <memory>

// SNOWPACK v11.08 headers - relative paths from phys/snowpack/
#include "meteoio/meteoio/MeteoIO.h"
#include "snowpack/snowpack/DataClasses.h"
#include "snowpack/snowpack/SnowpackConfig.h"
#include "snowpack/snowpack/snowpackCore/Snowpack.h"

// Global SNOWPACK instances - one per grid point
static std::map<std::pair<int, int>, std::unique_ptr<Snowpack>>
    snowpack_instances;
static std::map<std::pair<int, int>, SnowStation> snow_stations;
static std::unique_ptr<SnowpackConfig> global_config;
static bool config_initialized = false;

// Initialize SNOWPACK configuration
void initialize_snowpack_config() {
  if (config_initialized) return;
  
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Initializing SNOWPACK configuration...\n");

  // Create empty config first, then add keys
  mio::Config base_config;

  // Basic SNOWPACK configuration for WRF coupling
  base_config.addKey("CALCULATION_STEP_LENGTH", "Config", "900.0");  // 15 minutes
  base_config.addKey("SNOW_METAMORPHISM", "SnowpackAdvanced", "NIED");
  base_config.addKey("WATER_TRANSPORT", "SnowpackAdvanced", "BUCKET");
  base_config.addKey("CANOPY", "Snowpack", "false");
  base_config.addKey("SW_MODE", "Snowpack", "BOTH");
  base_config.addKey("HEIGHT_OF_METEO_VALUES", "Input", "2.0");
  base_config.addKey("HEIGHT_OF_WIND_VALUE", "Input", "10.0");
  base_config.addKey("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack", "false");
  base_config.addKey("SNOW_EROSION", "SnowpackAdvanced", "false");
  base_config.addKey("Alpine3D", "Snowpack", "false");
  base_config.addKey("RESEARCH", "SnowpackAdvanced", "false");
  base_config.addKey("METEO", "Input", "false");
  base_config.addKey("HARDNESS_PARAMETERIZATION", "SnowpackAdvanced", "MONTI");
  base_config.addKey("CHANGE_BC", "SnowpackAdvanced", "false");

  global_config = std::make_unique<SnowpackConfig>(base_config);
  config_initialized = true;
  
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Configuration initialized successfully!\n");
}

// C interface for Fortran binding
extern "C" {

/*
 * Real SNOWPACK physics interface called from WRF
 * This function must match the BIND(C) interface in module_sf_snowpack.F
 */
void snowpack_physics(double temp_air, double humidity, double wind_speed,
                      double wind_dir, double shortwave_in, double longwave_in,
                      double precipitation, double pressure, double height,
                      double* snow_swe, double* snow_depth,
                      double* surface_temp, double* heat_flux_sensible,
                      double* heat_flux_latent, double* albedo,
                      double* snow_coverage, double dt, int i_grid,
                      int j_grid) {
  
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: snowpack_physics() called for grid (%d,%d), temp=%.2fK\n", 
         i_grid, j_grid, temp_air);
  
  // Initialize configuration on first call
  initialize_snowpack_config();

  // Grid point identifier
  std::pair<int, int> grid_point(i_grid, j_grid);
  // Create SNOWPACK instance for this grid point if it doesn't exist
  if (snowpack_instances.find(grid_point) == snowpack_instances.end()) {
    snowpack_instances[grid_point] = std::make_unique<Snowpack>(*global_config);
    snow_stations[grid_point] = SnowStation();

    // Initialize snow station with basic setup
    SnowStation& xdata = snow_stations[grid_point];

    // Create minimal initialization data for SNOWPACK
    SN_SNOWSOIL_DATA ssdata;
    ssdata.meta.stationID =
        "WRF_" + std::to_string(i_grid) + "_" + std::to_string(j_grid);
    ssdata.meta.position.setXY(i_grid, j_grid,
                               0.0);  // Grid coordinates as position
    ssdata.Height = height;
    ssdata.nN = 1;       // Start with 1 node (ground)
    ssdata.nLayers = 0;  // No snow layers initially

    xdata.initialize(ssdata, 0);  // Initialize with sector 0

    printf("SNOWPACK-PHYSICS: Initialized grid point (%d,%d)\n", i_grid,
           j_grid);
  }

  // Get references to this grid point's SNOWPACK instance and data
  Snowpack& snowpack = *snowpack_instances[grid_point];
  SnowStation& xdata = snow_stations[grid_point];

  // Prepare meteorological data for SNOWPACK
  CurrentMeteo mdata(*global_config);
  mdata.date = mio::Date();   // Current timestep
  mdata.ta = temp_air;        // Air temperature [K]
  mdata.rh = humidity;        // Relative humidity [0-1]
  mdata.vw = wind_speed;      // Wind speed [m/s]
  mdata.dw = wind_dir;        // Wind direction [deg]
  mdata.iswr = shortwave_in;  // Incoming shortwave [W/m²]
  mdata.ea = longwave_in /
             (5.67e-8 * pow(temp_air, 4));  // Atmospheric emissivity from LW
  mdata.psum = precipitation * dt;          // Precipitation [mm]
  mdata.psum_ph = (temp_air > 273.15) ? 1.0 : 0.0;  // Phase: 0=snow, 1=rain

  // Update snow station state from WRF
  // Note: Real implementation would maintain multi-layer snowpack
  // For now, we'll use SNOWPACK's native state management

  // Boundary conditions
  BoundCond bdata;

  // Surface fluxes (output from SNOWPACK)
  SurfaceFluxes sdata;

  // Cumulative precipitation (maintained by SNOWPACK)
  static double cumu_precip = 0.0;

  try {
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: Calling SNOWPACK model for grid (%d,%d), precip=%.3f\n", 
           i_grid, j_grid, mdata.psum);
    
    // Call SNOWPACK model for this timestep
    snowpack.runSnowpackModel(mdata, xdata, cumu_precip, bdata, sdata);
    
    printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: SNOWPACK model completed successfully for grid (%d,%d)\n", 
           i_grid, j_grid);

    // Extract results for WRF
    *surface_temp =
        mdata.tss;  // Surface temperature from SNOWPACK (in CurrentMeteo)
    *snow_swe = xdata.swe;   // Snow water equivalent [mm] (direct member)
    *snow_depth = xdata.cH;  // Total calculated height [m] (direct member)

    // Energy fluxes from SNOWPACK
    *heat_flux_sensible = sdata.qs;  // Sensible heat flux [W/m²]
    *heat_flux_latent = sdata.ql;    // Latent heat flux [W/m²]

    // Surface properties
    *albedo = xdata.Albedo;                              // Snow albedo [-]
    *snow_coverage = (*snow_depth > 0.001) ? 1.0 : 0.0;  // Snow coverage [-]

  } catch (const std::exception& e) {
    // SNOWPACK physics failure - terminate with error
    printf("FATAL ERROR: SNOWPACK physics failed at grid(%d,%d): %s\n", i_grid,
           j_grid, e.what());
    printf(
        "SNOWPACK integration cannot continue. Check SNOWPACK configuration "
        "and data.\n");

    // Set error state and terminate
    *surface_temp = -999.0;  // Error indicator
    *snow_swe = -999.0;
    *snow_depth = -999.0;
    *heat_flux_sensible = -999.0;
    *heat_flux_latent = -999.0;
    *albedo = -999.0;
    *snow_coverage = -999.0;

    // Force program termination to prevent corrupted model state
    exit(1);
  }
  
  printf("SNOWPACK-DEBUG [C++/snowpack_wrf_bridge.cpp]: snowpack_physics() completed for grid (%d,%d)\n", 
         i_grid, j_grid);
}

}  // extern "C"
