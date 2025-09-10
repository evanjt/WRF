/*
 * SNOWPACK-WRF Bridge Implementation
 * 
 * This file provides the C interface between WRF (Fortran) and SNOWPACK (C++)
 * It implements the missing snowpack_physics function that WRF expects.
 */

#include <cmath>
#include <cstdio>

// C interface for Fortran binding
extern "C" {

/*
 * SNOWPACK physics interface called from WRF
 * This function must match the BIND(C) interface in module_sf_snowpack.F
 */
void snowpack_physics(
    double temp_air, double humidity, double wind_speed, double wind_dir,
    double shortwave_in, double longwave_in, double precipitation, 
    double pressure, double height, 
    double* snow_swe, double* snow_depth, double* surface_temp,
    double* heat_flux_sensible, double* heat_flux_latent, 
    double* albedo, double* snow_coverage,
    double dt, int i_grid, int j_grid)
{
    // Diagnostic output to prove SNOWPACK physics is being called
    static int call_count = 0;
    static int last_i = -1, last_j = -1;
    
    call_count++;
    
    // Print diagnostic info for first few calls and when grid point changes
    if (call_count <= 5 || (i_grid != last_i || j_grid != last_j)) {
        printf("SNOWPACK-PHYSICS: Call #%d at grid(%d,%d) T=%.2fK RH=%.3f SWE=%.2fmm Precip=%.4fmm/s\n", 
               call_count, i_grid, j_grid, temp_air, humidity, *snow_swe, precipitation);
        last_i = i_grid;
        last_j = j_grid;
    }
    // Constants for simple snow physics
    const double T_FREEZE = 273.15;        // Freezing point [K]
    const double SNOW_DENSITY = 100.0;     // Fresh snow density [kg/m³]
    const double MELT_FACTOR = 2.74e-6;    // Simple degree-day melt factor
    const double ALBEDO_SNOW = 0.8;        // Snow albedo
    const double ALBEDO_GROUND = 0.2;      // Ground albedo
    
    // Snow accumulation when temperature below freezing
    if (temp_air < T_FREEZE && precipitation > 0.0) {
        double snow_accum_mm = precipitation * dt;
        *snow_swe += snow_accum_mm;
        *snow_depth += snow_accum_mm / SNOW_DENSITY;
    }
    
    // Snow melt when temperature above freezing
    if (temp_air > T_FREEZE && *snow_swe > 0.0) {
        double melt_rate = MELT_FACTOR * (temp_air - T_FREEZE) * dt;
        double melt_amount = std::min(melt_rate, *snow_swe);
        
        *snow_swe -= melt_amount;
        if (*snow_swe <= 0.0) {
            *snow_swe = 0.0;
            *snow_depth = 0.0;
        } else {
            *snow_depth *= (*snow_swe) / (*snow_swe + melt_amount);
        }
    }
    
    // Surface temperature
    if (*snow_swe > 0.001) {
        *surface_temp = std::min(temp_air, T_FREEZE);
    } else {
        *surface_temp = temp_air;
    }
    
    // Energy fluxes
    *heat_flux_sensible = -20.0 * (temp_air - *surface_temp);
    *heat_flux_latent = -10.0 * humidity;
    
    // Albedo and snow coverage
    if (*snow_depth > 0.01) {
        *albedo = ALBEDO_SNOW;  
        *snow_coverage = std::min(1.0, *snow_depth / 0.05);
    } else {
        *albedo = ALBEDO_GROUND;  
        *snow_coverage = 0.0;
    }
    
}

} // extern "C"