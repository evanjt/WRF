/**
 * Pure C Bridge to SNOWPACK 
 * Avoids C++ complexity in WRF build system
 * Calls external SNOWPACK executable or simplified physics
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>  /* For getpid() */

/**
 * Simplified snow physics implementation in pure C
 * Based on basic energy balance and empirical relationships
 */
void snowpack_physics_simple(double temp_air, double humidity, double wind_speed, 
                            double shortwave, double longwave, double precipitation,
                            double pressure, double dt,
                            double *snow_swe, double *snow_depth, double *surface_temp,
                            double *heat_sensible, double *heat_latent, double *albedo) {
    
    const double SNOW_DENSITY = 300.0;  // kg/m3 - typical new snow density
    const double MELT_FACTOR = 0.003;   // Empirical melt factor
    const double FREEZE_TEMP = 273.15;  // K
    
    // Snow accumulation
    if (temp_air < FREEZE_TEMP && precipitation > 0.0) {
        double new_snow = precipitation * dt / 1000.0;  // Convert mm to m
        *snow_depth += new_snow;
        *snow_swe += precipitation * dt;  // mm
    }
    
    // Snow melt (simple degree-day approach)
    if (temp_air > FREEZE_TEMP && *snow_depth > 0.0) {
        double melt_rate = MELT_FACTOR * (temp_air - FREEZE_TEMP);  // m/s
        double melt_depth = melt_rate * dt;
        
        if (melt_depth > *snow_depth) {
            melt_depth = *snow_depth;
        }
        
        *snow_depth -= melt_depth;
        *snow_swe -= melt_depth * SNOW_DENSITY;  // Approximate SWE reduction
        
        if (*snow_depth < 0.001) {
            *snow_depth = 0.0;
            *snow_swe = 0.0;
        }
    }
    
    // Surface temperature (simple approach)
    if (*snow_depth > 0.01) {
        *surface_temp = fmin(temp_air, FREEZE_TEMP);  // Snow surface can't exceed 0°C
    } else {
        *surface_temp = temp_air - 2.0;  // Ground temperature offset
    }
    
    // Surface energy fluxes (simplified)
    double temp_diff = temp_air - *surface_temp;
    *heat_sensible = 10.0 * wind_speed * temp_diff;  // Simplified sensible heat
    *heat_latent = 50.0;  // Simple latent heat flux
    
    // Albedo
    if (*snow_depth > 0.01) {
        *albedo = 0.8 - 0.3 * fmax(0.0, (temp_air - FREEZE_TEMP) / 5.0);  // Aging snow
        *albedo = fmax(0.4, fmin(0.85, *albedo));
    } else {
        *albedo = 0.15;  // Soil albedo
    }
}

/**
 * External SNOWPACK executable interface
 * Calls standalone SNOWPACK model if available
 */
int call_external_snowpack(double temp_air, double humidity, double wind_speed, 
                          double shortwave, double longwave, double precipitation,
                          double *snow_swe, double *snow_depth, double *surface_temp,
                          double *heat_sensible, double *heat_latent, double *albedo) {
    
    char input_file[256], output_file[256];
    char command[512];
    FILE *fp_in, *fp_out;
    int result;
    
    // Create unique filenames for this process
    sprintf(input_file, "snowpack_input_%d.dat", getpid());
    sprintf(output_file, "snowpack_output_%d.dat", getpid());
    
    // Write input file
    fp_in = fopen(input_file, "w");
    if (!fp_in) return -1;
    
    fprintf(fp_in, "# SNOWPACK input\n");
    fprintf(fp_in, "%.2f  # Air temperature [K]\n", temp_air);
    fprintf(fp_in, "%.3f  # Relative humidity [0-1]\n", humidity);
    fprintf(fp_in, "%.2f  # Wind speed [m/s]\n", wind_speed);
    fprintf(fp_in, "%.1f  # Shortwave radiation [W/m2]\n", shortwave);
    fprintf(fp_in, "%.1f  # Longwave radiation [W/m2]\n", longwave);
    fprintf(fp_in, "%.3f  # Precipitation [mm]\n", precipitation);
    fprintf(fp_in, "%.2f  # Snow SWE [mm]\n", *snow_swe);
    fprintf(fp_in, "%.3f  # Snow depth [m]\n", *snow_depth);
    fclose(fp_in);
    
    // Call external SNOWPACK (if available)
    sprintf(command, "snowpack_simple %s %s 2>/dev/null", input_file, output_file);
    result = system(command);
    
    if (result == 0) {
        // Read results
        fp_out = fopen(output_file, "r");
        if (fp_out) {
            fscanf(fp_out, "%lf %lf %lf %lf %lf %lf", 
                   surface_temp, heat_sensible, heat_latent, 
                   albedo, snow_swe, snow_depth);
            fclose(fp_out);
            
            // Cleanup
            remove(input_file);
            remove(output_file);
            return 0;
        }
    }
    
    // Cleanup on failure
    remove(input_file);
    remove(output_file);
    return -1;
}

/**
 * Main interface called from Fortran
 * Falls back to simplified physics if external SNOWPACK not available
 */
void snowpack_physics(double temp_air, double humidity, double wind_speed, double wind_dir,
                     double shortwave, double longwave, double precipitation, 
                     double pressure, double height, double dt, int i_grid, int j_grid,
                     double *snow_swe, double *snow_depth, double *surface_temp,
                     double *heat_sensible, double *heat_latent, 
                     double *albedo, double *snow_coverage) {
    
    // Try external SNOWPACK first
    if (call_external_snowpack(temp_air, humidity, wind_speed, shortwave, longwave, 
                               precipitation, snow_swe, snow_depth, surface_temp,
                               heat_sensible, heat_latent, albedo) != 0) {
        
        // Fall back to simplified physics
        snowpack_physics_simple(temp_air, humidity, wind_speed, shortwave, longwave, 
                               precipitation, pressure, dt, snow_swe, snow_depth, 
                               surface_temp, heat_sensible, heat_latent, albedo);
    }
    
    // Calculate snow coverage
    *snow_coverage = (*snow_depth > 0.001) ? 1.0 : 0.0;
}