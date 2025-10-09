#ifndef SNOWPACK_BRIDGE_STRUCTS_H
#define SNOWPACK_BRIDGE_STRUCTS_H

// Data structures for organized parameter passing
struct MeteoInput {
    double temp_air;
    double humidity;
    double wind_speed;
    double wind_dir;
    double shortwave_in;
    double longwave_in;
    double precipitation;
    double pressure;
    double height;
    double dt;
    int i_grid;
    int j_grid;
    double wrf_lat;
    double wrf_lon;
};

struct SnowpackOutput {
    double snow_swe;
    double snow_depth;
    double surface_temp;
    double heat_flux_sensible;
    double heat_flux_latent;
    double albedo;
    double snow_coverage;
    double friction_velocity;
    double stability_param;
};

struct SnowpackLayerData {
    int n_layers;

    double layer_temp[100];
    double layer_thick[100];
    double layer_vol_ice[100];
    double layer_vol_water[100];
    double layer_vol_air[100];
    double layer_grain_radius[100];
    double layer_bond_radius[100];
    double layer_dendricity[100];
    double layer_sphericity[100];
};

struct BudgetData {
    double mass_precip;
    double mass_sublim;
    double mass_melt;
    double mass_swe;
    double mass_refreeze;
    double energy_lw_in;
    double energy_lw_out;
    double energy_sw_in;
    double energy_sw_out;
    double energy_sensible;
    double energy_latent;
    double energy_ground_flux;
    double energy_rain;
    double energy_total;
};

#endif // SNOWPACK_BRIDGE_STRUCTS_H
