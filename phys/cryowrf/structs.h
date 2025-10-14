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
  static constexpr int MAX_SOIL_LAYERS = 4;

  double snow_swe = 0.0;
  double snow_depth = 0.0;
  double surface_temp = 0.0;
  double heat_flux_sensible = 0.0;
  double heat_flux_latent = 0.0;
  double albedo = 0.0;
  double snow_coverage = 0.0;
  double friction_velocity = 0.0;
  double stability_param = 0.0;
  double roughness_mom = 0.0;
  double roughness_heat = 0.0;
  double latent_flux_kg = 0.0;
  double moisture_flux = 0.0;

  // Soil properties from SNOWPACK ground interface
  // CRYOWRF SOURCE: module_sf_snowpacklsm.F:364-366 updates soil moisture
  double soil_moisture_volumetric = 0.0;  // phiSoil from SNOWPACK [m³/m³]
  double soil_temperature = 273.15;       // ts0 from SNOWPACK [K]
  double soil_density = 0.0;              // SoilRho from SNOWPACK [kg/m³]
  double soil_conductivity = 0.0;         // SoilK from SNOWPACK
  double soil_heat_capacity = 0.0;        // SoilC from SNOWPACK
  double soil_moisture_liquid = 0.0;      // Calculated liquid fraction
  double soil_moisture_avail = 0.0;       // Calculated availability [0-1]
  double soil_moisture_total = 0.0;       // Calculated total water [mm]

  int soil_layer_count = 0;
  double soil_temp_layers[MAX_SOIL_LAYERS] = {0.0};
  double soil_moisture_vol_layers[MAX_SOIL_LAYERS] = {0.0};
  double soil_moisture_liq_layers[MAX_SOIL_LAYERS] = {0.0};
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
  double layer_cdot[100];
  double layer_meta[100];
  double layer_deposition_julian[100];
  double layer_graintype[100];
  double layer_marker[100];
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

#endif  // SNOWPACK_BRIDGE_STRUCTS_H
