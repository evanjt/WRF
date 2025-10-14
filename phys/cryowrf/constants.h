#ifndef SNOWPACK_CONSTANTS_H
#define SNOWPACK_CONSTANTS_H

#include <string>

class SnowpackConfig;

namespace SnowpackConstants {

struct Values {
  double t_crazy_min = 100.0;
  double t_crazy_max = 400.0;
  double precip_phase_threshold = 273.65;
  double surface_temp_guess = 400.0;
  double snow_density_fallback = 100.0;
  double snow_roughness_length = 0.002;
  double bare_roughness_length = 0.01;
  double snow_depth_roughness_threshold = 0.03;
  double measurement_height_floor = 1.0;
  double measurement_height_scale = 0.5;
  double snow_emissivity = 0.98;
  double soil_emissivity = 0.96;
  double default_soil_moisture_avail = 0.6;
  double default_soil_moisture_total = 80.0;
  double default_soil_moisture_vol = 0.3;
  double default_soil_moisture_liq = 0.3;
  double default_soil_temperature = 273.15;
  double default_soil_density = 1600.0;
  double default_soil_conductivity = 0.25;
  double default_soil_heat_capacity = 1200.0;
  std::string station_id_prefix = "WRF_GRID";
};

Values& get();
void load_from_config(const SnowpackConfig* config);

}  // namespace SnowpackConstants

extern "C" {
double snowpack_get_t_crazy_min();
double snowpack_get_t_crazy_max();
double snowpack_get_precip_phase_threshold();
double snowpack_get_snow_density_fallback();
double snowpack_get_snow_roughness_length();
double snowpack_get_bare_roughness_length();
double snowpack_get_snow_depth_roughness_threshold();
double snowpack_get_measurement_height_floor();
double snowpack_get_measurement_height_scale();
double snowpack_get_snow_emissivity();
double snowpack_get_soil_emissivity();
double snowpack_get_default_soil_moisture_avail();
double snowpack_get_default_soil_moisture_total();
double snowpack_get_default_soil_moisture_vol();
double snowpack_get_default_soil_moisture_liq();
double snowpack_get_default_soil_temperature();
double snowpack_get_default_soil_density();
double snowpack_get_default_soil_conductivity();
double snowpack_get_default_soil_heat_capacity();
}

#endif  // SNOWPACK_CONSTANTS_H
