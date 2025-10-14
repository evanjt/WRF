#include "constants.h"

#include "config.h"
#include "meteoio/meteoio/IOUtils.h"
#include "snowpack/snowpack/SnowpackConfig.h"

namespace SnowpackConstants {

Values& get() {
  static Values values;
  return values;
}

void load_from_config(const SnowpackConfig* config) {
  if (config == nullptr) {
    return;
  }

  auto& values = get();

  config->getValue("T_CRAZY_MIN", "SnowpackConstants", values.t_crazy_min,
                   mio::IOUtils::nothrow);
  config->getValue("T_CRAZY_MAX", "SnowpackConstants", values.t_crazy_max,
                   mio::IOUtils::nothrow);
  config->getValue("PRECIP_PHASE_THRESHOLD", "SnowpackConstants",
                   values.precip_phase_threshold, mio::IOUtils::nothrow);
  config->getValue("SURFACE_TEMP_GUESS", "SnowpackConstants",
                   values.surface_temp_guess, mio::IOUtils::nothrow);
  config->getValue("SNOW_DENSITY_FALLBACK", "SnowpackConstants",
                   values.snow_density_fallback, mio::IOUtils::nothrow);
  config->getValue("SNOW_ROUGHNESS_LENGTH", "SnowpackConstants",
                   values.snow_roughness_length, mio::IOUtils::nothrow);
  config->getValue("BARE_ROUGHNESS_LENGTH", "SnowpackConstants",
                   values.bare_roughness_length, mio::IOUtils::nothrow);
  config->getValue("SNOW_DEPTH_ROUGHNESS_THRESHOLD", "SnowpackConstants",
                   values.snow_depth_roughness_threshold,
                   mio::IOUtils::nothrow);
  config->getValue("MEASUREMENT_HEIGHT_FLOOR", "SnowpackConstants",
                   values.measurement_height_floor, mio::IOUtils::nothrow);
  config->getValue("MEASUREMENT_HEIGHT_SCALE", "SnowpackConstants",
                   values.measurement_height_scale, mio::IOUtils::nothrow);
  config->getValue("SNOW_EMISSIVITY", "SnowpackConstants",
                   values.snow_emissivity, mio::IOUtils::nothrow);
  config->getValue("SOIL_EMISSIVITY", "SnowpackConstants",
                   values.soil_emissivity, mio::IOUtils::nothrow);
  config->getValue("DEFAULT_SOIL_MOISTURE_AVAIL", "SnowpackConstants",
                   values.default_soil_moisture_avail, mio::IOUtils::nothrow);
  config->getValue("DEFAULT_SOIL_MOISTURE_TOTAL", "SnowpackConstants",
                   values.default_soil_moisture_total, mio::IOUtils::nothrow);
  config->getValue("DEFAULT_SOIL_MOISTURE_VOL", "SnowpackConstants",
                   values.default_soil_moisture_vol, mio::IOUtils::nothrow);
  config->getValue("DEFAULT_SOIL_MOISTURE_LIQ", "SnowpackConstants",
                   values.default_soil_moisture_liq, mio::IOUtils::nothrow);
  config->getValue("DEFAULT_SOIL_TEMPERATURE", "SnowpackConstants",
                   values.default_soil_temperature, mio::IOUtils::nothrow);
  config->getValue("DEFAULT_SOIL_DENSITY", "SnowpackConstants",
                   values.default_soil_density, mio::IOUtils::nothrow);
  config->getValue("DEFAULT_SOIL_CONDUCTIVITY", "SnowpackConstants",
                   values.default_soil_conductivity, mio::IOUtils::nothrow);
  config->getValue("DEFAULT_SOIL_HEAT_CAPACITY", "SnowpackConstants",
                   values.default_soil_heat_capacity, mio::IOUtils::nothrow);
  config->getValue("STATION_ID_PREFIX", "SnowpackConstants",
                   values.station_id_prefix, mio::IOUtils::nothrow);
}

}  // namespace SnowpackConstants

extern "C" {

double snowpack_get_t_crazy_min() {
  return SnowpackConstants::get().t_crazy_min;
}
double snowpack_get_t_crazy_max() {
  return SnowpackConstants::get().t_crazy_max;
}
double snowpack_get_precip_phase_threshold() {
  return SnowpackConstants::get().precip_phase_threshold;
}
double snowpack_get_snow_density_fallback() {
  return SnowpackConstants::get().snow_density_fallback;
}
double snowpack_get_snow_roughness_length() {
  return SnowpackConstants::get().snow_roughness_length;
}
double snowpack_get_bare_roughness_length() {
  return SnowpackConstants::get().bare_roughness_length;
}
double snowpack_get_snow_depth_roughness_threshold() {
  return SnowpackConstants::get().snow_depth_roughness_threshold;
}
double snowpack_get_measurement_height_floor() {
  return SnowpackConstants::get().measurement_height_floor;
}
double snowpack_get_measurement_height_scale() {
  return SnowpackConstants::get().measurement_height_scale;
}
double snowpack_get_snow_emissivity() {
  return SnowpackConstants::get().snow_emissivity;
}
double snowpack_get_soil_emissivity() {
  return SnowpackConstants::get().soil_emissivity;
}
double snowpack_get_default_soil_moisture_avail() {
  return SnowpackConstants::get().default_soil_moisture_avail;
}
double snowpack_get_default_soil_moisture_total() {
  return SnowpackConstants::get().default_soil_moisture_total;
}
double snowpack_get_default_soil_moisture_vol() {
  return SnowpackConstants::get().default_soil_moisture_vol;
}
double snowpack_get_default_soil_moisture_liq() {
  return SnowpackConstants::get().default_soil_moisture_liq;
}
double snowpack_get_default_soil_temperature() {
  return SnowpackConstants::get().default_soil_temperature;
}
double snowpack_get_default_soil_density() {
  return SnowpackConstants::get().default_soil_density;
}
double snowpack_get_default_soil_conductivity() {
  return SnowpackConstants::get().default_soil_conductivity;
}
double snowpack_get_default_soil_heat_capacity() {
  return SnowpackConstants::get().default_soil_heat_capacity;
}

}  // extern "C"
