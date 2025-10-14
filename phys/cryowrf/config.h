#ifndef CONFIG_H
#define CONFIG_H

#include <string>

#include "meteoio/meteoio/MeteoIO.h"

// Forward declarations to minimize dependencies
namespace mio {
class Config;
}

struct SnowpackConfig;
class SnowpackIO;
class Snowpack;

class SnowpackConfigManager {
 public:
  static mio::Config loadConfiguration(const std::string& ini_file_path);
  static void validateConfiguration(const mio::Config& cfg);
  static std::string getDefaultConfigPath();
};

#endif  // CONFIG_H
