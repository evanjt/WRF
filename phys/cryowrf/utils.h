#include <cmath>
#include <cstdio>
#include <algorithm>

#include "snowpack_bridge_structs.h"
#include "meteoio/meteoio/MeteoIO.h"

// Forward declarations to avoid circular dependency
namespace mio {
  class Date;
  class Config;
}

struct SnowStation;
struct SnowpackConfig;
class SnowpackIO;
class SnowpackBridge;
