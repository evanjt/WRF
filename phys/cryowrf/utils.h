#include <algorithm>
#include <cmath>
#include <cstdio>

#include "meteoio/meteoio/MeteoIO.h"
#include "structs.h"

// Forward declarations to avoid circular dependency
namespace mio {
class Date;
class Config;
}  // namespace mio

struct SnowStation;
struct SnowpackConfig;
class SnowpackIO;
class SnowpackBridge;
