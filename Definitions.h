#ifndef DEFINITIONS_CUH
#define DEFINITIONS_CUH 1

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <vector>

#include "KernelDefinitions.h"

using SensorHits = struct SensorHits;
using Hits = struct Hits;

const size_t EVENT_LEVEL_PARALLELISM = 4;
const size_t SENSOR_LEVEL_PARALLELISM = 4;
const size_t HIT_LEVEL_PARALLELISM = 4;

#endif
