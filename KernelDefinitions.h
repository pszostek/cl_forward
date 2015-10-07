
#ifndef KERNEL_DEFINITIONS_H
#define KERNEL_DEFINITIONS_H
// Used to prefer a device type over another one
#define DEVICE_CPU 0
#define DEVICE_GPU 1
#define DEVICE_ACCELERATOR 2
#define DEVICE_PREFERENCE DEVICE_GPU
#define DEVICE_NUMBER 0

const unsigned int NUMTHREADS_X = 64;
const unsigned int MAX_NUMTHREADS_Y = 16;
const unsigned int NUM_ATOMICS = 5;
#define USE_SHARED_FOR_HITS false
const unsigned int SH_HIT_MULT = 2;

const unsigned int MAX_TRACKS = 3000;
const unsigned int MAX_TRACK_SIZE = 24;

const double REQUIRED_UNIQUES = 0.6f;
const unsigned int MIN_HITS_TRACK = 3;
const float MAX_FLOAT = FLT_MAX;
const float MIN_FLOAT = -FLT_MAX;
const unsigned int MAX_SKIPPED_MODULES = 3;
const unsigned int TTF_MODULO = 2000;

const double PARAM_W = 3966.94f; // 0.050 / sqrt( 12. )
const double PARAM_MAXXSLOPE = 0.4f;
const double PARAM_MAXYSLOPE = 0.3f;
const double PARAM_MAXXSLOPE_CANDIDATES = 0.4f;

const double PARAM_TOLERANCE = 0.6f;
const double PARAM_TOLERANCE_CANDIDATES = 0.6f;

const double MAX_SCATTER = 0.000016f;
const unsigned int SENSOR_DATA_HITNUMS = 3;
#define RESULTS_FOLDER "results"

#define PRINT_SOLUTION true
#define PRINT_VERBOSE true
#define ASSERTS_ENABLED true

#if ASSERTS_ENABLED == true
#include <cassert>
#define ASSERT(EXPR) assert(EXPR);
#else
#define ASSERT(EXPR)
#endif

struct Sensor {
    unsigned int hitStart;
    unsigned int hitNums;
};

struct Hit {
    float x;
    float y;
    float z;
};

struct Track { // 4 + 24 * 4 = 100 B
    unsigned int hitsNum;
    unsigned int hits[MAX_TRACK_SIZE];
};

#endif