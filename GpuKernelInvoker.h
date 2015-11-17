#ifndef GPU_KERNEL_INVOKER
#define GPU_KERNEL_INVOKER 1

#include <iostream>
#include "Definitions.h"
#include "Tools.h"
#include "Event.h"
#include "Logger.h"

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <map>
#include <stdint.h>
#include <assert.h>

#ifdef __APPLE__
    #include "OpenCL/opencl.h"
#else
    #include <CL/cl.h>
#endif

void getMaxNumberOfHits(char*& input, int& maxHits);
void printOutSensorHits(int sensorNumber, int* prevs, int* nexts);
void printOutAllSensorHits(int* prevs, int* nexts);
void printInfo(int numberOfSensors, int numberOfHits);
void clChoosePlatform(cl_device_id*& devices, cl_platform_id& platform);
const char *getErrorString (cl_int error);
template <class T>
void clInitializeValue(cl_command_queue& commandQueue, cl_mem& param, size_t size, T value);

int invokeParallelSearch(
    const int startingEvent,
    const int eventsToProcess,
    const std::vector<const std::vector<uint8_t>* > & input,
    std::vector<std::vector<uint8_t> > & output);

struct EventBeginning {
  int numberOfSensors;
  int numberOfHits;
};

#endif
