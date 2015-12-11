
/**
 * Tools.h
 */

#ifndef TOOLS
#define TOOLS 1

#include <cstring>
#include <iostream>
#include <vector>
#include <sstream>
#include <map>
#include <cmath>
#include <stdint.h>
#include <stdexcept>

#ifdef __APPLE__
    #include "OpenCL/opencl.h"
#else
    #include <CL/cl.h>
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include "Logger.h"
#include "Event.h"
#include "Definitions.h"

#define clCheck(stmt) { \
  cl_int status = stmt; \
  if (status != CL_SUCCESS) { \
    std::cerr << "Error in function " << #stmt << std::endl; \
    std::cerr << "Error string: " << getErrorString(status) << std::endl; \
    exit(-1); \
  } \
}

template <class T>
std::string toString(T t){
    std::stringstream ss;
    std::string s;
    ss << t;
    ss >> s;
    return s;
}

enum class OutType {Binary, Text};

int convertClToString(const char *filename, std::string& s);
void setHPointersFromInput(uint8_t * input, size_t size,
  int& h_no_sensors, int& h_no_hits, int*& h_sensors_Zs, SensorHits& sensor_hits,
  unsigned int*& h_hit_IDs, Hits& hits);
void preorder_by_x(std::vector<const std::vector<uint8_t>* > & input);

// A non-efficient implementation that does what I need
void quicksort (float* a, float* b, float* c, unsigned int* d, int start, int end);
int divide (float* a, float* b, float* c, unsigned int* d, int first, int last);
template<typename T> void swap (T& a, T& b);
void printTrack(const Track& track, const Event& event, std::ofstream& outstream);
void writeTextTracks(const std::vector<Track>& tracks,
    const Event& event, std::ofstream& os);
void writeBinTracks(const std::vector<Track>& tracks, const Event& event, std::ofstream& os);
int findClosestModule(const int z, const std::map<int, int>& zhit_to_module);
std::map<std::string, float> calcResults(std::vector<float>& times);
void checkClError (const cl_int errcode_ret);




#endif
