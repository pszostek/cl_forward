#ifndef _DATA_FRAME_H_
#define _DATA_FRAME_H_

#include <vector>
#include <cstddef>
#include <cstdint>
#include "KernelDefinitions.h"


class DataFrame {
//private:
public:
    int number_of_sensors;
    int number_of_hits;
    int *sensor_Zs;
    unsigned int *hit_IDs;
    SensorHits sensor_hits;
    Hits hits;

    void set_h_pointers_from_input(uint8_t * input, size_t size);
    static void quicksort (float* a, float* b, float* c, unsigned int* d, int start, int end);
    template<typename T> static void swap (T& a, T& b);
    static int divide (float* a, float* b, float* c, unsigned int* d, int start, int end);
//public:
    DataFrame(uint8_t* input, size_t size);
    virtual ~DataFrame();

};


#endif
