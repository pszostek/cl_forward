#ifndef _DATA_FRAME_H_
#define _DATA_FRAME_H_

#include <vector>
#include <cstddef>
#include <cstdint>
#include "KernelDefinitions.h"


class Event {
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

    void fillCandidates(int* const hit_candidates,
        int* const hit_h2_candidates);

    void trackCreation(int* const sensor_data, int* const hit_candidates, int h0_index,
        bool* const hit_used, int* const hit_h2_candidates,
        std::vector<Track>& tracklets, std::vector<int>& tracks_to_follow);

    void trackForwarding(bool* const hit_used, int& tracks_insertPointer,
        int* const sensor_data, const unsigned int diff_ttf,
        std::vector<int>& tracks_to_follow, std::vector<int>& weak_tracks,
        const unsigned int prev_ttf, std::vector<Track>& tracklets,
        struct Track* const tracks);

    static float fitHitToTrack(const float tx, const float ty,
        const struct Hit* h0, const float h1_z, const struct Hit* h2);
//public:
    Event(uint8_t* input, size_t size);
    virtual ~Event();


    std::vector<Track> serialSearchByTriplets();

};


#endif
