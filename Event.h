#ifndef _DATA_FRAME_H_
#define _DATA_FRAME_H_

#include <vector>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <unordered_map>
#include <string>
#include "KernelDefinitions.h"

using CandidatesMap = std::unordered_map<int, std::pair<int, int>>;

class Event {
public:
    int number_of_sensors;
    int number_of_hits;
    int *sensor_Zs;
    unsigned int *hit_IDs;
    SensorHits sensor_hits;
    Hits hits;
    std::string filename;

    void set_h_pointers_from_input(uint8_t * input, size_t size);

    static void quicksort (float* a, float* b, float* c, unsigned int* d, int start, int end);

    template<typename T> static void swap (T& a, T& b);
    static int divide (float* a, float* b, float* c, unsigned int* d, int start, int end);

    Event(uint8_t* input, size_t size, std::string &fname);
    virtual ~Event();
};

void fillCandidates(const Event&, CandidatesMap& hit_candidates,
    CandidatesMap& hit_h2_candidates);

void trackCreation(const Event&, size_t cur_sensor, CandidatesMap& hit_candidates, int h0_index,
    std::vector<bool>& hit_used, CandidatesMap& hit_h2_candidates,
    std::vector<Track>& tracklets, std::vector<int>& tracks_to_follow);

void trackForwarding(const Event&, std::vector<bool>& hit_used,
    size_t cur_sensor,
    std::vector<int>& tracks_to_follow, std::vector<int>& weak_tracks,
    const unsigned int prev_ttf, std::vector<Track>& tracklets,
    std::vector<struct Track>& tracks);

static float fitHitToTrack(const float tx, const float ty,
    const struct Hit* h0, const float h1_z, const struct Hit* h2);

std::pair<float, float> findH2Boundaries(const Event& event, Hit h0, unsigned int cur_sensor, unsigned int second_sensor);

std::tuple<int, int, float> findBestFit(const Event& event, const Hit& h0,
                                        std::vector<bool>& hit_used, size_t cur_sensor,
                                        int first_h1, int last_h1);


std::vector<Track> serialSearchByTriplets(const Event&);

// OMP counterparts


void OMPFillCandidates(const Event&, std::pair<int, int> hit_candidates[],
    std::pair<int, int> hit_h2_candidates[]);

void OMPTrackCreation(const Event&, size_t cur_sensor, std::pair<int, int> hit_candidates[], int h0_index,
    std::vector<bool>& hit_used,std::pair<int, int> hit_h2_candidates[],
    std::vector<Track>& tracklets, std::vector<int>& tracks_to_follow);

void OMPTrackForwarding(const Event&, std::vector<bool>& hit_used,
    size_t cur_sensor,
    std::vector<int>& tracks_to_follow, std::vector<int>& weak_tracks,
    const unsigned int prev_ttf, std::vector<Track>& tracklets,
    std::vector<struct Track>& tracks);

static float OMPFitHitToTrack(const float tx, const float ty,
    const struct Hit* h0, const float h1_z, const struct Hit* h2);


std::pair<float, float> OMPFindH2Boundaries(const Event& event, Hit h0, unsigned int cur_sensor, unsigned int second_sensor);

std::tuple<int, int, float> OMPFindBestFit(const Event& event, const Hit& h0,
                                        const std::vector<bool>& hit_used, const size_t cur_sensor,
                                        int first_h1, int last_h1);

std::vector<Track> OMPSearchByTriplets(const Event&);


#endif
