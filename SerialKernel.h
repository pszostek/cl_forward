#ifndef SERIAL_KERNEL
#define SERIAL_KERNEL 1

float fitHitToTrack(const float tx, const float ty,
        const struct Hit* h0, const float h1_z, const struct Hit* h2);



void fillCandidates(int* const hit_candidates,
        int* const hit_h2_candidates, const int number_of_sensors,
        const SensorHits& sensors,
        const Hits& hits, const int* sensor_Zs);


void trackForwarding(const Hits& hits,
        bool* const hit_used, int& tracks_insertPointer,
        int& ttf_insertPointer, int& weaktracks_insertPointer,
        int* const sensor_data, const unsigned int diff_ttf,
        int* const tracks_to_follow, int* const weak_tracks,
        const unsigned int prev_ttf, struct Track* const tracklets,
        struct Track* const tracks, const int number_of_hits);



void trackCreation(const Hits& hits,
        int* const sensor_data, int* const hit_candidates, int h0_index,
        bool* const hit_used, int* const hit_h2_candidates,
        int& tracklets_insertPointer, int&  ttf_insertPointer,
        struct Track* const tracklets, int* const tracks_to_follow);



int serialSearchByTriplets(struct Track* const tracks, const uint8_t* input, size_t size); 



#endif // SERIAL_KERNEL
