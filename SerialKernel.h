#ifndef SERIAL_KERNEL
#define SERIAL_KERNEL 1

float fitHitToTrack(const float tx, const float ty,
        const struct Hit* h0, const float h1_z, const struct Hit* h2);



void fillCandidates(int* const hit_candidates,
        int* const hit_h2_candidates, const int number_of_sensors,
        const int* const sensor_hitStarts, const int* const sensor_hitNums,
        const float* const hit_Xs, const float* const hit_Ys,
        const float* const hit_Zs, const int* sensor_Zs);



void trackForwarding(const float* const hit_Xs,
        const float* const hit_Ys, const float* const hit_Zs,
        bool* const hit_used, int& tracks_insertPointer,
        int& ttf_insertPointer, int& weaktracks_insertPointer,
        int* const sensor_data, const unsigned int diff_ttf,
        int* const tracks_to_follow, int* const weak_tracks,
        const unsigned int prev_ttf, struct Track* const tracklets,
        struct Track* const tracks, const int number_of_hits);



void trackCreation(const float* const hit_Xs,
        const float* const hit_Ys, const float* const hit_Zs,
        int* const sensor_data, int* const hit_candidates, int h0_index,
        bool* const hit_used, int* const hit_h2_candidates,
        int& tracklets_insertPointer, int&  ttf_insertPointer,
        struct Track* const tracklets, int* const tracks_to_follow);



void serialSearchByTriplets(struct Track* const tracks, const uint8_t* input); 



#endif // SERIAL_KERNEL
