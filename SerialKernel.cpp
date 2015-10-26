#include <cstdint>
#include <cfloat>
#include <cmath>
#include <cassert>
#include "Logger.h"
#include "KernelDefinitions.h"
#include "SerialKernel.h"

/**
* @brief Fits hits to tracks.
* @details In case the tolerances constraints are met,
*          returns the chi2 weight of the track. Otherwise,
*          returns MAX_FLOAT.
*
* @param tx
* @param ty
* @param h0
* @param h1_z
* @param h2
* @return
*/

float fitHitToTrack(const float tx, const float ty,
        const struct Hit* h0, const float h1_z, const struct Hit* h2) {
    // tolerances
    const float dz = h2->z - h0->z;
    const float x_prediction = h0->x + tx * dz;
    const float dx = std::abs(x_prediction - h2->x);
    const bool tolx_condition = dx < PARAM_TOLERANCE;

    const float y_prediction = h0->y + ty * dz;
    const float dy = std::abs(y_prediction - h2->y);
    const bool toly_condition = dy < PARAM_TOLERANCE;

    // Scatter - Updated to last PrPixel
    const float scatterNum = (dx * dx) + (dy * dy);
    const float scatterDenom = 1.f / (h2->z - h1_z);
    const float scatter = scatterNum * scatterDenom * scatterDenom;

    const bool scatter_condition = scatter < MAX_SCATTER;
    const bool condition = tolx_condition && toly_condition && scatter_condition;

    return condition * scatter + !condition * MAX_FLOAT;
}

/**
* @brief Fills dev_hit_candidates.
*
* @param hit_candidates
* @param hit_h2_candidates
* @param number_of_sensors
* @param sensor_hitStarts
* @param sensor_hitNums
* @param hit_Xs
* @param hit_Ys
* @param hit_Zs
* @param sensor_Zs
*/

void fillCandidates(int* const hit_candidates,
        int* const hit_h2_candidates, const int number_of_sensors,
        const SensorHits& sensor_hits,
        const Hits& hits, const int* sensor_Zs) {

    //const int blockDim_product = get_local_size(0) * get_local_size(1);
    int cur_sensor = number_of_sensors - 1;
    while (cur_sensor >= 2) {
        const int second_sensor = cur_sensor - 2;

        const bool process_h1_candidates = cur_sensor >= 4;
        const bool process_h2_candidates = cur_sensor <= number_of_sensors - 3;

        // Sensor dependent calculations
        const int z_s0 = process_h2_candidates ? sensor_Zs[cur_sensor + 2] : 0;
        const int z_s2 = process_h2_candidates ? sensor_Zs[second_sensor] : 0;

        // Iterate in all hits in z0
        for (int h0_element=0; h0_element < sensor_hits.nums[cur_sensor]; ++h0_element) {
            assert(h0_element < sensor_hits.nums[cur_sensor]);
            const int h0_index = sensor_hits.starts[cur_sensor] + h0_element;
            struct Hit h0;
            h0.x = hits.Xs[h0_index];
            h0.z = hits.Zs[h0_index];
            const int hitstarts_s2 = sensor_hits.starts[second_sensor];
            const int hitnums_s2 = sensor_hits.nums[second_sensor];

            float xmin_h2, xmax_h2;
            if (process_h2_candidates) {
                // Note: Here, we take h0 as if it were h1, the rest
                // of the notation is fine.

                // Min and max possible x0s
                const float h_dist = std::abs(h0.z - z_s0);
                const float dxmax = PARAM_MAXXSLOPE_CANDIDATES * h_dist;
                const float x0_min = h0.x - dxmax;
                const float x0_max = h0.x + dxmax;

                // Min and max possible h1s for that h0
                float z2_tz = (((float) z_s2 - z_s0)) / (h0.z - z_s0);
                float x = x0_max + (h0.x - x0_max) * z2_tz;
                xmin_h2 = x - PARAM_TOLERANCE_CANDIDATES;

                x = x0_min + (h0.x - x0_min) * z2_tz;
                xmax_h2 = x + PARAM_TOLERANCE_CANDIDATES;
            }

            if (cur_sensor >= 4) {
                bool first_h1_found = false, last_h1_found = false;
                bool first_h2_found = false, last_h2_found = false;

                // Iterate in all hits in z1
                for (int h1_element=0; h1_element<hitnums_s2; ++h1_element) {
                    int h1_index = hitstarts_s2 + h1_element;
                    struct Hit h1;
                    h1.x = hits.Xs[h1_index];
                    h1.z = hits.Zs[h1_index];

                    if (process_h1_candidates && !last_h1_found) {
                        // Check if h0 and h1 are compatible
                        const float h_dist = std::abs(h1.z - h0.z);
                        const float dxmax = PARAM_MAXXSLOPE_CANDIDATES * h_dist;
                        const bool tol_condition = std::abs(h1.x - h0.x) < dxmax;

                        // Find the first one
                        if (!first_h1_found && tol_condition) {
                            ASSERT(2 * h0_index < 2 * (sensor_hits.starts[number_of_sensors-1] + sensor_hits.nums[number_of_sensors-1]))

                            hit_candidates[2 * h0_index] = h1_index;
                            first_h1_found = true;
                        }
                        // The last one, only if the first one has already been found
                        else if (first_h1_found && !tol_condition) {
                            ASSERT(2 * h0_index + 1 < 2 * (sensor_hit.starts[number_of_sensors-1] + sensor_hits.nums[number_of_sensors-1]))

                            hit_candidates[2 * h0_index + 1] = h1_index;
                            last_h1_found = true;
                        }
                    }

                    if (process_h2_candidates && !last_h2_found) {
                        if (!first_h2_found && h1.x > xmin_h2) {
                            ASSERT(2 * h0_index < 2 * (sensor_hits.starts[number_of_sensors-1] + sensor_hits.nums[number_of_sensors-1]))

                            hit_h2_candidates[2 * h0_index] = h1_index;
                            first_h2_found = true;
                        }
                        else if (first_h2_found && h1.x > xmax_h2) {
                            ASSERT(2 * h0_index + 1 < 2 * (sensor_hits.starts[number_of_sensors-1] + sensor_hits.nums[number_of_sensors-1]))

                            hit_h2_candidates[2 * h0_index + 1] = h1_index;
                            last_h2_found = true;
                        }
                    }

                    if ((!process_h1_candidates || last_h1_found) &&
                        (!process_h2_candidates || last_h2_found)) {
                        break;
                    }
                }

                // Note: If first is not found, then both should be -1
                // and there wouldn't be any iteration
                if (process_h1_candidates && first_h1_found && !last_h1_found) {
                    ASSERT(2 * h0_index + 1 < 2 * (sensor_hits.starts[number_of_sensors-1] + sensor_hits.nums[number_of_sensors-1]))

                    hit_candidates[2 * h0_index + 1] = hitstarts_s2 + hitnums_s2;
                }

                if (process_h2_candidates && first_h2_found && !last_h2_found) {
                    ASSERT(2 * h0_index + 1 < 2 * (sensor_hits.starts[number_of_sensors-1] + sensor_hits.nums[number_of_sensors-1]))

                    hit_h2_candidates[2 * h0_index + 1] = hitstarts_s2 + hitnums_s2;
                }
            }
        }

        --cur_sensor;
    }
}

/**
* @brief Performs the track forwarding.
*
* @param hit_Xs
* @param hit_Ys
* @param hit_Zs
* @param hit_used
* @param sensor_data
* @param diff_ttf
* @param tracks_to_follow
* @param weak_tracks
* @param prev_ttf
* @param tracklets
* @param tracks
* @param number_of_hits
*/
void trackForwarding(const Hits& hits,
        bool* const hit_used, int& tracks_insertPointer,
        int& ttf_insertPointer, int& weaktracks_insertPointer,
        int* const sensor_data, const unsigned int diff_ttf,
        int* const tracks_to_follow, int* const weak_tracks,
        const unsigned int prev_ttf, struct Track* const tracklets,
        struct Track* const tracks, const int number_of_hits) {

    for (unsigned int ttf_element=0; ttf_element<diff_ttf; ++ttf_element) {

        // These variables need to go here, shared memory and scope requirements
        float tx, ty, h1_z;
        unsigned int trackno, fulltrackno, skipped_modules;
        unsigned int best_hit_h2 = 0;
        struct Track t;
        struct Hit h0;
        // The logic is broken in two parts for shared memory loading
        DEBUG << "ttf_el: " << ttf_element << " diff_ttf " << diff_ttf << std::endl;

        // OA: tracks_to_follow is limited to TTF_MODULO elements.
        fulltrackno = tracks_to_follow[(prev_ttf + ttf_element) % TTF_MODULO];
        // OA: fulltrackno has not only the index into tracks_pointer encoded, but also:
        // bit 31: track_flag is set in the trackCreation function
        // bits 30-28: encodes skipped_modules
        const bool track_flag = (fulltrackno & 0x80000000) == 0x80000000;
        skipped_modules = (fulltrackno & 0x70000000) >> 28;
        trackno = fulltrackno & 0x0FFFFFFF;

        const struct Track* const track_pointer = track_flag ? tracklets : tracks;

        ASSERT(track_pointer==tracklets ? trackno < number_of_hits : true)
        ASSERT(track_pointer==tracks ? trackno < MAX_TRACKS : true)
        t = track_pointer[trackno];

        // Load last two hits in h0, h1
        const int t_hitsNum = t.hitsNum;
        ASSERT(t_hitsNum < MAX_TRACK_SIZE)
        const int h0_num = t.hits[t_hitsNum - 2];
        const int h1_num = t.hits[t_hitsNum - 1];

        ASSERT(h0_num < number_of_hits)
        h0.x = hits.Xs[h0_num];
        h0.y = hits.Ys[h0_num];
        h0.z = hits.Zs[h0_num];

        ASSERT(h1_num < number_of_hits)
        const float h1_x = hits.Xs[h1_num];
        const float h1_y = hits.Ys[h1_num];
        h1_z = hits.Zs[h1_num];

        // Track forwarding over t, for all hits in the next module
        // Line calculations
        const float td = 1.0f / (h1_z - h0.z);
        const float txn = (h1_x - h0.x);
        const float tyn = (h1_y - h0.y);
        tx = txn * td;
        ty = tyn * td;

        // Search for a best fit
        // Load shared elements

        // OA: once the shared mem goes away much of the logic here becomes redundant.
        // TODO: discuss this - I'm not entirely sure, but a simple loop k=0:hitNums
        // could be enugh
        // Iterate in the third list of hits
        // Tiled memory access on h2
        // Only load for get_local_id(1) == 0
        float best_fit = FLT_MAX;

        for (int k=0; k<sensor_data[SENSOR_DATA_HITNUMS + 2]; ++k) {
            const int h2_index = sensor_data[2] + k;
            struct Hit h2;
            h2.x = hits.Xs[h2_index];
            h2.y = hits.Ys[h2_index];
            h2.z = hits.Zs[h2_index];

            const float fit = fitHitToTrack(tx, ty, &h0, h1_z, &h2);
            const bool fit_is_better = fit < best_fit;

            best_fit = fit_is_better * fit + !fit_is_better * best_fit;
            best_hit_h2 = fit_is_better * h2_index + !fit_is_better * best_hit_h2;
        }


        // We have a best fit!
        // Fill in t, ONLY in case the best fit is acceptable
        if (best_fit != FLT_MAX) {
            // Mark h2 as used
            ASSERT(best_hit_h2 < number_of_hits)
            hit_used[best_hit_h2] = true;

            // Update the tracks to follow, we'll have to follow up
            // this track on the next iteration :)
            ASSERT(t.hitsNum < MAX_TRACK_SIZE)
            t.hits[t.hitsNum++] = best_hit_h2;

            // Update the track in the bag
            if (t.hitsNum <= 4) {
                ASSERT(t.hits[0] < number_of_hits)
                ASSERT(t.hits[1] < number_of_hits)
                ASSERT(t.hits[2] < number_of_hits)

                // Also mark the first three as used
                hit_used[t.hits[0]] = true;
                hit_used[t.hits[1]] = true;
                hit_used[t.hits[2]] = true;

                // If it is a track made out of less than or equal than 4 hits,
                // we have to allocate it in the tracks pointer
                // XXX OA: no more atomic_add needed
                //trackno = atomic_add(tracks_insertPointer, 1);
                trackno = tracks_insertPointer++;
            }

            // Copy the track into tracks
            ASSERT(trackno < number_of_hits)
            tracks[trackno] = t;

            // Add the tracks to the bag of tracks to_follow
            // XXX OA: no more atomic_add needed
            //const unsigned int ttfP = atomic_add(ttf_insertPointer, 1) % TTF_MODULO;
            tracks_to_follow[ttf_insertPointer++ % TTF_MODULO] = trackno;
        }
        // A track just skipped a module
        // We keep it for another round
        else if (skipped_modules <= MAX_SKIPPED_MODULES) {
            // Form the new mask
            trackno = ((skipped_modules + 1) << 28) | (fulltrackno & 0x8FFFFFFF);

            // Add the tracks to the bag of tracks to_follow
            //const unsigned int ttfP = atomic_add(ttf_insertPointer, 1) % TTF_MODULO;
            tracks_to_follow[ttf_insertPointer++ % TTF_MODULO] = trackno;
        }
        // If there are only three hits in this track,
        // mark it as "doubtful"
        else if (t.hitsNum == 3) {
            //const unsigned int weakP = atomic_add(weaktracks_insertPointer, 1);
            weak_tracks[weaktracks_insertPointer++] = trackno;
            ASSERT(weaktracks_insertPointer < number_of_hits)
        }
        // In the "else" case, we couldn't follow up the track,
        // so we won't be track following it anymore.
    }
}
/**
* @brief Track Creation
*
* @param hit_Xs
* @param hit_Ys
* @param hit_Zs
* @param sensor_data
* @param hit_candidates
* @param max_numhits_to_process
* @param sh_hit_x
* @param sh_hit_y
* @param sh_hit_z
* @param sh_hit_process
* @param hit_used
* @param hit_h2_candidates
* @param blockDim_sh_hit
* @param best_fits
* @param tracklets_insertPointer
* @param ttf_insertPointer
* @param tracklets
* @param tracks_to_follow
*/

void trackCreation(const Hits& hits,
        int* const sensor_data, int* const hit_candidates, int h0_index,
        bool* const hit_used, int* const hit_h2_candidates,
        int& tracklets_insertPointer, int&  ttf_insertPointer,
        struct Track* const tracklets, int* const tracks_to_follow) {

    DEBUG << "trackCreation: " << h0_index << std::endl;
    // Track creation starts
    unsigned int best_hit_h1 = 0;
    unsigned int best_hit_h2 = 0;
    struct Hit h0, h1;
    int first_h1;
    // PS: Removed, since never used
    // int last_h2;
    float dymax;

    unsigned int num_h1_to_process = 0;
    float best_fit = MAX_FLOAT;

    // We will repeat this for performance reasons
    h0.x = hits.Xs[h0_index];
    h0.y = hits.Ys[h0_index];
    h0.z = hits.Zs[h0_index];

    // Calculate new dymax
    const float s1_z = hits.Zs[sensor_data[1]];
    const float h_dist = std::abs(s1_z - h0.z);
    dymax = PARAM_MAXYSLOPE * h_dist;

    // Only iterate in the hits indicated by hit_candidates :)
    first_h1 = hit_candidates[2 * h0_index];
    const int last_h1 = hit_candidates[2 * h0_index + 1];
    num_h1_to_process = last_h1 - first_h1;
    // Only iterate max_numhits_to_process[0] iterations (with get_local_size(1) threads) :D :D :D
    for (unsigned int h1_element=0; h1_element < num_h1_to_process; ++h1_element) {
        int h1_index = first_h1 + h1_element;
        bool is_h1_used = hit_used[h1_index];
        float dz_inverted;

        if (!is_h1_used) {
            h1.x = hits.Xs[h1_index];
            h1.y = hits.Ys[h1_index];
            h1.z = hits.Zs[h1_index];

            dz_inverted = 1.f / (h1.z - h0.z);

            //int first_h2 = hit_h2_candidates[2 * h1_index];
            // PS: The following variable is removed, since it's never used
            //last_h2 = hit_h2_candidates[2 * h1_index + 1];
            //DEBUG << "first_h2 " << first_h2 << std::endl;

            // In case there be no h2 to process,
            // we can preemptively prevent further processing
            //inside_bounds &= first_h2 != -1;

            // Iterate in the third list of hits
            // Tiled memory access on h2
            for (int k=0; k<sensor_data[SENSOR_DATA_HITNUMS + 2]; ++k) {
                const int h2_index = sensor_data[2] + k;
                struct Hit h2;
                h2.x = hits.Xs[h2_index];
                h2.y = hits.Ys[h2_index];
                h2.z = hits.Zs[h2_index];

                // Predictions of x and y for this hit
                const float z2_tz = (h2.z - h0.z) * dz_inverted;
                const float x = h0.x + (h1.x - h0.x) * z2_tz;
                const float y = h0.y + (h1.y - h0.y) * z2_tz;
                const float dx = x - h2.x;
                const float dy = y - h2.y;
                if (std::abs(h1.y - h0.y) < dymax &&
                        std::abs(dx) < PARAM_TOLERANCE &&
                        std::abs(dy) < PARAM_TOLERANCE) {
                    // Calculate fit
                    const float scatterNum = (dx * dx) + (dy * dy);
                    const float scatterDenom = 1.f / (h2.z - h1.z);
                    const float scatter = scatterNum * scatterDenom * scatterDenom;
                    const bool condition = scatter < MAX_SCATTER;
                    const float fit = condition * scatter + !condition * MAX_FLOAT;

                    const bool fit_is_better = fit < best_fit;
                    best_fit = fit_is_better * fit + !fit_is_better * best_fit;
                    best_hit_h1 = fit_is_better * (h1_index) + !fit_is_better * best_hit_h1;
                    best_hit_h2 = fit_is_better * (h2_index) + !fit_is_better * best_hit_h2;
                }
            }
        }
    }

    if (best_fit < MAX_FLOAT) {
        // We have a best fit! - haven't we?
        // Fill in track information

        // Add the track to the bag of tracks
        const unsigned int trackP = tracklets_insertPointer++;
        tracklets[trackP].hitsNum = 3;
        tracklets[trackP].hits[0] = h0_index;
        tracklets[trackP].hits[1] = best_hit_h1;
        tracklets[trackP].hits[2] = best_hit_h2;

        // Add the tracks to the bag of tracks to_follow
        // Note: The first bit flag marks this is a tracklet (hitsNum == 3),
        // and hence it is stored in tracklets
        const unsigned int ttfP = ttf_insertPointer++ % TTF_MODULO;
        tracks_to_follow[ttfP] = 0x80000000 | trackP;
        DEBUG << "updated tracks_to_follow " << trackP << std::endl;
    }
}

/**
* @brief Track following algorithm, loosely based on Pr/PrPixel.
* @details It should be simplistic in its design, as is the Pixel VELO problem
*          Triplets are chosen based on a fit and forwarded using a typical track following algo.
*          Ghosts are inherently out of the equation, as the algorithm considers all possible
*          triplets and keeps the best. Upon following, if a hit is not found in the adjacent
*          module, the track[let] is considered complete.
*          Clones are removed based off a used-hit mechanism. A global array keeps track of
*          used hits when forming tracks consist of 4 or more hits.
*
*          The algorithm consists in two stages: Track following, and seeding. In each step [iteration],
*          the track following is performed first, hits are marked as used, and then the seeding is performed,
*          requiring the first two hits in the triplet to be unused.
*
* @param dev_tracks
* @param input The data of one event
* @param dev_tracks_to_follow
* @param dev_hit_used
* @param dev_atomicsStorage
* @param dev_tracklets
* @param dev_weak_tracks
* @param dev_event_offsets
* @param dev_hit_candidates
*/
int serialSearchByTriplets(struct Track* const tracks, const uint8_t* input, size_t size) {


    // Data initialization
    // Each event is treated with two blocks, one for each side.
    // OA: We treat one event at a time
    //const int event_number = get_group_id(0);
    // OA: the notion of blocks is not needed anymore.
    //const int events_under_process = get_num_groups(0);
    //const int tracks_offset = event_number * MAX_TRACKS;
    //const int blockDim_product = get_local_size(0) * get_local_size(1);

    // Pointers to data within the event


    const int number_of_sensors = *(uint32_t*)input; input += sizeof(uint32_t);
    const int number_of_hits = *(uint32_t*)input;    input += sizeof(uint32_t);
    SensorHits sensor_hits;
    const int* const sensor_Zs = (int*) input;  input += sizeof(int)*number_of_sensors;
    sensor_hits.starts =  (int*) input; input += sizeof(int)*number_of_sensors;
    sensor_hits.nums =  (int*) input; input += sizeof(int)*number_of_sensors;
    // PS: Removed since never used
    // const unsigned int* const hit_IDs =  (uint32_t*) input;
    input += sizeof(uint32_t)*number_of_hits;
    Hits hits;
    hits.Xs =  (float*) input; input += sizeof(float)*number_of_hits;
    hits.Ys =  (float*) input; input += sizeof(float)*number_of_hits;
    hits.Zs =  (float*) input;

    DEBUG << "number of sensors: " << number_of_sensors << std::endl;
    DEBUG << "number of hits: " << number_of_hits << std::endl;
    // Per event datatypes
    // OA: we pass a tracks pointer to be used in here.
    //__global struct Track* tracks = dev_tracks + tracks_offset;
    // OA: this is not needed in the serial version as atomic.
    //__global unsigned int* const tracks_insertPointer = (__global unsigned int*) dev_atomicsStorage + event_number;
    // OA: instead we allocate it locally and pass a reference.
    int tracks_insertPointer = 0;

    // Per side datatypes
    //const int hit_offset = dev_hit_offsets[event_number];
    bool* hit_used = new bool[number_of_hits];
    int* hit_candidates = new int[2 * number_of_hits];
    int* hit_h2_candidates = new int[2 * number_of_hits];
    for (int i = 0; i < 2*number_of_hits; ++i)
    {
        hit_used[i/2] = false;
        hit_candidates[i] = -1;
        hit_h2_candidates[i] = -1;
    }

    int* tracks_to_follow = new int[TTF_MODULO];
    int* weak_tracks = new int[number_of_hits];
    struct Track* const tracklets = new Track[number_of_hits];
    // TODO: replace by memset
    for (int i = 0; i < TTF_MODULO; ++i)
    {
        tracks_to_follow[i] = 0;
    }

    // Initialize variables according to event number and sensor side
    // Insert pointers (atomics)
    // OA: I don't think those are needed anymore
    //const int ip_shift = events_under_process + event_number * NUM_ATOMICS;
    //__global int* const weaktracks_insertPointer = (__global int*) dev_atomicsStorage + ip_shift + 1;
    //__global int* const tracklets_insertPointer = (__global int*) dev_atomicsStorage + ip_shift + 2;
    // OA: instead we declare local vars and pass references.
    int weaktracks_insertPointer = 0;
    int tracklets_insertPointer = 0;

    // OA: This was an atomic counter to be able to append to the global array tracks_to_follow while
    // processing several envents simultaneously - this will be needed later again probably
    //int* const ttf_insertPointer = (__global int*) dev_atomicsStorage + ip_shift + 3;
    int ttf_insertPointer = 0; // OA: instead we jsut keep a simple counter for book keeping
    //__global int* const sh_hit_lastPointer = (__global int*) dev_atomicsStorage + ip_shift + 4;
    //__global int* const max_numhits_to_process = (__global int*) dev_atomicsStorage + ip_shift + 5;

    // The fun begins
    // XXX OA: what's this for?
    // PS: Removed since never used
    // int sh_hit_process [NUMTHREADS_X];
    int sensor_data [6];

    // OA: This is used for coordinating between OpenCL threads within a group.
    // Not needed anymore
    //const int cond_sh_hit_mult = min((int) get_local_size(1), SH_HIT_MULT);
    //const int blockDim_sh_hit = NUMTHREADS_X * cond_sh_hit_mult;
    fillCandidates(hit_candidates, hit_h2_candidates, number_of_sensors,
            sensor_hits,
            hits, sensor_Zs);
    /*for (int i = 0; i < 2*number_of_hits; ++i)
        DEBUG << hit_h2_candidates[i] << ", ";
    DEBUG << std::endl;*/
    // Deal with odd or even in the same thread
    int cur_sensor = number_of_sensors - 1;

    // Prepare s1 and s2 for the first iteration
    unsigned int prev_ttf, last_ttf = 0;

    while (cur_sensor >= 4) {

        // OA: that's not needed anymore
        //barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

        // Iterate in sensors
        // Load in shared
        // XXX OA: I don't see the purpose of this. except for shared mem maybe but then
        // why is it not in an ifdef
        /*
        if (get_local_id(0) < 6 && get_local_id(1) == 0) {
            const int sensor_number = cur_sensor - (get_local_id(0) % 3) * 2;
            __global const int* const sensor_pointer = get_local_id(0) < 3 ? sensor_hitStarts : sensor_hitNums;

            sensor_data[get_local_id(0)] = sensor_pointer[sensor_number];
        }
        else if (get_local_id(0) == 7 && get_local_id(1) == 0) {
            max_numhits_to_process[0] = 0;
        }
        */
        // OA: this should do the same thing for the serial code as the above if branch:
        for (int i = 0; i < 6; ++i) {
            const int sensor_number = cur_sensor - (i % 3) * 2;
            const int* const sensor_pointer = i < 3 ? sensor_hits.starts : sensor_hits.nums;
            sensor_data[i] = sensor_pointer[sensor_number];
        }

        // We need this barrier if we are not using shared memory for the hits.
        // Removing shmem for hits removes the barriers in trackForwarding.
        // Otherwise the three statements from before could be executed before / after updating
        // the values inside trackForwarding
        prev_ttf = last_ttf;
        last_ttf = ttf_insertPointer;
        const unsigned int diff_ttf = last_ttf - prev_ttf;
        DEBUG << "diff_ttf: " << diff_ttf << std::endl;
        DEBUG << "prev_ttf: " << prev_ttf << std::endl;
        DEBUG << "last_ttf: " << last_ttf << std::endl;
        //barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

        // 2a. Track forwarding
        trackForwarding(hits, hit_used,
            tracks_insertPointer, ttf_insertPointer, weaktracks_insertPointer,
            sensor_data, diff_ttf, tracks_to_follow, weak_tracks, prev_ttf,
            tracklets, tracks, number_of_hits);

        // Iterate in all hits for current sensor
        // 2a. Seeding - Track creation

        // Pre-seeding
        // Get the hits we are going to iterate onto in sh_hit_process,
        // in groups of max NUMTHREADS_X

        DEBUG << "sensor_data[0]: " << sensor_data[0] << std::endl;
        for(int h0_index = sensor_data[0];
            h0_index < sensor_data[0] + sensor_data[SENSOR_DATA_HITNUMS];
            ++h0_index) {
            if (!hit_used[h0_index]) {
                trackCreation(hits, sensor_data,
                    hit_candidates, h0_index, hit_used, hit_h2_candidates,
                    tracklets_insertPointer, ttf_insertPointer, tracklets,
                    tracks_to_follow);
            }
        }

        cur_sensor -= 1;
    }

    prev_ttf = last_ttf;
    last_ttf = ttf_insertPointer;
    const unsigned int diff_ttf = last_ttf - prev_ttf;

    // Process the last bunch of track_to_follows
    for (unsigned int ttf_element = 0; ttf_element< diff_ttf; ++ttf_element) {

        const int fulltrackno = tracks_to_follow[(prev_ttf + ttf_element) % TTF_MODULO];
        const bool track_flag = (fulltrackno & 0x80000000) == 0x80000000;
        const int trackno = fulltrackno & 0x0FFFFFFF;

        // Here we are only interested in three-hit tracks,
        // to mark them as "doubtful"
        if (track_flag) {
            const unsigned int weakP = weaktracks_insertPointer++;
            ASSERT(weakP < number_of_hits)
            weak_tracks[weakP] = trackno;
        }
    }

    // Compute the three-hit tracks left
    const unsigned int weaktracks_total = weaktracks_insertPointer;
    for (unsigned int weaktrack_no = 0; weaktrack_no < weaktracks_total; ++weaktrack_no) {

        // Load the tracks from the tracklets
        const struct Track t = tracklets[weak_tracks[weaktrack_no]];

        // Store them in the tracks bag iff they
        // are made out of three unused hits
        if (!hit_used[t.hits[0]] &&
            !hit_used[t.hits[1]] &&
            !hit_used[t.hits[2]]) {
            const unsigned int trackno = tracks_insertPointer++;
            ASSERT(trackno < MAX_TRACKS)
            tracks[trackno] = t;
        }
    }
    delete[] hit_used;
    delete[] hit_candidates;
    delete[] hit_h2_candidates;
    delete[] tracks_to_follow;
    delete[] weak_tracks;
    delete[] tracklets;
    return tracks_insertPointer;
}
