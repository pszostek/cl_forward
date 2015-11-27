#include <cstdint>
#include <cfloat>
#include <cmath>
#include <cassert>
#include <utility>
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


std::pair<float, float> Event::findH2Boundaries(Hit h0, unsigned int cur_sensor, unsigned int second_sensor) {
        float xmin_h2, xmax_h2;
        const int z_s0 = sensor_Zs[cur_sensor + 2];
        const int z_s2 = sensor_Zs[second_sensor];

        // Note: Here, we take h0 as if it were h1, the rest
        // of the notation is fine.

        // Min and max possible x0s
        const float h_dist = std::abs(h0.z - z_s0);
        const float dxmax = PARAM_MAXXSLOPE_CANDIDATES * h_dist;
        const float x0_min = h0.x - dxmax;
        const float x0_max = h0.x + dxmax;

        // Min and max possible h1s for that h0
        float z2_tz = (((float) z_s2 - z_s0)) / (h0.z - z_s0);
        const float x_min = x0_max + (h0.x - x0_max) * z2_tz;
        xmin_h2 = x_min - PARAM_TOLERANCE_CANDIDATES;

        const float x_max = x0_min + (h0.x - x0_min) * z2_tz;
        xmax_h2 = x_max + PARAM_TOLERANCE_CANDIDATES;
        return std::make_pair(xmin_h2, xmax_h2);
}

float Event::fitHitToTrack(const float tx, const float ty,
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

// this guy should return <hit_candidates, hit_h2_candidates>
void Event::fillCandidates(CandidatesMap& hit_candidates,
        CandidatesMap& hit_h2_candidates) {
    //const int blockDim_product = get_local_size(0) * get_local_size(1);
    for(unsigned int cur_sensor = number_of_sensors - 1; cur_sensor >= 2; --cur_sensor) {
        const int second_sensor = cur_sensor - 2;

        const bool process_h1_candidates = cur_sensor >= 4;
        const bool process_h2_candidates = cur_sensor <= number_of_sensors - 3;

        // Sensor dependent calculations
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
            std::tie(xmin_h2, xmax_h2) = findH2Boundaries(h0, cur_sensor, second_sensor);

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

                            hit_candidates[h0_index].first = h1_index;
                            first_h1_found = true;
                        }
                        // The last one, only if the first one has already been found
                        else if (first_h1_found && !tol_condition) {
                            ASSERT(2 * h0_index + 1 < 2 * (sensor_hit.starts[number_of_sensors-1] + sensor_hits.nums[number_of_sensors-1]))

                            hit_candidates[h0_index].second = h1_index;
                            last_h1_found = true;
                        }
                    }

                    if (process_h2_candidates && !last_h2_found) {
                        if (!first_h2_found && h1.x > xmin_h2) {
                            ASSERT(2 * h0_index < 2 * (sensor_hits.starts[number_of_sensors-1] + sensor_hits.nums[number_of_sensors-1]))

                            hit_h2_candidates[h0_index].first = h1_index;
                            first_h2_found = true;
                        }
                        else if (first_h2_found && h1.x > xmax_h2) {
                            ASSERT(2 * h0_index + 1 < 2 * (sensor_hits.starts[number_of_sensors-1] + sensor_hits.nums[number_of_sensors-1]))

                            hit_h2_candidates[h0_index].second = h1_index;
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

                    hit_candidates[h0_index].second = hitstarts_s2 + hitnums_s2;
                }

                if (process_h2_candidates && first_h2_found && !last_h2_found) {
                    ASSERT(2 * h0_index + 1 < 2 * (sensor_hits.starts[number_of_sensors-1] + sensor_hits.nums[number_of_sensors-1]))

                    hit_h2_candidates[h0_index].second = hitstarts_s2 + hitnums_s2;
                }
            }
        }
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
void Event::trackForwarding(bool* const hit_used,
        int* const sensor_data,
        std::vector<int>& tracks_to_follow, std::vector<int>& weak_tracks,
        const unsigned int prev_ttf, std::vector<Track>& tracklets,
        std::vector<struct Track>& tracks) {

    std::vector<int> new_tracks_to_follow;
    for (auto it = tracks_to_follow.begin() + prev_ttf; it != tracks_to_follow.end(); ++it) {

        // These variables need to go here, shared memory and scope requirements
        float tx, ty, h1_z;
        unsigned int trackno, skipped_modules;
        unsigned int best_hit_h2 = 0;
        struct Track t;
        struct Hit h0;
        // The logic is broken in two parts for shared memory loading
        //DEBUG << "ttf_el: " << ttf_element << " diff_ttf " << diff_ttf << std::endl;

        // OA: tracks_to_follow is limited to TTF_MODULO elements.
        unsigned fulltrackno = *it;
        // OA: fulltrackno has not only the index into tracks_pointer encoded, but also:
        // bit 31: track_flag is set in the trackCreation function
        // bits 30-28: encodes skipped_modules
        const bool track_flag = (fulltrackno & 0x80000000) == 0x80000000;
        skipped_modules = (fulltrackno & 0x70000000) >> 28;
        trackno = fulltrackno & 0x0FFFFFFF;

        if (track_flag) {
            t = tracklets[trackno];
            ASSERT(trackno < number_of_hits);
        } else {
            t = tracks[trackno];
            ASSERT(trackno < MAX_TRACKS);
        }

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
        float best_fit = MAX_FLOAT;

        for (int k=0; k<sensor_data[SENSOR_DATA_HITNUMS + 2]; ++k) {
            const int h2_index = sensor_data[2] + k;
            struct Hit h2;
            h2.x = hits.Xs[h2_index];
            h2.y = hits.Ys[h2_index];
            h2.z = hits.Zs[h2_index];

            const float fit = fitHitToTrack(tx, ty, &h0, h1_z, &h2);

            if (fit < best_fit) {
                best_fit = fit;
                best_hit_h2 = h2_index;
            }
        }


        // We have a best fit!
        // Fill in t, ONLY in case the best fit is acceptable
        if (best_fit != MAX_FLOAT) {
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
                // we have to allocate it in the tracks vector
                // XXX OA: no more atomic_add needed
                trackno = tracks.size();
                // XXX PS: without this spell the whole thing doesnt work...
                tracks.push_back(Track());
            }

            // Copy the track into tracks
            ASSERT(trackno < number_of_hits)
            tracks[trackno] = t;

            // Add the tracks to the bag of tracks to_follow
            // XXX OA: no more atomic_add needed
            //const unsigned int ttfP = atomic_add(ttf_insertPointer, 1) % TTF_MODULO;
            new_tracks_to_follow.push_back(trackno);
        }
        // A track just skipped a module
        // We keep it for another round
        else if (skipped_modules <= MAX_SKIPPED_MODULES) {
            // Form the new mask
            trackno = ((skipped_modules + 1) << 28) | (fulltrackno & 0x8FFFFFFF);

            // Add the tracks to the bag of tracks to_follow
            //const unsigned int ttfP = atomic_add(ttf_insertPointer, 1) % TTF_MODULO;
            new_tracks_to_follow.push_back(trackno);
        }
        // If there are only three hits in this track,
        // mark it as "doubtful"
        else if (t.hitsNum == 3) {
            //const unsigned int weakP = atomic_add(weaktracks_insertPointer, 1);
            weak_tracks.push_back(trackno);
            ASSERT(weak_tracks.size() < number_of_hits)
        }
        // In the "else" case, we couldn't follow up the track,
        // so we won't be track following it anymore.
    }
    tracks_to_follow.insert(tracks_to_follow.end(),
                            new_tracks_to_follow.begin(),
                            new_tracks_to_follow.end());
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

void Event::trackCreation(int* const sensor_data, CandidatesMap& hit_candidates, int h0_index,
        bool* const hit_used, CandidatesMap& hit_h2_candidates,
        std::vector<Track>& tracklets, std::vector<int>& tracks_to_follow) {

    //DEBUG << "trackCreation: " << h0_index << std::endl;
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
    //if (h0_index == 249)
    //    DEBUG << "h0_index 249\n";

    // Calculate new dymax
    const float s1_z = hits.Zs[sensor_data[1]];
    const float h_dist = std::abs(s1_z - h0.z);
    dymax = PARAM_MAXYSLOPE * h_dist;

    // Only iterate in the hits indicated by hit_candidates :)
    first_h1 = hit_candidates[h0_index].first;
    const int last_h1 = hit_candidates[h0_index].second;
    num_h1_to_process = last_h1 - first_h1;
    // Only iterate max_numhits_to_process[0] iterations (with get_local_size(1) threads) :D :D :D
    for (unsigned int h1_element=0; h1_element < num_h1_to_process; ++h1_element) {
        int h1_index = first_h1 + h1_element;
        bool is_h1_used = hit_used[h1_index];
        //if (h0_index == 249)
        //    DEBUG << "h0_249 " << h1_index << " " << is_h1_used << std::endl;
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
                    if(scatter < MAX_SCATTER && scatter < best_fit) {
                        best_fit = scatter;
                        best_hit_h1 = h1_index;
                        best_hit_h2 = h2_index;
                    }
                }
            }
        }
    }

    if (best_fit < MAX_FLOAT) {
        // We have a best fit! - haven't we?
        // Fill in track information

        // Add the track to the bag of tracks
        // TODO: fix this in a clean way
        Track new_tracklet = {3, {static_cast<unsigned int>(h0_index), best_hit_h1, best_hit_h2}};
        if (h0_index == 249)
            DEBUG << h0_index << " " << best_hit_h1 << " " << best_hit_h2 << std::endl;
        tracklets.push_back(new_tracklet);
        unsigned int new_tracklet_idx = tracklets.size()-1;

        // Add the tracks to the bag of tracks to_follow
        // Note: The first bit flag marks this is a tracklet (hitsNum == 3),
        // and hence it is stored in tracklets
        tracks_to_follow.push_back(0x80000000 | new_tracklet_idx);
        //DEBUG << "updated tracks_to_follow " << new_tracklet_idx << std::endl;
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
std::vector<Track> Event::serialSearchByTriplets() {


    // Data initialization
    // Each event is treated with two blocks, one for each side.
    // OA: We treat one event at a time
    //const int event_number = get_group_id(0);
    // OA: the notion of blocks is not needed anymore.
    //const int events_under_process = get_num_groups(0);
    //const int tracks_offset = event_number * MAX_TRACKS;
    //const int blockDim_product = get_local_size(0) * get_local_size(1);

    // Pointers to data within the event


    //struct Track* const tracks = new Track[MAX_TRACKS];
    std::vector<struct Track> tracks;
    tracks.reserve(MAX_TRACKS);
    DEBUG << "number of sensors: " << number_of_sensors << std::endl;
    DEBUG << "number of hits: " << number_of_hits << std::endl;
    // Per event datatypes
    // OA: we pass a tracks pointer to be used in here.
    //__global struct Track* tracks = dev_tracks + tracks_offset;
    // OA: this is not needed in the serial version as atomic.
    //__global unsigned int* const tracks_insertPointer = (__global unsigned int*) dev_atomicsStorage + event_number;
    // OA: instead we allocate it locally and pass a reference.

    // Per side datatypes
    //const int hit_offset = dev_hit_offsets[event_number];
    bool* hit_used = new bool[number_of_hits];
    for(int i = 0; i < number_of_hits; ++i) {
        hit_used[i] = false;
    }
    CandidatesMap hit_candidates;
    CandidatesMap hit_h2_candidates;

    //std::vector<int> tracks_to_follow(TTF_MODULO, 0);
    std::vector<int> tracks_to_follow;
    //int* tracks_to_follow = new int[TTF_MODULO];

    std::vector<int> weak_tracks;
    weak_tracks.reserve(number_of_hits);

    std::vector<Track> tracklets;
    tracklets.reserve(number_of_hits);

    // Initialize variables according to event number and sensor side
    // Insert pointers (atomics)
    // OA: I don't think those are needed anymore
    //const int ip_shift = events_under_process + event_number * NUM_ATOMICS;
    //__global int* const weaktracks_insertPointer = (__global int*) dev_atomicsStorage + ip_shift + 1;
    //__global int* const tracklets_insertPointer = (__global int*) dev_atomicsStorage + ip_shift + 2;

    // OA: This was an atomic counter to be able to append to the global array tracks_to_follow while
    // processing several envents simultaneously - this will be needed later again probably
    //int* const ttf_insertPointer = (__global int*) dev_atomicsStorage + ip_shift + 3;
    //int ttf_insertPointer = 0; // OA: instead we jsut keep a simple counter for book keeping
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
    fillCandidates(hit_candidates, hit_h2_candidates);
    /*for (int i = 0; i < 2*number_of_hits; ++i)
        DEBUG << hit_h2_candidates[i] << ", ";
    DEBUG << std::endl;*/
    // Deal with odd or even in the same thread
    // Prepare s1 and s2 for the first iteration
    unsigned int prev_ttf, last_ttf = 0;

    for (unsigned int cur_sensor = number_of_sensors-1; cur_sensor >= 4; --cur_sensor) {

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

        sensor_data[0] = sensor_hits.starts[cur_sensor];  // never used in the code
        sensor_data[1] = sensor_hits.starts[cur_sensor-2];
        sensor_data[2] = sensor_hits.starts[cur_sensor-4];
        sensor_data[3] = sensor_hits.nums[cur_sensor];  // never used in the code
        sensor_data[4] = sensor_hits.nums[cur_sensor-2];  // never used in the code
        sensor_data[5] = sensor_hits.nums[cur_sensor-4];

        // We need this barrier if we are not using shared memory for the hits.
        // Removing shmem for hits removes the barriers in trackForwarding.
        // Otherwise the three statements from before could be executed before / after updating
        // the values inside trackForwarding
        prev_ttf = last_ttf;
        last_ttf = tracks_to_follow.size();
        //DEBUG << "diff_ttf: " << diff_ttf << std::endl;
        //DEBUG << "prev_ttf: " << prev_ttf << std::endl;
        //DEBUG << "last_ttf: " << last_ttf << std::endl;
        //barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

        // 2a. Track forwarding
        trackForwarding(hit_used,
            sensor_data, tracks_to_follow, weak_tracks, prev_ttf,
            tracklets, tracks);

        // Iterate in all hits for current sensor
        // 2a. Seeding - Track creation

        // Pre-seeding
        // Get the hits we are going to iterate onto in sh_hit_process,
        // in groups of max NUMTHREADS_X

        //DEBUG << "sensor_data[0]: " << sensor_data[0] << std::endl;
        for(int h0_index = sensor_data[0];
            h0_index < sensor_data[0] + sensor_data[SENSOR_DATA_HITNUMS];
            ++h0_index) {
            if (!hit_used[h0_index]) {
                trackCreation(sensor_data,
                    hit_candidates, h0_index, hit_used, hit_h2_candidates,
                    tracklets,
                    tracks_to_follow);
            }
        }
    }

    // Process the last bunch of track_to_follows
    for (unsigned int ttf_element = last_ttf; ttf_element< tracks_to_follow.size(); ++ttf_element) {

        const int fulltrackno = tracks_to_follow[ttf_element];
        const bool track_flag = (fulltrackno & 0x80000000) == 0x80000000;
        const int trackno = fulltrackno & 0x0FFFFFFF;

        // Here we are only interested in three-hit tracks,
        // to mark them as "doubtful"
        if (track_flag) {
            weak_tracks.push_back(trackno);
            ASSERT(weak_tracks.size() < number_of_hits);
        }
    }

    // Compute the three-hit tracks left
    for (const auto& weak_track: weak_tracks) {

        // Load the tracks from the tracklets
        const struct Track t = tracklets[weak_track];

        // Store them in the tracks bag iff they
        // are made out of three unused hits
        if (!hit_used[t.hits[0]] &&
            !hit_used[t.hits[1]] &&
            !hit_used[t.hits[2]]) {
            ASSERT(trackno < MAX_TRACKS)
            tracks.push_back(t);
        }
    }

    delete[] hit_used;
    return tracks;
}
