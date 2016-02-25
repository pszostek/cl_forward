#include <cstdint>
#include <cfloat>
#include <cmath>
#include <cassert>
#include <utility>
#include "Logger.h"
#include "KernelDefinitions.h"
#include "SerialKernel.h"
#include "Definitions.h"


/* This methods takes a h0 event.hits and looks for the best h1 and h2, somehow.
    first_h1 and last_h1 specify the iteration range, i.e. the function looks in event.hits[first_h1, last_h1]
    */

std::tuple<int, int, float> OMPFindBestFit(const Event& event,
            const Hit& h0, const std::vector<bool>& hit_used,
            const size_t cur_sensor, int first_h1, int last_h1) {

    unsigned int best_hit_h1 = 0;
    unsigned int best_hit_h2 = 0;
    float best_fit = MAX_FLOAT;
    struct Hit h1;

    int* __restrict__ sensor_starts = event.sensor_hits.starts;
    int* __restrict__ sensor_nums = event.sensor_hits.nums;

    // TODO: PS: event.hits.{Xs,Ys,Zs} should be aligned to 64 bytes. They are set
    //           in Event.cpp:5
    float* __restrict__ hits_xs = event.hits.Xs; 
    float* __restrict__ hits_ys = event.hits.Ys; 
    float* __restrict__ hits_zs = event.hits.Zs; 

    // Calculate new dymax
    const float s1_z = hits_zs[sensor_starts[cur_sensor-2]];
    const float h_dist = std::abs(s1_z - h0.z);
    const float dymax = PARAM_MAXYSLOPE * h_dist;

    const unsigned HIT_LEVEL_PARALLELISM_SQRT = sqrt(HIT_LEVEL_PARALLELISM);

    #pragma omp parallel for num_threads(HIT_LEVEL_PARALLELISM_SQRT) schedule(static)
    for (unsigned int h1_index=first_h1; h1_index < last_h1; ++h1_index) {
        bool is_h1_used = hit_used[h1_index];

        if (is_h1_used)
            continue;

        h1.x = hits_xs[h1_index];
        h1.y = hits_ys[h1_index];
        h1.z = hits_zs[h1_index];

        float dz_inverted = 1.f / (h1.z - h0.z);

        // In case there be no h2 to process,
        // we can preemptively prevent further processing
        //inside_bounds &= first_h2 != -1;

        // Iterate in the third list of event.hits
        // PS: this loop *GETS VECTORIZED*
        #pragma omp parallel for simd num_threads(HIT_LEVEL_PARALLELISM_SQRT) schedule(static)
        for(size_t h2_element=0; h2_element < sensor_nums[cur_sensor-4]; ++h2_element) {
            size_t h2_index = h2_element + sensor_starts[cur_sensor-4];
            struct Hit h2;
            h2.x = hits_xs[h2_index];
            h2.y = hits_ys[h2_index];
            h2.z = hits_zs[h2_index];

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
    return std::make_tuple(best_hit_h1, best_hit_h2, best_fit);
}


// PS: this function *GETS VECTORIZED*
#pragma omp declare simd notinbranch
void OMPFindH2Boundaries(const int* __restrict__ sensor_Zs,
    const float h0_x, const float h0_z, const unsigned int cur_sensor,
    const unsigned int second_sensor, float* __restrict__ xmin_h2_ptr,
    float* __restrict__ xmax_h2_ptr) {
        const int z_s0 = sensor_Zs[cur_sensor + 2];
        const int z_s2 = sensor_Zs[second_sensor];


        // Note: Here, we take h0 as if it were h1, the rest
        // of the notation is fine.

        // Min and max possible x0s
        const float h_dist = std::abs(h0_z - z_s0);
        const float dxmax = PARAM_MAXXSLOPE_CANDIDATES * h_dist;
        const float x0_min = h0_x - dxmax;
        const float x0_max = h0_x + dxmax;

        // Min and max possible h1s for that h0
        float z2_tz = (((float) z_s2 - z_s0)) / (h0_z - z_s0);
        const float x_min = x0_max + (h0_x - x0_max) * z2_tz;
        *xmin_h2_ptr = x_min - PARAM_TOLERANCE_CANDIDATES;

        const float x_max = x0_min + (h0_x - x0_min) * z2_tz;
        *xmax_h2_ptr = x_max + PARAM_TOLERANCE_CANDIDATES;
}


/**
* @brief Fits hits to tracks.
* @details In case the tolerances constraints are met,
*          returns the chi2 weight of the track. Otherwise,
*          returns MAX_FLOAT.
*/

#pragma omp declare simd notinbranch
float OMPFitHitToTrack(const float tx, const float ty,
        const struct Hit* __restrict__ h0, const float h1_z, const struct Hit* __restrict__ h2) {

    // tolerances
    const float dz = h2->z - h0->z;
    const float x_prediction = h0->x + tx * dz;
    const float dx = std::abs(x_prediction - h2->x);

    const float y_prediction = h0->y + ty * dz;
    const float dy = std::abs(y_prediction - h2->y);

    const float scatterNum = (dx * dx) + (dy * dy);
    const float scatterDenom = 1.f / (h2->z - h1_z);
    const float scatter = scatterNum * scatterDenom * scatterDenom;

    const bool scatter_condition = scatter < MAX_SCATTER;

    if (dx < PARAM_TOLERANCE &&
        dy < PARAM_TOLERANCE &&
        scatter < MAX_SCATTER)
        return scatter;
    else
        return MAX_FLOAT;
}


/**
* @brief Fills dev_hit_candidates.
*
*/

void OMPFillCandidates(const Event& event, std::pair<int, int> hit_candidates[],
        std::pair<int, int> hit_h2_candidates[]) {
    /*
     * This loop used to iterate to cur_sensors >= 2, but then a check was made whether
     * cur_sensor is greater or equal to four.
     */
    #pragma omp parallel num_threads(SENSOR_LEVEL_PARALLELISM) shared(hit_candidates, hit_h2_candidates)
    {
        #pragma omp for schedule(static)
        for(unsigned int cur_sensor = event.number_of_sensors - 1; cur_sensor >= 4; --cur_sensor) {
            const int second_sensor = cur_sensor - 2;
            const bool process_h2_candidates = cur_sensor <= event.number_of_sensors - 3;

            std::pair<float, float> h2_boundaries[event.sensor_hits.nums[cur_sensor]];

            // PS: this loop *GETS VECTORIZED*
            #pragma omp simd
            for (int h0_element=0; h0_element < event.sensor_hits.nums[cur_sensor]; ++h0_element) {
                const int h0_index = event.sensor_hits.starts[cur_sensor] + h0_element;
                // TODO: PS: these accesses have to be aligned to 64 bytes
                const float h0_x = event.hits.Xs[h0_index];
                const float h0_z = event.hits.Zs[h0_index];
                float xmin_h2, xmax_h2;
                OMPFindH2Boundaries(event.sensor_Zs, h0_x, h0_z, cur_sensor, second_sensor, &xmin_h2, &xmax_h2);
                h2_boundaries[h0_element].first = xmin_h2;
                h2_boundaries[h0_element].second = xmax_h2;

            }

            // Sensor dependent calculations
            // Iterate in all hits in z0
            #pragma omp parallel num_threads(HIT_LEVEL_PARALLELISM)
            for (int h0_element=0; h0_element < event.sensor_hits.nums[cur_sensor]; ++h0_element) {
                assert(h0_element < event.sensor_hits.nums[cur_sensor]);
                const int h0_index = event.sensor_hits.starts[cur_sensor] + h0_element;
                struct Hit h0;
                h0.x = event.hits.Xs[h0_index];
                h0.z = event.hits.Zs[h0_index];
                const int hitstarts_s2 = event.sensor_hits.starts[second_sensor];
                const int hitnums_s2 = event.sensor_hits.nums[second_sensor];

                float xmin_h2, xmax_h2;
                std::tie(xmin_h2, xmax_h2) =  h2_boundaries[h0_element];

                bool first_h1_found = false, last_h1_found = false;
                bool first_h2_found = false, last_h2_found = false;

                // Iterate in all hits in z1
                for (int h1_element=0; h1_element<hitnums_s2; ++h1_element) {
                    int h1_index = hitstarts_s2 + h1_element;
                    struct Hit h1;
                    h1.x = event.hits.Xs[h1_index];
                    h1.z = event.hits.Zs[h1_index];

                    if (!last_h1_found) {
                        // Check if h0 and h1 are compatible
                        const float h_dist = std::abs(h1.z - h0.z);
                        const float dxmax = PARAM_MAXXSLOPE_CANDIDATES * h_dist;
                        const bool tol_condition = std::abs(h1.x - h0.x) < dxmax;

                        // Find the first one
                        if (!first_h1_found && tol_condition) {
                            ASSERT(2 * h0_index < 2 * (event.sensor_hits.starts[number_of_sensors-1] + event.sensor_hits.nums[number_of_sensors-1]))

                            hit_candidates[h0_index].first = h1_index;
                            first_h1_found = true;
                        }
                        // The last one, only if the first one has already been found
                        else if (first_h1_found && !tol_condition) {
                            ASSERT(2 * h0_index + 1 < 2 * (event.sensor_hits.starts[number_of_sensors-1] + event.sensor_hits.nums[number_of_sensors-1]))

                            hit_candidates[h0_index].second = h1_index;
                            last_h1_found = true;
                        }
                    }

                    if (process_h2_candidates && !last_h2_found) {
                        if (!first_h2_found && h1.x > xmin_h2) {
                            ASSERT(2 * h0_index < 2 * (event.sensor_hits.starts[number_of_sensors-1] + event.sensor_hits.nums[number_of_sensors-1]))

                            hit_h2_candidates[h0_index].first = h1_index;
                            first_h2_found = true;
                        }
                        else if (first_h2_found && h1.x > xmax_h2) {
                            ASSERT(2 * h0_index + 1 < 2 * (event.sensor_hits.starts[number_of_sensors-1] + event.sensor_hits.nums[number_of_sensors-1]))

                            hit_h2_candidates[h0_index].second = h1_index;
                            last_h2_found = true;
                        }
                    }

                    if (last_h1_found &&
                        (!process_h2_candidates || last_h2_found)) {
                        break;
                    }
                }

                // Note: If first is not found, then both should be -1
                // and there wouldn't be any iteration
                if (first_h1_found && !last_h1_found) {
                    ASSERT(2 * h0_index + 1 < 2 * (event.sensor_hits.starts[number_of_sensors-1] + event.sensor_hits.nums[number_of_sensors-1]))

                    hit_candidates[h0_index].second = hitstarts_s2 + hitnums_s2;
                }

                if (process_h2_candidates && first_h2_found && !last_h2_found) {
                    ASSERT(2 * h0_index + 1 < 2 * (event.sensor_hits.starts[number_of_sensors-1] + event.sensor_hits.nums[number_of_sensors-1]))

                    hit_h2_candidates[h0_index].second = hitstarts_s2 + hitnums_s2;
                }
            }
        }
    } // omp parallel
}

/**
* @brief Performs the track forwarding.
*
*/

void OMPTrackForwarding(const Event& event, std::vector<bool>& hit_used,
        const size_t cur_sensor,
        std::vector<int>& tracks_to_follow, std::vector<int>& weak_tracks,
        const unsigned int prev_ttf, std::vector<Track>& tracklets,
        std::vector<struct Track>& tracks) {

    std::vector<int> new_tracks_to_follow;

    #pragma omp parallel for num_threads(HIT_LEVEL_PARALLELISM) schedule(static)
    for (unsigned ttf_index = prev_ttf; ttf_index < tracks_to_follow.size() ; ++ttf_index) {

        unsigned int trackno, skipped_modules;
        unsigned int best_hit_h2 = 0;
        struct Track t;
        struct Hit h0;

        unsigned fulltrackno = tracks_to_follow[ttf_index];
        // OA: fulltrackno has not only the index into tracks_pointer encoded, but also:
        // bit 31: track_flag is set in the OMPTrackCreation function
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
        h0.x = event.hits.Xs[h0_num];
        h0.y = event.hits.Ys[h0_num];
        h0.z = event.hits.Zs[h0_num];

        ASSERT(h1_num < number_of_hits)
        const float h1_x = event.hits.Xs[h1_num];
        const float h1_y = event.hits.Ys[h1_num];
        const float h1_z = event.hits.Zs[h1_num];

        // Track forwarding over t, for all hits in the next module
        // Line calculations
        const float td = 1.0f / (h1_z - h0.z);
        const float txn = (h1_x - h0.x);
        const float tyn = (h1_y - h0.y);
        const float tx = txn * td;
        const float ty = tyn * td;

        // Search for a best fit
        // Load shared elements

        // OA: once the shared mem goes away much of the logic here becomes redundant.
        // TODO: discuss this - I'm not entirely sure, but a simple loop k=0:hitNums
        // could be enugh
        float best_fit = MAX_FLOAT;

        for (size_t h2_element=0; h2_element < event.sensor_hits.nums[cur_sensor-4]; ++h2_element) {
            size_t h2_index = h2_element + event.sensor_hits.starts[cur_sensor-4];
            struct Hit h2 = {event.hits.Xs[h2_index], event.hits.Ys[h2_index], event.hits.Zs[h2_index]};

            const float fit = OMPFitHitToTrack(tx, ty, &h0, h1_z, &h2);

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
            #pragma omp critical
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

                #pragma omp critical
                {
                    trackno = tracks.size();
                    // XXX PS: without this spell the whole thing doesnt work..
                    tracks.emplace_back();
                }
            }

            // Copy the track into tracks
            ASSERT(trackno < number_of_hits)
            #pragma omp critical
            tracks[trackno] = t;

            // Add the tracks to the bag of tracks to_follow
            #pragma omp critical
            new_tracks_to_follow.push_back(trackno);
        }
        // A track just skipped a module
        // We keep it for another round
        else if (skipped_modules <= MAX_SKIPPED_MODULES) {
            // Form the new mask
            trackno = ((skipped_modules + 1) << 28) | (fulltrackno & 0x8FFFFFFF);

            // Add the tracks to the bag of tracks to_follow
            //const unsigned int ttfP = atomic_add(ttf_insertPointer, 1) % TTF_MODULO;
            #pragma omp critical
            new_tracks_to_follow.push_back(trackno);
        }
        // If there are only three hits in this track,
        // mark it as "doubtful"
        else if (t.hitsNum == 3) {
            //const unsigned int weakP = atomic_add(weaktracks_insertPointer, 1);
            #pragma omp critical
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
*/

void OMPTrackCreation(const Event& event, const size_t cur_sensor, std::pair<int, int> hit_candidates[], int h0_index,
        std::vector<bool>& hit_used, std::pair<int, int> hit_h2_candidates[],
        std::vector<Track>& tracklets, std::vector<int>& tracks_to_follow) {

    // Track creation starts
    struct Hit h0 = {event.hits.Xs[h0_index], event.hits.Ys[h0_index], event.hits.Zs[h0_index]};
    unsigned int best_hit_h1, best_hit_h2;
    float best_fit;


    // Only iterate in the hits indicated by hit_candidates :)
    const int first_h1 = hit_candidates[h0_index].first;
    const int last_h1 = hit_candidates[h0_index].second;

    // Only iterate max_numhits_to_process[0] iterations (with get_local_size(1) threads) :D :D :D
    std::tie(best_hit_h1, best_hit_h2, best_fit) = OMPFindBestFit(event, h0, hit_used, cur_sensor, first_h1, last_h1);

    if (best_fit < MAX_FLOAT) {
        // We have a best fit! - haven't we?
        // Fill in track information

        // Add the track to the bag of tracks
        Track new_tracklet = {3, {static_cast<unsigned int>(h0_index), best_hit_h1, best_hit_h2}};
        tracklets.push_back(new_tracklet);
        unsigned int new_tracklet_idx = tracklets.size()-1;

        // Add the tracks to the bag of tracks to_follow
        // Note: The first bit flag marks this is a tracklet (hitsNum == 3),
        // and hence it is stored in tracklets
        tracks_to_follow.push_back(0x80000000 | new_tracklet_idx);
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
*/
std::vector<Track> OMPSearchByTriplets(const Event& event) {

    std::vector<struct Track> tracks;
    tracks.reserve(MAX_TRACKS);
    DEBUG << "number of sensors: " << event.number_of_sensors << std::endl;
    DEBUG << "number of hits: " << event.number_of_hits << std::endl;

    std::vector<bool> hit_used = std::vector<bool>(event.number_of_hits, false);
    std::pair<int, int> hit_candidates[event.number_of_hits];
    std::pair<int, int> hit_h2_candidates[event.number_of_hits];

    std::vector<int> tracks_to_follow;

    std::vector<int> weak_tracks;
    weak_tracks.reserve(event.number_of_hits);

    std::vector<Track> tracklets;
    tracklets.reserve(event.number_of_hits);


    OMPFillCandidates(event, hit_candidates, hit_h2_candidates);
    unsigned int prev_ttf, last_ttf = 0;

    for (unsigned int cur_sensor = event.number_of_sensors-1; cur_sensor >= 4; --cur_sensor) {

        prev_ttf = last_ttf;
        last_ttf = tracks_to_follow.size();

        OMPTrackForwarding(event, hit_used,
            cur_sensor, tracks_to_follow, weak_tracks, prev_ttf,
            tracklets, tracks);

        for(int h0_element = 0;
            h0_element < event.sensor_hits.nums[cur_sensor];
            ++h0_element) {
            int h0_index = h0_element + event.sensor_hits.starts[cur_sensor];
            if (!hit_used[h0_index]) {
                OMPTrackCreation(event, cur_sensor,
                    hit_candidates, h0_index, hit_used, hit_h2_candidates,
                    tracklets,
                    tracks_to_follow);
            }
        }
    }

    // Process the last bunch of track_to_follows
    #pragma omp parallel num_threads(SENSOR_LEVEL_PARALLELISM*HIT_LEVEL_PARALLELISM)
    {
        #pragma omp for
        for (unsigned int ttf_element = last_ttf; ttf_element< tracks_to_follow.size(); ++ttf_element) {

            const int fulltrackno = tracks_to_follow[ttf_element];
            const bool is_track = fulltrackno & 0x80000000;
            const int trackno = fulltrackno & 0x0FFFFFFF;

            // Here we are only interested in three-hit tracks,
            // to mark them as "doubtful"
            if (is_track) {
                #pragma omp critical
                {
                    weak_tracks.push_back(trackno);
                }
                ASSERT(weak_tracks.size() < number_of_hits);
            }
        } // omp for

        // Compute the three-hit tracks left
        #pragma omp for
        for (unsigned weak_track_idx=0; weak_track_idx<weak_tracks.size(); ++weak_track_idx) {
            auto& weak_track = weak_tracks[weak_track_idx];
            // Load the tracks from the tracklets
            const struct Track t = tracklets[weak_track];

            // Store them in the tracks bag iff they
            // are made out of three unused hits
            if (!hit_used[t.hits[0]] &&
                !hit_used[t.hits[1]] &&
                !hit_used[t.hits[2]]) {
                ASSERT(trackno < MAX_TRACKS)
                #pragma omp critical
                {
                    tracks.push_back(t);
                }
            }
        } // omp for
    } // omp parallel
    return tracks;
}
