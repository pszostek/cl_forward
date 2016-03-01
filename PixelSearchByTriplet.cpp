
#include "Event.h"
#include "PixelSearchByTriplet.h"
#include "SerialKernel.h"


// extern int*   h_no_sensors;
// extern int*   h_no_hits;
// extern int*   h_sensor_Zs;
// extern unsigned int* h_hit_IDs;

#ifdef TIMING_ENABLED
#include <omp.h>
#include <cstdint>
struct Timing {
    uint32_t cycles_low, cycles_high;
    uint64_t cycles;

    Timing() : cycles(0) {}

    void start() {
        asm volatile (
        "CPUID\n\t"/*serialize*/
        "RDTSC\n\t"/*read the clock*/
        "mov %%edx, %0\n\t"
        "mov %%eax, %1\n\t"
        : "=r" (cycles_high), "=r" (cycles_low)
        :: "%rax", "%rbx", "%rcx", "%rdx");
    }
    void stop() {
        uint32_t cycles_high1, cycles_low1;
        asm volatile (
            "RDTSCP\n\t"/*read the clock*/
            "mov %%edx, %0\n\t"
            "mov %%eax, %1\n\t"
            "CPUID\n\t"
            : "=r" (cycles_high1), "=r" (cycles_low1)
            :: "%rax", "%rbx", "%rcx", "%rdx");
        cycles += ((uint64_t(cycles_high1) << 32) + uint64_t(cycles_low1))
            - ((uint64_t(cycles_high) << 32) + uint64_t(cycles_low));
    }
} timing;
#endif


std::string get_output_filename(std::string input_filename, OutType output_type) {
    std::size_t last_slash_idx = input_filename.find_last_of('/');
    std::string basename = input_filename.substr(last_slash_idx+1);
    std::string basename_wo_extension = basename.substr(0, basename.size()-4);
    std::string output_filename;

    if (output_type == OutType::Text) {
        output_filename = std::string(RESULTS_FOLDER) + "/" + basename_wo_extension + "_serial_txt.out";
    } else if (output_type == OutType::Binary) {
        output_filename = std::string(RESULTS_FOLDER) + "/" + basename_wo_extension + "_serial_bin.out";
    }
    return output_filename;
}

int independent_execute(
    const std::vector<std::vector<uint8_t> > & input,
    std::vector<std::vector<uint8_t> > & output,
    std::vector<std::string> &filenames,
    ExecMode mode, OutType outtype) {

  std::vector<const std::vector<uint8_t>* > converted_input;
  converted_input.resize(input.size());

  for (unsigned int i=0; i<input.size(); ++i) {
    converted_input[i] = &(input[i]);
  }

  std::cout << std::fixed << std::setprecision(2);

  // Order input hits by X
  preorder_by_x(converted_input);

  if (mode == ExecMode::Serial){
      return cpuPixelSearchByTripletSerialRun(converted_input, output, filenames, outtype);
#ifdef WITH_OPENCL
  } else if (mode == ExecMode::OpenCl) {
      return gpuPixelSearchByTripletInvocation(converted_input, output);
#endif
  } else if (mode == ExecMode::OpenMP) {
      return cpuPixelSearchByTripletOpenMPRun(converted_input, output, filenames, outtype);
  } else {
      DEBUG << "not yet implemented";
      return 0;
  }
}

void independent_post_execute(const std::vector<std::vector<uint8_t> > & output) {
    // DEBUG << "post_execute invoked" << std::endl;
    DEBUG << std::endl << "Size of output: " << output.size() << " entries" << std::endl;
}

#ifdef WITH_OPENCL
int gpuPixelSearchByTriplet(
    const std::vector<const std::vector<uint8_t>* > & input,
    std::vector<std::vector<uint8_t> > & output) {

  // Silent execution
  std::cout << std::fixed << std::setprecision(2);
  logger::ll.verbosityLevel = 0;
  return gpuPixelSearchByTripletInvocation(input, output);
}

/**
 * Common entrypoint for Gaudi and non-Gaudi
 * @param input
 * @param output
 */
int gpuPixelSearchByTripletInvocation(
    const std::vector<const std::vector<uint8_t>* > & input,
    std::vector<std::vector<uint8_t> > & output) {
  DEBUG << "Invoking gpuPixelSearchByTriplet with " << input.size() << " events" << std::endl;

  // Define how many blocks / threads we need to deal with numberOfEvents
  // Each execution will return a different output
  output.resize(input.size());

  // Execute maximum n number of events every time
  const int max_events_to_process_per_kernel = 16000;

  for (unsigned int i=0; i<input.size(); i+=max_events_to_process_per_kernel){
    int events_to_process = input.size() - i;
    if (events_to_process > max_events_to_process_per_kernel)
      events_to_process = max_events_to_process_per_kernel;

    invokeParallelSearch(i, events_to_process, input, output);
  }

  return 0;
}
#endif

/**
  * Common entrypoint for Gaudi and non-Gaudi
  * @param input
  * @param output
  */
int cpuPixelSearchByTripletSerialRun(
        const std::vector<const std::vector<uint8_t>* > & input,
        std::vector<std::vector<uint8_t> > & output, std::vector<std::string> &filenames, OutType outtype) {
    DEBUG << "executing cpuPixelSearchByTriplet with " << input.size() << " events" << std::endl;
    output.resize(input.size());

    std::vector<std::vector<Track>> event_tracks;
    std::vector<Event> events;

    event_tracks.reserve(input.size());
    events.reserve(input.size());

    for (unsigned int input_index = 0; input_index < input.size(); ++input_index) {
        DEBUG << "Processing data frame " << input_index << std::endl;

        const std::vector<uint8_t>* event_input = input[input_index];
        Event event((uint8_t*) &(*event_input)[0], event_input->size(), filenames[input_index]);
#ifdef TIMING_ENABLED
        timing.start();
#endif
        std::vector<Track> tracks = serialSearchByTriplets(event);
#ifdef  TIMING_ENABLED
        timing.stop();
#endif
        DEBUG << "Done. Found " << tracks.size() <<" tracks." << std::endl;
        events.push_back(std::move(event));
        event_tracks.push_back(std::move(tracks));
    }
#ifdef TIMING_ENABLED
    #include <iostream>
     INFO << "%%% TIME: " << timing.cycles << " cycles." << std::endl;
#endif
    for(size_t output_idx=0; output_idx<event_tracks.size(); ++output_idx) { 
        // Print to output file with event no.
        const Event event = events.at(output_idx);
        const std::vector<Track> tracks = event_tracks.at(output_idx);

        if (outtype == OutType::Text) {
            std::string fileName = get_output_filename(event.filename, OutType::Text);
            std::ofstream outfile(fileName, std::ios::out);
            DEBUG << "writing to: " << fileName << std::endl;
            unsigned int track_idx = 0;
            writeTextTracks(tracks, event, outfile);

            outfile.close();
        } else if (outtype == OutType::Binary) {
            std::string fileName = get_output_filename(event.filename, OutType::Binary);
            std::ofstream outfile(fileName, std::ios::out | std::ios::binary);
            DEBUG << "writing to: " << fileName << std::endl;
            writeBinTracks(tracks, event, outfile);
            outfile.close();
        }
    }

    return 0;
}

int cpuPixelSearchByTripletOpenMPRun(
        const std::vector<const std::vector<uint8_t>* > & input,
        std::vector<std::vector<uint8_t> > & output, std::vector<std::string> &filenames, OutType outtype) {
    DEBUG << "executing cpuPixelSearchByTripletOpenMPRun with " << input.size() << " events" << std::endl;
    output.resize(input.size());

    // std::vector<std::vector<Track>> event_tracks;
    // std::vector<Event> events;
    std::vector<Track>* event_tracks[input.size()];
    Event* events[input.size()];
    double runtime;

    // event_tracks.reserve(input.size());
    // events.reserve(input.size());

#if(EVENT_LEVEL_PARALLELISM != 0)
    #pragma omp parallel shared(event_tracks, events, input) num_threads(EVENT_LEVEL_PARALLELISM)
    {
#endif
#ifdef TIMING_ENABLED
#if(EVENT_LEVEL_PARALLELISM != 0)
    #pragma omp single
#endif
    {
    runtime = -omp_get_wtime();
    timing.start();
    }
#endif // TIMING_ENABLED
#if(EVENT_LEVEL_PARALLELISM != 0)
        #pragma omp for schedule(dynamic,1)
#endif
        for (unsigned int input_index = 0; input_index < input.size(); ++input_index) {
            DEBUG << "Processing data frame " << input_index << std::endl;

            const std::vector<uint8_t>* event_input = input[input_index];
            Event* event = new Event((uint8_t*) &(*event_input)[0], event_input->size(), filenames[input_index]);

            std::vector<Track> tracks = OMPSearchByTriplets(*event);
            DEBUG << "Done. Found " << tracks.size() <<" tracks." << std::endl;
            // events.push_back(std::move(event));
            // event_tracks.push_back(std::move(tracks));
            events[input_index] = event;
            event_tracks[input_index] = new std::vector<Track>(std::move(tracks));
        }
#ifdef  TIMING_ENABLED
#if(EVENT_LEVEL_PARALLELISM != 0)
    #pragma omp single
#endif
    {
    runtime += omp_get_wtime();
    timing.stop();
    }
#endif
#if(EVENT_LEVEL_PARALLELISM != 0)
    } // omp parallel
#endif
#ifdef TIMING_ENABLED
    #include <iostream>
    std::cout << "%%% TIME: " << runtime << " seconds, "<< timing.cycles << " cycles." << std::endl;
#endif
    for(size_t output_idx=0; output_idx<input.size(); ++output_idx) { 
        // Print to output file with event no.
        const Event* event = events[output_idx];
        // const std::vector<Track> tracks = event_tracks.at(output_idx);
        const std::vector<Track>* tracks = event_tracks[(output_idx)];

        if (outtype == OutType::Text) {
            std::string fileName = get_output_filename(event->filename, OutType::Text);
            std::ofstream outfile(fileName, std::ios::out);
            DEBUG << "writing to: " << fileName << std::endl;
            unsigned int track_idx = 0;
            writeTextTracks(*tracks, *event, outfile);

            outfile.close();
        } else if (outtype == OutType::Binary) {
            std::string fileName = get_output_filename(event->filename, OutType::Binary);
            std::ofstream outfile(fileName, std::ios::out | std::ios::binary);
            DEBUG << "writing to: " << fileName << std::endl;
            writeBinTracks(*tracks, *event, outfile);
            outfile.close();
        }
        delete event;
        delete tracks;
    }

    return 0;
}
