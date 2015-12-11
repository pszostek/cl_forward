
#include "Event.h"
#include "PixelSearchByTriplet.h"
#include "SerialKernel.h"

// extern int*   h_no_sensors;
// extern int*   h_no_hits;
// extern int*   h_sensor_Zs;
// extern unsigned int* h_hit_IDs;

#ifdef TIMING_ENABLED
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
  logger::ll.verbosityLevel = 3;

  // Order input hits by X
  preorder_by_x(converted_input);

  if (mode == ExecMode::OpenCl) {
      return gpuPixelSearchByTripletInvocation(converted_input, output);
  } else if (mode == ExecMode::Serial){
      return cpuPixelSearchByTripletSerialRun(converted_input, output, filenames, outtype);
  } else {
      DEBUG << "not yet implemented";
      return 0;
  }
}

void independent_post_execute(const std::vector<std::vector<uint8_t> > & output) {
    // DEBUG << "post_execute invoked" << std::endl;
    DEBUG << std::endl << "Size of output: " << output.size() << " entries" << std::endl;
}

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

    for (unsigned int input_index = 0; input_index < input.size(); ++input_index) {
        DEBUG << "Processing data frame " << input_index << std::endl;

        const std::vector<uint8_t>* event_input = input[input_index];
        Event event((uint8_t*) &(*event_input)[0], event_input->size(), filenames[input_index]);
#ifdef TIMING_ENABLED
        timing.start();
#endif
        auto tracks = event.serialSearchByTriplets();
#ifdef  TIMING_ENABLED
        timing.stop();
#endif
        DEBUG << "Done. Found " << tracks.size() <<" tracks." << std::endl;

      // Calculate z to sensor map
        std::map<int, int> zhit_to_module;
        // TODO: Why this guy was called another time here?
        // setHPointersFromInput((uint8_t*) &(*(input[input_index]))[0], input[input_index]->size(),
        //     number_of_sensors, number_of_hits, sensor_Zs, sensor_hits,
        //     hit_IDs, hits);

        // map to convert from z of hit to module
        for(int j=0; j<event.number_of_sensors; ++j){
            const int z = event.sensor_Zs[j];
            zhit_to_module[z] = j;
        }
      // Some hits z may not correspond to a sensor's,
      // but be close enough
        for(int j=0; j<event.number_of_hits; ++j){
            const int z = (int) event.hits.Zs[j];
            if (zhit_to_module.find(z) == zhit_to_module.end()){
                const int sensor = findClosestModule(z, zhit_to_module);
                zhit_to_module[z] = sensor;
            }
        }

        char *c_inputfname = new char[event.filename.size()];
        std::strcpy(c_inputfname, event.filename.c_str());
        const char *c_inputbase = basename(c_inputfname);
        std::string input_basename(c_inputbase);
        // Print to output file with event no.
        if (outtype == OutType::Text) {
            std::string fileName = std::string(RESULTS_FOLDER) + std::string("/") + input_basename.substr(0,input_basename.size() - 4) + std::string("_serial_txt.out");
            std::ofstream outfile(fileName, std::ios::out);
            DEBUG << "writing to: " << fileName << std::endl;
            unsigned int track_idx = 0;
            for(auto track: tracks) {
                outfile << "Track #" << track_idx << ", length " << (int) track.hitsNum << std::endl;
                printTrack(track, zhit_to_module, event, outfile);
                ++track_idx;
            }
            outfile.close();
        } else if (outtype == OutType::Binary) {
            std::string fileName = std::string(RESULTS_FOLDER) + std::string("/") + input_basename.substr(0,input_basename.size() - 4) + std::string("_serial_bin.out");
            std::ofstream outfile(fileName, std::ios::out | std::ios::binary);
            DEBUG << "writing to: " << fileName << std::endl;
            writeBinTracks(tracks, event, outfile);
            outfile.close();
        }
    }
#ifdef TIMING_ENABLED
    #include <iostream>
    std::cout << "%%% TIME: " << timing.cycles << " cycles." << std::endl;
#endif

    return 0;
}
