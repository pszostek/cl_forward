
#include "PixelSearchByTriplet.h"
#include "SerialKernel.h"

extern int*   h_no_sensors;
extern int*   h_no_hits;
extern int*   h_sensor_Zs;
extern int*   h_sensor_hitStarts;
extern int*   h_sensor_hitNums;
extern unsigned int* h_hit_IDs;


int independent_execute(
    const std::vector<std::vector<uint8_t> > & input,
    std::vector<std::vector<uint8_t> > & output,
    ExecMode mode) {

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
      return cpuPixelSearchByTripletSerialRun(converted_input, output);
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
        std::vector<std::vector<uint8_t> > & output) {
    DEBUG << "executing cpuPixelSearchByTriplet with " << input.size() << " events" << std::endl;

    Hits hits;
    SensorHits sensor_hits;
    // Define how many blocks / threads we need to deal with numberOfEvents
    // Each execution will return a different output
    output.resize(input.size());

    for (unsigned int i = 0; i < input.size(); ++i) {
        int numTracks = 0;
        DEBUG << "Processing event " << i << std::endl;
        Track *tracks = new Track[MAX_TRACKS];

        const std::vector<uint8_t>* event_input = input[i];
        setHPointersFromInput((uint8_t*) &(*event_input)[0], event_input->size(), sensor_hits, hits);

        numTracks = serialSearchByTriplets(tracks,(uint8_t*) &(*event_input)[0]);
        DEBUG << "Done." << std::endl;
        DEBUG << "Found " << numTracks << " tracks." << std::endl;




      // Calculate z to sensor map
      std::map<int, int> zhit_to_module;
      setHPointersFromInput((uint8_t*) &(*(input[i]))[0], input[i]->size(), sensor_hits, hits);
      int number_of_sensors = *h_no_sensors;
        // map to convert from z of hit to module
      for(int j=0; j<number_of_sensors; ++j){
          const int z = h_sensor_Zs[j];
          zhit_to_module[z] = j;
      }
      // Some hits z may not correspond to a sensor's,
      // but be close enough
      for(int j=0; j<*h_no_hits; ++j){
          const int z = (int) hits.Zs[j];
          if (zhit_to_module.find(z) == zhit_to_module.end()){
            const int sensor = findClosestModule(z, zhit_to_module);
            zhit_to_module[z] = sensor;
          }
      }


        // Print to output file with event no.
        std::ofstream outfile (std::string(RESULTS_FOLDER) + std::string("/") + toString(i) + std::string("_serial.out"));
        DEBUG << "writing to: " << std::string(RESULTS_FOLDER) + std::string("/") + toString(i) + std::string("_serial.out") << std::endl;
        for(int j=0; j<numTracks; ++j){
            printTrack(tracks, j, zhit_to_module, hits, outfile);
        }
        outfile.close();

        delete[] tracks;
    }

    return 0;
}
