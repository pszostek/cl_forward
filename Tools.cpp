
#include "Tools.h"
#include "Definitions.h"

/*int*   h_no_sensors;
int*   h_no_hits;
int*   h_sensor_Zs;
int*   h_sensor_hitStarts;
int*   h_sensor_hitNums;
unsigned int* h_hit_IDs;
float* h_hit_Xs;
float* h_hit_Ys;
float* h_hit_Zs;*/

/* convert the kernel file into a string */
int convertClToString(const char *filename, std::string& s)
{
  char*  str;
  std::fstream f(filename, (std::fstream::in | std::fstream::binary));

  if (f.is_open()) {
    size_t size;
    size_t fileSize;
    f.seekg(0, std::fstream::end);
    size = fileSize = (size_t)f.tellg();
    f.seekg(0, std::fstream::beg);
    str = new char[size+1];

    if (!str) {
      f.close();
      return 0;
    }

    f.read(str, fileSize);
    f.close();
    str[size] = '\0';
    s = str;
    delete[] str;
    return 0;
  }

  std::cout << "Error: failed to open file\n:" << filename << std::endl;
  return -1;
}

void preorder_by_x(std::vector<const std::vector<uint8_t>* > & input) {
  // Order *all* the input vectors by h_hit_Xs natural order
  // per sensor
  const int number_of_input_files = input.size();
  //const std::vector<uint8_t>* startingEvent_input = input[0];
  //this call is just to set h_no_sensors and number_of_sensors
  //setHPointersFromInput((uint8_t*) &(*startingEvent_input)[0], startingEvent_input->size(),
  //  number_of_sensors, number_of_hits, h_sensors_Zs, sensor_hits, hits);
#if(EVENT_LEVEL_PARALLELISM != 0)
  #pragma omp parallel for num_threads(EVENT_LEVEL_PARALLELISM)
#endif
  for (int i=0; i < number_of_input_files; ++i) {
    int number_of_sensors, number_of_hits;
    int *sensor_Zs;
    unsigned int *hit_IDs;
    SensorHits sensor_hits;
    Hits hits;

    int acc_hitnums = 0;
    const std::vector<uint8_t>* event_input = input[i];
    setHPointersFromInput((uint8_t*) &(*event_input)[0], event_input->size(),
      number_of_sensors, number_of_hits, sensor_Zs, sensor_hits,
      hit_IDs, hits);

#if(SENSOR_LEVEL_PARALLELISM != 0)
    #pragma omp parallel for num_threads(SENSOR_LEVEL_PARALLELISM) shared(acc_hitnums)
#endif
    for (int j=0; j<number_of_sensors; j++) {
      const int hitnums = sensor_hits.nums[j];
      quicksort(hits.Xs, hits.Ys, hits.Zs, hit_IDs, acc_hitnums, acc_hitnums + hitnums - 1);
#if(SENSOR_LEVEL_PARALLELISM != 0)
      #pragma omp critical
#endif
      acc_hitnums += hitnums;
    }
  }
}

void setHPointersFromInput(uint8_t * input, size_t size,
  int& number_of_sensors, int& number_of_hits, int*& sensor_Zs, SensorHits& sensor_hits,
  unsigned int*& hit_IDs, Hits& hits) {
  uint8_t * end = input + size;

  number_of_sensors   = *((int32_t*)input); input += sizeof(int32_t);
  number_of_hits      = *((int32_t*)input); input += sizeof(int32_t);
  sensor_Zs           = (int32_t*)input; input += sizeof(int32_t) * number_of_sensors;
  sensor_hits.starts = (int32_t*)input; input += sizeof(int32_t) * number_of_sensors;
  sensor_hits.nums   = (int32_t*)input; input += sizeof(int32_t) * number_of_sensors;
  hit_IDs             = (uint32_t*)input; input += sizeof(uint32_t) * number_of_hits;
  hits.Xs             = (float*)  input; input += sizeof(float)   * number_of_hits;
  hits.Ys             = (float*)  input; input += sizeof(float)   * number_of_hits;
  hits.Zs             = (float*)  input; input += sizeof(float)   * number_of_hits;
  ASSERT (input == end);
}

std::map<std::string, float> calcResults(std::vector<float>& times){
    // sqrt ( E( (X - m)2) )
    std::map<std::string, float> results;
    float deviation = 0.0f, variance = 0.0f, mean = 0.0f, min = MAX_FLOAT, max = 0.0f;

    for(auto it = times.begin(); it != times.end(); it++){
        const float seconds = (*it);
        mean += seconds;
        variance += seconds * seconds;

        if (seconds < min) min = seconds;
        if (seconds > max) max = seconds;
    }

    mean /= times.size();
    variance = (variance / times.size()) - (mean * mean);
    deviation = std::sqrt(variance);

    results["variance"] = variance;
    results["deviation"] = deviation;
    results["mean"] = mean;
    results["min"] = min;
    results["max"] = max;

    return results;
}

#ifdef WITH_OPENCL
void checkClError(const cl_int errcode_ret) {
  // CHECK_OPENCL_ERROR(errcode_ret, "Error ");
  if (errcode_ret != CL_SUCCESS) {
    std::cerr << "Error " << errcode_ret << std::endl;
    exit(-1);
  }
}
#endif


std::map<int, int> map_z_coordinate_to_sensors(const Event& event) {
    std::map<int, int> zhit_to_module;
    for(int j=0; j<event.number_of_sensors; ++j){
        const int z = event.sensor_Zs[j];
        zhit_to_module[z] = j;
    }
  // some hits z may not correspond to a sensor's,
  // but be close enough
    for(int j=0; j<event.number_of_hits; ++j){
        const int z = (int) event.hits.Zs[j];
        if (zhit_to_module.find(z) == zhit_to_module.end()){
            const int sensor = findClosestModule(z, zhit_to_module);
            zhit_to_module[z] = sensor;
        }
    }
    return zhit_to_module;
}

/**
 * Prints tracks
 * Track #n, length <length>:
 *  <ID> module <module>, x <x>, y <y>, z <z>
 *
 * @param tracks
 * @param trackNumber
 */
void printTrack(const Track& track, const Event& event, std::ofstream& outstream){
    std::map<int, int> zhit_to_module = map_z_coordinate_to_sensors(event);

    for(unsigned int hit_idx = 0; hit_idx < track.hitsNum; ++hit_idx){
        const int hitNumber = track.hits[hit_idx];
        const unsigned int id = event.hit_IDs[hitNumber];
        const float x = event.hits.Xs[hitNumber];
        const float y = event.hits.Ys[hitNumber];
        const float z = event.hits.Zs[hitNumber];
        const int module = zhit_to_module.at((int) z);

        outstream << " " << std::setw(8) << id << " (" << hitNumber << ")"
        << " module " << std::setw(2) << module
        << ", x " << std::setw(6) << x
        << ", y " << std::setw(6) << y
        << ", z " << std::setw(6) << z << std::endl;
    }
    outstream << std::endl;
}

void writeTextTracks(const std::vector<Track>& tracks,
    const Event& event, std::ofstream& os) {

    size_t track_idx = 0;
    for(auto track: tracks) {
        os << "Track #" << track_idx << ", length " << (int) track.hitsNum << std::endl;
        printTrack(track, event, os);
        ++track_idx;
    }
}


/**
 * Write tracks in binary format.
 * We write both the event information as well as the tracks themselves. This will make
 * post-processing easier
 */
 void writeBinTracks(const std::vector<Track>& tracks, const Event& event, std::ofstream& os) {

     // first get position and value of infile data size
     std::ifstream infile(event.filename.c_str(), std::ifstream::binary);

     uint32_t funcNameLen;
     uint32_t dataSize;

     infile.read((char*) &funcNameLen, sizeof(uint32_t));
     infile.seekg(4+funcNameLen);
     infile.read((char*) &dataSize, sizeof(uint32_t));
     infile.seekg(0);
     // copy data to output file
     os << infile.rdbuf();
     infile.close();

     // write track info while keeping track of extra bytes written
     uint32_t extra_bytes = 0;
     int ntracks = static_cast<int>(tracks.size());
     os.write((char*) &ntracks, sizeof(ntracks));
     extra_bytes += sizeof(ntracks);
     for (auto track: tracks) {
         os.write((char*) &(track.hitsNum), sizeof(track.hitsNum));
         extra_bytes += sizeof(track.hitsNum);
          for(unsigned int hit_idx = 0; hit_idx < track.hitsNum; ++hit_idx){
            const int hitNumber = track.hits[hit_idx];
            const unsigned int id = event.hit_IDs[hitNumber];
            os.write((char*) &id, sizeof(id));
            extra_bytes += sizeof(id);
        }
     }
     // update file size in header
     dataSize += extra_bytes;
     os.seekp(4+funcNameLen);
     os.write((char*) &dataSize, sizeof(uint32_t));
 }

/**
 * The z of the hit may not correspond to any z in the sensors.
 * @param  z
 * @param  zhit_to_module
 * @return                sensor number
 */
int findClosestModule(const int z, const std::map<int, int>& zhit_to_module){
  auto it = zhit_to_module.find(z);
  if (it != zhit_to_module.end())
    return it->second;

  int error = 0;
  while(true){
    error++;
    const int lowerAttempt = z - error;
    const int higherAttempt = z + error;

    auto it_lowerAttempt = zhit_to_module.find(lowerAttempt);
    if (it_lowerAttempt != zhit_to_module.end()){
      return it_lowerAttempt->second;
    }

    auto it_higherAttempt = zhit_to_module.find(higherAttempt);
    if (it_higherAttempt != zhit_to_module.end()){
      return it_higherAttempt->second;
    }
  }
}

void quicksort (float* a, float* b, float* c, unsigned int* d, int start, int end) {
#if(SENSOR_LEVEL_PARALLELISM != 0)
    #pragma omp parallel
    {
        #pragma omp single
        {
#endif
            if (start < end) {
                const int pivot = divide(a, b, c, d, start, end);
#if(SENSOR_LEVEL_PARALLELISM != 0)
                #pragma omp task
#endif
                    quicksort(a, b, c, d, start, pivot - 1);
#if(SENSOR_LEVEL_PARALLELISM != 0)
                #pragma omp task
#endif
                    quicksort(a, b, c, d, pivot + 1, end);
            }
#if(SENSOR_LEVEL_PARALLELISM != 0)
        }
    }
#endif
}

int divide (float* a, float* b, float* c, unsigned int* d, int start, int end) {
    int left;
    int right;
    float pivot;

    pivot = a[start];
    left = start;
    right = end;

    while (left < right) {
        while (a[right] > pivot) {
            right--;
        }

        while ((left < right) && (a[left] <= pivot)) {
            left++;
        }

        if (left < right) {
            swap(a[left], a[right]);
            swap(b[left], b[right]);
            swap(c[left], c[right]);
            swap(d[left], d[right]);
        }
    }

    swap(a[right], a[start]);
    swap(b[right], b[start]);
    swap(c[right], c[start]);
    swap(d[right], d[start]);

    return right;
}

template<typename T>
void swap (T& a, T& b) {
    T temp = a;
    a = b;
    b = temp;
}
