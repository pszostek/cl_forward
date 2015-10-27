#include "GpuKernelInvoker.h"
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS

int    h_no_sensors;
int    h_no_hits;
int*   h_sensor_Zs;
unsigned int* h_hit_IDs;
SensorHits sensor_hits;
Hits hits;

int invokeParallelSearch(
    const int startingEvent,
    const int eventsToProcess,
    const std::vector<const std::vector<uint8_t>* > & input,
    std::vector<std::vector<uint8_t> > & output) {
  cl_int errcode_ret;
  const std::vector<uint8_t>* startingEvent_input = input[startingEvent];
  setHPointersFromInput((uint8_t*) &(*startingEvent_input)[0], startingEvent_input->size(),
    h_no_sensors, h_no_hits, h_sensor_Zs, sensor_hits, h_hit_IDs, hits);

  // Startup settings
  size_t global_work_size[2] = { (size_t) NUMTHREADS_X * eventsToProcess, 2 };
  size_t local_work_size[2] = { (size_t) NUMTHREADS_X, 2 };
  cl_uint work_dim = 2;

  // Choose platform according to the macro DEVICE_PREFERENCE
  cl_device_id* devices;
  cl_platform_id platform = NULL;
  clChoosePlatform(devices, platform);

  // Step 3: Create context
  cl_context context = clCreateContext(NULL, 1, devices, NULL, NULL, &errcode_ret); checkClError(errcode_ret);

  // Step 4: Creating command queue associate with the context
  cl_command_queue commandQueue = clCreateCommandQueue(context, devices[DEVICE_NUMBER], CL_QUEUE_PROFILING_ENABLE, NULL);

  // Step 5: Create program object - KernelDefinitions.h + Kernel.cl
  std::string definitions_str, kernel_str, source_str;
  clCheck(convertClToString("KernelDefinitions.h", definitions_str));
  clCheck(convertClToString("Kernel.cl", kernel_str));
  source_str = definitions_str + kernel_str;
  const char* source = source_str.c_str();
  size_t sourceSize[] = { source_str.size() };
  cl_program program = clCreateProgramWithSource(context, 1, &source, sourceSize, NULL);

  // Step 6: Build program
  const char* buildOptions = "";
  // const char* buildOptions = "-cl-nv-maxrregcount=32";
  // const char* buildOptions = "-g -s /home/dcampora/nfs/projects/gpu/tf_opencl/KernelDefinitions.cl -s /home/dcampora/nfs/projects/gpu/tf_opencl/Kernel.cl";
  cl_int status = clBuildProgram(program, 1, devices, buildOptions, NULL, NULL);

  if (status != CL_SUCCESS) {
    std::cerr << "Error string: " << getErrorString(status) << std::endl;

    if (status == CL_BUILD_PROGRAM_FAILURE) {
      size_t log_size;
      clGetProgramBuildInfo(program, devices[DEVICE_NUMBER], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
      char* log = (char *) malloc(log_size);
      clGetProgramBuildInfo(program, devices[DEVICE_NUMBER], CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
      std::cerr << "Build log: " << std::endl << log << std::endl;
    }

    exit(-1);
  }

  size_t size;
  clCheck(clGetProgramBuildInfo(program, devices[DEVICE_NUMBER], CL_PROGRAM_BUILD_LOG , 0, NULL, &size));

  // Step 7: Memory

  // Allocate memory
  // Prepare event offset and hit offset
  std::vector<int> event_offsets;
  std::vector<int> hit_offsets;
  int acc_size = 0, acc_hits = 0;
  for (int i=0; i<eventsToProcess; ++i) {
    EventBeginning* event = (EventBeginning*) &(*(input[startingEvent + i]))[0];
    const int event_size = input[startingEvent + i]->size();

    event_offsets.push_back(acc_size);
    hit_offsets.push_back(acc_hits);

    acc_size += event_size;
    acc_hits += event->numberOfHits;
  }

  // Allocate CPU buffers
  const int atomic_space = NUM_ATOMICS + 1;
  int* atomics = (int*) malloc(eventsToProcess * atomic_space * sizeof(int));
  int* hit_candidates = (int*) malloc(2 * acc_hits * sizeof(int));

  // Allocate GPU buffers
  cl_mem dev_tracks = clCreateBuffer(context, CL_MEM_READ_WRITE, eventsToProcess * MAX_TRACKS * sizeof(Track), NULL, &errcode_ret); checkClError(errcode_ret);
  cl_mem dev_tracklets = clCreateBuffer(context, CL_MEM_READ_WRITE, acc_hits * sizeof(Track), NULL, &errcode_ret); checkClError(errcode_ret);
  cl_mem dev_weak_tracks = clCreateBuffer(context, CL_MEM_READ_WRITE, acc_hits * sizeof(cl_int), NULL, &errcode_ret); checkClError(errcode_ret);
  cl_mem dev_tracks_to_follow = clCreateBuffer(context, CL_MEM_READ_WRITE, eventsToProcess * TTF_MODULO * sizeof(cl_int), NULL, &errcode_ret); checkClError(errcode_ret);
  cl_mem dev_atomicsStorage = clCreateBuffer(context, CL_MEM_READ_WRITE, eventsToProcess * atomic_space * sizeof(cl_int), NULL, &errcode_ret); checkClError(errcode_ret);
  cl_mem dev_event_offsets = clCreateBuffer(context, CL_MEM_READ_ONLY, event_offsets.size() * sizeof(cl_int), NULL, &errcode_ret); checkClError(errcode_ret);
  cl_mem dev_hit_offsets = clCreateBuffer(context, CL_MEM_READ_ONLY, hit_offsets.size() * sizeof(cl_int), NULL, &errcode_ret); checkClError(errcode_ret);
  cl_mem dev_hit_used = clCreateBuffer(context, CL_MEM_READ_WRITE, acc_hits * sizeof(cl_bool), NULL, &errcode_ret); checkClError(errcode_ret);
  cl_mem dev_input = clCreateBuffer(context, CL_MEM_READ_ONLY, acc_size * sizeof(cl_char), NULL, &errcode_ret); checkClError(errcode_ret);
  cl_mem dev_best_fits = clCreateBuffer(context, CL_MEM_READ_WRITE, eventsToProcess * NUMTHREADS_X * MAX_NUMTHREADS_Y * sizeof(cl_float), NULL, &errcode_ret); checkClError(errcode_ret);
  cl_mem dev_hit_candidates = clCreateBuffer(context, CL_MEM_READ_WRITE, 2 * acc_hits * sizeof(cl_int), NULL, &errcode_ret); checkClError(errcode_ret);
  cl_mem dev_hit_h2_candidates = clCreateBuffer(context, CL_MEM_READ_WRITE, 2 * acc_hits * sizeof(cl_int), NULL, &errcode_ret); checkClError(errcode_ret);

  clCheck(clEnqueueWriteBuffer(commandQueue, dev_event_offsets, CL_TRUE, 0, event_offsets.size() * sizeof(cl_int), &event_offsets[0], 0, NULL, NULL));
  clCheck(clEnqueueWriteBuffer(commandQueue, dev_hit_offsets, CL_TRUE, 0, hit_offsets.size() * sizeof(cl_int), &hit_offsets[0], 0, NULL, NULL));

  acc_size = 0;
  for (int i=0; i<eventsToProcess; ++i) {
    clCheck(clEnqueueWriteBuffer(commandQueue, dev_input, CL_TRUE, acc_size, input[startingEvent + i]->size(), &(*(input[startingEvent + i]))[0], 0, NULL, NULL));
    acc_size += input[startingEvent + i]->size();
  }

  // Step 8: Create kernel object
  cl_kernel kernel = clCreateKernel(program, "clSearchByTriplets", NULL);

  // Step 9: Sets Kernel arguments
  clCheck(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &dev_tracks));
  clCheck(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &dev_input));
  clCheck(clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &dev_tracks_to_follow));
  clCheck(clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &dev_hit_used));
  clCheck(clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &dev_atomicsStorage));
  clCheck(clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &dev_tracklets));
  clCheck(clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &dev_weak_tracks));
  clCheck(clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &dev_event_offsets));
  clCheck(clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &dev_hit_offsets));
  clCheck(clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &dev_best_fits));
  clCheck(clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &dev_hit_candidates));
  clCheck(clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &dev_hit_h2_candidates));

  // Adding timing
  // Timing calculation
  // XXX: reset to 4!
  unsigned int niterations = 1;
  unsigned int nexperiments = 1;

  std::vector<std::vector<float>> time_values {nexperiments};
  std::vector<std::map<std::string, float>> mresults {nexperiments};

  // Get and log the OpenCL device name
  char deviceName [1024];
  clCheck(clGetDeviceInfo(devices[DEVICE_NUMBER], CL_DEVICE_NAME, 1024, deviceName, NULL));
  DEBUG << "Invoking kernels on your " << deviceName << std::endl;

  for (auto i=0; i<nexperiments; ++i) {
    // Update the number of threads in Y if more than 1 experiment
    if (nexperiments!=1) {
      global_work_size[1] = i+1;
      local_work_size[1] = i+1;

      DEBUG << i+1 << ": " << std::flush;
    }

    for (auto j=0; j<niterations; ++j) {
      // Initialize values to zero
      clInitializeValue<cl_bool>(commandQueue, dev_hit_used, acc_hits, false);
      clInitializeValue<cl_int>(commandQueue, dev_atomicsStorage, eventsToProcess * atomic_space, 0);
      clInitializeValue<cl_int>(commandQueue, dev_hit_candidates, 2 * acc_hits, -1);
      clInitializeValue<cl_int>(commandQueue, dev_hit_h2_candidates, 2 * acc_hits, -1);

      // Just for debugging
      clInitializeValue<cl_char>(commandQueue, dev_tracks, eventsToProcess * MAX_TRACKS * sizeof(Track), 0);
      clInitializeValue<cl_char>(commandQueue, dev_tracklets, acc_hits * sizeof(Track), 0);
      clInitializeValue<cl_int>(commandQueue, dev_tracks_to_follow, eventsToProcess * TTF_MODULO, 0);
      clCheck(clFinish(commandQueue));

      cl_event kernelEvent;

      clCheck(clEnqueueNDRangeKernel(commandQueue, kernel, work_dim, NULL, global_work_size, local_work_size, 0, NULL, &kernelEvent));
      // clCheck(clFinish(commandQueue));
      clCheck(clWaitForEvents(1 , &kernelEvent));

      // Start and end of event
      unsigned long tstart = 0;
      unsigned long tend = 0;
      clGetEventProfilingInfo(kernelEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong) , &tstart, NULL);
      clGetEventProfilingInfo(kernelEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &tend, NULL);
      clReleaseEvent(kernelEvent);

      // Compute the duration in nanoseconds
      unsigned long tduration = tend - tstart;

      // DEBUG << "Execution time (ms): " << tduration / 1000000.0 << std::endl;
      time_values[i].push_back(tduration / 1000000.0f);

      DEBUG << "." << std::flush;
    }
    DEBUG << std::endl;
  }

  // Step 11: Get results
  if (PRINT_SOLUTION) DEBUG << "Number of tracks found per event:" << std::endl << " ";
  clCheck(clEnqueueReadBuffer(commandQueue, dev_atomicsStorage, CL_TRUE, 0, eventsToProcess * atomic_space * sizeof(int), atomics, 0, NULL, NULL));
  for (int i=0; i<eventsToProcess; ++i){
    const int numberOfTracks = atomics[i];
    if (PRINT_SOLUTION) DEBUG << numberOfTracks << ", ";

    output[startingEvent + i].resize(numberOfTracks * sizeof(Track));
    if (numberOfTracks > 0) {
      clCheck(clEnqueueReadBuffer(commandQueue, dev_tracks, CL_TRUE, i * MAX_TRACKS * sizeof(Track), numberOfTracks * sizeof(Track), &(output[startingEvent + i])[0], 0, NULL, NULL));
    }
  }
  if (PRINT_SOLUTION) DEBUG << std::endl;

  if (PRINT_VERBOSE) {
    // Print solution of all events processed, to results
    for (int i=0; i<eventsToProcess; ++i) {

      // Calculate z to sensor map
      std::map<int, int> zhit_to_module;
      setHPointersFromInput((uint8_t*) &(*(input[startingEvent + i]))[0], input[startingEvent + i]->size(),
        h_no_sensors, h_no_hits, h_sensor_Zs, sensor_hits, h_hit_IDs, hits);
      if (logger::ll.verbosityLevel > 0){
        // map to convert from z of hit to module
        for(int j=0; j<h_no_sensors; ++j){
          const int z = h_sensor_Zs[j];
          zhit_to_module[z] = j;
        }
        // Some hits z may not correspond to a sensor's,
        // but be close enough
        for(int j=0; j<h_no_hits; ++j){
          const int z = (int) hits.Zs[j];
          if (zhit_to_module.find(z) == zhit_to_module.end()){
            const int sensor = findClosestModule(z, zhit_to_module);
            zhit_to_module[z] = sensor;
          }
        }
      }

      // Print to output file with event no.
      const int numberOfTracks = output[i].size() / sizeof(Track);
      Track* tracks_in_solution = (Track*) &(output[startingEvent + i])[0];
      std::ofstream outfile (std::string(RESULTS_FOLDER) + std::string("/") + toString(i) + std::string(".out"));
      DEBUG << "writing to: " << std::string(RESULTS_FOLDER) + std::string("/") + toString(i) + std::string(".out") << std::endl;
      for(int j=0; j<numberOfTracks; ++j){
        // TODO: PS: I commented out this code, because it would be too much
        //           hassle to make it work, but maybe this would be useful.
        //   outstream << "Track #" << j << ", length " << (int) tracks[j].hitsNum << std::endl;
        //printTrack(tracks_in_solution, j, zhit_to_module, hits, h_hit_IDs, outfile);
      }
      outfile.close();
    }
  }

  DEBUG << std::endl << "Time averages:" << std::endl;
  for (auto i=0; i<nexperiments; ++i){
    mresults[i] = calcResults(time_values[i]);
    DEBUG << " nthreads (" << NUMTHREADS_X << ", " << (nexperiments==1 ? local_work_size[1] : i+1) <<  "): " << mresults[i]["mean"]
      << " ms (std dev " << mresults[i]["deviation"] << ")" << std::endl;
  }

  // Step 12: Clean the resources
  clCheck(clReleaseKernel(kernel));
  clCheck(clReleaseProgram(program));
  clCheck(clReleaseCommandQueue(commandQueue));
  clCheck(clReleaseContext(context));

  clCheck(clReleaseMemObject(dev_tracks));
  clCheck(clReleaseMemObject(dev_tracklets));
  clCheck(clReleaseMemObject(dev_weak_tracks));
  clCheck(clReleaseMemObject(dev_tracks_to_follow));
  clCheck(clReleaseMemObject(dev_atomicsStorage));
  clCheck(clReleaseMemObject(dev_event_offsets));
  clCheck(clReleaseMemObject(dev_hit_offsets));
  clCheck(clReleaseMemObject(dev_hit_used));
  clCheck(clReleaseMemObject(dev_input));
  clCheck(clReleaseMemObject(dev_best_fits));
  clCheck(clReleaseMemObject(dev_hit_candidates));
  clCheck(clReleaseMemObject(dev_hit_h2_candidates));

  free(atomics);
  free(devices);
  free(hit_candidates);

  return 0;
}

const char *getErrorString (cl_int error) {
switch(error){
    // run-time and JIT compiler errors
    case 0: return "CL_SUCCESS";
    case -1: return "CL_DEVICE_NOT_FOUND";
    case -2: return "CL_DEVICE_NOT_AVAILABLE";
    case -3: return "CL_COMPILER_NOT_AVAILABLE";
    case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case -5: return "CL_OUT_OF_RESOURCES";
    case -6: return "CL_OUT_OF_HOST_MEMORY";
    case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case -8: return "CL_MEM_COPY_OVERLAP";
    case -9: return "CL_IMAGE_FORMAT_MISMATCH";
    case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case -11: return "CL_BUILD_PROGRAM_FAILURE";
    case -12: return "CL_MAP_FAILURE";
    case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
    case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
    case -15: return "CL_COMPILE_PROGRAM_FAILURE";
    case -16: return "CL_LINKER_NOT_AVAILABLE";
    case -17: return "CL_LINK_PROGRAM_FAILURE";
    case -18: return "CL_DEVICE_PARTITION_FAILED";
    case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

    // compile-time errors
    case -30: return "CL_INVALID_VALUE";
    case -31: return "CL_INVALID_DEVICE_TYPE";
    case -32: return "CL_INVALID_PLATFORM";
    case -33: return "CL_INVALID_DEVICE";
    case -34: return "CL_INVALID_CONTEXT";
    case -35: return "CL_INVALID_QUEUE_PROPERTIES";
    case -36: return "CL_INVALID_COMMAND_QUEUE";
    case -37: return "CL_INVALID_HOST_PTR";
    case -38: return "CL_INVALID_MEM_OBJECT";
    case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case -40: return "CL_INVALID_IMAGE_SIZE";
    case -41: return "CL_INVALID_SAMPLER";
    case -42: return "CL_INVALID_BINARY";
    case -43: return "CL_INVALID_BUILD_OPTIONS";
    case -44: return "CL_INVALID_PROGRAM";
    case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
    case -46: return "CL_INVALID_KERNEL_NAME";
    case -47: return "CL_INVALID_KERNEL_DEFINITION";
    case -48: return "CL_INVALID_KERNEL";
    case -49: return "CL_INVALID_ARG_INDEX";
    case -50: return "CL_INVALID_ARG_VALUE";
    case -51: return "CL_INVALID_ARG_SIZE";
    case -52: return "CL_INVALID_KERNEL_ARGS";
    case -53: return "CL_INVALID_WORK_DIMENSION";
    case -54: return "CL_INVALID_WORK_GROUP_SIZE";
    case -55: return "CL_INVALID_WORK_ITEM_SIZE";
    case -56: return "CL_INVALID_GLOBAL_OFFSET";
    case -57: return "CL_INVALID_EVENT_WAIT_LIST";
    case -58: return "CL_INVALID_EVENT";
    case -59: return "CL_INVALID_OPERATION";
    case -60: return "CL_INVALID_GL_OBJECT";
    case -61: return "CL_INVALID_BUFFER_SIZE";
    case -62: return "CL_INVALID_MIP_LEVEL";
    case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
    case -64: return "CL_INVALID_PROPERTY";
    case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
    case -66: return "CL_INVALID_COMPILER_OPTIONS";
    case -67: return "CL_INVALID_LINKER_OPTIONS";
    case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";

    // extension errors
    case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
    case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
    case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
    case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
    case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
    case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
    default: return "Unknown OpenCL error";
    }
}

void clChoosePlatform(cl_device_id*& devices, cl_platform_id& platform) {
  // Choose the first available platform
  // PS: removed, since never used
  // cl_platform_id* clPlatformIDs;
  cl_uint numPlatforms;
  clCheck(clGetPlatformIDs(0, NULL, &numPlatforms));
  if(numPlatforms > 0)
  {
    cl_platform_id* platforms = (cl_platform_id*) malloc(numPlatforms * sizeof(cl_platform_id));
    clCheck(clGetPlatformIDs(numPlatforms, platforms, NULL));
    platform = platforms[0];
    free(platforms);
  }

  // Choose a device from the platform according to DEVICE_PREFERENCE
  cl_uint numCpus = 0;
  cl_uint numGpus = 0;
  cl_uint numAccelerators = 0;
  clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 0, NULL, &numCpus);
  clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &numGpus);
  clGetDeviceIDs(platform, CL_DEVICE_TYPE_ACCELERATOR, 0, NULL, &numAccelerators);
  devices = (cl_device_id*) malloc(numAccelerators * sizeof(cl_device_id));

  DEBUG << std::endl << "Devices available: " << std::endl
    << "CPU: " << numCpus << std::endl
    << "GPU: " << numGpus << std::endl
    << "Accelerators: " << numAccelerators << std::endl;

  if (DEVICE_PREFERENCE == DEVICE_CPU && numCpus > 0) {
    DEBUG << "Choosing CPU" << std::endl;
    clCheck(clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, numCpus, devices, NULL));
  }
  else if (DEVICE_PREFERENCE == DEVICE_GPU && numGpus > 0) {
    DEBUG << "Choosing GPU" << std::endl;
    clCheck(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numGpus, devices, NULL));
  }
  else if (DEVICE_PREFERENCE == DEVICE_ACCELERATOR && numAccelerators > 0) {
    DEBUG << "Choosing accelerator" << std::endl;
    clCheck(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ACCELERATOR, numAccelerators, devices, NULL));
  }
  else {
    // We couldn't match the preference.
    // Let's try the first device that appears available.
    cl_uint numDevices = 0;
    clCheck(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices));
    if (numDevices > 0) {
      DEBUG << "Preference device couldn't be met" << std::endl
            << "Choosing an available OpenCL capable device" << std::endl;
      clCheck(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, numDevices, devices, NULL));
    }
    else {
      DEBUG << "No OpenCL capable device detected" << std::endl
            << "Check the drivers, OpenCL runtime or ICDs are available" << std::endl;
      exit(-1);
    }
  }
  DEBUG << std::endl;
}


template <class T>
void clInitializeValue(cl_command_queue& commandQueue, cl_mem& param, size_t size, T value) {
    T* temp;
    if (value == 0) temp = (T*) calloc(size, sizeof(T));
    else {
        temp = (T*) malloc(size * sizeof(T));
        for (int i=0; i<size; ++i) temp[i] = value;
    }

    clCheck(clEnqueueWriteBuffer(commandQueue, param, CL_TRUE, 0, size * sizeof(T), temp, 0, NULL, NULL));
    free(temp);
}


void printOutAllSensorHits(int* prevs, int* nexts){
  DEBUG << "All valid sensor hits: " << std::endl;
  for(int i=0; i<h_no_sensors; ++i){
    for(int j=0; j < sensor_hits.nums[i]; ++j){
      int hit = sensor_hits.starts[i] + j;

      if(nexts[hit] != -1){
        DEBUG << hit << ", " << nexts[hit] << std::endl;
      }
    }
  }
}

void printOutSensorHits(int sensorNumber, int* prevs, int* nexts){
  for(int i=0; i < sensor_hits.nums[sensorNumber]; ++i){
    int hstart = sensor_hits.starts[sensorNumber];

    DEBUG << hstart + i << ": " << prevs[hstart + i] << ", " << nexts[hstart + i] << std::endl;
  }
}

void printInfo(int numberOfSensors, int numberOfHits, const Hits& hits) {
  numberOfSensors = numberOfSensors>52 ? 52 : numberOfSensors; // 42 would fit better

  DEBUG << "Read info:" << std::endl
    << " no sensors: " << h_no_sensors << std::endl
    << " no hits: " << h_no_hits << std::endl
    << numberOfSensors << " sensors: " << std::endl;

  for (int i=0; i<numberOfSensors; ++i){
    DEBUG << " Zs: " << h_sensor_Zs[i] << std::endl
      << " hitStarts: " << sensor_hits.starts[i] << std::endl
      << " hitNums: " << sensor_hits.nums[i] << std::endl << std::endl;
  }

  DEBUG << numberOfHits << " hits: " << std::endl;

  for (int i=0; i<numberOfHits; ++i){
    DEBUG << " hit_id: " << h_hit_IDs[i] << std::endl
      << " hit_X: " << hits.Xs[i] << std::endl
      << " hit_Y: " << hits.Ys[i] << std::endl
      << " hit_Z: " << hits.Zs[i] << std::endl << std::endl;
  }
}

void getMaxNumberOfHits(char*& input, int& maxHits){
  int* l_no_sensors = (int*) &input[0];
  int* l_no_hits = (int*) (l_no_sensors + 1);
  int* l_sensor_Zs = (int*) (l_no_hits + 1);
  int* l_sensor_hitStarts = (int*) (l_sensor_Zs + l_no_sensors[0]);
  int* l_sensor_hitNums = (int*) (l_sensor_hitStarts + l_no_sensors[0]);

  maxHits = 0;
  for(int i=0; i<l_no_sensors[0]; ++i){
    if(l_sensor_hitNums[i] > maxHits)
      maxHits = l_sensor_hitNums[i];
  }
}
