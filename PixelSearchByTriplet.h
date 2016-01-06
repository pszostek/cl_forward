
#ifndef PIXELSEARCHBYTRIPLET
#define PIXELSEARCHBYTRIPLET 1

#include "FileStdLogger.h"
#include "Tools.h"

#ifdef WITH_OPENCL
#include "GpuKernelInvoker.h"
#endif

#include "Logger.h"

#include <stdint.h>

enum class ExecMode {Serial, TBB, OpenCl};

int independent_execute(
    const std::vector<std::vector<uint8_t> > & input,
    std::vector<std::vector<uint8_t> > & output,
    std::vector<std::string> &filenames,
    ExecMode mode, OutType outtype);

void independent_post_execute(const std::vector<std::vector<uint8_t> > & output);

#ifdef WITH_OPENCL
int gpuPixelSearchByTriplet(
    const std::vector<const std::vector<uint8_t>* > & input,
    std::vector<std::vector<uint8_t> > & output);

/**
 * Common entrypoint for Gaudi and non-Gaudi
 * @param input
 * @param output
 */
int gpuPixelSearchByTripletInvocation(
    const std::vector<const std::vector<uint8_t>* > & input,
    std::vector<std::vector<uint8_t> > & output);
#endif

/**
 * Common entrypoint for Gaudi and non-Gaudi
 * @param input
 * @param output
 */
int cpuPixelSearchByTripletSerialRun(
        const std::vector<const std::vector<uint8_t>* > & input,
        std::vector<std::vector<uint8_t> > & output,
        std::vector<std::string> &filenames, OutType outtype);

#endif
