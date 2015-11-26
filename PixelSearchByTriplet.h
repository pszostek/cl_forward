
#ifndef PIXELSEARCHBYTRIPLET
#define PIXELSEARCHBYTRIPLET 1

#include "FileStdLogger.h"
#include "Tools.h"
#include "GpuKernelInvoker.h"
#include "Logger.h"

#include <stdint.h>

enum class ExecMode {Serial, TBB, OpenCl};

int independent_execute(
    const std::vector<std::vector<uint8_t> > & input,
    std::vector<std::vector<uint8_t> > & output,
    ExecMode mode, OutType outtype);

void independent_post_execute(const std::vector<std::vector<uint8_t> > & output);

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

/**
 * Common entrypoint for Gaudi and non-Gaudi
 * @param input
 * @param output
 */
int cpuPixelSearchByTripletSerialRun(
        const std::vector<const std::vector<uint8_t>* > & input,
        std::vector<std::vector<uint8_t> > & output, OutType outtype);

#endif
