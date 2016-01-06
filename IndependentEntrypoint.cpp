
/**
 *      Autocontained cross-platform independent entrypoint
 *      for Gaudi offloaded algorithms.
 *
 *      Intended for testing and debugging, completely decoupled
 *      from the Gaudi framework.
 *
 *      author  -   Daniel Campora
 *      email   -   dcampora@cern.ch
 *
 *      June, 2014
 *      CERN
 */

#include <iostream>
#include <string>
#include <cstring>
#include <exception>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "PixelSearchByTriplet.h"



void printUsage(char* argv[]){
    std::cerr << "Usage: "
        << argv[0]
        << " -serial|-ocl"
        << " -bin|-tex"
        << " <comma separated input filenames>"
        << std::endl;
}

/**
 * Generic StrException launcher
 */
class StrException : public std::exception
{
public:
    std::string s;
    StrException(std::string ss) : s(ss) {}
    ~StrException() throw () {} // Updated
    const char* what() const throw() { return s.c_str(); }
};

/**
 * Checks file existence
 * @param  name
 * @return
 */
bool fileExists (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

/**
 * Reads some data from an input file, following the
 * specified format of choice
 *
 * Format expected in file:
 *
 * int funcNameLen
 * char* funcName
 * int dataSize
 * char* data
 */
void readFileIntoVector(std::string filename, std::vector<unsigned char> & output){
    // Check if file exists
    if (!fileExists(filename)){
        throw StrException("Error: File " + filename + " does not exist.");
    }

    std::ifstream infile (filename.c_str(), std::ifstream::binary);

    // get size of file
    // PS: removed, since never used
    // infile.seekg(0, std::ifstream::end);
    // int size = infile.tellg();
    // infile.seekg(0);

    // Read format expected:
    //  int funcNameLen
    //  char* funcName
    //  int dataSize
    //  char* data
    int funcNameLen;
    int dataSize;
    std::vector<char> funcName;

    char* pFuncNameLen = (char*) &funcNameLen;
    char* pDataSize = (char*) &dataSize;
    infile.read(pFuncNameLen, sizeof(int));

    funcName.resize(funcNameLen);
    infile.read(&(funcName[0]), funcNameLen);
    infile.read(pDataSize, sizeof(int));

    // read content of infile with a vector
    output.resize(dataSize);
    infile.read ((char*) &(output[0]), dataSize);
    infile.close();
}

/**
 * This is if the function is called on its own
 * (ie. non-gaudi execution)
 *
 * In that case, the file input is expected.
 * As a convention, multiple files would be specified
 * with comma-separated values
 *
 * @param  argc
 * @param  argv
 * @return
 */
int main(int argc, char *argv[])
{
    std::string filename_str, mode_opt, out_opt;

    // PS: removed, since never used
    // int fileNumber = 1;
    std::string delimiter = ",";
    std::vector<std::vector<unsigned char> > input;
    std::vector<std::string> filenames;
    // Get params (getopt independent)
    if (argc != 4){
        printUsage(argv);
        return 0;
    }

    mode_opt = std::string(argv[1]);
    std::transform(mode_opt.begin(), mode_opt.end(), mode_opt.begin(), ::tolower);
    out_opt = std::string(argv[2]);
    filename_str = std::string(argv[3]);

    OutType outtype;
    if (out_opt.compare("-bin") == 0) {
        outtype = OutType::Binary;
    } else if (out_opt.compare("-tex") == 0) {
        outtype = OutType::Text;
    } else {
        std::cout << "Output type " << out_opt << " not known." << std::endl;
    }

    // Check how many files were specified and
    // call the entrypoint with the suggested format
    if(filename_str.empty()){
        std::cerr << "No filename specified" << std::endl;
        printUsage(argv);
        return -1;
    }

    size_t numberOfOcurrences = std::count(filename_str.begin(), filename_str.end(), ',') + 1;
    input.resize(numberOfOcurrences);
    filenames.resize(numberOfOcurrences);

    size_t posFound = filename_str.find(delimiter);
    if (posFound != std::string::npos){
        int input_index = 0;
        size_t prevFound = 0;
        while(prevFound != std::string::npos){
            filenames[input_index] = filename_str.substr(prevFound, posFound-prevFound);
            if (posFound == std::string::npos){
                readFileIntoVector(filenames[input_index], input[input_index]);
                prevFound = posFound;
            }
            else {
                readFileIntoVector(filenames[input_index], input[input_index]);
                prevFound = posFound + 1;
                posFound = filename_str.find(delimiter, posFound + 1);
            }
            input_index++;
        }
    }
    else {
        filenames[0] = filename_str;
        readFileIntoVector(filename_str, input[0]);
    }

    // Print out first byte from formatter->inputPointer
    std::cout << input.size() << " files read" << std::endl;

    // Call offloaded algo
    std::vector<std::vector<unsigned char> > output;
    if (mode_opt.compare("-serial") == 0) {
        independent_execute(input, output, filenames, ExecMode::Serial, outtype);
    } else if (mode_opt.compare("-ocl") == 0) {
#ifdef WITH_OPENCL
        independent_execute(input, output, filenames, ExecMode::OpenCl, outtype);
#else
        std::cout << "Please recompile with -DBUILD_OPENCL=ON" << std::endl;
#endif
    } else {
        std::cout << "Execution mode " << mode_opt << " not yet supported" << std::endl;
    }

    // Post execution entrypoint
    independent_post_execute(output);

    return 0;
}
