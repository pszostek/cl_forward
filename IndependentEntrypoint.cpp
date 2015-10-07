
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
    std::string filename, mode_opt;
    // PS: removed, since never used
    // int fileNumber = 1;
    std::string delimiter = ",";
    std::vector<std::vector<unsigned char> > input;

    // Get params (getopt independent)
    if (argc != 3){
        printUsage(argv);
        return 0;
    }

    mode_opt = std::string(argv[1]);
    filename = std::string(argv[2]);

    // Check how many files were specified and
    // call the entrypoint with the suggested format
    if(filename.empty()){
        std::cerr << "No filename specified" << std::endl;
        printUsage(argv);
        return -1;
    }

    size_t numberOfOcurrences = std::count(filename.begin(), filename.end(), ',') + 1;
    input.resize(numberOfOcurrences);

    size_t posFound = filename.find(delimiter);
    if (posFound != std::string::npos){
        int input_index = 0;
        size_t prevFound = 0;
        while(prevFound != std::string::npos){
            if (posFound == std::string::npos){
                readFileIntoVector(filename.substr(prevFound, posFound-prevFound), input[input_index]);
                prevFound = posFound;
            }
            else {
                readFileIntoVector(filename.substr(prevFound, posFound-prevFound), input[input_index]);
                prevFound = posFound + 1;
                posFound = filename.find(delimiter, posFound + 1);
            }
            input_index++;
        }
    }
    else {
        readFileIntoVector(filename, input[0]);
    }

    // Print out first byte from formatter->inputPointer
    std::cout << input.size() << " files read" << std::endl;

    // Call offloaded algo
    std::vector<std::vector<unsigned char> > output;
    if (mode_opt.compare("-ocl") == 0) {
        independent_execute(input, output, ExecMode::OpenCl);
    } else if (mode_opt.compare("-serial") == 0) {
        independent_execute(input, output, ExecMode::Serial);
    } else {
        std::cout << "Execution mode " << mode_opt << " not yet supported" << std::endl;
    }

    // Post execution entrypoint
    independent_post_execute(output);

    return 0;
}
