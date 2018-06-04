#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <future>
#include <memory>
#include <stdlib.h>
#include <functional>
#include "Utils.h"
#include "FastaStructs.h"
#include "MultithreadedFastaReader.h"
#include "md5.h"
#include "LSH.h"

using namespace fer::zesoi::bioinfo;

int main(int argc, char** argv) {
  // Input file path
  const auto inFile = argv[1];
  // Thread pool size
  const auto numIOThreads = atoi(argv[2]);
  std::cout << "Reading input file: " << inFile
            << ", size " << Utils::getFileSizeGB(inFile) << " GB" 
            << std::endl;
  std::cout << "Using " << numIOThreads << " io threads" << std::endl;

  // K-gram
  const auto K = atoi(argv[3]);
  // Window size
  const auto W = atoi(argv[4]);
  // HashType
  const auto hashTypeInt = atoi(argv[5]);
  const auto hashType = LSH::intToHashType(hashTypeInt);
  // MOD value
  // atoi(argv[5])
  auto parseCb = [K, W, hashType](const FastaHeader& fh) {
    const auto& features = LSH::getMinHash(
      fh.getData(),
      K, 
      W, 
      hashType); 
    return FeaturizedSample(fh.getID(), "", features); 
  }; 
  auto fileReader = make_unique<MultithreadedFastaReader<FeaturizedSample>>(
    inFile, 
    numIOThreads, 
    parseCb
  );
  auto futureResults = fileReader->read();
 
  return 0;
}
