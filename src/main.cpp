#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <future>
#include <memory>
#include <stdlib.h>
#include <functional>
#include "utils/Utils.h"
#include "fasta/FastaStructs.h"
#include "fasta/MultithreadedFastaReader.h"
#include "hash/md5.h"
#include "hash/LSH.h"

using namespace fer::zesoi::bioinfo;

int main(int argc, char** argv) {
  // Input file path
  const auto inFile = argv[1];
  // Thread pool size
  const auto numThreads = atoi(argv[2]);
  std::cout << "Reading input file: " << inFile
            << ", size " << Utils::getFileSizeGB(inFile) << " GB" 
            << std::endl;
  std::cout << "Using " << numThreads << " threads" << std::endl;

  // K-gram
  const auto K = atoi(argv[3]);
  // Window size
  const auto W = atoi(argv[4]);
  // HashType
  const auto hashTypeInt = atoi(argv[5]);
  const auto hashType = LSH::intToHashType(hashTypeInt);
  // MOD value
  // atoi(argv[5])

  // Thread pool
  auto threadPool = std::make_shared<ThreadPool>(numThreads);
  // Define how do we convert protein data to feature vector
  auto parseCb = [K, W, hashType](const FastaHeader& fh) {
    const auto& features = LSH::getMinHash(
      fh.getData(),
      K, 
      W, 
      hashType); 
    return FeaturizedSample(fh.getID(), "", features); 
  }; 

  // **** Data processing pipeline **** 
  // Kick in multithreaded file reader and vectorization
  auto fileReader = make_unique<MultithreadedFastaReader<FeaturizedSample>>(
    inFile, 
    threadPool,
    numThreads,
    parseCb
  );
  auto futureResults = fileReader->read();
  // Collect data from the threads
//  for (const auto future : futureResults) {
 // } 

  threadPool->stop();
 
  return 0;
}
