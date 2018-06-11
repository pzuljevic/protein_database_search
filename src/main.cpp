#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <future>
#include <memory>
#include <mutex>
#include <stdlib.h>
#include <set>
#include <unordered_map>
#include <functional>
#include "utils/Utils.h"
#include "fasta/FastaStructs.h"
#include "fasta/MultithreadedFastaReader.h"
#include "hash/md5.h"
#include "hash/LSH.h"
#include "structs/Structs.h"
#include "algorithm/Clustering.h"

using namespace fer::zesoi::bioinfo;

int main(int argc, char** argv) {
  if (argc != 10) {
    std::cout << "Run: ./binary "
              << "<inFile> <numThreads> <K> <W> <hashType> <MOD> "
              << "<numLengthClusters> <looseThreshold> <tightThreshold>" 
              << std::endl;
  }
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
  const auto MOD = atoi(argv[6]);
  const auto numLengthClusters = atoi(argv[7]);
  
  const auto looseThreshold = atof(argv[8]);
  const auto tightThreshold = atof(argv[9]);

  // Thread pool
  auto threadPool = std::make_shared<ThreadPool>(numThreads);
  // Define how do we convert protein data to feature vector
  auto parseCb = [K, W, hashType, MOD](const FastaHeader& fh) {
    const auto& features = LSH::getMinHash(
      fh.getData(),
      K, 
      W, 
      hashType,
      MOD); 
    return std::make_shared<DataSample>(
      fh.getHeader(),
      features, 
      fh.getData().size()
    ); 
  }; 

  // **** Data processing pipeline **** 
  // Kick in multithreaded file reader and vectorization
  auto fileReader = make_unique<MultithreadedFastaReader<DataSamplePtr>>(
    inFile, 
    threadPool,
    numThreads,
    parseCb
  );
  auto futureResults = fileReader->read();

  // Collect data from the threads 
  std::vector<DataSamplePtr> samples; 
  std::mutex lck;
  for (auto& future : *futureResults) {
    threadPool->enqueue([&future, &samples, &lck] {
      const auto& result = future.get();
      std::lock_guard<std::mutex> guard(lck); 
      samples.insert(samples.end(), result.begin(), result.end()); 
    });
  } 
  threadPool->stop();
  threadPool.reset();
  std::cout << "Done loading " << samples.size() 
            << " samples into memory" << std::endl;

  // **** Clustering pipeline **** 
  // 1) Kick in clustering per length 
  const auto& lengthClusters = LengthClustering::run(
    samples, 
    numLengthClusters,
    25, // tight step
    50 // loose step
  );
  // 2) Divide data to threads for tight clustering
  auto jaccardClustering = JaccardClustering(samples, tightThreshold, looseThreshold);
  threadPool = std::make_shared<ThreadPool>(numThreads);
  for (const auto& it : lengthClusters) {
    threadPool->enqueue([&it, &jaccardClustering, &lck] {      
      {
        std::lock_guard<std::mutex> guard(lck); 
        std::cout << "Length bucket id: " <<  it.first 
                  << " total number of samples:" 
                  << it.second.size() << std::endl;
      }

      std::unordered_map<int64_t, std::vector<int64_t>> tightClusters; 
      std::unordered_map<int64_t, std::vector<int64_t>> looseClusters; 
      jaccardClustering.run(
        it.second,     // vector::sample IDs vector
        tightClusters, // map::<tight id, <list of sample IDs>>
        looseClusters  // map::<loose id, <lost of tight IDs>>
      );

      {
        std::lock_guard<std::mutex> guard(lck); 
        std::cout << "Total tightClusters: " << tightClusters.size() 
                  << ", total loose clusters:" << looseClusters.size() 
                  << " among " << it.second.size() 
                  << " samples" << std::endl;
      }
    });   
  }
  threadPool->stop();
  
  for (auto sample : samples) {
    sample->print();
  }
  return 0;
}
