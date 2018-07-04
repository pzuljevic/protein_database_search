#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <future>
#include <memory>
#include <mutex>
#include <stdlib.h>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <functional>
#include <math.h>
#include "utils/Utils.h"
#include "fasta/FastaStructs.h"
#include "fasta/MultithreadedFastaReader.h"
#include "hash/md5.h"
#include "hash/LSH.h"
#include "structs/Structs.h"
#include "algorithm/Clustering.h"

#define OFFSET_LC 1000
#define READ_HEADER_SIZE 11

using namespace fer::zesoi::bioinfo;

int main(int argc, char** argv) {
  if (argc < 11) {
    std::cout << "Run: ./binary "
              << "<inFile> <numThreads> <K> <W> <hashType> <modulo> "
              << "<numLengthClusters> <looseThreshold> <tightThreshold> " 
              << "<outFile> <queryInFile>"
              << std::endl;
    exit(1);
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
  std::cout << "K=" << K << " ";
  // Window size
  const auto W = atoi(argv[4]);
  std::cout << "W=" << W << " ";
  // HashType
  const auto hashTypeInt = atoi(argv[5]);
  const auto hashType = LSH::intToHashType(hashTypeInt);
  // modulo value
  const auto modulo = atoi(argv[6]);
  std::cout << "modulo=" << modulo << std::endl;
  const auto numLengthClusters = atoi(argv[7]);
  
  const auto looseThreshold = atof(argv[8]);
  const auto tightThreshold = atof(argv[9]);
  const auto outFile = argv[10];
  const auto queryInFile = argc == 12? argv[11] : "";

  // Thread pool
  auto threadPool = std::make_shared<ThreadPool>(numThreads);
  // Define how do we convert protein data to feature vector
  // std::vector<int64_t> a = LSH::getCoeff(W);
  // std::vector<int64_t> b = LSH::getCoeff(W);
  auto parseCb = [K, W, hashType, modulo](const FastaHeader& fh) {
   const auto& features = LSH::getMinHash(
     fh.getData(), K, W, hashType, modulo); 
    return std::make_shared<DataSample>(
      fh.getHeader(),
      features,
      fh.getData(),
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
    numLengthClusters
  );
  // 2) Divide data to threads for tight clustering
  auto jaccardClustering = JaccardClustering(
    samples, 
    tightThreshold, 
    looseThreshold);
  threadPool = std::make_shared<ThreadPool>(numLengthClusters);
  int64_t offset = 0;

  std::unordered_map<int64_t, std::vector<int64_t>> tightClustersGlobal; 
  std::unordered_map<int64_t, std::vector<int64_t>> looseClustersGlobal; 

  std::map<std::string, std::vector<int64_t>> representClustersMap;
  for (const auto& it : lengthClusters) {
    threadPool->enqueue([
        &it, 
        &jaccardClustering, 
        &tightClustersGlobal,
        &looseClustersGlobal,
        &lck, 
        offset, 
        &samples,
        &representClustersMap] {      
      {
        std::lock_guard<std::mutex> guard(lck); 
        std::cout << "Length bucket id: " <<  it.first 
                  << " total number of samples:" 
                  << it.second.size() << std::endl;
      }

      std::unordered_map<int64_t, std::vector<int64_t>> tightClusters; 
      std::unordered_map<int64_t, std::vector<int64_t>> looseClusters; 
      jaccardClustering.run(
        offset,
        it.second,     // vector::sample IDs vector
        tightClusters, // map::<tight id, <list of sample IDs>>
        looseClusters  // map::<loose id, <lost of tight IDs>>
      );

      {
        std::lock_guard<std::mutex> guard(lck); 
        for (const auto it : looseClusters) {
          int64_t id = tightClusters[it.second[0]][0];
          const auto& tag = samples[id]->getHeader().substr(0, READ_HEADER_SIZE);
          representClustersMap[tag] = samples[id]->getFeatures(); 
          looseClustersGlobal[it.first] = it.second;
        }
        for (const auto it : tightClusters) {
          tightClustersGlobal[it.first] = it.second;
        }
        std::cout << "Total tightClusters: " << tightClusters.size() 
                  << ", total loose clusters:" << looseClusters.size() 
                  << " among " << it.second.size() 
                  << " samples" << std::endl;
      }
    });   
    offset += OFFSET_LC;
  }
  threadPool->stop();
  std::set<int64_t> tc; 
  std::set<int64_t> lc; 
  std::map<int64_t, std::set<std::string>> clusters;

  for (auto sample : samples) {
    for (auto c : sample->getTightClusters()) { 
    tc.insert(c + sample->getOffset());
    }
    for (auto c : sample->getLooseClusters()) { 
      lc.insert(c + sample->getOffset());
      clusters[c + sample->getOffset()].insert(
        sample->getHeader().substr(0, READ_HEADER_SIZE)); 
    }
  }

  std::ofstream outputFile;
  outputFile.open(outFile);
  for (const auto it1 : representClustersMap) {
    outputFile << it1.first << " ";
  }
  outputFile << std::endl;
  auto distanceMetricFactory = DistanceMetricFactory(); 
  const auto distanceMetric = 
  distanceMetricFactory.getMetric(DistanceMetricType::JACCARD);
  for (const auto it1 : representClustersMap) {
    for (const auto it2 : representClustersMap) {
      const auto d = 1 - distanceMetric->getSimilarity(
        it1.second, 
        it2.second);
      outputFile << d << " "; 
    }
    outputFile << std::endl;
  } 
  outputFile.close();
  std::cout << "Total number of loose clusters: " << lc.size() 
            << ", total number of tight clusters: " << tc.size() 
            << ", total number of samples: " << samples.size()
            << std::endl;
  if (argc == 12) {
    std::cout << "Reading query input file: " << queryInFile << std::endl;
    auto poaThreadPool = std::make_shared<ThreadPool>(1);
    auto poaFileReader = make_unique<MultithreadedFastaReader<DataSamplePtr>>(
      queryInFile, 
      poaThreadPool,
      1,
      parseCb
    );
    auto poaFutureResults = poaFileReader->read();
    // Collect data from the threads 
    std::vector<DataSamplePtr> poaSamples; 
    for (auto& future : *poaFutureResults) {
      poaThreadPool->enqueue([&future, &poaSamples, &lck] {
        const auto& result = future.get();
        std::lock_guard<std::mutex> guard(lck); 
        poaSamples.insert(poaSamples.end(), result.begin(), result.end()); 
      });
    } 

    poaThreadPool->stop();

    // For online phase read each sample from the input file and run POA against
    // it on loose clusters with other script.
    // This will generate log files for each sample for query input file.
    const int64_t numTC = 20; 
    int64_t bestLooseId = 0; 
    for (const auto& sample : poaSamples) {
      double maxSim = 0.0;
      for (auto& it : looseClustersGlobal) {
        const auto looseID = it.first; 
        for (auto& tightID : looseClustersGlobal[looseID]) {
          const auto d = distanceMetric->getSimilarity(
            samples[tightClustersGlobal[tightID][0]]->getFeatures(), 
            sample->getFeatures()
          );
          if (d > maxSim) {
            maxSim = d;
            bestLooseId = looseID;
          }
        }
      }
      std::ofstream queryOutputFile;
      queryOutputFile.open(sample->getHeader().substr(0, 15) + ".poa");
      queryOutputFile << ">" << sample->getHeader() << std::endl;
      queryOutputFile << sample->getData() << std::endl;
      int cntTC = 0;
      for (auto& tightID : looseClustersGlobal[bestLooseId]) {
        if (cntTC++ > numTC) break;
        queryOutputFile << ">" << samples[tightClustersGlobal[tightID][0]]->getHeader() << std::endl;
        queryOutputFile << samples[tightClustersGlobal[tightID][0]]->getData() << std::endl;
      }
      queryOutputFile.close();
    }
  }
  return 0;
}
