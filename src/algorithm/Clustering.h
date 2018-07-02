#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <algorithm>
#include <unordered_map>
#include <math.h>
#include <chrono>
#include "../utils/DistanceMetric.h"
#include "../structs/Structs.h"

namespace fer {
namespace zesoi {
namespace bioinfo {

class LengthClustering {
 public:
  // Returns a map where the key is cluster ID and 
  // the value is the list of sample IDs
  static std::unordered_map<int64_t, std::vector<int64_t>> run(
      std::vector<DataSamplePtr>& samples, 
      int64_t numBuckets) {
    std::unordered_map<int64_t, std::vector<int64_t>> lengthClusters;
    for (int64_t i = 0; i < (int64_t)samples.size(); ++i) {
      const auto looseID = (samples[i]->getLength() / 50) % numBuckets;
      lengthClusters[looseID].push_back(i);
    } 
    return lengthClusters;
  }

 protected:
  LengthClustering() {}
}; 

class JaccardClustering { 
 public:
  JaccardClustering(
      std::vector<DataSamplePtr>& samples, 
      double tightThreshold, 
      double looseThreshold) : 
    samples_(samples),
    tightThreshold_(tightThreshold),
    looseThreshold_(looseThreshold) {
    std::cout << "Initialized Jaccard based clustering with thresholds: "
              << "tight=" << tightThreshold
              << " loose=" << looseThreshold
              << std::endl;
  }

  // Will set loose and tight clusters
  void run(
      int64_t offset,
      std::vector<int64_t> sampleIDs,
      std::unordered_map<int64_t, std::vector<int64_t>>& 
      tightClusters,   // Contains the lists of sample IDs
      std::unordered_map<int64_t, std::vector<int64_t>>& 
      looseClusters) { // Contains the lists od tight cluster IDs

    // Use Jaccard distance
    auto distanceMetricFactory = DistanceMetricFactory(); 
    const auto distanceMetric = 
    distanceMetricFactory.getMetric(DistanceMetricType::JACCARD);
    int k = 0;
    int lcThresholdUpdateCnt = 0;
    int previousSizeLooseCluster = 0;
    int tcThresholdUpdateCnt = 0;
    int previousSizeTightCluster = 0;

    for (const auto i : sampleIDs) {
      auto start = std::chrono::high_resolution_clock::now();
      if (++k % 1000 == 0) 
        std::cout << sampleIDs.size() << " status: " 
                  <<  (double) k / sampleIDs.size() 
                  << " " << looseClusters.size() 
                  << " tight: " << tightClusters.size() << std::endl;
      if (looseClusters.size() - previousSizeLooseCluster > sqrt(++lcThresholdUpdateCnt) 
          && lcThresholdUpdateCnt > 1000) {
        std::cout << "Relaxing threshold " << sampleIDs.size()
                  << " " << looseThreshold_ << " -> " 
                  << looseThreshold_ * 0.9 << std::endl;
        looseThreshold_ *= 0.9; 
        previousSizeLooseCluster = looseClusters.size();
        lcThresholdUpdateCnt = 0;
      } 
      std::vector<int64_t> looseIDs;
      auto bestLooseID = 0; 
      auto maxLooseSim = 0;
      int br = 0;
      for (auto& it : looseClusters) {
        if (br++ > 10) break;
        const auto looseID = it.first; 
        const auto tightID = it.second[0];
        const auto d = distanceMetric->getSimilarity(
          samples_[tightClusters[tightID][0]]->getFeatures(), 
          samples_[i]->getFeatures()
        ); 
        if (d > maxLooseSim && d > looseThreshold_) {
          maxLooseSim = d;
          bestLooseID = looseID;
          looseIDs.push_back(bestLooseID);
        }
      }       
      auto t1 = std::chrono::high_resolution_clock::now();

      auto bestTightID = -1; 
      auto maxTightSim = 0;
      if (looseIDs.size() == 0) {
        // Create loose cluster
        const auto looseID = looseClusters.size();
        looseIDs.push_back(looseID);
        bestLooseID = looseID;
        maxLooseSim = 1.00;
        // Create tight cluster
        const auto tightID = tightClusters.size(); 
        looseClusters[looseID].push_back(tightID);
        tightClusters[tightID].push_back(i);
        bestTightID = tightID;
        maxTightSim = 1.00;
      }

      // Find the best tight and loose cluster
      int cnt2 = 0;
      for (auto& looseID : looseIDs) {
        cnt2++;
        int cnt = 0;
        for (auto& tightID : looseClusters[looseID]) {
          cnt++;
          const auto d = distanceMetric->getSimilarity(
            samples_[tightClusters[tightID][0]]->getFeatures(), 
            samples_[i]->getFeatures()
          ); 
          if (d > maxTightSim) {
            maxTightSim = d;
            bestTightID = tightID;
          }
          if (cnt > 10) break;
        }
        if (cnt2 > 10) break;
      }
      auto t2 = std::chrono::high_resolution_clock::now();

      // Tight cluster is not found
      if (bestTightID == -1 || maxTightSim < tightThreshold_) {
        // Create tight cluster
        const auto tightID = tightClusters.size(); 
        looseClusters[bestLooseID].push_back(tightID);
        tightClusters[tightID].push_back(i);
        bestTightID = tightID;
      } 

      samples_[i]->addLooseCluster(bestLooseID);
      samples_[i]->addTightCluster(bestTightID);
      samples_[i]->setOffset(offset);

      if (looseClusters[bestLooseID].size() - previousSizeTightCluster > sqrt(tcThresholdUpdateCnt) 
          && tcThresholdUpdateCnt > 100) {
        tightThreshold_ *= 0.9;
        tcThresholdUpdateCnt = 0;
        previousSizeTightCluster = looseClusters[bestLooseID].size();
        std::cout << "Relaxing threshold " << sampleIDs.size()
                  << " " << tightThreshold_ << " -> " 
                  << tightThreshold_ * 0.9 << std::endl;
      }

      auto finish = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = finish - start;
      std::chrono::duration<double> t1e = t1 - start;
      std::chrono::duration<double> t2e = t2 - t1;

      if (k % 1000 == 0 && elapsed.count() > 0.1)  
        std::cout << "Size: " << sampleIDs.size() 
                  << " t:" << elapsed.count() 
                  << " time, t1=" << t1e.count() 
                  << ", t2=" << t2e.count() 
                  << " loose: " << looseClusters.size() 
                  << " tight: " << tightClusters.size() 
                  << " total: " << sampleIDs.size() << std::endl;
    } 
    std::cout << sampleIDs.size() << " done " << std::endl;
  }

 protected:
  std::vector<DataSamplePtr> samples_; 
  double tightThreshold_;
  double looseThreshold_;
};

}
}
}
#endif
