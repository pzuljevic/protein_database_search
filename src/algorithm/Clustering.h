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

    for (const auto i : sampleIDs) {
      if (++k % 10000 == 0) { 
        std::cout << sampleIDs.size() << " processing status: " 
                  <<  (double) k / sampleIDs.size() 
                  << " loose cluster count" << looseClusters.size() 
                  << " tight cluster count: " << tightClusters.size() << std::endl; 
      }
      if (looseClusters.size() - previousSizeLooseCluster > sqrt(++lcThresholdUpdateCnt) 
          && lcThresholdUpdateCnt > sqrt(samples_.size())) {
        std::cout << "Relaxing threshold " << sampleIDs.size()
                  << " " << looseThreshold_ << " -> " 
                  << looseThreshold_ * 0.9 << std::endl;
        looseThreshold_ *= 0.8; 
        previousSizeLooseCluster = looseClusters.size();
        lcThresholdUpdateCnt = 0;
      } 
      std::vector<int64_t> looseIDs;
      auto bestLooseID = 0; 
      auto maxLooseSim = 0;
      int cntLC = 0;
      for (auto& it : looseClusters) {
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

      auto bestTightID = -1; 
      double maxTightSim = 0;
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

      cntLC = 0;
      for (auto& looseID : looseIDs) {
        cntLC++;
        int cntTC = 0;
        for (auto& tightID : looseClusters[looseID]) {
          cntTC++;
          const auto d = distanceMetric->getSimilarity(
            samples_[tightClusters[tightID][0]]->getFeatures(), 
            samples_[i]->getFeatures()
          ); 
          if (d > maxTightSim) {
            maxTightSim = d;
            bestTightID = tightID;
          }
        }
      }

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
    } 
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
