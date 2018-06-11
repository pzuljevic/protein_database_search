#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <algorithm>
#include <unordered_map>
#include <math.h>
#include "../utils/DistanceMetric.h"
#include "../structs/Structs.h"

namespace fer {
namespace zesoi {
namespace bioinfo {

#define TIGHT_STEP 25
#define LOOSE_STEP 50
#define MAX_L 25*1000
#define MAX_TIGHT_BUCKET MAX_L / TIGHT_STEP
#define MAX_LOOSE_BUCKET MAX_L / LOOSE_STEP

class LengthClustering {
 public:
  static int64_t getThightCluster(
      int64_t l, 
      int64_t tStep = TIGHT_STEP,
      int64_t maxTightBucket = MAX_TIGHT_BUCKET) {
    return std::min(l / tStep, maxTightBucket);
  }

  static int64_t getLooseCluster(
      int64_t l, 
      int64_t lStep = LOOSE_STEP,
      int64_t maxLooseBucket = MAX_LOOSE_BUCKET) {
    return std::min(l / lStep, maxLooseBucket);
  }
  
  // Returns a map where the key is cluster ID and 
  // the value is the list of sample IDs
  static std::unordered_map<int64_t, std::vector<int64_t>> run(
      std::vector<DataSamplePtr>& samples, 
      int64_t numBuckets,
      int64_t tightStep = TIGHT_STEP, 
      int64_t looseStep = LOOSE_STEP) { 
    std::unordered_map<int64_t, std::vector<int64_t>> lengthClusters;
    for (int64_t i = 0; i < samples.size(); ++i) {
      const auto tightID = 
        LengthClustering::getThightCluster(samples[i]->getLength());
      const auto looseID = 
        LengthClustering::getLooseCluster(samples[i]->getLength());
      lengthClusters[tightID % numBuckets].push_back(i);
      if ((looseID % numBuckets) != (tightID % numBuckets)) {
      //  lengthClusters[looseID % numBuckets].push_back(i);
      }
      samples[i]->setLengthClusterTag(tightID, looseID);
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
    for (const auto i : sampleIDs) {
      std::vector<int64_t> looseIDs;
      auto bestLooseID = 0; 
      auto maxLooseSim = 0;
      for (auto& it : looseClusters) {
        const auto looseID = it.first; 
        const auto tightID = it.second[0];
        const auto d = distanceMetric->getSimilarity(
          samples_[tightClusters[tightID][0]]->getFeatures(), 
          samples_[i]->getFeatures()
        ); 
        if (d > looseThreshold_) {
          looseIDs.push_back(looseID);
          if (d > maxLooseSim) {
            maxLooseSim = d;
            bestLooseID = looseID;
          }
        }
      }       

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
      for (auto& looseID : looseIDs) {
        for (auto& tightID : looseClusters[looseID]) {
          const auto d = distanceMetric->getSimilarity(
            samples_[tightClusters[tightID][0]]->getFeatures(), 
            samples_[i]->getFeatures()
          ); 
          if (d > tightThreshold_ && d > maxTightSim) {
            maxTightSim = d;
            bestTightID = tightID;
          }
        }
      }

      // Tight cluster is not found
      if (bestTightID == -1) {
        // Create tight cluster
        const auto tightID = tightClusters.size(); 
        looseClusters[bestLooseID].push_back(tightID);
        tightClusters[tightID].push_back(i);
        bestTightID = tightID;
      } 

      samples_[i]->addLooseCluster(bestLooseID);
      samples_[i]->addTightCluster(bestTightID);
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
