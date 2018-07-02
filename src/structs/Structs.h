#ifndef STRUCTS_H
#define STRUCTS_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <mutex>
#include <memory>

namespace fer {
namespace zesoi {
namespace bioinfo {

struct LengthClusterTag {
  int64_t tightID;
  int64_t looseID;
};

class DataSample {
 public:
  DataSample(
    const std::string& header, 
    const std::vector<int64_t>& features,
    int64_t length) 
    : id_((int64_t)this), header_(header), features_(features), length_(length) {} 
 
  int64_t getID() const {
    return id_;
  }

  std::string getHeader() const {
    return header_;
  }

  int64_t getLength() const {
    return length_;
  }
  
  const std::vector<int64_t>& getFeatures() const {
    return features_;
  }
  
  void setLengthClusterTag(int64_t tightID, int64_t looseID) {
    lct_.tightID = tightID;
    lct_.looseID = looseID;
  }
  LengthClusterTag getLengthClusterTag() const {
    return lct_;
  }

  void addLooseCluster(int64_t id) {
    std::lock_guard<std::mutex> guard(lck_); 
    looseClusters_.insert(id);
  }
  void eraseLooseCluster(int64_t id) {
    std::lock_guard<std::mutex> guard(lck_); 
    looseClusters_.erase(id);
  }
  void addTightCluster(int64_t id) {
    std::lock_guard<std::mutex> guard(lck_); 
    tightClusters_.insert(id);
  }

  std::unordered_set<int64_t> getLooseClusters() {
    return looseClusters_;
  } 

  std::unordered_set<int64_t> getTightClusters() {
    return tightClusters_;
  } 
  
  void setOffset(int64_t offset) {
    offset_ = offset;
  }
 
  int64_t getOffset() const {
    return offset_;
  } 

  void print() const {
    std::cout << "ID: " << id_ << " " 
              << "Header: " << header_ << " " 
              << "Length: " << length_ << " " 
              << "LengthClusterTag (tight/loose): " 
              << lct_.tightID << " " << lct_.looseID << " " 
              << "Features size: " << features_.size() 
              << " Features:  "; 
    for (const auto f : features_) {
      std::cout << f << " ";
    }
    std::cout << std::endl;
    std::cout << "Tight clusters: "; 
    for (const auto t: tightClusters_) {
      std::cout << t << " ";
    }
    std::cout << std::endl;
    std::cout << "Loose clusters: "; 
    for (const auto l: looseClusters_) {
      std::cout << l << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }

 private:
  int64_t id_;
  std::string header_;
  std::vector<int64_t> features_; 
  int64_t length_;
  LengthClusterTag lct_;
  std::unordered_set<int64_t> tightClusters_; 
  std::unordered_set<int64_t> looseClusters_; 
  std::mutex lck_; 
  int64_t offset_{0};
};

typedef std::shared_ptr<DataSample> DataSamplePtr;

}
}
}
#endif
