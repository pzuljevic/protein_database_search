#ifndef DISTANCE_METRIC_H
#define DISTANCE_METRIC_H

#include <string>
#include <vector>
#include <map>
#include <memory>

namespace fer {
namespace zesoi {
namespace bioinfo {

enum DistanceMetricType {
  JACCARD,
};

class IDistanceMetric {
 public:
  virtual ~IDistanceMetric() {};
  
  virtual double getSimilarity(
      const std::vector<int64_t>& seq1, 
      const std::vector<int64_t>& seq2) const = 0;
};

class JaccardDistance : public IDistanceMetric {
 public: 
  JaccardDistance() = default;  
  ~JaccardDistance() = default;

  double getSimilarity(
      const std::vector<int64_t>& seq1, 
      const std::vector<int64_t>& seq2) const override;
};

class DistanceMetricFactory {
 public:
  DistanceMetricFactory();

  const IDistanceMetric* getMetric(DistanceMetricType metricType); 
  
 private:
  IDistanceMetric* jaccardDistance_;
};

}
}
}
#endif
