#include <string>
#include <vector>
#include <map>
#include <memory>

namespace fer {
namespace zesoi {
namespace bioinfo {

enum DistanceMetricType {
  JACCARD,
  EUCLIDEAN,
  PLAIN,
};

class IDistanceMetric {
 public:
  virtual ~IDistanceMetric() {};
  
  // DEPRECATED
  virtual double getDistance(
      const std::string& seq1, const std::string& seq2) const = 0;

  virtual double getDistance(
      const std::vector<int64_t>& seq1, 
      const std::vector<int64_t>& seq2) const = 0;
};

class JaccardDistance : public IDistanceMetric {
 public: 
  JaccardDistance() = default;  
  ~JaccardDistance() = default;

  double getDistance(
      const std::string& seq1, const std::string& seq2) const override;

  double getDistance(
      const std::vector<int64_t>& seq1, 
      const std::vector<int64_t>& seq2) const override;
};

class EuclideanDistance : public IDistanceMetric {
 public: 
  EuclideanDistance() = default;  
  ~EuclideanDistance() = default;

  double getDistance(
      const std::string& seq1, const std::string& seq2) const override;

  double getDistance(
      const std::vector<int64_t>& seq1, 
      const std::vector<int64_t>& seq2) const override;
};

class PlainDistance : public IDistanceMetric {
 public: 
  PlainDistance() = default;  
  ~PlainDistance() = default;

  double getDistance(
      const std::string& seq1, const std::string& seq2) const override;

  double getDistance(
      const std::vector<int64_t>& seq1, 
      const std::vector<int64_t>& seq2) const override;
};

class DistanceMetricFactory {
 public:
  DistanceMetricFactory();

  const IDistanceMetric* getMetric(DistanceMetricType metricType); 
  
 private:
  IDistanceMetric* jaccardDistance_;
  IDistanceMetric* euclideanDistance_;
  IDistanceMetric* plainDistance_;
};

}
}
}
