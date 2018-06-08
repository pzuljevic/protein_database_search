#include "DistanceMetric.h"

#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <math.h>
#include <vector>
#include <algorithm>

using namespace fer::zesoi::bioinfo;

// DEPRECATED
double JaccardDistance::getDistance(
    const std::string& seq1, const std::string& seq2) const {
  std::unordered_map<char, int64_t> d1;
  for (const auto c : seq1) {
    d1[c]++;
  }

  std::unordered_map<char, int64_t> d2;
  for (const auto c : seq2) {
    d2[c]++;
  }

  int64_t intersectCnt = 0;
  for (const auto& c : d1) {
    intersectCnt += std::min(c.second, d2[c.first]);
  }
  int64_t unionCnt = seq1.size() + seq2.size() - intersectCnt;
  if (unionCnt == 0) {
    return 0;
  }
  return (double) intersectCnt / unionCnt;
}

double JaccardDistance::getDistance(
    const std::vector<int64_t>& seq1, 
    const std::vector<int64_t>& seq2) const {
  int64_t l1 = seq1.size(), l2 = seq2.size();
  int64_t intersectionCnt = 0; 
  for (size_t i = 0, j = 0; i < seq1.size() && i < seq2.size(); ) {
    if (seq1[i] == seq2[j]) {
      intersectionCnt++;
      ++i; 
      ++j;
    } else if (seq1[i] > seq2[j]) {
       j++; 
    } else {
       i++;
    } 
  }
  int64_t unionCnt = l1 + l2 - intersectionCnt;
  if (unionCnt == 0) {
    return 0;
  }
  return (double) intersectionCnt / unionCnt;
}

double EuclideanDistance::getDistance(
    const std::string& seq1, const std::string& seq2) const {
  throw std::runtime_error("Not implemented");
}

double EuclideanDistance::getDistance(
    const std::vector<int64_t>& seq1, 
    const std::vector<int64_t>& seq2) const {
  double dist = 0;
  for (size_t i = 0; i < seq1.size(); ++i) {
    const auto delta = std::abs(seq1[i] - seq2[i]);
    dist += delta * delta;
  }
  return sqrt(dist);
}

double PlainDistance::getDistance(
    const std::string& seq1, const std::string& seq2) const {
  double dist = (double) (seq1.size() - seq2.size()) / seq2.size();
  if (dist > 0) return dist;
  return -dist;
}

double PlainDistance::getDistance(
    const std::vector<int64_t>& seq1, 
    const std::vector<int64_t>& seq2) const {
  double dist = (double) (seq1[0] - seq2[0]) / seq2[0];
  if (dist > 0) return dist;
  return -dist;
}

DistanceMetricFactory::DistanceMetricFactory() : 
  jaccardDistance_(new JaccardDistance()),
  euclideanDistance_(new EuclideanDistance()),
  plainDistance_(new PlainDistance()) {
}

const IDistanceMetric* 
DistanceMetricFactory::getMetric(DistanceMetricType metricType) {
  switch (metricType) {
    case DistanceMetricType::JACCARD:
      return jaccardDistance_;
    case DistanceMetricType::EUCLIDEAN:
      return euclideanDistance_;
    case DistanceMetricType::PLAIN:
      return plainDistance_;
    default:
      throw std::runtime_error("Unknown distance metric"); 
  }
} 
