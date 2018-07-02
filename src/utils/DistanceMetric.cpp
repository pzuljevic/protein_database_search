#include "DistanceMetric.h"

#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <math.h>
#include <vector>
#include <algorithm>

using namespace fer::zesoi::bioinfo;

static double sum = 0;
static int n = 0;

double JaccardDistance::getSimilarity(
    const std::vector<int64_t>& seq1, 
    const std::vector<int64_t>& seq2) const {
  int64_t l1 = seq1.size(), l2 = seq2.size();
  int64_t intersectionCnt = 0; 
  for (size_t i = 0, j = 0; i < seq1.size() && j < seq2.size(); ) {
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
  sum += (double) intersectionCnt / unionCnt;
  n++;
  if (n == 1000) std::cout << sum / n << std::endl;
  if (n == 10000) std::cout << "D " << sum / n << std::endl;

  return (double) intersectionCnt / unionCnt;
}

DistanceMetricFactory::DistanceMetricFactory() : 
  jaccardDistance_(new JaccardDistance()) {
}

const IDistanceMetric* 
DistanceMetricFactory::getMetric(DistanceMetricType metricType) {
  switch (metricType) {
    case DistanceMetricType::JACCARD:
      return jaccardDistance_;
    default:
      throw std::runtime_error("Unknown distance metric"); 
  }
} 
