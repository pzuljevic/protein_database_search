#include "LSH.h"

#include <set>
#include <cstdint>
#include "md5.h"

using namespace fer::zesoi::bioinfo;

std::function<int64_t(const std::string&, int64_t, int64_t)> 
LSH::getHashFunc(
    HashType hashType, 
    int MOD) {
  switch (hashType) {
    case HashType::MD5:
    default:
      return [](const std::string& str, int64_t pos, int64_t L) {
        return getKMerValueMD5(str, pos, L);
      };
  } 
}

std::vector<int64_t> LSH::getMinHash(
    const std::string& seq, 
    int K, 
    int W,
    HashType hashType, 
    int MOD) {
  std::set<int64_t> features;
  auto hashEval = getHashFunc();
  for (int64_t i = 0; i < (int64_t)(seq.size() - W); i += W) {
    int64_t maxIndex = i;
    int64_t maxValue = INT64_MIN;
    for (int j = 0; j < W; ++j) {
      const auto val = hashEval(seq, i+j, K);
      if (val > maxValue) {
        maxIndex = i + j;
        maxValue = val;
      }
    }
    features.insert(maxValue);
  } 
  return std::vector<int64_t>(features.begin(), features.begin());
}

std::vector<int64_t> LSH::getMaxHash(
    const std::string& seq, 
    int K, 
    int W,
    HashType hashType, 
    int MOD) {
  std::set<int64_t> features;
  auto hashEval = getHashFunc();
  for (int64_t i = 0; i < (int64_t)(seq.size() - W); i += W) {
    int64_t minIndex = i;
    int64_t minValue = INT64_MAX;
    for (int j = 0; j < W; ++j) {
      const auto val = hashEval(seq, i+j, K);
      if (val < minValue) {
        minIndex = i + j;
        minValue = val;
      }
    }
    features.insert(minValue);
  } 
  return std::vector<int64_t>(features.begin(), features.begin());
}

int64_t LSH::getKMerValuePlain(
    const std::string& seq, int64_t offset, int64_t K) {
  int64_t val = 0;
  for (int64_t i = offset; i < offset + K && i < seq.size(); ++i) {
    val *= ALPHABET_SIZE;
    val += (int)seq[i] - 'A' + 1; 
  }
  return val;
}

int64_t LSH::getKMerValueMOD(
    const std::string& seq, int64_t offset, int64_t K, int64_t MOD) {
  return getKMerValuePlain(seq, offset, K) % MOD;
}

int64_t LSH::getKMerValueMD5(
    const std::string& seq, int64_t offset, int64_t K) {
  return md5ToInt(md5(seq.substr(offset, K)).c_str());
}

int64_t LSH::md5ToInt(const char* fingerprint) {
  int v1, v2, v3, v4;
  sscanf(&fingerprint[0], "%x", &v1);
  sscanf(&fingerprint[8], "%x", &v2);
  sscanf(&fingerprint[16], "%x", &v3);
  sscanf(&fingerprint[24], "%x", &v4);
  return v1 ^ v2 ^ v3 ^v4;
}
