#include "LSH.h"

#include <set>
#include <cstdint>
#include "md5.h"

using namespace fer::zesoi::bioinfo;

HashType LSH::intToHashType(int type) {
  switch (type) {
   case 0: 
     std::cout << "'MD5' is being used for hashing" << std::endl;
     return HashType::MD5;
   case 1:
     std::cout << "'PLAIN' is being used for hashing" << std::endl;
     return HashType::PLAIN;
   case 2:
     std::cout << "'PLAIN_MOD' is being used for hashing" << std::endl;
     return HashType::PLAIN_MOD;
   default:
     std::cout << "'MD5' is being used for hashing" << std::endl;
     return HashType::MD5;
  }
};

std::function<int64_t(const std::string&, int64_t, int64_t)> 
LSH::getHashFunc(
    HashType hashType, 
    int MOD) {
  switch (hashType) {
    case HashType::MD5:
      return [](const std::string& str, int64_t pos, int64_t L) {
        return getKMerValueMD5(str, pos, L);
      };
    case HashType::PLAIN:
      return [](const std::string& str, int64_t pos, int64_t L) {
        return getKMerValuePlain(str, pos, L);
      };
    case HashType::PLAIN_MOD:
      return [MOD](const std::string& str, int64_t pos, int64_t L) {
        return getKMerValuePlainMOD(str, pos, L, MOD);
      };
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
  auto hashEval = getHashFunc(hashType, MOD);
  for (int64_t i = 0; i < (int64_t)seq.size(); i += W) {
    int64_t minIndex = i;
    int64_t minValue = INT64_MAX;
    for (int j = 0; j < W && (i + j) < seq.size(); ++j) {
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

std::vector<int64_t> LSH::getMaxHash(
    const std::string& seq, 
    int K, 
    int W,
    HashType hashType, 
    int MOD) {
  std::set<int64_t> features;
  auto hashEval = getHashFunc(hashType, MOD);
  for (int64_t i = 0; i < (int64_t)seq.size(); i += W) {
    int64_t maxIndex = i;
    int64_t maxValue = INT64_MIN;
    for (int j = 0; j < W && (i + j) < seq.size(); ++j) {
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

int64_t LSH::getKMerValuePlain(
    const std::string& seq, int64_t offset, int64_t K) {
  int64_t val = 0;
  for (int64_t i = offset; i < offset + K && i < seq.size(); ++i) {
    val *= ALPHABET_SIZE;
    val += (int)seq[i] - 'A' + 1; 
  }
  return val;
}

int64_t LSH::getKMerValuePlainMOD(
    const std::string& seq, int64_t offset, int64_t K, int64_t MOD) {
  return getKMerValuePlain(seq, offset, K) % MOD;
}

int64_t LSH::getKMerValueMD5(
    const std::string& seq, int64_t offset, int64_t K) {
  return md5ToInt(md5(seq.substr(offset, K)).c_str());
}

int64_t LSH::md5ToInt(const char* fingerprint) {
  size_t len = strlen(fingerprint);
  size_t offset = len < 16 ? 0 : len-16;
  unsigned long long hash_tail = strtoull(fingerprint + offset, NULL, 16);
  return hash_tail % INT64_MAX;
}
