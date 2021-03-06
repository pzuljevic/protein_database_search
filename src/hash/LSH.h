#ifndef LSH_H
#define LSH_H

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <future>
#include <memory>
#include <stdlib.h>

namespace fer {
namespace zesoi {
namespace bioinfo {

#define ALPHABET_SIZE 24
#define MAX_KMER_VALUE 4586471424 // 24^7

enum HashType {
  MD5,
  PLAIN,
  PLAIN_MOD,
};

class LSH {
 public:

  static HashType intToHashType(int type);

  static std::function<int64_t(const std::string&, int64_t, int64_t)> getHashFunc(
    HashType hashType = HashType::MD5, 
    int64_t MOD = 0);

  static std::vector<int64_t> getCoeff(int64_t K);

  static std::vector<int64_t> getKMinHash(
    const std::string& seq, 
    int64_t K, 
    const std::vector<int64_t>& a,
    const std::vector<int64_t>& b
  );

  static std::vector<int64_t> getMinHash(
    const std::string& seq, 
    int64_t K, 
    int64_t W,
    HashType hashType = HashType::MD5, 
    int64_t MOD = 0
  );

  static std::vector<int64_t> getMinHashW(
    const std::string& seq, 
    int64_t K, 
    int64_t W,
    HashType hashType = HashType::MD5, 
    int64_t MOD = 0
  );

  static std::vector<int64_t> getMaxHash(
    const std::string& seq, 
    int64_t K, 
    int64_t W,
    HashType hashType = HashType::MD5, 
    int64_t MOD = 0
  );

  static int64_t getKMerValuePlain(
    const std::string& seq, int64_t offset, int64_t K);

  static int64_t getKMerValuePlainMOD(
    const std::string& seq, int64_t offset, int64_t K, int64_t MOD);

  static int64_t getKMerValueMD5(
    const std::string& seq, int64_t offset, int64_t K);

  static int64_t md5ToInt(const char* fingerprint);

 protected:
  LSH() {}
};

}
}
}
#endif
