#ifndef STRUCTS_H
#define STRUCTS_H

#include <string>
#include <vector>

namespace fer {
namespace zesoi {
namespace bioinfo {

class DataSample {
 public:
  DataSample(
    int64_t id, 
    const std::string& header, 
    const std::vector<int64_t>& features,
    int64_t length) 
    : id_(id), header_(header), features_(features), length_(length) {}
 
  int64_t getLength() const {
    return length_;
  }

 private:
  int64_t id_;
  std::string header_;
  std::vector<int64_t> features_; 
  int64_t length_;
};

}
}
}
#endif
