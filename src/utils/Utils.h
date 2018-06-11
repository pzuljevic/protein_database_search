#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <functional>
#include <fstream>
#include <memory>

#define BYTES_TO_GIGABYTES 1024.0 * 1024.0 * 1024.0

namespace fer {
namespace zesoi {
namespace bioinfo {

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
/*
template<typename T, typename... Args>
std::shared_ptr<T> make_shared(Args&&... args) {
  return std::shared_ptr<T>(new T(std::forward<Args>(args)...));
}*/

class Utils {
 public:
  static uint64_t getFileSizeB(const char* filename);
  static double getFileSizeGB(const char* filename);
  static double convertBtoGB(uint64_t bytes);

 protected:
  Utils() {}
};

}
}
}

#endif
