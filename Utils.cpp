#include "Utils.h"

#include <sys/stat.h>
#include <iostream>
#include <fstream>

using namespace fer::zesoi::bioinfo;

uint64_t Utils::getFileSizeB(const char* filename) {
  struct stat stat_buf;
  int rc = stat(filename, &stat_buf);
  return rc == 0 ? stat_buf.st_size : 0;
} 

double Utils::getFileSizeGB(const char* filename) {
  return convertBtoGB(getFileSizeB(filename));
}

double Utils::convertBtoGB(uint64_t bytes) {
  return bytes / (BYTES_TO_GIGABYTES); 
}
