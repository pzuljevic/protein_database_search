#include "FastaStructs.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

using namespace fer::zesoi::bioinfo;

FastaHeader::FastaHeader(
    int64_t id, const std::string& header, const std::string& data)
  : id_(id), header_(header), data_(data) { 
}

void FastaHeader::toStr() const {
  std::cout << "ID: " << id_ << std::endl;
  std::cout << "Header: " << header_ << std::endl;
  std::cout << "Data: " << data_ << std::endl;
}

const std::string FastaHeader::getHeader() const {
  return header_;
}

const std::string FastaHeader::getData() const {
  return data_;
}

const int64_t FastaHeader::getID() const {
  return id_;
}

std::vector<DataSamplePtr> FastaParser::parseFile(
    const char* filePath, 
    int64_t startSeek, 
    int64_t endSeek,
    std::function<DataSamplePtr(const FastaHeader& fh)> parseCb) {
  auto start = std::chrono::high_resolution_clock::now();
  std::vector<DataSamplePtr> data;
  std::ifstream input(filePath);
  if(!input.good()){
    std::cerr << "Error opening the file ..." << std::endl;
    return data;
  }

  std::string line, name, content;
  int64_t id = 0;
  input.seekg(startSeek, input.beg);
  while (std::getline(input, line).good()) {
    if (line.empty() || line[0] == '>') { // Identifier marker
      if (!name.empty() ) { // Print out what we read from the last entry
          data.push_back(parseCb(FastaHeader(id++, name, content)));
          name.clear();
      }
      if (input.tellg() > endSeek) {
        break;
      }
      if (!line.empty()) {
          name = line.substr(1);
      }
      content.clear();
    } else if (!name.empty()) {
      if (line.find(' ') != std::string::npos) { // Invalid sequence--no spaces allowed
        name.clear();
        content.clear();
      } else {
        content += line;
      }
    } 
  }
  if (!name.empty()) {
    data.push_back(parseCb(FastaHeader(id++, name, content)));
  }
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  processingTimeSec_ = elapsed.count(); 
  return data;
}

double FastaParser::getProcessingTimeSec() const {
  return processingTimeSec_;
} 
