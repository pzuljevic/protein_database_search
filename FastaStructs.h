#ifndef FASTA_STRUCTS_H
#define FASTA_STRUCTS_H

#include <string>
#include <vector>
#include <functional>

namespace fer {
namespace zesoi {
namespace bioinfo {

// >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
class FastaHeader {
 public:
  FastaHeader(
    int64_t id, const std::string& header, const std::string& data);
  ~FastaHeader() = default;

  void toStr() const;
  const std::string getHeader() const;
  const int64_t getID() const;
  const std::string getData() const;

 private:
  int64_t id_;
  std::string header_;
  std::string data_;
};

class FeaturizedSample {
 public:
  FeaturizedSample(
    int64_t id, 
    const std::string& header, 
    const std::vector<int64_t>& features) 
    : id_(id), header_(header), features_(features) {}

 private:
  int64_t id_;
  std::string header_;
  std::vector<int64_t> features_; 
};

class FastaParser {
 public:
  FastaParser() {}

  std::vector<FeaturizedSample> parseFile(
    const char* filePath, 
    uint64_t startSeek, 
    uint64_t endSeek,
    std::function<FeaturizedSample(const FastaHeader& fh)> parseCb
  );

  double getProcessingTimeSec() const; 

 private:
  double processingTimeSec_;
};

}
}
}

#endif
