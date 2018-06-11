#ifndef FASTA_STRUCTS_H
#define FASTA_STRUCTS_H

#include <string>
#include <vector>
#include <functional>
#include "../structs/Structs.h"

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

class FastaParser {
 public:
  FastaParser() {}

  std::vector<DataSamplePtr> parseFile(
    const char* filePath, 
    int64_t startSeek, 
    int64_t endSeek,
    std::function<DataSamplePtr(const FastaHeader& fh)> parseCb
  );

  double getProcessingTimeSec() const; 

 private:
  double processingTimeSec_;
};

}
}
}

#endif
