#ifndef MULTITHREADED_FILE_READER
#define MULTITHREADED_FILE_READER

#include <string>
#include <vector>
#include <functional>
#include <memory>
#include <future>
#include <mutex>
#include <fstream>
#include "../threading/SimplePoolExecutor.h"
#include "../fasta/FastaStructs.h"
#include "../utils/Utils.h"

namespace fer {
namespace zesoi {
namespace bioinfo {

template <class ParseResult>
class MultithreadedFastaReader {
 typedef std::vector<ParseResult> ParseResults;
 typedef std::future<std::vector<ParseResult>> FutureParseResults;

 public:
  MultithreadedFastaReader (
    const char* inFile, 
    std::shared_ptr<ThreadPool>& threadPool,
    size_t numAllocatedThreads,
    std::function<ParseResult(const FastaHeader& fh)> parseCb) 
    : inFile_(inFile),
      threadPool_(threadPool), 
      numAllocatedThreads_(numAllocatedThreads),
      parseCb_(parseCb) {
  }

  inline std::unique_ptr<std::vector<FutureParseResults>> read() {
    auto results = make_unique<std::vector<FutureParseResults>>();
    uint64_t prevOffsetB = 0;
    uint64_t nextOffsetB = 0;
    uint64_t deltaOffsetB = Utils::getFileSizeB(inFile_) / numAllocatedThreads_; 
    // Schedule file read workers
    for (int i = 0; i < numAllocatedThreads_; ++i) {
      prevOffsetB = nextOffsetB;
      nextOffsetB += deltaOffsetB;
      std::cout << "Enqueuing reading file: worker=" << i + 1 
                << "/" << numAllocatedThreads_ << std::endl 
                << " start_seek_gigabytes=" << Utils::convertBtoGB(prevOffsetB)
                << "\tend_seek_gigabytes=" << Utils::convertBtoGB(nextOffsetB)
                << "\tdelta_seek_gigabytes=" << Utils::convertBtoGB(deltaOffsetB)
                << std::endl
                << " start_seek_bytes=" << prevOffsetB 
                << "\tend_seek_bytes=" << nextOffsetB 
                << "\tdelta_seek_bytes=" << deltaOffsetB
                << std::endl;
      auto resultFuture = threadPool_->enqueue(
        [this, i](uint64_t startSeek, uint64_t endSeek) { 
          auto fp = FastaParser();
          auto result = fp.parseFile(inFile_, startSeek, endSeek, parseCb_);
          {
            std::unique_lock<std::mutex> lock(stderrMutex_);
            std::cout << "Processing data for worker: " << i + 1 
                      << " from: "<< startSeek 
                      << " to: " << endSeek 
                      << " took " << fp.getProcessingTimeSec() << " seconds, " 
                      << " read " << result.size() << " samples" << std::endl;
          }
          return result;
        }, 
        prevOffsetB, 
        nextOffsetB
      );
      results->push_back(std::move(resultFuture));
    }
    return std::move(results);
  }
  
 private:
  const char* inFile_;
  std::shared_ptr<ThreadPool> threadPool_;
  const size_t numAllocatedThreads_;
  std::function<ParseResult(const FastaHeader& fh)> parseCb_; 
  std::mutex stderrMutex_;
}; 

}
}
}

#endif
