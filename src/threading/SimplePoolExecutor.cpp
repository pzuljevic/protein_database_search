#include "SimplePoolExecutor.h"

#include <iostream>

using namespace fer::zesoi::bioinfo;

ThreadPool::ThreadPool(size_t numThreads) : stopped_(false) {
  for (size_t i = 0; i < numThreads ;++i) {
    workers_.emplace_back([this] {
      for(;;) {
        std::function<void()> task;
        {
          std::unique_lock<std::mutex> lock(this->queueMutex_);
          this->condition_.wait(lock, [this] { 
            return this->stopped_ || !this->tasks_.empty(); 
          });
          if (this->stopped_ && this->tasks_.empty()) {
              return; 
          }
          task = std::move(this->tasks_.front());
          this->tasks_.pop();
        }
        task();
      }
    });
  }
}

void ThreadPool::stop() {
  if (stopped_ ) {
   return;
  }
  std::cout << "Stopping thread pool executor..." << std::endl; 
  {
    std::unique_lock<std::mutex> lock(queueMutex_);
    stopped_ = true;
  }
  condition_.notify_all();
  for(std::thread &worker: workers_) {
    worker.join();
  }
}

ThreadPool::~ThreadPool() {
  stop();
}
