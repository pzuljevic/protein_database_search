#ifndef THREAD_POOL_H
#define THREAD_POOL_H

// Reference: https://github.com/progschj/ThreadPool/blob/master/ThreadPool.h
#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

namespace fer {
namespace zesoi {
namespace bioinfo {

class ThreadPool {
 public:
  ThreadPool(size_t numThreads);
  template<class F, class... Args>
  auto enqueue(F&& f, Args&&... args) 
      -> std::future<typename std::result_of<F(Args...)>::type>;
  ~ThreadPool();
  void stop();

 private:
  // need to keep track of threads so we can join them
  std::vector<std::thread> workers_;
  // the task queue
  std::queue<std::function<void()>> tasks_;
  
  // synchronization
  std::mutex queueMutex_;
  std::condition_variable condition_;
  bool stopped_;
};
 
// add new work item to the pool
template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args) 
    -> std::future<typename std::result_of<F(Args...)>::type> {
  using return_type = typename std::result_of<F(Args...)>::type;
  auto task = std::make_shared<std::packaged_task<return_type()>>(
    std::bind(std::forward<F>(f), std::forward<Args>(args)...)
  );
      
  std::future<return_type> res = task->get_future();
  {
    std::unique_lock<std::mutex> lock(queueMutex_);
    // don't allow enqueueing after stopped_ping the pool
    if (stopped_) {
      throw std::runtime_error("enqueue on stopped_ped ThreadPool");
    }
    tasks_.emplace([task](){ (*task)(); });
  }
  condition_.notify_one();
  return res;
}

}
}
}
#endif
