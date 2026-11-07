#ifndef COLLIE_PARALLEL_EXECUTOR_HH
#define COLLIE_PARALLEL_EXECUTOR_HH

#include "RunContext.hh"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace collie {

template <typename Function>
void parallelFor(std::size_t count, Function&& function, std::size_t requestedWorkers = 0) {
  if (count == 0)
    return;

  const std::size_t configured =
      requestedWorkers == 0 ? RunContext::global().workers() : requestedWorkers;
  const std::size_t workers =
      std::min(count, std::min(RunContext::maximumWorkers(), std::max<std::size_t>(1, configured)));
  if (workers == 1) {
    for (std::size_t index = 0; index < count; ++index)
      function(index);
    return;
  }

  std::atomic<std::size_t> next{0};
  std::atomic<bool> failed{false};
  std::exception_ptr error;
  std::mutex errorMutex;
  std::vector<std::thread> threads;
  threads.reserve(workers);

  try {
    for (std::size_t worker = 0; worker < workers; ++worker) {
      threads.emplace_back([&]() {
        while (!failed.load(std::memory_order_relaxed)) {
          const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
          if (index >= count)
            break;
          try {
            function(index);
          } catch (...) {
            {
              std::lock_guard<std::mutex> lock(errorMutex);
              if (!error)
                error = std::current_exception();
            }
            failed.store(true, std::memory_order_relaxed);
          }
        }
      });
    }
  } catch (...) {
    failed.store(true, std::memory_order_relaxed);
    for (std::thread& thread : threads)
      thread.join();
    throw;
  }

  for (std::thread& thread : threads)
    thread.join();
  if (error)
    std::rethrow_exception(error);
}

} // namespace collie

#endif
