#ifndef COLLIE_RUN_CONTEXT_HH
#define COLLIE_RUN_CONTEXT_HH

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <string>

namespace collie {

/**
 * Reproducible execution settings shared by stochastic Collie algorithms.
 *
 * The global context reads COLLIE_RANDOM_SEED and COLLIE_NUM_THREADS once at
 * startup. Applications can override it explicitly before constructing any
 * calculator. Derived seeds depend only on the master seed, a stream name and
 * an index, so results do not depend on thread scheduling.
 */
class RunContext {
public:
  explicit RunContext(std::uint64_t seed = 0, std::size_t workers = 0);

  std::uint64_t seed() const noexcept {
    return seed_;
  }
  std::size_t workers() const noexcept {
    return workers_;
  }
  static constexpr std::size_t maximumWorkers() noexcept {
    return 256;
  }

  void setSeed(std::uint64_t seed) noexcept;
  void setWorkers(std::size_t workers) noexcept;

  std::uint64_t seedFor(const std::string& stream, std::uint64_t index = 0) const noexcept;
  std::uint64_t nextSeed(const std::string& stream) noexcept;

  std::string manifestJson(const std::string& algorithm = "") const;

  static RunContext& global();

private:
  static std::uint64_t makeDefaultSeed();
  static std::size_t defaultWorkers() noexcept;
  static std::size_t normalizedWorkers(std::size_t workers) noexcept;

  std::uint64_t seed_;
  std::size_t workers_;
  std::atomic<std::uint64_t> sequence_;
};

void configureGlobalRun(std::uint64_t seed, std::size_t workers = 0) noexcept;

} // namespace collie

#endif
