#include "RunContext.hh"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <thread>

namespace {

std::uint64_t splitmix64(std::uint64_t value) noexcept {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

std::uint64_t hashStream(const std::string& stream) noexcept {
  std::uint64_t hash = 14695981039346656037ULL;
  for (unsigned char value : stream) {
    hash ^= value;
    hash *= 1099511628211ULL;
  }
  return hash;
}

std::uint64_t readSeedFromEnvironment() {
  const char* value = std::getenv("COLLIE_RANDOM_SEED");
  if (value == nullptr || *value == '\0')
    return 0;

  char* end = nullptr;
  const unsigned long long parsed = std::strtoull(value, &end, 10);
  return end != value && *end == '\0' ? parsed : 0;
}

std::size_t readWorkersFromEnvironment() {
  const char* value = std::getenv("COLLIE_NUM_THREADS");
  if (value == nullptr || *value == '\0' || *value == '-')
    return 0;

  char* end = nullptr;
  errno = 0;
  const unsigned long long parsed = std::strtoull(value, &end, 10);
  if (errno == ERANGE || end == value || *end != '\0' ||
      parsed > std::numeric_limits<std::size_t>::max() ||
      parsed > collie::RunContext::maximumWorkers()) {
    return 0;
  }
  return static_cast<std::size_t>(parsed);
}

std::string jsonEscape(const std::string& value) {
  std::ostringstream escaped;
  for (const char c : value) {
    switch (c) {
      case '\\':
        escaped << "\\\\";
        break;
      case '"':
        escaped << "\\\"";
        break;
      case '\n':
        escaped << "\\n";
        break;
      case '\r':
        escaped << "\\r";
        break;
      case '\t':
        escaped << "\\t";
        break;
      default:
        escaped << c;
        break;
    }
  }
  return escaped.str();
}

} // namespace

namespace collie {

RunContext::RunContext(std::uint64_t seed, std::size_t workers)
    : seed_(seed == 0 ? makeDefaultSeed() : seed),
      workers_(normalizedWorkers(workers)),
      sequence_(0) {}

void RunContext::setSeed(std::uint64_t seed) noexcept {
  seed_ = seed == 0 ? 1 : seed;
  sequence_.store(0, std::memory_order_relaxed);
}

void RunContext::setWorkers(std::size_t workers) noexcept {
  workers_ = normalizedWorkers(workers);
}

std::uint64_t RunContext::seedFor(const std::string& stream, std::uint64_t index) const noexcept {
  return splitmix64(seed_ ^ splitmix64(hashStream(stream)) ^ splitmix64(index));
}

std::uint64_t RunContext::nextSeed(const std::string& stream) noexcept {
  return seedFor(stream, sequence_.fetch_add(1, std::memory_order_relaxed));
}

std::string RunContext::manifestJson(const std::string& algorithm) const {
  std::ostringstream output;
  output << "{\"schema\":1,\"master_seed\":" << seed_ << ",\"workers\":" << workers_
         << ",\"algorithm\":\"" << jsonEscape(algorithm) << "\""
         << ",\"seed_derivation\":\"splitmix64-indexed\""
         << ",\"uniform_engine\":\"std::mt19937_64\""
         << ",\"distributions\":[\"std::normal_distribution\","
            "\"std::poisson_distribution\"]}";
  return output.str();
}

RunContext& RunContext::global() {
  static RunContext context(readSeedFromEnvironment(), readWorkersFromEnvironment());
  return context;
}

std::uint64_t RunContext::makeDefaultSeed() {
  std::random_device source;
  const auto clock = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  return splitmix64(static_cast<std::uint64_t>(clock) ^
                    (static_cast<std::uint64_t>(source()) << 32U) ^
                    static_cast<std::uint64_t>(source()));
}

std::size_t RunContext::defaultWorkers() noexcept {
  return std::min(maximumWorkers(), std::max<std::size_t>(1, std::thread::hardware_concurrency()));
}

std::size_t RunContext::normalizedWorkers(std::size_t workers) noexcept {
  return workers == 0 ? defaultWorkers() : std::min(workers, maximumWorkers());
}

void configureGlobalRun(std::uint64_t seed, std::size_t workers) noexcept {
  RunContext::global().setSeed(seed);
  if (workers != 0)
    RunContext::global().setWorkers(workers);
}

} // namespace collie
