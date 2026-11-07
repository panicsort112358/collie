#ifndef COLLIE_RANDOM_ENGINE_HH
#define COLLIE_RANDOM_ENGINE_HH

#include "RunContext.hh"

#include <cmath>
#include <cstdint>
#include <random>

namespace collie {

/**
 * Small standard-library random adapters used by the legacy calculators.
 *
 * Explicit seeds are supported for indexed, schedule-independent parallel
 * work. The default constructors remain convenient for serial legacy paths
 * and draw from the process RunContext.
 */
class GaussianRandom {
public:
  explicit GaussianRandom(std::uint64_t seed = RunContext::global().nextSeed("legacy.gaussian"))
      : engine_(seed),
        distribution_(0.0, 1.0) {}

  double fire() {
    return distribution_(engine_);
  }

  void setSeed(std::uint64_t seed) {
    engine_.seed(seed);
    distribution_.reset();
  }

private:
  std::mt19937_64 engine_;
  std::normal_distribution<double> distribution_;
};

class PoissonRandom {
public:
  explicit PoissonRandom(std::uint64_t seed = RunContext::global().nextSeed("legacy.poisson"))
      : engine_(seed) {}

  double fire(double mean) {
    if (!(mean > 0.0) || !std::isfinite(mean))
      return 0.0;
    return static_cast<double>(std::poisson_distribution<std::uint64_t>(mean)(engine_));
  }

  void setSeed(std::uint64_t seed) {
    engine_.seed(seed);
  }

private:
  std::mt19937_64 engine_;
};

} // namespace collie

#endif
