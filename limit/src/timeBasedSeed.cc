#include "timeBasedSeed.hh"
#include "RunContext.hh"

unsigned int timeBasedSeed() {
  const std::uint64_t seed = collie::RunContext::global().nextSeed("legacy-compatibility");
  const unsigned int narrowed = static_cast<unsigned int>(seed ^ (seed >> 32U));
  return narrowed == 0 ? 1U : narrowed;
}
