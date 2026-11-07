#include <CLcompute.hh>

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <limits>

const int CLcompute::LEVEL_NONE = 200;
const int CLcompute::LEVEL_TEST = 100;
const int CLcompute::LEVEL_COMBO3 = 90;
const int CLcompute::LEVEL_COMBO2 = 80;
const int CLcompute::LEVEL_COMBO = 60;
const int CLcompute::LEVEL_VERYVERYFAST = 40;
const int CLcompute::LEVEL_VERYFAST = 20;
const int CLcompute::LEVEL_FAST = 10;
const int CLcompute::LEVEL_STANDARD = 0;
const int CLcompute::LEVEL_FINE = -10;
const int CLcompute::LEVEL_VERYFINE = -20;
const int CLcompute::LEVEL_VERYVERYFINE = -40;

bool CLcompute::parseIterations(const char* options, int& iterations) {
  if (options == nullptr)
    return false;
  const char* setting = std::strstr(options, "its=");
  if (setting == nullptr)
    return false;

  char* end = nullptr;
  errno = 0;
  const long long requested = std::strtoll(setting + 4, &end, 10);
  if (errno == ERANGE || end == setting + 4 || requested < 0 ||
      requested > MAX_MONTE_CARLO_ITERATIONS) {
    return false;
  }
  iterations = static_cast<int>(requested);
  return true;
}
