#ifndef COLLIE_TOOLS_COLLIECLI_HH
#define COLLIE_TOOLS_COLLIECLI_HH

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace collie {
namespace cli {

struct Point {
  int var1 = 0;
  int var2 = 0;
  int var3 = 0;

  bool operator<(const Point& other) const {
    if (var1 != other.var1)
      return var1 < other.var1;
    if (var2 != other.var2)
      return var2 < other.var2;
    return var3 < other.var3;
  }

  bool operator==(const Point& other) const {
    return var1 == other.var1 && var2 == other.var2 && var3 == other.var3;
  }
};

struct PointSelection {
  bool all = true;
  int dimensions = 0;
  Point point;
};

struct InputConfig {
  std::vector<std::string> files;
  PointSelection selection;
};

enum class LimitMethod { Fast, Systematics, Profile, Profile2 };

struct RunConfig {
  std::uint64_t seed = 1;
  std::size_t workers = 0;
  bool seedExplicit = false;
  bool workersExplicit = false;
  bool verbose = false;
  bool force = false;
};

struct InspectConfig {
  std::vector<std::string> files;
  bool details = false;
};

struct ValidateConfig {
  std::vector<std::string> files;
  bool strict = false;
};

struct LimitConfig {
  InputConfig input;
  RunConfig run;
  std::string output;
  LimitMethod method = LimitMethod::Profile2;
  int precision = 1;
  int iterations = 0;
  int sigma = 0;
  double confidenceLevel = 0.95;
  double accuracy = 0.001;
  double searchStart = 1.0;
  double signalExclusion = 0.005;
  bool calculateExpected = true;
  bool calculateObserved = true;
  bool useHistogramStatistics = false;
  bool fitSignal = false;
};

enum class DataHypothesis { Observed, BackgroundOnly, SignalPlusBackground };

struct FitConfig {
  InputConfig input;
  RunConfig run;
  std::string output;
  int precision = 1;
  double lowerBound = 0.0;
  double upperBound = 10.0;
  double granularity = 0.1;
  bool scan = false;
  bool makePlots = false;
  int pseudoExperiments = 0;
  double pseudoSignal = 0.0;
  DataHypothesis hypothesis = DataHypothesis::Observed;
};

struct CombineConfig {
  RunConfig run;
  std::string output;
  std::vector<std::string> inputs;
  bool allowDuplicatePoints = false;
};

int run(int argc, char* argv[]);

} // namespace cli
} // namespace collie

#endif
