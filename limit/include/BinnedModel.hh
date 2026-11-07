#ifndef COLLIE_BINNED_MODEL_HH
#define COLLIE_BINNED_MODEL_HH

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace collie {

struct NuisanceEffect {
  std::string name;
  std::vector<double> downFraction;
  std::vector<double> upFraction;
};

struct Sample {
  std::string name;
  bool signal = false;
  std::size_t signalIndex = 0;
  std::vector<double> nominal;
  std::vector<NuisanceEffect> effects;
};

struct Channel {
  std::string name;
  std::vector<double> observed;
  std::vector<Sample> samples;
};

/**
 * Plain C++ representation of a binned statistical model.
 *
 * This type deliberately contains no ROOT objects. ROOT files are adapters at
 * the repository boundary; likelihood evaluation and validation can therefore
 * be tested without dictionaries or file ownership.
 */
class BinnedModel {
public:
  std::vector<Channel> channels;

  std::vector<std::string> validate() const;
  std::size_t binCount() const noexcept;

  std::vector<double> expected(const Channel& channel,
                               const std::vector<double>& signalStrengths = {},
                               const std::map<std::string, double>& nuisances = {}) const;

  double poissonDeviance(const std::vector<double>& signalStrengths = {},
                         const std::map<std::string, double>& nuisances = {}) const;
};

} // namespace collie

#endif
