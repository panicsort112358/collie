#include "BinnedModel.hh"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <sstream>

namespace collie {

std::vector<std::string> BinnedModel::validate() const {
  std::vector<std::string> errors;
  if (channels.empty())
    errors.emplace_back("model contains no channels");
  std::set<std::string> channelNames;

  for (const Channel& channel : channels) {
    if (channel.name.empty())
      errors.emplace_back("channel name is empty");
    else if (!channelNames.insert(channel.name).second)
      errors.emplace_back("channel name '" + channel.name + "' is duplicated");
    if (channel.observed.empty()) {
      errors.emplace_back("channel '" + channel.name + "' contains no bins");
      continue;
    }
    if (channel.samples.empty()) {
      errors.emplace_back("channel '" + channel.name + "' contains no samples");
    }
    for (std::size_t bin = 0; bin < channel.observed.size(); ++bin) {
      if (!std::isfinite(channel.observed[bin]) || channel.observed[bin] < 0.0) {
        errors.emplace_back("channel '" + channel.name + "' has invalid observed content");
        break;
      }
    }
    std::set<std::string> sampleNames;
    for (const Sample& sample : channel.samples) {
      if (sample.name.empty()) {
        errors.emplace_back("channel '" + channel.name + "' contains a sample with an empty name");
      } else if (!sampleNames.insert(sample.name).second) {
        errors.emplace_back("sample name '" + sample.name + "' is duplicated in channel '" +
                            channel.name + "'");
      }
      if (sample.nominal.size() != channel.observed.size()) {
        errors.emplace_back("sample '" + sample.name +
                            "' has a different bin count from its channel");
      }
      for (const double value : sample.nominal) {
        if (!std::isfinite(value) || value < 0.0) {
          errors.emplace_back("sample '" + sample.name + "' has invalid nominal content");
          break;
        }
      }
      std::set<std::string> effectNames;
      for (const NuisanceEffect& effect : sample.effects) {
        if (effect.name.empty()) {
          errors.emplace_back("sample '" + sample.name +
                              "' contains a nuisance with an empty name");
        } else if (!effectNames.insert(effect.name).second) {
          errors.emplace_back("nuisance '" + effect.name + "' is duplicated in sample '" +
                              sample.name + "'");
        }
        if (effect.downFraction.size() != sample.nominal.size() ||
            effect.upFraction.size() != sample.nominal.size()) {
          errors.emplace_back("nuisance '" + effect.name +
                              "' has a different bin count from sample '" + sample.name + "'");
        }
        const auto invalidFraction = [](double value) { return !std::isfinite(value); };
        if (std::any_of(effect.downFraction.begin(), effect.downFraction.end(), invalidFraction) ||
            std::any_of(effect.upFraction.begin(), effect.upFraction.end(), invalidFraction)) {
          errors.emplace_back("nuisance '" + effect.name + "' has a non-finite fractional effect");
        }
      }
    }
  }
  return errors;
}

std::size_t BinnedModel::binCount() const noexcept {
  return std::accumulate(
      channels.begin(), channels.end(), std::size_t{0},
      [](std::size_t total, const Channel& channel) { return total + channel.observed.size(); });
}

std::vector<double> BinnedModel::expected(const Channel& channel,
                                          const std::vector<double>& signalStrengths,
                                          const std::map<std::string, double>& nuisances) const {
  std::vector<double> values(channel.observed.size(), 0.0);
  for (const Sample& sample : channel.samples) {
    const double strength = sample.signal && sample.signalIndex < signalStrengths.size()
                                ? signalStrengths[sample.signalIndex]
                                : 1.0;

    for (std::size_t bin = 0; bin < channel.observed.size() && bin < sample.nominal.size(); ++bin) {
      double fractionalShift = 0.0;
      for (const NuisanceEffect& effect : sample.effects) {
        const auto nuisance = nuisances.find(effect.name);
        if (nuisance == nuisances.end())
          continue;
        if (bin >= effect.downFraction.size() || bin >= effect.upFraction.size()) {
          continue;
        }
        const double theta = nuisance->second;
        const double fraction = theta >= 0.0 ? effect.upFraction[bin] : effect.downFraction[bin];
        fractionalShift += std::abs(theta) * fraction;
      }
      values[bin] += sample.nominal[bin] * strength * std::max(0.0, 1.0 + fractionalShift);
    }
  }
  return values;
}

double BinnedModel::poissonDeviance(const std::vector<double>& signalStrengths,
                                    const std::map<std::string, double>& nuisances) const {
  double result = 0.0;
  for (const Channel& channel : channels) {
    const std::vector<double> prediction = expected(channel, signalStrengths, nuisances);
    for (std::size_t bin = 0; bin < channel.observed.size(); ++bin) {
      const double observed = channel.observed[bin];
      const double mean = std::max(prediction[bin], 1e-12);
      result += mean - observed;
      if (observed > 0.0)
        result += observed * std::log(observed / mean);
    }
  }
  for (const auto& nuisance : nuisances) {
    result += 0.5 * nuisance.second * nuisance.second;
  }
  return 2.0 * result;
}

} // namespace collie
