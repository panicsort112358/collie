#include "SigBkgdDistAdapter.hh"

#include "CollieDistribution.hh"
#include "SigBkgdDist.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <utility>
#include <vector>

namespace collie {
namespace {

std::size_t distributionBinCount(const CollieDistribution& distribution) {
  return static_cast<std::size_t>(distribution.getNXbins()) *
         static_cast<std::size_t>(std::max(1, distribution.getNYbins()));
}

std::vector<double> binContents(const CollieDistribution& distribution, double scale) {
  std::vector<double> contents;
  contents.reserve(distributionBinCount(distribution));
  const int yBins = std::max(1, distribution.getNYbins());
  for (int y = 0; y < yBins; ++y) {
    for (int x = 0; x < distribution.getNXbins(); ++x) {
      contents.push_back(distribution.getEfficiency(x, y) * scale);
    }
  }
  return contents;
}

std::vector<NuisanceEffect> nuisanceEffects(const CollieDistribution& distribution,
                                            const SigBkgdDist& combined,
                                            const std::vector<double>& nominal, double scale) {
  std::vector<NuisanceEffect> effects;
  std::vector<double> fluctuations(static_cast<std::size_t>(combined.getNsyst()), 0.0);

  for (int systematic = 0; systematic < combined.getNsyst(); ++systematic) {
    const std::string name = combined.getSystName(systematic);
    if (!distribution.hasSystematic(name))
      continue;

    NuisanceEffect effect;
    effect.name = name;
    effect.downFraction.reserve(nominal.size());
    effect.upFraction.reserve(nominal.size());

    for (int direction = -1; direction <= 1; direction += 2) {
      fluctuations[static_cast<std::size_t>(systematic)] = static_cast<double>(direction);
      std::vector<double>& fractions = direction < 0 ? effect.downFraction : effect.upFraction;

      std::size_t bin = 0;
      const int yBins = std::max(1, distribution.getNYbins());
      for (int y = 0; y < yBins; ++y) {
        for (int x = 0; x < distribution.getNXbins(); ++x, ++bin) {
          const double varied = distribution.getEfficiencyVaried(x, y, fluctuations.data()) * scale;
          fractions.push_back(nominal[bin] == 0.0 ? 0.0 : varied / nominal[bin] - 1.0);
        }
      }
      fluctuations[static_cast<std::size_t>(systematic)] = 0.0;
    }

    const int localSystematic = distribution.getSystIndex(name);
    const int yBins = std::max(1, distribution.getNYbins());
    const TH1* positive =
        yBins > 1 ? static_cast<const TH1*>(distribution.getPositiveSyst2D(localSystematic))
                  : static_cast<const TH1*>(distribution.getPositiveSyst(localSystematic));
    const TH1* negative =
        yBins > 1 ? static_cast<const TH1*>(distribution.getNegativeSyst2D(localSystematic))
                  : static_cast<const TH1*>(distribution.getNegativeSyst(localSystematic));
    const bool matchingShape = positive != nullptr && negative != nullptr &&
                               positive->GetNbinsX() == distribution.getNXbins() &&
                               negative->GetNbinsX() == distribution.getNXbins() &&
                               positive->GetNbinsY() == yBins && negative->GetNbinsY() == yBins;
    std::size_t bin = 0;
    for (int y = 0; y < yBins; ++y) {
      for (int x = 0; x < distribution.getNXbins(); ++x, ++bin) {
        const double rawUp = matchingShape ? (yBins > 1 ? positive->GetBinContent(x + 1, y + 1)
                                                        : positive->GetBinContent(x + 1))
                                           : std::numeric_limits<double>::quiet_NaN();
        const double rawDown = matchingShape ? (yBins > 1 ? negative->GetBinContent(x + 1, y + 1)
                                                          : negative->GetBinContent(x + 1))
                                             : std::numeric_limits<double>::quiet_NaN();
        if (!std::isfinite(rawUp))
          effect.upFraction[bin] = rawUp;
        if (!std::isfinite(rawDown))
          effect.downFraction[bin] = rawDown;
      }
    }
    effects.push_back(std::move(effect));
  }
  return effects;
}

Sample makeSample(const CollieDistribution& distribution, const SigBkgdDist& combined,
                  const std::string& name, double scale, bool signal, std::size_t signalIndex) {
  Sample sample;
  sample.name = name;
  sample.signal = signal;
  sample.signalIndex = signalIndex;
  sample.nominal = binContents(distribution, scale);
  sample.effects = nuisanceEffects(distribution, combined, sample.nominal, scale);
  return sample;
}

BinnedModel makeAggregateModel(const SigBkgdDist& distribution, const std::string& channelName) {
  Channel channel;
  channel.name = channelName;
  channel.observed.reserve(distribution.nbins());

  Sample signal;
  signal.name = "signal";
  signal.signal = true;
  signal.nominal.reserve(distribution.nbins());

  Sample background;
  background.name = "background";
  background.nominal.reserve(distribution.nbins());

  for (int bin = 0; bin < distribution.nbins(); ++bin) {
    channel.observed.push_back(distribution.data(bin));
    signal.nominal.push_back(distribution.signal(bin));
    background.nominal.push_back(distribution.bkgd(bin));
  }
  channel.samples.push_back(std::move(signal));
  channel.samples.push_back(std::move(background));

  BinnedModel model;
  model.channels.push_back(std::move(channel));
  return model;
}

} // namespace

BinnedModel makeBinnedModel(const SigBkgdDist& distribution,
                            const std::string& fallbackChannelName) {
  const std::vector<std::string> channelNames = distribution.getChannelNames();
  std::size_t nativeBins = 0;
  for (const std::string& channelName : channelNames) {
    const CollieDistribution* data = distribution.getDataDist(channelName);
    if (data != nullptr)
      nativeBins += distributionBinCount(*data);
  }

  if (channelNames.empty() || nativeBins != static_cast<std::size_t>(distribution.nbins())) {
    return makeAggregateModel(distribution, fallbackChannelName);
  }

  std::map<std::string, std::size_t> signalIndices;
  for (int index = 0; index < distribution.getNsignals(); ++index) {
    signalIndices[distribution.getSigName(index)] = static_cast<std::size_t>(index);
  }

  BinnedModel model;
  std::size_t offset = 0;
  for (const std::string& channelName : channelNames) {
    const CollieDistribution* data = distribution.getDataDist(channelName);
    Channel channel;
    channel.name = channelName;
    const std::size_t bins = distributionBinCount(*data);
    channel.observed.reserve(bins);
    for (std::size_t bin = 0; bin < bins; ++bin) {
      channel.observed.push_back(distribution.data(static_cast<int>(offset + bin)));
    }

    for (unsigned int index = 0; index < distribution.getNsigDist(channelName); ++index) {
      const CollieDistribution* signal = distribution.getSigDist(channelName, index);
      if (signal == nullptr)
        continue;
      const std::string rawName = signal->GetName();
      const std::map<std::string, std::size_t>::const_iterator knownIndex =
          signalIndices.find(rawName);
      const std::size_t signalIndex =
          knownIndex == signalIndices.end() ? index : knownIndex->second;
      channel.samples.push_back(makeSample(*signal, distribution, "signal:" + rawName,
                                           distribution.getSignalScale(channelName, rawName), true,
                                           signalIndex));
    }

    for (unsigned int index = 0; index < distribution.getNbkgdDist(channelName); ++index) {
      const CollieDistribution* background = distribution.getBkgdDist(channelName, index);
      if (background == nullptr)
        continue;
      channel.samples.push_back(
          makeSample(*background, distribution, "background:" + std::string(background->GetName()),
                     distribution.getBackgroundScale(static_cast<int>(index)), false, 0));
    }
    offset += bins;
    model.channels.push_back(std::move(channel));
  }
  return model;
}

} // namespace collie
