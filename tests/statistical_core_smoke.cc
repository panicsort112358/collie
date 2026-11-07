#include "BinnedModel.hh"
#include "ParallelExecutor.hh"
#include "RunContext.hh"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
  collie::Sample signal;
  signal.name = "signal";
  signal.signal = true;
  signal.nominal = {2.0, 2.0};
  collie::NuisanceEffect calibration;
  calibration.name = "calibration";
  calibration.downFraction = {-0.1, -0.1};
  calibration.upFraction = {0.2, 0.2};
  signal.effects = {calibration};

  collie::Sample background;
  background.name = "background";
  background.nominal = {5.0, 6.0};

  collie::Channel channel;
  channel.name = "counting";
  channel.observed = {7.0, 8.0};
  channel.samples = {signal, background};

  collie::BinnedModel model;
  model.channels = {channel};
  if (!model.validate().empty() || model.binCount() != 2) {
    std::cerr << "Synthetic binned model failed validation\n";
    return 1;
  }
  if (std::abs(model.poissonDeviance({1.0})) > 1e-12) {
    std::cerr << "Asimov model did not recover its injected strength\n";
    return 2;
  }
  const std::vector<double> shifted =
      model.expected(model.channels.front(), {1.0}, {{"calibration", 1.0}});
  if (shifted.size() != 2 || std::abs(shifted[0] - 7.4) > 1e-12 ||
      std::abs(shifted[1] - 8.4) > 1e-12) {
    std::cerr << "Synthetic nuisance interpolation regressed\n";
    return 6;
  }

  const collie::RunContext first(12345, 2);
  const collie::RunContext second(12345, 8);
  if (first.seedFor("pseudoexperiment", 17) != second.seedFor("pseudoexperiment", 17)) {
    std::cerr << "Derived seed depends on worker count\n";
    return 3;
  }
  const collie::RunContext capped(12345, collie::RunContext::maximumWorkers() + 1);
  if (capped.workers() != collie::RunContext::maximumWorkers()) {
    std::cerr << "Worker cap was not enforced\n";
    return 7;
  }

  std::vector<std::uint64_t> serial(64);
  std::vector<std::uint64_t> parallel(64);
  collie::parallelFor(
      serial.size(),
      [&](std::size_t index) { serial[index] = first.seedFor("pseudoexperiment", index); }, 1);
  collie::parallelFor(
      parallel.size(),
      [&](std::size_t index) { parallel[index] = first.seedFor("pseudoexperiment", index); }, 4);
  if (serial != parallel) {
    std::cerr << "Parallel execution changed deterministic random streams\n";
    return 4;
  }

  bool propagated = false;
  try {
    collie::parallelFor(
        32,
        [](std::size_t index) {
          if (index == 7)
            throw std::runtime_error("expected");
        },
        4);
  } catch (const std::runtime_error&) {
    propagated = true;
  }
  if (!propagated) {
    std::cerr << "Parallel worker exception was not propagated\n";
    return 8;
  }

  return 0;
}
