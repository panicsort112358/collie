#include "CLfast.hh"
#include "CLpoint.hh"
#include "RunContext.hh"
#include "SigBkgdDist.hh"

#include <array>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

namespace {

std::vector<double> fingerprint(const CLpoint& point) {
  return {point.cls_obs,
          point.clb_obs,
          point.clsb_obs,
          point.clb_med,
          point.clsb_med,
          point.cls_med,
          point.cls_med_p1s,
          point.cls_med_p2s,
          point.cls_med_m1s,
          point.cls_med_m2s,
          point.llrobs,
          point.llrb,
          point.llrb_p1s,
          point.llrb_p2s,
          point.llrb_m1s,
          point.llrb_m2s,
          point.llrsb,
          point.llrsb_p1s,
          point.llrsb_p2s,
          point.llrsb_m1s,
          point.llrsb_m2s,
          point.clb_sb,
          point.clb_sb_p1s,
          point.clb_sb_p2s,
          point.clb_sb_m1s,
          point.clb_sb_m2s,
          static_cast<double>(point.nTrials_obs),
          static_cast<double>(point.nTrials_exp)};
}

CLpoint calculate(SigBkgdDist& model, std::size_t workers) {
  collie::configureGlobalRun(0x5eed1234ULL, workers);
  std::unique_ptr<CLfast> calculator(new CLfast);
  calculator->configure("its=2000");
  CLpoint result(125);
  if (!calculator->calculateCLs(model, result, CLcompute::LEVEL_TEST)) {
    std::cerr << "CLfast calculation failed\n";
    std::exit(2);
  }
  if (calculator->getNMCiter() < 2000 || calculator->getNMCiter() > 4000) {
    std::cerr << "Configured adaptive trial count was not honored\n";
    std::exit(3);
  }
  return result;
}

} // namespace

int main() {
  std::array<double, 5> signal{};
  std::array<double, 5> background{};
  std::array<double, 5> observed{};
  SigBkgdDist model("synthetic", 5, 0.0, 5.0, signal.data(), background.data(), observed.data(),
                    125);

  signal = {0.3, 1.4, 2.2, 0.8, 3.1};
  background = {8.0, 4.5, 12.0, 1.7, 6.2};
  observed = {8.0, 7.0, 15.0, 1.0, 10.0};

  const CLpoint serial = calculate(model, 1);
  const CLpoint parallel = calculate(model, 4);
  if (fingerprint(serial) != fingerprint(parallel)) {
    std::cerr << "CLfast result changed with worker count\n";
    return 1;
  }
  return 0;
}
