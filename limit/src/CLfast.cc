#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "CLfast.hh"
#include "ParallelExecutor.hh"
#include "RunContext.hh"

namespace {

struct TrialResult {
  double background;
  double signalBackground;
};

double llrWeight(double background, double signal) {
  if (background < 1e-6)
    return std::log(1.0 + signal * 1e6);
  return std::log(1.0 + signal / background);
}

double calculateLLR(const std::vector<double>& background, const std::vector<double>& signal,
                    const std::vector<double>& data) {
  double total = 0.0;
  for (std::size_t bin = 0; bin < background.size(); ++bin) {
    if (background[bin] == 0.0 && signal[bin] == 0.0)
      continue;
    total += -signal[bin] + data[bin] * llrWeight(background[bin], signal[bin]);
  }
  return -2.0 * total;
}

double poisson(std::mt19937_64& engine, double mean) {
  if (!(mean > 0.0) || !std::isfinite(mean))
    return 0.0;
  std::poisson_distribution<unsigned long long> distribution(mean);
  return static_cast<double>(distribution(engine));
}

} // namespace

template <class T> inline const T& MAXof(const T& a, const T& b) {
  return (a > b) ? a : b;
}

CLfast::CLfast() {
  iterations_ = 0;
  callIndex_ = 0;
  nmcdone_ = 0;
  useHistoStats_ = false;
  medianExpected_ = true;
  setNoviceFlag(true);
}

bool CLfast::calculateCLs(const SigBkgdDist& sbd, CLpoint& CLs, int effort) {
  int effortIterations = 25000;
  if (effort == LEVEL_TEST)
    effortIterations = 100;
  else if (effort == LEVEL_VERYVERYFAST)
    effortIterations = 5000;
  else if (effort == LEVEL_VERYFAST)
    effortIterations = 15000;
  else if (effort == LEVEL_FAST)
    effortIterations = 25000;
  else if (effort == LEVEL_STANDARD)
    effortIterations = 50000;
  else if (effort == LEVEL_FINE)
    effortIterations = 100000;
  else if (effort == LEVEL_VERYFINE)
    effortIterations = 200000;
  else if (effort == LEVEL_VERYVERYFINE)
    effortIterations = 500000;

  return doCLs(sbd, CLs, iterations_ > 0 ? iterations_ : effortIterations);
}

bool CLfast::doCLs(const SigBkgdDist& sbd, CLpoint& CLs, int its) {
  const int nbins = sbd.nbins();
  if (nbins <= 0) {
    nmcdone_ = 0;
    return false;
  }

  const std::vector<double> background(sbd.bkgd(), sbd.bkgd() + nbins);
  const std::vector<double> signal(sbd.signal(), sbd.signal() + nbins);
  const std::vector<double> observed(sbd.data(), sbd.data() + nbins);
  const auto validMean = [](double value) { return std::isfinite(value) && value >= 0.0; };
  if (!std::all_of(background.begin(), background.end(), validMean) ||
      !std::all_of(signal.begin(), signal.end(), validMean) ||
      !std::all_of(observed.begin(), observed.end(),
                   [](double value) { return std::isfinite(value); })) {
    nmcdone_ = 0;
    return false;
  }

  std::vector<double> backgroundExpected(background);
  std::vector<double> signalBackgroundExpected(background);
  for (int bin = 0; bin < nbins; ++bin)
    signalBackgroundExpected[bin] += signal[bin];

  const double llr_d = calculateLLR(background, signal, observed);
  const double llr_b_expected = calculateLLR(background, signal, backgroundExpected);
  const double llr_sb_expected = calculateLLR(background, signal, signalBackgroundExpected);

  std::vector<double> weights(static_cast<std::size_t>(nbins), 0.0);
  for (int bin = 0; bin < nbins; ++bin) {
    if (background[bin] != 0.0 || signal[bin] != 0.0)
      weights[bin] = llrWeight(background[bin], signal[bin]);
  }

  int numer = 0;
  int denom = 0;
  int denom_bexp = 0;

  CLs.llrobs = llr_d;
  CLs.llrb = llr_b_expected;
  CLs.llrsb = llr_sb_expected;
  const collie::RunContext& runContext = collie::RunContext::global();
  const std::string randomStream = "CLfast.poisson." + std::to_string(callIndex_++);
  const std::size_t maxTrials = static_cast<std::size_t>(its) * 2;
  const std::size_t keepCapacity = std::min(MaxKeep_, maxTrials);
  btrial_.resize(keepCapacity);
  sbtrial_.resize(keepCapacity);
  const std::size_t batchCapacity = std::min<std::size_t>(8192, maxTrials);
  std::vector<TrialResult> trials(batchCapacity);

  int nmc = 0;
  while (static_cast<std::size_t>(nmc) < maxTrials && (denom < its / 2.0 || nmc < its)) {
    const std::size_t firstTrial = static_cast<std::size_t>(nmc);
    const std::size_t batchSize = std::min(batchCapacity, maxTrials - firstTrial);

    collie::parallelFor(batchSize, [&](std::size_t offset) {
      const std::size_t trialIndex = firstTrial + offset;
      std::mt19937_64 engine(runContext.seedFor(randomStream, trialIndex));
      double backgroundTotal = 0.0;
      double signalBackgroundTotal = 0.0;

      for (int bin = 0; bin < nbins; ++bin) {
        if (background[bin] == 0.0 && signal[bin] == 0.0)
          continue;
        const double backgroundDraw = poisson(engine, background[bin]);
        const double signalDraw = poisson(engine, signal[bin]);
        backgroundTotal += -signal[bin] + backgroundDraw * weights[bin];
        signalBackgroundTotal += -signal[bin] + (backgroundDraw + signalDraw) * weights[bin];
      }

      trials[offset].background = -2.0 * backgroundTotal;
      trials[offset].signalBackground = -2.0 * signalBackgroundTotal;
    });

    for (std::size_t offset = 0; offset < batchSize && (denom < its / 2.0 || nmc < its);
         ++offset, ++nmc) {
      const double llr_b = trials[offset].background;
      const double llr_sb = trials[offset].signalBackground;

      if (llr_d <= llr_b) {
        denom++; // background "less than" data
        if (llr_d <= llr_sb)
          numer++; // s+b "less than" data
      }
      if (llr_b_expected <= llr_b)
        denom_bexp++; // background "less than" med expected b
      if (static_cast<std::size_t>(nmc) < MaxKeep_) {
        btrial_[nmc] = llr_b;
        sbtrial_[nmc] = llr_sb;
      }
    }
  }

  nmcdone_ = nmc;
  CLs.clb_obs = (double(denom) / double(nmcdone_));
  CLs.clsb_obs = (double(numer) / double(nmcdone_));
  CLs.clb_med = (double(denom_bexp) / double(nmcdone_));

  double cls_obs = 1.0;
  if (denom > 0)
    cls_obs = 1.0 - double(numer) / double(denom); // fraction of trials where data < s+b
  if (sbd.totSignal() < 1e-5)
    cls_obs = 0.0;
  CLs.cls_obs = cls_obs;

  const int lim = static_cast<int>(std::min<std::size_t>(nmcdone_, MaxKeep_));

  CLs.nTrials_obs = nmcdone_;
  CLs.nTrials_exp = lim;

  std::sort(btrial_.begin(), btrial_.begin() + lim);
  std::sort(sbtrial_.begin(), sbtrial_.begin() + lim);

  CLs.cls_med = -1.0;
  CLs.cls_med_m1s = -1.0;
  CLs.cls_med_m2s = -1.0;
  CLs.cls_med_p1s = -1.0;
  CLs.cls_med_p2s = -1.0;
  const int medianDenominator = std::max(1, lim / 2);

  int m1s_pos = std::min(lim - 1, int(lim * (1.0 - 0.3173 / 2.0) + 0.5) + 1);
  int p1s_pos = std::min(lim - 1, int(lim * (0.3173 / 2.0) + 0.5) + 1);
  int m2s_pos = std::min(lim - 1, int(lim * (1.0 - 0.0455 / 2.0) + 0.5) + 1);
  int p2s_pos = std::min(lim - 1, int(lim * (0.0455 / 2.0) + 0.5) + 1);

  CLs.llrsb_m2s = sbtrial_[m2s_pos];
  CLs.llrsb_m1s = sbtrial_[m1s_pos];
  CLs.llrsb_p1s = sbtrial_[p1s_pos];
  CLs.llrsb_p2s = sbtrial_[p2s_pos];

  CLs.llrb_m2s = btrial_[m2s_pos];
  CLs.llrb_m1s = btrial_[m1s_pos];
  CLs.llrb_p1s = btrial_[p1s_pos];
  CLs.llrb_p2s = btrial_[p2s_pos];

  double val = btrial_[lim / 2];
  for (int i = 0; i < lim; i++) {
    if (sbtrial_[i] >= val && CLs.cls_med < 0) {
      CLs.clsb_med = double(lim - i) / double(lim);
      CLs.cls_med = MAXof(0.0, double(i) / medianDenominator - 1.0);
    }
    if (sbtrial_[i] >= btrial_[m1s_pos] && CLs.cls_med_m1s < 0) {
      CLs.cls_med_m1s = 1.0 - double(lim - i) / double(lim - m1s_pos);
      CLs.clsb_med_m1s = double(lim - i) / double(lim);
    }
    if (sbtrial_[i] >= btrial_[m2s_pos] && CLs.cls_med_m2s < 0) {
      CLs.cls_med_m2s = 1.0 - double(lim - i) / double(lim - m2s_pos);
      CLs.clsb_med_m2s = double(lim - i) / double(lim);
    }
    if (sbtrial_[i] >= btrial_[p1s_pos] && CLs.cls_med_p1s < 0) {
      CLs.cls_med_p1s = 1.0 - double(lim - i) / double(lim - p1s_pos);
      CLs.clsb_med_p1s = double(lim - i) / double(lim);
    }
    if (sbtrial_[i] >= btrial_[p2s_pos] && CLs.cls_med_p2s < 0) {
      CLs.cls_med_p2s = 1.0 - double(lim - i) / double(lim - p2s_pos);
      CLs.clsb_med_p2s = double(lim - i) / double(lim);
    }
  }
  if (CLs.cls_med == -1.0)
    CLs.cls_med = 1.0;
  if (CLs.cls_med_m1s == -1.0)
    CLs.cls_med_m1s = 1.0;
  if (CLs.cls_med_m2s == -1.0)
    CLs.cls_med_m2s = 1.0;
  if (CLs.cls_med_p1s == -1.0)
    CLs.cls_med_p1s = 1.0;
  if (CLs.cls_med_p2s == -1.0)
    CLs.cls_med_p2s = 1.0;

  if (CLs.cls_med > 1.0 && CLs.cls_med < 1.001)
    CLs.cls_med = 1.0;
  if (CLs.cls_med < 0.0)
    CLs.cls_med = 0.0;
  CLs.clsb_med = (1.0 - CLs.cls_med) / 2.0;

  if (!medianExpected_) {
    val = llr_b_expected;
    CLs.cls_med = -1.0;
    for (int i = 0; i < lim; i++) {
      if (sbtrial_[i] >= val && CLs.cls_med < 0) {
        CLs.clsb_med = 1.0 - double(i) / double(lim);
        CLs.cls_med = CLs.clb_med > 0.0 ? MAXof(0.0, 1 - CLs.clsb_med / CLs.clb_med) : 1.0;
      }
    }
  }

  CLs.clb_sb = -1.0;
  CLs.clb_sb_p1s = -1.0;
  CLs.clb_sb_p2s = -1.0;
  CLs.clb_sb_m1s = -1.0;
  CLs.clb_sb_m2s = -1.0;
  for (int i = 0; i < lim; i++) {
    if (btrial_[i] > CLs.llrsb && CLs.clb_sb < 0)
      CLs.clb_sb = 1.0 - double(i) / double(lim);
    if (btrial_[i] > CLs.llrsb_p1s && CLs.clb_sb_p1s < 0)
      CLs.clb_sb_p1s = 1.0 - double(i) / double(lim);
    if (btrial_[i] > CLs.llrsb_p2s && CLs.clb_sb_p2s < 0)
      CLs.clb_sb_p2s = 1.0 - double(i) / double(lim);
    if (btrial_[i] > CLs.llrsb_m1s && CLs.clb_sb_m1s < 0)
      CLs.clb_sb_m1s = 1.0 - double(i) / double(lim);
    if (btrial_[i] > CLs.llrsb_m2s && CLs.clb_sb_m2s < 0)
      CLs.clb_sb_m2s = 1.0 - double(i) / double(lim);
  }
  if (CLs.clb_sb == -1.0)
    CLs.clb_sb = 0;
  if (CLs.clb_sb_m1s == -1.0)
    CLs.clb_sb_m1s = 0;
  if (CLs.clb_sb_m2s == -1.0)
    CLs.clb_sb_m2s = 0;
  if (CLs.clb_sb_p1s == -1.0)
    CLs.clb_sb_p1s = 0;
  if (CLs.clb_sb_p2s == -1.0)
    CLs.clb_sb_p2s = 0;

  return true;
}

void CLfast::configure(const char* options) {
  parseIterations(options, iterations_);
}

TH1D* CLfast::getLLRdist_sb(const char* title, int bins, double min, double max) {
  if (nmcdone_ <= 0) {
    if (max < min) {
      min = 0.0;
      max = 1.0;
    }
    return new TH1D(title, title, bins, min, max);
  }

  if (max < min) {
    min = sbtrial_[0];
    max = sbtrial_[std::min<std::size_t>(nmcdone_, MaxKeep_) - 1];
    std::cout << " SB MIN: " << min << " MAX " << max << std::endl;
  }
  TH1D* out = new TH1D(title, title, bins, min, max);
  const int lim = static_cast<int>(std::min<std::size_t>(nmcdone_, MaxKeep_));

  for (int i = 0; i < lim; i++)
    out->Fill(sbtrial_[i]);

  return out;
}

TH1D* CLfast::getLLRdist_b(const char* title, int bins, double min, double max) {
  if (nmcdone_ <= 0) {
    if (max < min) {
      min = 0.0;
      max = 1.0;
    }
    return new TH1D(title, title, bins, min, max);
  }

  if (max < min) {
    min = btrial_[0];
    max = btrial_[std::min<std::size_t>(nmcdone_, MaxKeep_) - 1];
    std::cout << " B MIN: " << min << " MAX " << max << std::endl;
  }
  TH1D* out = new TH1D(title, title, bins, min, max);
  const int lim = static_cast<int>(std::min<std::size_t>(nmcdone_, MaxKeep_));

  for (int i = 0; i < lim; i++)
    out->Fill(btrial_[i]);

  return out;
}

void CLfast::noviceWarning() {
  if (useHistoStats_) {
    printf("\n*****************************************************************\n");
    printf(
        "CLfast Warning:  You're operating in novice mode and you're using\nstatistical uncertainties from your input histograms.  Are you\nsure you've set the values properly in the histograms?\n");
    printf("*****************************************************************\n\n");
  } else {
    printf("\n*****************************************************************\n");
    printf(
        "CLfast Warning:  The histogram statistical uncertainty flag is turned\noff by default.  This calculation does not include statistical uncertainties.\n");
    printf("*****************************************************************\n\n");
  }
}
