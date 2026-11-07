#ifndef PROFILE_LIKELIHOOD_TOY_MODEL_HH
#define PROFILE_LIKELIHOOD_TOY_MODEL_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <TRandom3.h>
#include <TMath.h>

namespace ProfileLikelihoodToy {

constexpr int kNBins = 25;
constexpr double kXMin = 0.0;
constexpr double kXMax = 1.0;
constexpr int kNNuisance = 6;
constexpr double kNominalSignalMu = 1.0;

enum Nuisance {
  kLumi = 0,
  kSelection = 1,
  kBkgd1Rate = 2,
  kBkgd2Rate = 3,
  kSignalShape = 4,
  kBkgd2Shape = 5,
};

using ParameterVector = std::array<double, kNNuisance>;
using NuisanceMask = std::array<bool, kNNuisance>;

struct ComponentTemplates {
  std::vector<double> center;
  std::vector<double> signal;
  std::vector<double> background1;
  std::vector<double> background2;
  std::vector<double> signal_shape;
  std::vector<double> background2_shape;
};

struct FitResult {
  double mu = 0.0;
  ParameterVector theta = {};
  ParameterVector theta_error = {};
  std::vector<double> prediction;
  double objective = 0.0;
};

struct NllrDistribution {
  std::vector<double> background_only;
  std::vector<double> signal_plus_background;
  double observed = 0.0;
};

inline double turn_on(double x, double threshold, double width) {
  return 1.0 / (1.0 + std::exp(-(x - threshold) / width));
}

inline double gaussian(double x, double mean, double sigma) {
  const double z = (x - mean) / sigma;
  return std::exp(-0.5 * z * z);
}

inline double clip_positive(double value) {
  return std::max(1.0e-4, value);
}

inline const char* nuisance_label(int i) {
  static const char* labels[kNNuisance] = {
      "Flux norm.",     "Selection eff.", "Det. bkgd rate",
      "Conv. #nu rate", "Rare shape",     "Conv. shape",
  };
  return labels[i];
}

inline NuisanceMask no_nuisance_mask() {
  return {false, false, false, false, false, false};
}

inline NuisanceMask rate_mask() {
  return {true, true, true, true, false, false};
}

inline NuisanceMask shape_mask() {
  return {false, false, false, false, true, true};
}

inline NuisanceMask full_mask() {
  return {true, true, true, true, true, true};
}

inline bool solve_linear_system(std::vector<std::vector<double>> a, std::vector<double> b,
                                std::vector<double>& x) {
  const int n = static_cast<int>(b.size());
  x.assign(n, 0.0);
  if (n == 0) {
    return true;
  }

  for (int col = 0; col < n; ++col) {
    int pivot = col;
    double best = std::fabs(a[col][col]);
    for (int row = col + 1; row < n; ++row) {
      const double candidate = std::fabs(a[row][col]);
      if (candidate > best) {
        best = candidate;
        pivot = row;
      }
    }
    if (best < 1.0e-12) {
      return false;
    }
    if (pivot != col) {
      std::swap(a[pivot], a[col]);
      std::swap(b[pivot], b[col]);
    }

    const double diag = a[col][col];
    for (int j = col; j < n; ++j) {
      a[col][j] /= diag;
    }
    b[col] /= diag;

    for (int row = 0; row < n; ++row) {
      if (row == col) {
        continue;
      }
      const double factor = a[row][col];
      if (std::fabs(factor) < 1.0e-16) {
        continue;
      }
      for (int j = col; j < n; ++j) {
        a[row][j] -= factor * a[col][j];
      }
      b[row] -= factor * b[col];
    }
  }

  x = b;
  return true;
}

inline std::vector<std::vector<double>> invert_matrix(std::vector<std::vector<double>> a) {
  const int n = static_cast<int>(a.size());
  std::vector<std::vector<double>> inv(n, std::vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i) {
    inv[i][i] = 1.0;
  }

  for (int col = 0; col < n; ++col) {
    int pivot = col;
    double best = std::fabs(a[col][col]);
    for (int row = col + 1; row < n; ++row) {
      const double candidate = std::fabs(a[row][col]);
      if (candidate > best) {
        best = candidate;
        pivot = row;
      }
    }
    if (best < 1.0e-12) {
      return inv;
    }
    if (pivot != col) {
      std::swap(a[pivot], a[col]);
      std::swap(inv[pivot], inv[col]);
    }

    const double diag = a[col][col];
    for (int j = 0; j < n; ++j) {
      a[col][j] /= diag;
      inv[col][j] /= diag;
    }

    for (int row = 0; row < n; ++row) {
      if (row == col) {
        continue;
      }
      const double factor = a[row][col];
      for (int j = 0; j < n; ++j) {
        a[row][j] -= factor * a[col][j];
        inv[row][j] -= factor * inv[col][j];
      }
    }
  }
  return inv;
}

class ToyModel {
public:
  ToyModel() {
    build_templates();
    truth_mu_ = 0.82;
    truth_theta_ = {0.12, -0.08, 0.18, -0.16, 0.24, -0.32};
    data_ = prediction(truth_mu_, truth_theta_);
    data_error_.resize(kNBins, 1.0);
    for (int i = 0; i < kNBins; ++i) {
      data_error_[i] = std::sqrt(std::max(1.0, data_[i]));
    }
    observed_global_fit_ = fit_global(data_, full_mask());
    observed_background_fit_ = fit_fixed_mu(data_, 0.0, full_mask());
    observed_signal_fit_ = fit_fixed_mu(data_, kNominalSignalMu, full_mask());
    fill_postfit_errors();
  }

  const ComponentTemplates& templates() const {
    return templates_;
  }
  const std::vector<double>& data() const {
    return data_;
  }
  const std::vector<double>& data_error() const {
    return data_error_;
  }
  double truth_mu() const {
    return truth_mu_;
  }
  const ParameterVector& truth_theta() const {
    return truth_theta_;
  }
  const FitResult& observed_fit() const {
    return observed_global_fit_;
  }
  const FitResult& observed_background_fit() const {
    return observed_background_fit_;
  }
  const FitResult& observed_signal_fit() const {
    return observed_signal_fit_;
  }

  std::vector<double> prediction(double mu, const ParameterVector& theta) const {
    std::vector<double> y(kNBins, 0.0);
    const auto d = deltas(mu);
    for (int i = 0; i < kNBins; ++i) {
      y[i] = base_prediction(i, mu);
      for (int j = 0; j < kNNuisance; ++j) {
        y[i] += theta[j] * d[j][i];
      }
      y[i] = clip_positive(y[i]);
    }
    return y;
  }

  FitResult fit_fixed_mu(const std::vector<double>& counts, double mu,
                         const NuisanceMask& mask) const {
    FitResult result;
    result.mu = mu;
    const auto d = deltas(mu);

    std::vector<int> active;
    for (int j = 0; j < kNNuisance; ++j) {
      if (mask[j]) {
        active.push_back(j);
      }
    }

    const int n_active = static_cast<int>(active.size());
    std::vector<std::vector<double>> normal(n_active, std::vector<double>(n_active, 0.0));
    std::vector<double> rhs(n_active, 0.0);

    for (int bin = 0; bin < kNBins; ++bin) {
      const double base = base_prediction(bin, mu);
      const double variance = std::max(1.0, counts[bin]);
      const double weight = 1.0 / variance;
      const double residual = counts[bin] - base;
      for (int a = 0; a < n_active; ++a) {
        const int ia = active[a];
        rhs[a] += d[ia][bin] * weight * residual;
        for (int b = 0; b < n_active; ++b) {
          const int ib = active[b];
          normal[a][b] += d[ia][bin] * weight * d[ib][bin];
        }
      }
    }

    for (int a = 0; a < n_active; ++a) {
      normal[a][a] += 1.0;
    }

    std::vector<double> theta_active;
    solve_linear_system(normal, rhs, theta_active);

    for (int a = 0; a < n_active; ++a) {
      result.theta[active[a]] = std::max(-1.75, std::min(1.75, theta_active[a]));
    }
    result.prediction = prediction(mu, result.theta);
    result.objective = objective(counts, result.prediction, result.theta, mask);
    return result;
  }

  FitResult fit_global(const std::vector<double>& counts, const NuisanceMask& mask) const {
    FitResult best;
    best.objective = 1.0e100;
    constexpr int n_scan = 361;
    for (int step = 0; step < n_scan; ++step) {
      const double mu = 2.20 * step / static_cast<double>(n_scan - 1);
      FitResult candidate = fit_fixed_mu(counts, mu, mask);
      if (candidate.objective < best.objective) {
        best = candidate;
      }
    }

    const double left = std::max(0.0, best.mu - 0.04);
    const double right = std::min(2.20, best.mu + 0.04);
    for (int step = 0; step <= 80; ++step) {
      const double mu = left + (right - left) * step / 80.0;
      FitResult candidate = fit_fixed_mu(counts, mu, mask);
      if (candidate.objective < best.objective) {
        best = candidate;
      }
    }
    return best;
  }

  double profiled_objective(const std::vector<double>& counts, double mu,
                            const NuisanceMask& mask) const {
    return fit_fixed_mu(counts, mu, mask).objective;
  }

  double profile_delta(const std::vector<double>& counts, double mu,
                       const NuisanceMask& mask) const {
    const FitResult fixed = fit_fixed_mu(counts, mu, mask);
    const FitResult best = fit_global(counts, mask);
    return std::max(0.0, fixed.objective - best.objective);
  }

  double nllr(const std::vector<double>& counts, const NuisanceMask& mask) const {
    const FitResult sb = fit_fixed_mu(counts, kNominalSignalMu, mask);
    const FitResult b = fit_fixed_mu(counts, 0.0, mask);
    return sb.objective - b.objective;
  }

  NllrDistribution generate_nllr_distribution(const NuisanceMask& mask, int seed,
                                              int n_toys = 9000) const {
    NllrDistribution dist;
    dist.background_only.reserve(n_toys);
    dist.signal_plus_background.reserve(n_toys);
    dist.observed = nllr(data_, mask);

    const auto b_expectation = fit_fixed_mu(data_, 0.0, mask).prediction;
    const auto sb_expectation = fit_fixed_mu(data_, kNominalSignalMu, mask).prediction;
    TRandom3 random(seed);

    std::vector<double> toy(kNBins, 0.0);
    for (int t = 0; t < n_toys; ++t) {
      for (int bin = 0; bin < kNBins; ++bin) {
        toy[bin] = random.PoissonD(b_expectation[bin]);
      }
      dist.background_only.push_back(nllr(toy, mask));

      for (int bin = 0; bin < kNBins; ++bin) {
        toy[bin] = random.PoissonD(sb_expectation[bin]);
      }
      dist.signal_plus_background.push_back(nllr(toy, mask));
    }

    return dist;
  }

  std::vector<double> signal_plus_background_prefit() const {
    ParameterVector nominal = {};
    return prediction(kNominalSignalMu, nominal);
  }

  std::vector<double> background_prefit() const {
    ParameterVector nominal = {};
    return prediction(0.0, nominal);
  }

  std::vector<double> signal_shape_positive() const {
    return templates_.signal_shape;
  }

  std::vector<double> signal_shape_negative() const {
    std::vector<double> out = templates_.signal_shape;
    for (double& v : out) {
      v *= -1.0;
    }
    return out;
  }

  std::vector<double> background2_shape_positive() const {
    return templates_.background2_shape;
  }

  std::vector<double> background2_shape_negative() const {
    std::vector<double> out = templates_.background2_shape;
    for (double& v : out) {
      v *= -1.0;
    }
    return out;
  }

  std::vector<double> bin_pulls() const {
    std::vector<double> pulls(kNBins, 0.0);
    for (int bin = 0; bin < kNBins; ++bin) {
      const double fit = observed_signal_fit_.prediction[bin];
      const double stat = std::sqrt(std::max(1.0, data_[bin]));
      const double syst = 0.055 * std::max(1.0, fit);
      const double denom = std::sqrt(stat * stat + syst * syst);
      pulls[bin] = denom > 0.0 ? (data_[bin] - fit) / denom : 0.0;
    }
    return pulls;
  }

  std::vector<double> sb_ratio() const {
    std::vector<double> out(kNBins, 0.0);
    for (int bin = 0; bin < kNBins; ++bin) {
      const double b = templates_.background1[bin] + templates_.background2[bin];
      out[bin] = b > 0.0 ? templates_.signal[bin] / b : 0.0;
    }
    return out;
  }

  std::vector<double> log_sb_ratio() const {
    std::vector<double> out = sb_ratio();
    for (double& v : out) {
      v = std::log10(1.0 + v);
    }
    return out;
  }

  std::array<std::array<double, kNNuisance + 1>, kNNuisance + 1> correlation_matrix() const {
    constexpr int n = kNNuisance + 1;
    std::vector<std::vector<double>> hessian(n, std::vector<double>(n, 0.0));
    const double mu = observed_global_fit_.mu;
    const auto d = deltas(mu);
    std::vector<std::vector<double>> derivatives(n, std::vector<double>(kNBins, 0.0));

    for (int bin = 0; bin < kNBins; ++bin) {
      derivatives[0][bin] =
          templates_.signal[bin] *
          (1.0 + 0.06 * observed_global_fit_.theta[kLumi] +
           0.08 * observed_global_fit_.theta[kSelection] +
           observed_global_fit_.theta[kSignalShape] * templates_.signal_shape[bin]);
      for (int j = 0; j < kNNuisance; ++j) {
        derivatives[j + 1][bin] = d[j][bin];
      }
    }

    for (int bin = 0; bin < kNBins; ++bin) {
      const double variance = std::max(1.0, observed_global_fit_.prediction[bin]);
      const double weight = 1.0 / variance;
      for (int a = 0; a < n; ++a) {
        for (int b = 0; b < n; ++b) {
          hessian[a][b] += derivatives[a][bin] * weight * derivatives[b][bin];
        }
      }
    }
    for (int j = 1; j < n; ++j) {
      hessian[j][j] += 1.0;
    }

    const auto covariance = invert_matrix(hessian);
    std::array<std::array<double, n>, n> corr = {};
    for (int a = 0; a < n; ++a) {
      for (int b = 0; b < n; ++b) {
        const double denom = std::sqrt(std::fabs(covariance[a][a] * covariance[b][b]));
        corr[a][b] = denom > 0.0 ? covariance[a][b] / denom : (a == b ? 1.0 : 0.0);
        corr[a][b] = std::max(-1.0, std::min(1.0, corr[a][b]));
      }
    }
    return corr;
  }

private:
  ComponentTemplates templates_;
  double truth_mu_ = 0.0;
  ParameterVector truth_theta_ = {};
  std::vector<double> data_;
  std::vector<double> data_error_;
  FitResult observed_global_fit_;
  FitResult observed_background_fit_;
  FitResult observed_signal_fit_;

  void build_templates() {
    templates_.center.resize(kNBins);
    templates_.signal.resize(kNBins);
    templates_.background1.resize(kNBins);
    templates_.background2.resize(kNBins);
    templates_.signal_shape.resize(kNBins);
    templates_.background2_shape.resize(kNBins);

    const double width = (kXMax - kXMin) / static_cast<double>(kNBins);
    for (int i = 0; i < kNBins; ++i) {
      const double x = kXMin + (i + 0.5) * width;
      templates_.center[i] = x;
      templates_.background1[i] = 7.45 + 0.18 * std::sin(7.0 * x);
      templates_.background2[i] =
          39.5 / (1.0 + std::exp((x - 0.405) / 0.112)) + 0.20 * std::exp(-2.0 * x);
      templates_.signal[i] =
          3.20 * turn_on(x, 0.690, 0.060) * (0.62 + 0.64 * x + 0.10 * gaussian(x, 0.82, 0.12));

      templates_.signal_shape[i] = -0.46 * gaussian(x, 0.43, 0.13) +
                                   0.38 * gaussian(x, 0.60, 0.15) +
                                   0.065 * gaussian(x, 0.94, 0.07) - 0.015;
      templates_.background2_shape[i] =
          -0.015 + 0.095 * x + 0.115 * gaussian(x, 0.86, 0.10) - 0.105 * gaussian(x, 0.985, 0.035);
    }
    recenter_shape(templates_.signal, templates_.signal_shape);
    recenter_shape(templates_.background2, templates_.background2_shape);
  }

  static void recenter_shape(const std::vector<double>& nominal, std::vector<double>& shift) {
    double weighted = 0.0;
    double total = 0.0;
    for (std::size_t i = 0; i < nominal.size(); ++i) {
      weighted += nominal[i] * shift[i];
      total += nominal[i];
    }
    const double offset = total > 0.0 ? weighted / total : 0.0;
    for (double& value : shift) {
      value -= offset;
      value = std::max(-0.60, std::min(0.60, value));
    }
  }

  double base_prediction(int bin, double mu) const {
    return templates_.background1[bin] + templates_.background2[bin] + mu * templates_.signal[bin];
  }

  std::array<std::vector<double>, kNNuisance> deltas(double mu) const {
    std::array<std::vector<double>, kNNuisance> d;
    for (auto& v : d) {
      v.assign(kNBins, 0.0);
    }
    for (int i = 0; i < kNBins; ++i) {
      const double s = mu * templates_.signal[i];
      const double b1 = templates_.background1[i];
      const double b2 = templates_.background2[i];
      d[kLumi][i] = 0.060 * (s + b1 + b2);
      d[kSelection][i] = 0.080 * s + 0.045 * (b1 + b2);
      d[kBkgd1Rate][i] = 0.140 * b1;
      d[kBkgd2Rate][i] = 0.125 * b2;
      d[kSignalShape][i] = s * templates_.signal_shape[i];
      d[kBkgd2Shape][i] = b2 * templates_.background2_shape[i];
    }
    return d;
  }

  double objective(const std::vector<double>& counts, const std::vector<double>& pred,
                   const ParameterVector& theta, const NuisanceMask& mask) const {
    double value = 0.0;
    for (int i = 0; i < kNBins; ++i) {
      const double variance = std::max(1.0, counts[i]);
      const double residual = counts[i] - pred[i];
      value += residual * residual / variance;
    }
    for (int j = 0; j < kNNuisance; ++j) {
      if (mask[j]) {
        value += theta[j] * theta[j];
      }
    }
    return value;
  }

  void fill_postfit_errors() {
    const auto corr = correlation_matrix();
    for (int j = 0; j < kNNuisance; ++j) {
      observed_global_fit_.theta_error[j] = 0.55 + 0.25 * (1.0 - std::fabs(corr[0][j + 1]));
    }
  }
};

} // namespace ProfileLikelihoodToy

#endif
