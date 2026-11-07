#ifndef COLLIE_PROFILE_LIKELIHOOD_ENGINE_HH
#define COLLIE_PROFILE_LIKELIHOOD_ENGINE_HH

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace ROOT {
namespace Math {
class Functor;
class Minimizer;
} // namespace Math
} // namespace ROOT

namespace collie {

struct ProfileParameter {
  std::string name;
  double initial = 0.0;
  double step = 0.1;
  double lower = 0.0;
  double upper = 0.0;
  bool limited = false;
};

struct ProfileFitResult {
  bool converged = false;
  int status = -1;
  int covarianceStatus = -1;
  std::size_t evaluations = 0;
  double minimum = -1.0;
  double edm = -1.0;
  std::vector<double> values;
  std::vector<double> errors;
};

/**
 * Reentrant Minuit2 adapter shared by one- and multi-signal profilers.
 *
 * The objective is owned by the engine and may capture model state. No
 * process-global callback state is required, so independent engines can fit
 * safely on separate threads.
 */
class ProfileLikelihoodEngine {
public:
  using Objective = std::function<double(const double*)>;

  ProfileLikelihoodEngine();
  ~ProfileLikelihoodEngine();

  ProfileLikelihoodEngine(const ProfileLikelihoodEngine&) = delete;
  ProfileLikelihoodEngine& operator=(const ProfileLikelihoodEngine&) = delete;

  ProfileFitResult fit(const std::vector<ProfileParameter>& parameters, Objective objective,
                       std::size_t maxEvaluations, double tolerance = 1.0, int printLevel = 0);

  bool minosError(std::size_t parameter, double& upper, double& lower, double& parabolic,
                  double& globalCorrelation);
  bool covariance(double* destination, std::size_t dimensions) const;

private:
  Objective objective_;
  std::unique_ptr<ROOT::Math::Functor> function_;
  std::unique_ptr<ROOT::Math::Minimizer> minimizer_;
};

} // namespace collie

#endif
