#include "ProfileLikelihoodEngine.hh"

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace collie {

ProfileLikelihoodEngine::ProfileLikelihoodEngine() = default;
ProfileLikelihoodEngine::~ProfileLikelihoodEngine() = default;

ProfileFitResult ProfileLikelihoodEngine::fit(const std::vector<ProfileParameter>& parameters,
                                              Objective objective, std::size_t maxEvaluations,
                                              double tolerance, int printLevel) {
  if (parameters.empty()) {
    throw std::invalid_argument("profile fit requires at least one parameter");
  }

  // Combined runs Migrad first and falls back to Simplex, matching the
  // resilience of the legacy MINUIT "MINIMIZE" command.
  minimizer_.reset(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined"));
  if (!minimizer_) {
    throw std::runtime_error("ROOT Minuit2 minimizer is unavailable");
  }

  objective_ = std::move(objective);
  function_ = std::make_unique<ROOT::Math::Functor>(
      [this](const double* values) { return objective_(values); },
      static_cast<unsigned int>(parameters.size()));

  minimizer_->SetFunction(*function_);
  minimizer_->SetMaxFunctionCalls(static_cast<unsigned int>(maxEvaluations));
  minimizer_->SetMaxIterations(static_cast<unsigned int>(maxEvaluations));
  minimizer_->SetTolerance(tolerance);
  minimizer_->SetPrintLevel(printLevel);
  minimizer_->SetStrategy(1);
  minimizer_->SetErrorDef(1.0);

  for (std::size_t index = 0; index < parameters.size(); ++index) {
    const ProfileParameter& parameter = parameters[index];
    const bool accepted =
        parameter.limited
            ? minimizer_->SetLimitedVariable(static_cast<unsigned int>(index), parameter.name,
                                             parameter.initial, parameter.step, parameter.lower,
                                             parameter.upper)
            : minimizer_->SetVariable(static_cast<unsigned int>(index), parameter.name,
                                      parameter.initial, parameter.step);
    if (!accepted) {
      throw std::runtime_error("Minuit2 rejected parameter '" + parameter.name + "'");
    }
  }

  const bool minimized = minimizer_->Minimize();
  ProfileFitResult result;
  result.status = minimizer_->Status();
  result.converged = minimized && result.status == 0;
  result.covarianceStatus = minimizer_->CovMatrixStatus();
  result.evaluations = minimizer_->NCalls();
  result.minimum = minimizer_->MinValue();
  result.edm = minimizer_->Edm();

  const double* values = minimizer_->X();
  const double* errors = minimizer_->Errors();
  result.values.assign(values, values + parameters.size());
  result.errors.assign(errors, errors + parameters.size());
  return result;
}

bool ProfileLikelihoodEngine::minosError(std::size_t parameter, double& upper, double& lower,
                                         double& parabolic, double& globalCorrelation) {
  if (!minimizer_ || parameter >= minimizer_->NDim())
    return false;

  double minosLower = 0.0;
  double minosUpper = 0.0;
  const bool ok =
      minimizer_->GetMinosError(static_cast<unsigned int>(parameter), minosLower, minosUpper);
  lower = minosLower;
  upper = minosUpper;
  parabolic = minimizer_->Errors()[parameter];
  globalCorrelation = minimizer_->GlobalCC(static_cast<unsigned int>(parameter));
  return ok;
}

bool ProfileLikelihoodEngine::covariance(double* destination, std::size_t dimensions) const {
  if (!minimizer_ || destination == nullptr || dimensions > minimizer_->NDim()) {
    return false;
  }
  for (std::size_t row = 0; row < dimensions; ++row) {
    for (std::size_t column = 0; column < dimensions; ++column) {
      destination[row * dimensions + column] =
          minimizer_->CovMatrix(static_cast<unsigned int>(row), static_cast<unsigned int>(column));
    }
  }
  return true;
}

} // namespace collie
