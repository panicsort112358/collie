#include "ProfileLikelihoodEngine.hh"

#include <cmath>
#include <iostream>

int main() {
  collie::ProfileParameter strength;
  strength.name = "mu";
  strength.initial = 0.5;
  strength.step = 0.1;
  strength.lower = 0.0;
  strength.upper = 4.0;
  strength.limited = true;

  collie::ProfileLikelihoodEngine profiler;
  const collie::ProfileFitResult fit = profiler.fit(
      {strength},
      [](const double* parameters) {
        const double displacement = parameters[0] - 2.0;
        return displacement * displacement;
      },
      1000, 0.01);
  if (!fit.converged || fit.values.size() != 1 || std::abs(fit.values.front() - 2.0) > 1e-4) {
    std::cerr << "Reentrant Minuit2 smoke fit did not converge\n";
    return 1;
  }

  return 0;
}
