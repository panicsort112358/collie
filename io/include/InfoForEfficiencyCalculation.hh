#ifndef InfoForEfficiencyCalculation_hh_included
#define InfoForEfficiencyCalculation_hh_included

#include <TObjString.h>
#include <TNamed.h>
#include <TArrayD.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>

/** \brief constants and cached values per s, bin to calculate effeciency

    This class represents the numbers needed to calculate the partial
    efficiency (used as fval) of one distribution for one bin, given a
    value of one of the s parameters.

    Two of the needed numbers are the assymetric sensitivities sigmaP and
    sigmaN.

    The other needed numbers are remembered from the last base point
    (that is, the last point for which multiple s-parameters had
    changed).  For optimal speed when only one parameter has change,
    we want the sum of the fval values for all the OTHER s-parameters.
    For the ability to speed up the case where two parameters change,
    we also keep the base fval for this s-parameter.
*/

struct InfoForEfficiencyCalculation {
  double sigmaP;
  double sigmaN;
  bool assym;
  double exclusionSum;
  double baseDeltaEfficiency;

  ClassDef(InfoForEfficiencyCalculation, 1)
};

#endif
