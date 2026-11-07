#ifndef ThreeSigmaEvidence_hh_included
#define ThreeSigmaEvidence_hh_included

#include "SigBkgdDist.hh"
#include "CLcompute.hh"

/** \brief General cross-section limit seeking code.

    This code seeks the factors by which the signal must be scaled for
    a specified confidence level to be obtained.  These factors can be
    converted directly into cross-section limits.

    \warning The cross-section limit seeking process can be quite slow!  Several
    confidence level calculation cycles must be performed.  Expect the process
    to require around 10x as long as a single CLcompute::LEVEL_STANDARD computation.
 */
class ThreeSigmaEvidence {
public:
  /// constructor
  ThreeSigmaEvidence() {
    p_cl = NULL;
    m_accuracy = 0.00008;
    m_precision = 1;
    m_sigma = 0;
    m_verbose = false;
    m_cllevel = 0.99865; // =1 - (1-0.9973)/2
    m_seed = 0.5;
    m_calcExp = true;
    m_calcObs = true;
    m_expFact = -1;
    m_obsFact = -1;
  }

  /// set the CLcompute class to use when seeking for a limit
  void setup(CLcompute* computer) {
    p_cl = computer;
  }

  /// set the accuracy required in the calculation
  void setAccuracy(double acc) {
    m_accuracy = acc;
  }
  /// get the accuracy required in the calculation
  double getAccuracy() const {
    return m_accuracy;
  }

  /// set the precision required in the calculation
  inline void setPrecision(int pre) {
    m_precision = pre;
  }
  /// get the precision required in the calculation
  inline int getPrecision() const {
    return m_precision;
  }

  /// set the N-sigma required in the calculation
  inline void setNSigma(int sig) {
    m_sigma = sig;
  }
  /// get the sigma required in the calculation
  inline int getNSigma() const {
    return m_sigma;
  }

  /// turn on/off expected or observed limit calculations (default=true for both)
  inline void calculateExpected(bool calc) {
    m_calcExp = calc;
  }
  inline void calculateObserved(bool calc) {
    m_calcObs = calc;
  }

  /// set the low threshold for limit searching
  void setSearchSeed(double t) {
    m_seed = t;
  }
  /// get the low threshold for limit searching
  double getSearchSeed() const {
    return m_seed;
  }

  /// get the verbose setting
  bool getVerbose() const {
    return m_verbose;
  }
  /// set whether the limit-seeking code produces status text as it seeks
  void setVerbose(bool verby = true) {
    m_verbose = verby;
  }

  /// get the confidence level limit to seek for
  double getCLlevel() const {
    return m_cllevel;
  }
  /** \brief set the confidence level limit to seek for
      \warning Limits can be difficult to obtain above 99% CL
  */
  void setCLlevel(double cl) {
    m_cllevel = cl;
    double acc = (1.0f - m_cllevel) / 50.0f;
    if (m_accuracy > acc)
      m_accuracy = acc;
  }
  /// perform the cross-section limit seeking
  bool calculate(const SigBkgdDist& dist, CLpoint& CLs);
  // dump statistics and results
  void print();

private:
  double m_accuracy;
  double m_cllevel;
  double m_seed;
  double m_expFact;
  double m_obsFact;

  int m_precision;
  int m_sigma;
  CLcompute* p_cl;
  bool m_verbose;
  bool m_calcExp;
  bool m_calcObs;
};

#endif // ThreeSigmaEvidence_hh_included
