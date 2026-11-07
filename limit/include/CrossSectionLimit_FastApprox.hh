#ifndef CrossSectionLimit_FA_hh_included
#define CrossSectionLimit_FA_hh_included

#include "SigBkgdDist.hh"
#include "CLcompute.hh"

/** \brief Fast approximation to CLfit2 cross-section limit seeking code.

    This code seeks the factors by which the signal must be scaled for
    a specified confidence level to be obtained.  These factors can be
    converted directly into cross-section limits.  This calculation
    approximates the results you may obtain using the CLfit2 algorithm,
    but in much less time.  The code calculates all LLR values, as well
    as +/-1,2 sigma limits in addition to exp/obs limits.

    \warning This calculation is a close approximation, but only an
    approximation.  The expected limits should be accurate to within
    ~1-2% on average, while the observed limits should be accurate to
    withing ~4-5% on average.
 */

class CrossSectionLimit_FastApprox {
public:
  /// constructor
  CrossSectionLimit_FastApprox() {
    m_accuracy = 0.001;
    m_precision = 1;
    m_seed = 1.00;
    m_cllevel = 0.95;

    m_verbose = false;

    m_fitLLRobs = 0;
    m_expFact = -1;
    m_obsFact = -1;
    m_dLLRnomExp = -1;
    m_dLLRfitExp = -1;
    m_dLLRfitObs = -1;

    m_fitBOnlyLLR.clear();
    m_fitBOnlyDLLR.clear();
  }

  /// set the accuracy required in the calculation
  inline void setAccuracy(double acc) {
    m_accuracy = acc;
  }
  /// get the accuracy required in the calculation
  inline double getAccuracy() const {
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

  /// get the verbose setting
  inline bool getVerbose() const {
    return m_verbose;
  }
  /// set whether the limit-seeking code produces status text as it seeks
  inline void setVerbose(bool verby = true) {
    m_verbose = verby;
  }

  /// get the confidence level limit to seek for
  inline double getCLlevel() const {
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
  void calculateFast(const SigBkgdDist& dist, CLpoint& CLs);
  void calculateFit(const SigBkgdDist& dist, CLpoint& CLs);

  CLpoint myCLpoint;

  double m_accuracy;
  double m_cllevel;
  double m_seed;

  int m_precision;
  bool m_verbose;

  double m_expFact;
  double m_obsFact;
  double m_fitLLRobs;
  double m_dLLRnomExp;
  double m_dLLRfitExp;
  double m_dLLRfitObs;

  vector<double> m_fitBOnlyLLR;
  vector<double> m_fitBOnlyDLLR;
};

#endif // CrossSectionLimit_FA_hh_included
