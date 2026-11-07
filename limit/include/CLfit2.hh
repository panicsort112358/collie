#ifndef CLfit2_hh_included
#define CLfit2_hh_included

#include "SigBkgdDist.hh"
#include "CLcompute.hh"

/** Confidence level calculation by Monte Carlo trials. */
/**
   This class calculates confidence levels via integration of Poisson
   log-likelihood ratios for background-only (NULL) hypotheses and
   signal plus background (TEST) hypotheses.  The marginalization due
   to uncertainties on nuisance parameters is ameliorated by fitting
   outcomes to the data (or pseudodata).

   In this class, the full range of bins is fit **TWICE**, once each for
   the TEST and NULL hypotheses.  The test statistic is formed from
   the log likelihood ratio of these two minimized likelihoods.
**/
class CLfit2 : public CLcompute {
public:
  /// constructor
  CLfit2();

  /** Defined options include:
      "its=(number of Monte Carlo trials to throw)"
  */
  virtual void configure(const char* options);
  virtual bool calculateCLs(const SigBkgdDist& sbd, CLpoint& CLs, int effort);

  /// Calculate statistical uncertainty using contents of histograms.  Default = false
  //    Only turn this on if you're sure your histograms contain the correct
  //    statistical errors per bin!!
  void useHistoStats(bool use) {
    useHistoStats_ = use;
  }

  /// Use the Asimov / Wilks / Wald approximation to speed things up!
  void useAWWapproximation(bool use) {
    doAWW_ = use;
  }

  // Get the total number of MC iterations performed
  inline int getNMCiter() const {
    return nmcdone_;
  }

  //Obtain histograms of S+B and B-only LLR distributions
  TH1D* getLLRdist_sb(const char* title, int bins, double min, double max);
  TH1D* getLLRdist_b(const char* title, int bins, double min, double max);

  void noviceWarning();

private:
  bool doCLs(const SigBkgdDist& sbd, CLpoint& CLs, int its);

  bool doCLsAWW(const SigBkgdDist& sbd, CLpoint& CLs);

  static const int MaxKeep_ = (4 * 1024 * 1024) / sizeof(double) / 2; // allocate 4 MB to this task

  double btrial_[MaxKeep_];
  double sbtrial_[MaxKeep_];

  collie::PoissonRandom rp_;

  int iterations_;

  bool useHistoStats_;

  bool doAWW_;

  int nmcdone_;

  int lun_;
};

#endif
