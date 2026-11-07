#ifndef CLfit_hh_included
#define CLfit_hh_included

#include "SigBkgdDist.hh"
#include "CLcompute.hh"

/** Confidence level calculation by Monte Carlo trials. */
/**
   This class calculates confidence levels via integration of Poisson
   log-likelihood ratios for background-only (NULL) hypotheses and
   signal plus background (TEST) hypotheses.  The marginalization due
   to uncertainties on nuisance parameters is ameliorated by fitting
   outcomes to the data (or pseudodata).

   In this class, the full range of bins is fit **ONCE** for the
   background-only hypothesis.  The fit is restricted to sideband
   regions by ignoring bins with log(1.0+s/b)>llrCut, where llrCut is
   a parameter defined by the user (default = 0.005).  This class runs
   much faster than CLfit2 but can be biased by the choice of llrCut.
   Final variable distributions with no clear sidebands will benefit
   little from this formulation and will see more improvement with CLfit2.

   Users can choose to fit based on the TEST or NULL hypotheses, while
   the default is NULL.  When choosing the TEST hypothesis, the llrCut
   chosen will be ignored.
**/

class CLfit : public CLcompute {
public:
  /// constructor
  CLfit();

  /** Defined options include:
      "its=(number of Monte Carlo trials to throw)"
  */
  virtual void configure(const char* options);
  virtual bool calculateCLs(const SigBkgdDist& sbd, CLpoint& CLs, int effort);

  //Include signal in fit?  Default is false.
  inline void fitSignal(bool fitsig) {
    fitSignal_ = fitsig;
  }

  //ignore all bins with log(1.0 + s/b)>llrCut.  Default = 0.005
  inline void logSigExclusion(double llrCut) {
    sigLLR_ = llrCut;
  }

  // Calculate statistical uncertainty using contents of histograms.  Default = false
  //    Only turn this on if you're sure your histograms contain the correct
  //    statistical errors per bin!!
  void useHistoStats(bool use) {
    useHistoStats_ = use;
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

  static const int MaxKeep_ = (4 * 1024 * 1024) / sizeof(double) / 2; // allocate 5 MB to this task

  double btrial_[MaxKeep_];
  double sbtrial_[MaxKeep_];

  collie::PoissonRandom rp_;

  int iterations_;

  int nmcdone_;

  bool useHistoStats_;

  int lun_;

  bool fitSignal_;

  double sigLLR_;
};

#endif
