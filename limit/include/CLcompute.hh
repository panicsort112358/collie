#ifndef CLcompute_hh_included
#define CLcompute_hh_included

#include "SigBkgdDist.hh"
#include "CLpoint.hh"
#include <TMath.h>

/** Abstract class for a confidence-level calculator. */
class CLcompute {
public:
  CLcompute()
      : medianExpected_(true),
        noviceFlag_(true) {}
  virtual ~CLcompute() {}

  ///used in the scanning approach
  static const int LEVEL_COMBO;
  ///used in the scanning approach
  static const int LEVEL_NONE;
  ///used in the scanning approach
  static const int LEVEL_COMBO3;
  ///used in the scanning approach
  static const int LEVEL_COMBO2;
  ///used for testing only
  static const int LEVEL_TEST;
  /// Each CL calculation should take O(5 ms)
  static const int LEVEL_VERYVERYFAST;
  /// Each CL calculation should take O(10 ms)
  static const int LEVEL_VERYFAST;
  /// Each CL calculation should take O(25 ms)
  static const int LEVEL_FAST;
  /// Each CL calculation may take O(50 ms)
  static const int LEVEL_STANDARD;
  /// Each CL calculation may take O(100 ms)
  static const int LEVEL_FINE;
  /// Each CL calculation may take O(200 ms)
  static const int LEVEL_VERYFINE;
  /// Each CL calculation may take O(500 ms)
  static const int LEVEL_VERYVERYFINE;
  /// Largest safe override: calculators may internally evaluate twice this many trials.
  static constexpr int MAX_MONTE_CARLO_ITERATIONS = 100000000;

  /** Accessor methods for parent LLR distributions **/
  virtual TH1D* getLLRdist_sb(const char* title, int bins, double min, double max) = 0;
  virtual TH1D* getLLRdist_b(const char* title, int bins, double min, double max) = 0;

  /** \brief Configure the confidence-level calculator.  (implementation dependent) */
  virtual void configure(const char* options) {}
  /** \brief Calculate confidence levels at a default level of effort,
      accuracy, and execution time. */
  inline bool calculateCLs(const SigBkgdDist& sbd, CLpoint& CLs) {
    return calculateCLs(sbd, CLs, LEVEL_STANDARD);
  }
  /** \brief Calculate confidence levels at the specified level of effort,
      accuracy, and execution time. */
  virtual bool calculateCLs(const SigBkgdDist& sbd, CLpoint& CLs, int level) = 0;

  //to use or not use statistical uncertainties in input histograms
  virtual void useHistoStats(bool use) = 0;

  //Flags and warnings for novice users, default = true
  virtual void noviceWarning() {}
  inline bool getNoviceFlag() const {
    return noviceFlag_;
  }
  inline void setNoviceFlag(bool flag) {
    noviceFlag_ = flag;
  }

  //Do we want to calculate the median expected limit (appropriate for gaussian problems)
  //  or do we want the data==bkgd limit (appropriate for poisson problems)
  //  Default = true
  inline bool getMedianExpected() const {
    return medianExpected_;
  }
  inline void setMedianExpected(bool flag) {
    medianExpected_ = flag;
  }

protected:
  static bool parseIterations(const char* options, int& iterations);

  bool medianExpected_;

private:
  bool noviceFlag_;
};

#endif
