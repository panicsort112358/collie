#ifndef CrossSectionCalc2D_hh_included
#define CrossSectionCalc2D_hh_included

#include "SigBkgdDist.hh"
#include "CLcompute.hh"

/** \brief General cross-section calculation code

     This code scans the magnitude of the input cross section for the
     signal to determine the value which agrees best with data. These
     results can be converted into a central (most-likely) value and
     an associated width.

    \warning The cross-section calculation process can be quite slow!  Many
    confidence level calculation cycles must be performed.  Expect the process
    to require around 10x as long as a single CLcompute::LEVEL_STANDARD computation.
 */
class CrossSectionCalc2D {
public:
  /// constructor
  CrossSectionCalc2D() {
    p_cl = NULL;
    m_granularity = 0.1;
    m_precision = 1;
    m_runPE = false;
    m_verbose = false;
    m_seed = 0.0;
    m_maximum = 10.0;
    m_sbhypo = false;
    m_bohypo = false;
    m_fitValues = {};
  }

  /// set the CLcompute class to use when seeking for a limit
  inline void setup(CLcompute* computer) {
    p_cl = computer;
  }

  /// set the granularity of the xsec steps
  inline void setGranularity(double gr) {
    m_granularity = gr;
  }
  /// get the granularity of the xsec steps
  inline double getGranularity() const {
    return m_granularity;
  }

  /// set the precision required in the calculation
  inline void setPrecision(int pre) {
    m_precision = pre;
  }
  /// get the precision required in the calculation
  inline int getPrecision() const {
    return m_precision;
  }

  /// set the low threshold for xsec searching
  void setLowerBound(double t) {
    m_seed = t;
  }
  /// get the low threshold for limit searching
  double getLowerBound() const {
    return m_seed;
  }

  /// set the upper threshold for xsec searching
  void setUpperBound(double t) {
    m_maximum = t;
  }
  /// get the upper threshold for limit searching
  double getUpperBound() const {
    return m_maximum;
  }

  /// get the verbose setting
  inline bool getVerbose() const {
    return m_verbose;
  }
  /// set whether the limit-seeking code produces status text as it seeks
  inline void setVerbose(bool verby = true) {
    m_verbose = verby;
  }

  inline void runPseudoExp(bool runPE) {
    m_runPE = runPE;
  }

  inline void setSBHypothesis(bool sb) {
    m_sbhypo = sb;
    m_bohypo = false;
  }
  inline void setBOnlyHypothesis(bool sb) {
    m_bohypo = sb;
    m_sbhypo = false;
  }

  /// perform the cross-section limit seeking
  bool calculate(const SigBkgdDist& dist, CLpoint& CLs, bool doPlots = true);

  // perform a series of fits to pseudo-data for determining significances
  bool testFitPE(const SigBkgdDist& dist, double sigVal1, double sigVal2, int nPE);

  // generate a 2D grid from which to extract correlation info
  TH2D* get2DContour(const SigBkgdDist& dist, double xMin, double xMax, int binsx, double yMin,
                     double yMax, int binsy);

  // generate fitted correlation plots for signal rate vs one other variable
  //  TH1D* fitCorrPlot(const char* fname, histo vars);
  // generate fixed parameter correlation plots for signal rate vs one other variable
  //  TH1D* fixedCorrPlot(const char* fname, histo vars);

  // dump statistics and results
  void print();

  // dump error matrix from the fit
  void printEMAT();

private:
  CLpoint calcValue(SigBkgdDist dist, double sf);

  CLcompute* p_cl;
  double m_granularity;
  double m_seed;
  double m_maximum;
  int m_precision;
  bool m_verbose;
  bool m_runPE;
  bool m_sbhypo;
  bool m_bohypo;

  struct {
    double data;

    double signalI;
    double bkgdI;

    double signalF_sb;
    double bkgdF_sb;
    double chi2_sb;
    int status_sb;

    double sigScale1;
    double sigScaleErrTot1;
    double sigScaleErrStat1;
    double sigScaleErrSyst1;

    double sigScale2;
    double sigScaleErrTot2;
    double sigScaleErrStat2;
    double sigScaleErrSyst2;

    double signalF_b;
    double bkgdF_b;
    double chi2_b;
    int status_b;

    vector<double> params;
    vector<double> sigmas;
    vector<string> names;

    vector<vector<double>> emat;

  } m_fitValues;
};

#endif // CrossSectionCalc2D_hh_included
