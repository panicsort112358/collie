#ifndef CrossSectionCalc_hh_included
#define CrossSectionCalc_hh_included

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
class CrossSectionCalc {
public:
  /// constructor
  CrossSectionCalc() {
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
  bool calculate(const SigBkgdDist& dist, CLpoint& CLs, bool makePlots = true);

  /// just do a fit to the specified hypothesis and record values via ROOT
  bool generateFitValues(const SigBkgdDist& dist, bool doEMAT = false);

  // perform a series of fits to pseudo-data for determining significances
  bool testFitPE(const SigBkgdDist& dist, double sigVal, int nPE);

  // generate fitted correlation plots for signal rate vs one other variable
  //  TH1D* fitCorrPlot(const char* fname, histo vars);
  // generate fixed parameter correlation plots for signal rate vs one other variable
  //  TH1D* fixedCorrPlot(const char* fname, histo vars);

  // dump statistics and results
  void print();

  // dump error matrix from the fit
  void printEMAT();

  /*
  void get2DContour(const SigBkgdDist& dist, double errVal, int idx1, int idx2, int nPts, double* xVals, double* yVals);
  void get2DContour(const SigBkgdDist& dist, double errVal, const char* syst1, const char* syst2, int nPts, double* xVals, double* yVals);
  void get2DContour(const SigBkgdDist& dist, double errVal, const char* syst1, const char* syst2, int nPts, TH2D* hist);

  void get1DContour(const SigBkgdDist& dist, int idx1, int nPts, double xMin, double xMax, double* xVals);
  void get1DContour(const SigBkgdDist& dist, const char* syst1, int nPts, double xMin, double xMax, double* xVals);
  */

  TH1D* get1DContour(const SigBkgdDist& dist, double xMin, double xMax, int bins);
  TH2D* get2DContour(const SigBkgdDist& dist, string syst, double xMin, double xMax, int bins,
                     double yMin, double yMax, int binsy);

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
    double sigScale;
    double sigScaleErrTot;
    double sigScaleErrP;
    double sigScaleErrM;
    double sigScaleErrStatP;
    double sigScaleErrStatM;
    double sigScaleErrStat;
    double sigScaleErrSyst;
    double sigScaleErrSystP;
    double sigScaleErrSystM;

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

#endif // CrossSectionCalc_hh_included
