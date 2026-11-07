#ifndef SigBkgdDist_HH_INCLUDED
#define SigBkgdDist_HH_INCLUDED

#include "CollieDistribution.hh"
#include "RandomEngine.hh"
#include "TrialPoint.hh"
#include "TH1D.h"
#include "TObjString.h"
#include "TObjArray.h"
#include <vector>
#include <map>

static const double ln10 = 2.302585092994046;

/** \brief Signal, background, and data together as distributions. */
class SigBkgdDist {
public:
  /// create an empty SigBkgdDist
  SigBkgdDist(int var1, int var2 = 0, int var3 = 0);
  /// construct a SigBkgdDist
  SigBkgdDist(string chan, int bins, double min, double max, int var1, int var2 = 0, int var3 = 0);
  /// assemble a SigBkgdDist from arrays.  These arrays will not be deleted when the SigBkgdDist is destroyed
  SigBkgdDist(string chan, int bins, double min, double max, double* sigarray, double* bkgdarray,
              double* dataarray, int var1, int var2 = 0, int var3 = 0);
  /// copy constructor
  SigBkgdDist(const SigBkgdDist& sbd); // Warning - at present, the copy ctor
                                       // does not have the usual C++ semantics
  /// destructor
  ~SigBkgdDist();

  ///  initiate basic setup of distributions
  void setup(string chan, int bins, double min, double max, int var1, int var2, int var3);

  /// set absolute signal scale, non-multiplicative
  void setSignalScale(double f);
  /// set absolute signal scale, non-multiplicative
  void setSignalScale(string sig, double f);
  void setSignalScale(int sigIndex, double f);
  /// set absolute signal scale, non-multiplicative
  void setSignalScaleChan(string chan, double f);
  /// set absolute signal scale, non-multiplicative
  void setSignalScale(string chan, string sig, double f);
  /// set absolute bkgd scale, non-multiplicative
  void setBackgroundScale(double f);
  /// set absolute data scale, non-multiplicative
  void setDataScale(double f);

  /// scale all signals by a factor common to all channels
  void scaleSignal(double f);
  /// scale one signal by a factor, common to all channels
  void scaleSignal(string sig, double f,
                   std::vector<double> const& fluctMap = std::vector<double>());

  /// scale all signals in one channel by a factor
  void scaleSignalChan(string sig, double f,
                       std::vector<double> const& fluctMap = std::vector<double>());
  /// scale one signal by a factor for one channel
  void scaleSignal(string chan, string sig, double f,
                   std::vector<double> const& fluctMap = std::vector<double>());

  /// scale the background by a factor
  void scaleBackground(double f);
  /// scale one background by a factor
  void scaleBackground(uint idx, double f,
                       std::vector<double> const& fluctMap = std::vector<double>());
  /// scale the data by a factor
  void scaleData(double f);

  /// get the signal scale factor
  inline double getSignalScale() const {
    return gl_cSigScale;
  }
  /// get the signal scale factor for one signal in one channel
  double getSignalScale(string chan, string sig) const;

  /// get the background scale factor
  inline double getBackgroundScale() const {
    return m_bkgdScale[0];
  }
  double getBackgroundScale(int i) const;
  /// get the data scale factor
  inline double getDataScale() const {
    return m_dataScale;
  }

  /// access the (unscaled) background distribution
  inline double& bkgd(int i) {
    return m_bkgd[i];
  }
  /// access the (scaled) background distribution
  inline double bkgd(int i) const {
    return m_bkgd[i];
  }
  /// access the (unscaled) background distribution
  inline const double* bkgd() const {
    return m_bkgd;
  }

  /// access the (unscaled) data distribution
  inline double& data(int i) {
    return m_data[i];
  }
  /// access the (scaled) data distribution
  inline double data(int i) const {
    return m_data[i];
  }
  /// access the (unscaled) data distribution
  inline const double* data() const {
    return m_data;
  }

  /// access the (unscaled) signal distribution
  inline double& signal(int i) {
    return m_signal[i];
  }
  /// access the (scaled) signal distribution
  inline double signal(int i) const {
    return m_signal[i];
  }
  /// access the (unscaled) signal distribution
  inline const double* signal() const {
    return m_signal;
  }

  /// access the total statistical uncertainty in a given bin
  inline double getSignalStatErr(int i) const {
    return m_signalErr[i];
  }
  inline double getBkgdStatErr(int i) const {
    return m_bkgdErr[i];
  }

  /// get the N of bins
  inline int nbins() const {
    return n_bins;
  }
  /// get the min of the distribution
  inline double min() const {
    return m_min;
  }
  /// get the max of the distribution
  inline double max() const {
    return m_max;
  }
  /// get first independent variable (often a mass)
  inline int var1() const {
    return m_var1;
  }
  /// get first independent variable (if used)
  inline int var2() const {
    return m_var2;
  }
  /// get first independent variable (if used)
  inline int var3() const {
    return m_var3;
  }

  /// add up all signal
  double totSignal() const;
  /// add up all background
  double totBkgd() const;
  /// add up all data
  double totData() const;

  // fill data with S+B prediction
  void fillDataSB(bool sb = true);
  // fill data with B-only prediction
  void fillDataBOnly(bool bo = true);

  /// append another distribution to this distribution (ignoring the meaning of
  ///  max and min and just appending bins to make a longer array)
  int append(const SigBkgdDist& asbd);

  /** \brief calculate the Log-Likelihood Ratio for this distribution set.
      The LRR is defined as <i> -2 * sum (-s_i + d_i * log(1+s_i/b_i)) </i>.
  */
  double calculateDeltaLLR() const;
  double calculateDeltaLLRobs() const;
  double calculateLLR() const;
  double calculateLLR(int i) const;
  double calculateLLRmax() const;

  double calculateLLRw() const;

  double calculateLLRb() const;

  double calculateLLRsb() const;

  double calculateFmax() const;

  double calculateF(int i) const;

  //  double calculateChi2LLR(bool signal=false, double sigLLR=0) const;
  double calculateChi2LLR(bool signal, double sigLLR) const;

  //  double calculateChi2LLR_deletion(bool signal, vector<int> delbins) const;

  //MS Sept10
  double calculateChi2LLRderivative(int ipar, bool signal, double sigLLR, double* fluctMap);

  //  USED TO BE: inline double calculateSystDiff(const double* syst) const;
  double calculateSystDiff(std::vector<double> const& syst) const;

  //create new random values for systematics...
  void varySystematics();

  //zero fluctuations
  void clearSystFluctValues(bool zeroSyst);

  //get random values for systematics...
  // USED TO BE: const double* getFluctValues() const { return gl_systRand;}
  const double* getSystFluctValues() const {
    return &gl_systRand[0];
  }
  double getSystFluctValue(unsigned int i) const {
    return (i < 0 || i >= gl_systNames.size()) ? 0 : gl_systRand[i];
  }
  double getSystFluctValue(string name);

  //set fluctuations to specified values
  bool setSystFluctValue(unsigned int i, double value);
  bool setSystFluctValue(string name, double value);

  //toggle the constraint for a specified nuisance parameter (ie, float this parameter)
  bool getFloatFlag(unsigned int i) const {
    return (i < 0 || i >= gl_systNames.size()) ? 0 : gl_systFloat[i];
  }
  bool getFloatFlag(string name);
  bool setFloatFlag(unsigned int i, bool floatit);
  bool setFloatFlag(string name, bool floatit);

  //toggle the fit/no fit decision for a specified nuisance parameter (ie, do/don't fit this parameter)
  bool getSystFitFlag(unsigned int i) const {
    return (i < 0 || i >= gl_systNames.size()) ? 0 : gl_systFit[i];
  }
  bool getSystFitFlag(string name);
  bool setSystFitFlag(unsigned int i, bool fitit);
  bool setSystFitFlag(string name, bool fitit);

  //toggle the fit/no fit decision for a specified background (ie, do/don't fit this background)
  bool getBkgdFitFlag(unsigned int i) const {
    return (i < 0 || i >= gl_bkgdNames.size()) ? 0 : gl_bkgdFit[i];
  }
  bool getBkgdFitFlag(string name);
  bool setBkgdFitFlag(unsigned int i, bool fitit);
  bool setBkgdFitFlag(string name, bool fitit);

  //toggle the fit/no fit decision for a specified signal (ie, do/don't fit this signal)
  bool getSigFitFlag(unsigned int i) const {
    return (i < 0 || i >= gl_sigNames.size()) ? 0 : gl_sigFit[i];
  }
  bool getSigFitFlag(string name);
  bool setSigFitFlag(unsigned int i, bool fitit);
  bool setSigFitFlag(string name, bool fitit);

  //get number of total systematics, backgrounds, signals...
  inline int getNsyst() const {
    return gl_systNames.size();
  }
  inline int getNbkgds() const {
    return gl_bkgdNames.size();
  }
  inline int getNsignals() const {
    return gl_sigNames.size();
  }

  //Get the number of signal or background dists for a given channel
  unsigned int getNsigDist(string chanName) const;
  unsigned int getNbkgdDist(string chanName) const;

  //Get the names of all channels...
  vector<string> getChannelNames() const;

  //Get the name of the specified channel...
  string getChannelName(uint chan) const;

  //get name of systematic, background, signal
  inline string getSystName(int i) const {
    return gl_systNames[i];
  }
  inline string getBkgdName(int i) const {
    return gl_bkgdNames[i];
  }
  inline string getSigName(int i) const {
    return gl_sigNames[i];
  }

  //get index of systematic, background, signal
  int getSystIndex(string name);
  //  int getSigIndex(string name);
  //  int getBkgdIndex(string name);

  //fill arrays with fluctuated values, use someone else's map
  inline void fluctuate(std::vector<double> const& fluctMap, TrialPoint& lastMinuitTrialPoint) {
    fillArrays(false, true, fluctMap, lastMinuitTrialPoint);
    return;
  }

  inline void fluctuate(std::vector<double> const& fluctMap) {
    TrialPoint dummy;
    fillArrays(false, true, fluctMap, dummy);
    return;
  }

  //Read in fluctuations from a histogram
  void fluctuate(TH1D* fluctMap);

  //Read in fluctuations from a previous fit output file
  void fluctuate(TFile* fitFile, string fitType = "NULL");

  //fill arrays with fluctuated values, use our map
  inline void fluctuate() {
    varySystematics();
    TrialPoint dummy;
    fillArrays(false, true, gl_systRand, dummy);
    return;
  }

  //turn on/off signal fluctuation...default=true
  void fluctuateSignal(bool value) {
    m_sigFluct = value;
  }

  //is this SBD fluctuated?
  inline bool isVaried() const {
    return m_varied;
  }

  //is this SBD initialized?
  inline bool isInit() const {
    return m_init;
  }

  //fill all arrays with baseline model values
  void setBaselineModel(bool copy = false);

  //fill all arrays with adjusted model values
  void zeroModel();

  /// Calculate statistical uncertainty using contents of histograms.  Default = false
  void useHistoStats(bool use) {
    m_useStat = use;
  }
  inline bool usingHistoStats() const {
    return m_useStat;
  }

  //Add a signal, background, or data distribution to this collection
  void addSigDist(string chan, string sig, CollieDistribution* d);
  void addBkgdDist(string s, CollieDistribution* d);
  void addDataDist(string s, CollieDistribution* d);

  //Lookup for a given distribution within this collection
  const CollieDistribution* getSigDist(string chan, string sig) const;
  const CollieDistribution* getSigDist(string chan, unsigned int i) const;
  const CollieDistribution* getBkgdDist(string chan, unsigned int i) const;
  const CollieDistribution* getDataDist(string chan) const;

  /// rebin in log_10(signal/background)
  int rebin_sob(double min_l10_sob, double max_l10_sob, int bins);
  /// rebin in log_10(s/b), but shrink to fit actual min/max of s/b dist
  int rebin_sob_adaptive(double min_l10_sob, double max_l10_sob, int bins);
  /// rebin in log_10(s/b), but shrink to fit only 99% of signal
  int rebin_sob_minimized(double min_l10_sob, double max_l10_sob, int bins);

  /// compute expectation values in the arrays, so fit can be evaluated
  inline void fillArrays(bool fillData = false, bool varySyst = false,
                         std::vector<double> const& fluctMap = std::vector<double>()) {
    TrialPoint trialPoint;
    fillArrays(fillData, varySyst, fluctMap, trialPoint);
  }
  void fillArrays(bool fillData, bool varySyst, std::vector<double> const& fluctMap,
                  TrialPoint& lastMinuitTrialPoint);

  inline void setSystVar(double scl, double off) {
    m_adjSyst = true;
    m_systScale = scl;
    m_systOffset = off;
  }
  inline void clearSystVar() {
    m_adjSyst = false;
    m_systScale = 1;
    m_systOffset = 0;
  }
  void setSystCentralValue(uint i, double par);
  void setSystCentralValue(string name, double par);

  //Exclude a channel in the fit calculations...
  void excludeChannel(string chan);
  void excludeChannel(int chan);
  string getExcludedChannel() const {
    return gl_exclChan;
  }

  //Exclude bins from the fit calculations...
  void excludeBins(vector<uint> delBins);
  void clearExludedBins();

  //Calculate average background uncertainties
  void calculateBkgdUncertainty(double& total, double& average, double& sigWeighted);

  //verbosity switch
  inline void setVerbose(bool val) {
    m_verbose = val;
  }

  /// Add some ROOT accessors....
  TH1D* getSigHistogram();
  TH1D* getBkgdHistogram();
  TH1D* getDataHistogram();
  void drawPlots();
  //  void generateFitHistos(vector<double> signalSF, vector<double> fitValues, TH2* errMat, string type);
  void generateFitHistos(TH1* signalSF, TH1* fitValues, TH2* errMat, string type);

private:
  int getSOBbin(double sob);

  void generateFitMap(TH1* fitValuesIn, TH2* errMatIn, string type);

  void fillCondensedMap(bool copy = false);
  void sobMinMax(double& min, double& max);
  void sobMinMax_minimized(double& min, double& max);
  void fillRebinnedArrays(bool fillData, bool varySyst, const double* fluctMap);

  // functions to implement partial likelihood caching optimization
  void setNewPointSignalExpectations(const double* fluctMap, double* sig);
  void setNewPointBkgdExpectations(const double* fluctMap, double* sig);
  void setPerturbedSignalExpectations(std::vector<double> const& fluctMap,
                                      std::vector<double>::const_iterator onlyNonBaseS,
                                      double* sig);

  void setPerturbedBkgdExpectations(std::vector<double> const& fluctMap,
                                    std::vector<double>::const_iterator onlyNonBaseS, double* sig);

  void prepareSignalExclusionSums();

  void prepareBkgdExclusionSums();

  void fillSignalArrays(std::vector<double> const& fluctMap);

  collie::GaussianRandom m_randgaus;

  double m_systScale;
  double m_systOffset;
  double m_min, m_max;

  map<string, map<string, double>>::iterator iterSclO;
  map<string, double>::iterator iterSclI;
  map<string, map<string, double>> m_sigScale;
  bool m_sigScalesVary;
  double gl_cSigScale;

  vector<double> m_bkgdScale;
  double m_dataScale;

  int n_trueBins, n_bins;
  std::vector<uint> m_delBins;
  std::vector<uint> m_delBinsLookup;
  uint n_cacheData;
  int m_var1, m_var2, m_var3;

  string m_chan;

  bool m_adjSyst;
  bool m_useStat;
  bool m_deleteArrays;
  bool m_rebinned;
  bool m_varied;
  bool m_init;
  bool m_condensed;
  bool m_sigFluct;
  bool m_verbose;
  bool m_sbhypo;
  bool m_bohypo;

  int* m_rebinBins;

  double* m_signal;
  double* m_bkgd;

  double* m_signalErr;
  double* m_bkgdErr;

  //  double** d_signal;
  //  double** d_bkgd;

  double* m_signalParent;
  double* m_bkgdParent;
  double* m_data;

  TH1D* m_fitValueHist;
  TH2D* m_errMatHist;

  std::vector<double> gl_systRand; // USED TO BE: double* gl_systRand;

  bool* gl_systFloat;
  bool* gl_systFit;
  bool* gl_bkgdFit;
  bool* gl_sigFit;
  double* gl_systCenter;

  int gl_addedBins, gl_thisBin, gl_thisIdx;
  bool gl_fitAllSources;

  map<string, double> gl_fluctMap;
  map<string, int> gl_bkgdMap;
  map<string, int> gl_sigMap;
  map<string, double> gl_systCuts;
  map<string, double>::iterator fluctIter;

  vector<string> gl_systNames;
  vector<string> gl_bkgdNames;
  vector<string> gl_sigNames;
  string gl_exclChan;

  map<string, map<string, CollieDistribution*>>::iterator iterS;
  map<string, vector<CollieDistribution*>>::iterator iterB;
  map<string, CollieDistribution*>::iterator iterD;
  vector<CollieDistribution*>::iterator iterV;

  map<string, map<string, CollieDistribution*>> m_SigDist;
  map<string, vector<CollieDistribution*>> m_BkgdDist;
  map<string, CollieDistribution*> m_DataDist;

  bool m_signalExclusionSumsReady;
  bool m_bkgdExclusionSumsReady;
};

#endif // SigBkgdDist_HH_INCLUDED
