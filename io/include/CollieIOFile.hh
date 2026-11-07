#ifndef CollieIOFile_hh_included
#define CollieIOFile_hh_included

#include "CollieHistogram.hh"
#include "CollieQuickKEYS.hh"
#include "CollieDistribution.hh"
#include "CollieMasspoint.hh"
#include "CollieChannel.hh"
#include "CollieXsec.hh"
#include <iostream>
#include <list>

using namespace std;

class CollieIOFile {

  static constexpr double UNDEF = -1.2345e10;
  static constexpr double HMAX = 100000;

public:
  // Constructors, destructor
  CollieIOFile();
  CollieIOFile(string outfile, string chanName);
  ~CollieIOFile();

  //Novice flag switch
  // Certain methods and features are off limits to novice users.
  // Switching off this flag implicitly implies that users know
  // what they're doing.
  // flag = true; means the user is a novice (default)
  // flag = false; means the user is not a novice
  inline void setNoviceFlag(bool flag) {
    noviceFlag_ = flag;
    if (flag == true)
      systOvrd_ = false;
  }

  //Switch to override default limitations on systematics
  //  Distributions of shape systematics will be flattened in certain
  //  instances when this flag is set to false.  If the novice flag
  //  (see above) is set, this flag is to false.
  // flag = true;  All shape systematics are taken at face value.
  // flag = false; Shape systematics may be flattened.
  inline void setSystematicsOverride(bool flag) {
    systOvrd_ = flag;
    if (noviceFlag_ == true)
      systOvrd_ = false;
  }

  //Initialization method
  void initFile(string outfile, string chanName);

  // Write the collieIO output file and save the Final Variable inspection file
  void storeFile();

  // Given a properly initialized channel, specify the background model via a vector
  //   of background names.  Specify the names of model parameters in varNames.
  void createChannel(vector<string> bkgdNames);
  void createChannel(vector<string> bkgdNames, vector<string> sigNames);
  void createChannel(vector<string> bkgdNames, vector<string> sigNames, vector<string> varNames);

  // Specify the intended parameters of the input distributions.  All input distrubutions will be
  //   forced into this pre-specified template.  Any unrecoverable binning issues will be rejected.
  inline void setInputHist(double min, double max, int nbins) {
    histMinX_ = min;
    histMaxX_ = max;
    histBinsX_ = nbins;
  }
  inline void setInputHist2D(double minX, double maxX, int nbinsX, double minY, double maxY,
                             int nbinsY) {
    histMinX_ = minX;
    histMaxX_ = maxX;
    histBinsX_ = nbinsX;
    histMinY_ = minY;
    histMaxY_ = maxY;
    histBinsY_ = nbinsY;
  }

  // Optionally rebin all input distributions by the specified parameter
  inline void setRebin(int rebinX, int rebinY = 1) {
    rebinX_ = rebinX;
    rebinY_ = rebinY;
  }

  // Optionally cut a square region of the input distributions..
  inline void setCutoffs(double low, double high) {
    cutLowX_ = low;
    cutHighX_ = high;
  }
  inline void setCutoffs2D(double lowX, double highX, double lowY, double highY) {
    cutLowX_ = lowX;
    cutHighX_ = highX;
    cutLowY_ = lowY;
    cutHighY_ = highY;
  }

  // Optionally smooth input distributions.  Smoothing parameters are passed when creating
  //   mass points, set smoothing parameters to -1 to disable smoothing for that distribution
  inline void setSmooth(bool smth) {
    smooth_ = smth;
  }

  // Smoothing input information on statistics of lowest stats channel
  inline double getHistNorm() const {
    return histNorm_;
  }
  inline void setHistNorm(int norm) {
    histNorm_ = norm;
  }

  // verbosity flag
  inline void setVerbosity(bool verby) {
    verb_ = verby;
  }

  // Generate a set of input distributions corresponding to a
  //   specified model parameter (eg, mass) in 1-D or 2-D final variables
  //   Alpha parameters are smoothing parameters and are ignored if smoothing is
  //   disabled.  The order of the backgrounds in the vector<TH1D*> specifies the
  //   background indices used throughout.
  //   Model points can be parameterized in up to three dimensions.

  void createMassPoint(int par, TH1D* data, TH1D* sig, double alphaS, vector<TH1D*>& bkgd,
                       vector<double>& alphaB);
  void createMassPoint(int parX, int parY, TH1D* data, TH1D* sig, double alphaS,
                       vector<TH1D*>& bkgd, vector<double>& alphaB);
  void createMassPoint(int parX, int parY, int parZ, TH1D* data, TH1D* sig, double alphaS,
                       vector<TH1D*>& bkgd, vector<double>& alphaB);

  void createMassPoint(int par, TH1D* data, vector<TH1D*>& sig, vector<double>& alphaS,
                       vector<TH1D*>& bkgd, vector<double>& alphaB);
  void createMassPoint(int parX, int parY, TH1D* data, vector<TH1D*>& sig, vector<double>& alphaS,
                       vector<TH1D*>& bkgd, vector<double>& alphaB);
  void createMassPoint(int parX, int parY, int parZ, TH1D* data, vector<TH1D*>& sig,
                       vector<double>& alphaS, vector<TH1D*>& bkgd, vector<double>& alphaB);

  void createMassPoint2D(int par, TH2D* data, TH2D* sig, double alphaS, vector<TH2D*>& bkgd,
                         vector<double>& alphaB);
  void createMassPoint2D(int parX, int parY, TH2D* data, TH2D* sig, double alphaS,
                         vector<TH2D*>& bkgd, vector<double>& alphaB);
  void createMassPoint2D(int parX, int parY, int parZ, TH2D* data, TH2D* sig, double alphaS,
                         vector<TH2D*>& bkgd, vector<double>& alphaB);

  void createMassPoint2D(int par, TH2D* data, vector<TH2D*>& sig, vector<double>& alphaS,
                         vector<TH2D*>& bkgd, vector<double>& alphaB);
  void createMassPoint2D(int parX, int parY, TH2D* data, vector<TH2D*>& sig, vector<double>& alphaS,
                         vector<TH2D*>& bkgd, vector<double>& alphaB);
  void createMassPoint2D(int parX, int parY, int parZ, TH2D* data, vector<TH2D*>& sig,
                         vector<double>& alphaS, vector<TH2D*>& bkgd, vector<double>& alphaB);

  // Same as last four createMassPoint methods, but allow data to be stored as an
  //   indexed event list.  The event list can be indexed by a series of doubles,
  //   integers, or both.
  void createMassPoint_EvtList(int par, int nVarD, int nVarI, TH1D* sig, double alphaS,
                               vector<TH1D*>& bkgd, vector<double>& alphaB);
  void createMassPoint_EvtList(int parX, int parY, int nVarD, int nVarI, TH1D* sig, double alphaS,
                               vector<TH1D*>& bkgd, vector<double>& alphaB);
  void createMassPoint_EvtList(int parX, int parY, int parZ, int nVarD, int nVarI, TH1D* sig,
                               double alphaS, vector<TH1D*>& bkgd, vector<double>& alphaB);

  void createMassPoint2D_EvtList(int par, int nVarD, int nVarI, TH2D* sig, double alphaS,
                                 vector<TH2D*>& bkgd, vector<double>& alphaB);
  void createMassPoint2D_EvtList(int parX, int parY, int nVarD, int nVarI, TH2D* sig, double alphaS,
                                 vector<TH2D*>& bkgd, vector<double>& alphaB);
  void createMassPoint2D_EvtList(int parX, int parY, int parZ, int nVarD, int nVarI, TH2D* sig,
                                 double alphaS, vector<TH2D*>& bkgd, vector<double>& alphaB);

  void createMassPoint_EvtList(int par, int nVarD, int nVarI, vector<TH1D*>& sig,
                               vector<double>& alphaS, vector<TH1D*>& bkgd, vector<double>& alphaB);
  void createMassPoint_EvtList(int parX, int parY, int nVarD, int nVarI, vector<TH1D*>& sig,
                               vector<double>& alphaS, vector<TH1D*>& bkgd, vector<double>& alphaB);
  void createMassPoint_EvtList(int parX, int parY, int parZ, int nVarD, int nVarI,
                               vector<TH1D*>& sig, vector<double>& alphaS, vector<TH1D*>& bkgd,
                               vector<double>& alphaB);

  void createMassPoint2D_EvtList(int par, int nVarD, int nVarI, vector<TH2D*>& sig,
                                 vector<double>& alphaS, vector<TH2D*>& bkgd,
                                 vector<double>& alphaB);
  void createMassPoint2D_EvtList(int parX, int parY, int nVarD, int nVarI, vector<TH2D*>& sig,
                                 vector<double>& alphaS, vector<TH2D*>& bkgd,
                                 vector<double>& alphaB);
  void createMassPoint2D_EvtList(int parX, int parY, int parZ, int nVarD, int nVarI,
                                 vector<TH2D*>& sig, vector<double>& alphaS, vector<TH2D*>& bkgd,
                                 vector<double>& alphaB);

  /// if possible, access the event list (read-only) (NULL if not available)
  const CollieEventList* getDataEventList(int massX, int massY = -1) const;
  /// if possible, access the event list in read/write mode (NULL if not available)
  CollieEventList* getDataEventListMutable(int massX, int massY = -1);

  /// if possible, store data event list
  void finalizeDataEventList(int massX, int massY = -1);

  // Override the data distribution specified via createMass point.
  //   Useful when generating interpolated mass points.
  void setDataDist(TH1D* data, int parX, int parY = -1, int parZ = -1);
  void setDataDist2D(TH2D* data, int parX, int parY = -1, int parZ = -1);

  // Given a set of input mass points, generate a grid (1D or 2D) based
  //   The grid parameter specifies the granularity (step size)
  void interpolateMassGrid(int dataFlag, int grid, int massX_start, int massX_end,
                           int massY_start = -1, int massY_end = -1);

  // Create a signal or background systematic with a constant uncertainty for all bins
  //  of the input distribution.  Must specify mass point and background index, if appropriate.
  //  Positive fluctuations are given in the fractional deviation in the positive direction
  //  posFluct = [ ((+1sigma)-Nominal)/Nominal ]
  //  Negative fluctuations are given in the fractional deviation in the negative direction
  //  negFluct = [ (Nominal - (-1sigma))/Nominal   **Collie handles the sign internally**
  void createFlatSigSystematic(int idx, string syst, double posFluct, double negFluct, int parX,
                               int parY = -1, int parZ = -1);
  void createFlatSigSystematic(string syst, double posFluct, double negFluct, int parX,
                               int parY = -1, int parZ = -1);
  void createFlatBkgdSystematic(int bkgdIndex, string syst, double posFluct, double negFluct,
                                int parX, int parY = -1, int parZ = -1);

  void createFlatSigSystematic2D(int idx, string syst, double pos, double neg, int parX,
                                 int parY = -1, int parZ = -1);
  void createFlatSigSystematic2D(string syst, double pos, double neg, int parX, int parY = -1,
                                 int parZ = -1);
  void createFlatBkgdSystematic2D(int bkgdIndex, string syst, double pos, double neg, int parX,
                                  int parY = -1, int parZ = -1);

  // Create a signal or background systematic with per-bin fluctuations specified by the per-bins
  //  values in the input histograms. Must specify mass point and background index, if appropriate.
  //  Positive fluctuations are given in the fractional deviation in the positive direction
  //  posFluct = [ ((+1sigma)-Nominal)/Nominal ] (per-bin)
  //  Negative fluctuations are given in the fractional deviation in the negative direction
  //  negFluct = [ (Nominal - (-1sigma))/Nominal (per-bin)  **Collie handles the sign internally**
  void createSigSystematic(string syst, TH1D* posFluct, TH1D* negFluct, int parX, int parY = -1,
                           int parZ = -1);
  void createSigSystematic(int sigIndex, string syst, TH1D* posFluct, TH1D* negFluct, int parX,
                           int parY = -1, int parZ = -1);
  void createBkgdSystematic(int bkgdIndex, string syst, TH1D* posFluct, TH1D* negFluct, int parX,
                            int parY = -1, int parZ = -1);
  void createSigSystematic2D(string syst, TH2D* posFluct, TH2D* negFluct, int parX, int parY = -1,
                             int parZ = -1);
  void createSigSystematic2D(int sigIndex, string syst, TH2D* posFluct, TH2D* negFluct, int parX,
                             int parY = -1, int parZ = -1);
  void createBkgdSystematic2D(int bkgdIndex, string syst, TH2D* posFluct, TH2D* negFluct, int parX,
                              int parY = -1, int parZ = -1);

  // Create a signal or background systematic derived from the final variable shapes obtained after
  //  varying the nuisance parameter in question.
  //  posShape = histogram obtained via +1sigma deviation
  //  posShape = histogram obtained via -1sigma deviation
  //
  //  The "scaleFactor" parameter allows users to change the amplitude of the shape variation
  //  (eg, 0.5 means 50% of original)
  //
  //  The "norm" parameter allows users to remove any normalization change in the shape variation
  //  relative to the nominal
  //
  //  The "flatten" parameter will convert the alternative shape change into a flat systematic
  //  (same value for all bins) using the overall change in normalization.  By construction, the "norm"
  //  parameter will be set to false if the "flatten" parameter is set to true.
  //
  //  The "testShape" parameter will determine whether the systematic shapes will be evaluated using bins
  //  with roughly equivalent statistical uncertainty, thus eliminating the potential for wild fluctuations
  //  due to low statistics populations.

  void createShapeSigSystematic(int idx, string syst, TH1D* posShape, TH1D* negShape, int parX,
                                int parY = -1, int parZ = -1, double scaleFactor = 1.0,
                                bool norm = false, int flatten = 0, int testShape = 1);
  void createShapeSigSystematic(string syst, TH1D* posShape, TH1D* negShape, int parX,
                                int parY = -1, int parZ = -1, double scaleFactor = 1.0,
                                bool norm = false, int flatten = 0, int testShape = 1);

  void createShapeBkgdSystematic(int bkgdIndex, string syst, TH1D* posShape, TH1D* negShape,
                                 int parX, int parY = -1, int parZ = -1, double scaleFactor = 1.0,
                                 bool norm = false, int flatten = 0, int testShape = 1);

  void createShapeSigSystematic2D(string syst, TH2D* posShape, TH2D* negShape, int parX,
                                  int parY = -1, int parZ = -1, double scaleFactor = 1.0,
                                  bool norm = false, int flatten = 0, int testShape = 1);

  void createShapeSigSystematic2D(int idx, string syst, TH2D* posShape, TH2D* negShape, int parX,
                                  int parY = -1, int parZ = -1, double scaleFactor = 1.0,
                                  bool norm = false, int flatten = 0, int testShape = 1);

  void createShapeBkgdSystematic2D(int bkgdIndex, string syst, TH2D* posShape, TH2D* negShape,
                                   int parX, int parY = -1, int parZ = -1, double scaleFactor = 1.0,
                                   bool norm = false, int flatten = 0, int testShape = 1);

  //  void setPoissonBkgd(int bkgdIndex, double bkgdNorm, double bkgdErrP, double bkgdErrN, int massX, int massY=-1);

  // Double (ie, don't constrain) the specified systematic in all fits
  void setBkgdFloatFlag(int bkgdIndex, string syst, bool doubleIt, int parX, int parY = -1,
                        int parZ = -1);
  void setSigFloatFlag(string syst, bool doubleIt, int parX, int parY = -1, int parZ = -1);

  // Specify that this systematic should have a log-normal PDF (rather than Gaussian)
  void setLogNormalFlag(string syst, bool doLogNorm, int parX, int parY = -1, int parZ = -1);

  //  Interpolation routines
  void interpolate2D(int axis, CollieHistogram2d* h1, double p1, double xsec1,
                     CollieHistogram2d* h2, double p2, double xsec2, CollieHistogram2d* hf,
                     double pf, double xsecF);
  void interpolate(CollieHistogram* h1, double p1, double xsec1, CollieHistogram* h2, double p2,
                   double xsec2, CollieHistogram* hf, double pf, double xsecF);

  void interpolate(TH1D* h1, double p1, double xsec1, TH1D* h2, double p2, double xsec2, TH1D* hf,
                   double pf, double xsecF);

  bool getSmoothedROOT(double alpha, const TH1D* in, TH1D* out);
  bool getSmoothedROOT(double alpha, const CollieHistogram* in, TH1D* out);
  //  bool getSmoothedROOT2D(double alpha, TH2D* in, TH2D& out);
  //  bool getSmoothedROOT2D(double alpha, CollieHistogram& in, TH2D& out);

  bool getSmoothedCollie(double alpha, const TH1D* in, CollieHistogram* out);
  bool getSmoothedCollie(double alpha, const CollieHistogram* in, CollieHistogram* out);
  bool getSmoothedCollie2D(double alpha, const TH2D* in, CollieHistogram2d* out);
  bool getSmoothedCollie2D(double alpha, const CollieHistogram2d* in, CollieHistogram2d* out);

  bool getGroupedBinCollie(double alpha, const TH1D* in, CollieHistogram* out);
  bool getGroupedBinCollie(double alpha, const CollieHistogram* in, CollieHistogram* out);
  bool getGroupedBinCollie2D(double alpha, const TH2D* in, CollieHistogram2d* out);
  bool getGroupedBinCollie2D(double alpha, const CollieHistogram2d* in, CollieHistogram2d* out);

  bool getROOT(const CollieHistogram* in, TH1D* out);
  bool getROOT2D(const CollieHistogram2d* in, TH2D* out);
  bool getCollie(const TH1D* in, CollieHistogram* out);
  bool getCollie2D(const TH2D* in, CollieHistogram2d* out);

  // utility routines
  inline void doCD() {
    if (channel_)
      channel_->doCD();
  }

public: //getter methods
  const double getCutLowX() const {
    return cutLowX_;
  }
  const double getCutHighX() const {
    return cutHighX_;
  }
  const double getCutLowY() const {
    return cutLowY_;
  }
  const double getCutHighY() const {
    return cutHighY_;
  }
  const int getRebinX() const {
    return rebinX_;
  }
  const int getRebinY() const {
    return rebinY_;
  }
  const double getHistMinX() const {
    return histMinX_;
  }
  const double getHistMaxX() const {
    return histMaxX_;
  }
  const int getHistBinsX() const {
    return histBinsX_;
  }
  const double getHistMinY() const {
    return histMinY_;
  }
  const double getHistMaxY() const {
    return histMaxY_;
  }
  const int getHistBinsY() const {
    return histBinsY_;
  }
  int getBkgdIndex(const string name);
  int getSigIndex(const string name);
  const int getNBkgd() const {
    return nBkgd_;
  }
  const int getNSig() const {
    return nSig_;
  }

  //Rebinning Methods
  //  This method provides a tool to generate a 1D histogram binning that
  //  is safe WRT background fluctuations.  The design forces a minimum
  //  number of events in each bin with an additional requirement
  //  that the statistical uncertainty on the bin contents meets a minimum
  //  requirement.  The number of output bins and histogram range will match
  //  the number of bins specified in the cfile->setInputHist(...) method.
  //  Bins will not be split, thus if the number of input bins is smaller
  //  than the number of desired output bins, you will get empty bins.
  //
  //  There are a few input parameters:
  //
  //  TH1D* bkgd: The total background desired for rebinning purposes.
  //  TH1D* sig:  The total signal desired for rebinning purposes.
  //
  //  double highS: The fraction of total bins in the high S/B region.
  //                This region will require only a minimum number of
  //                total events per bin and a minimum statistical uncertainty.
  //
  //  double midS:  The fraction of total bins in the middle S/B region.
  //                This region has a linearly increasing bin content
  //                requirement.
  //
  //  Note: The remaining bins (Ndesired - Nhi - Nmid) will include the
  //        total remaining background.  The algorithm will choose either
  //        a linear or quadratically increasing bin content, depending on
  //        the number of events remaining.
  //
  //  double minSB: The minimum number of signal+background events in a
  //                given bin.
  //
  //  double minB:  The minimum number of background events in a given bin.
  //
  //  double minStatSB: The minimum fractional statistical uncertainty on
  //                    the sum of signal+bkgd in a given bin.
  //
  //  double minStatB:  The minimum fractional statistical uncertainty on
  //                    the total background in a given bin.
  int generateBinMap(TH1D* bkgd, TH1D* signal, string style, double highS = 0.20,
                     double midS = 0.40, double minSB = 0.1001, double minB = 1e-4,
                     double minStatSB = 0.201, double minStatB = 0.251);

  void setSystPriors(TH1D* inHist);

private:
  void createSigSystematic_Intl(int sigIndex, string syst, TH1D* posFluct, TH1D* negFluct, int parX,
                                int parY = -1, int parZ = -1);
  void createBkgdSystematic_Intl(int bkgdIndex, string syst, TH1D* posFluct, TH1D* negFluct,
                                 int parX, int parY = -1, int parZ = -1);
  void createSigSystematic2D_Intl(int sigIndex, string syst, TH2D* posFluct, TH2D* negFluct,
                                  int parX, int parY = -1, int parZ = -1);
  void createBkgdSystematic2D_Intl(int bkgdIndex, string syst, TH2D* posFluct, TH2D* negFluct,
                                   int parX, int parY = -1, int parZ = -1);

  bool check(bool twoD = false);
  void applyCuts(TH1* hist);
  void applyCuts(CollieHistogram* hist);
  bool checkROOTHisto(TH1* hist, bool norm);
  bool checkForShape(CollieHistogram* ref, TH1* pos, TH1* neg, string syst = "");

  /// Given a nominal and fluctuated histogram, return an output histogram
  /// describing the fractional uncertainties per bin in global bins of equal
  /// normalization.  The prob argument describes the ratio of the maximum bin
  /// to minimum bin in the rebinned distribution.  The norm argument determines
  /// if the fractional difference will include normalization differences. The
  /// output maintains the same number of bins as the input.
  void EqualProbDifference(TH1* nominal, TH1* fluctuation, TH1* output, string syst,
                           bool norm = false, double prob = 4.0);
  void EqualProbDifference(CollieHistogram* nominal, TH1* fluctuation, TH1* output, string syst,
                           bool norm = false, double prob = 4.0);

  bool testShapeSystematics(CollieHistogram* nom, TH1* fluct, TH1* out, bool norm, string syst,
                            double cutoff = 0.035);

  void checkRates(CollieHistogram* ref, TH1* pos, TH1* neg, string syst = "");
  void checkBins(TH1D* data, vector<TH1D*>& sig, vector<TH1D*>& bkgd);
  //  void checkBins(TH1D* data, TH1D* sig, vector<TH1D*>& bkgd);
  void checkBins2D(TH2D* data, vector<TH2D*>& sig, vector<TH2D*>& bkgd);

  double getSystErr(double nv, double sv, double nerr, double serr);
  double totErr(CollieHistogram* hist, int i, int j);
  double totErr(TH1* hist, int i, int j);

  void report();

  void fillDist(CollieHistogram* hist, CollieDistribution* dist);

  void interpolateMassList(int dataFlag, int grid, int m1, int m2, list<int> l1, int xMass = 0,
                           int yMass = 0);
  void oneStepPass(int grid, int massLo, int massHi, list<int> masses, int massX = 0,
                   int massY = 0);

  inline int mIndex(int pX, int pY = -1, int pZ = -1) const {
    return pX + (pY < 0 ? 0 : (int)(pY * 10000)) + (pZ < 0 ? 0 : (int)(pZ * 10000000));
  }
  inline void mLookUp(const int idx, int& pX, int& pY, int& pZ) {
    pX = idx % 10000;
    pY = int(idx / 10000) % 1000;
    pZ = int(idx / 10000000);
  }

  bool logMassPoint(int pX, int pY = -1, int pZ = -1);
  bool checkMassPoint(int pX, int pY = -1, int pZ = -1);

  void generateSmoothedHistos();

  TH1D* h1FMT(TH1D* hist);
  // TH2D* h2FMT(TH2D* hist);

  int generateMVABinMap(TH1D* bkgd, TH1D* signal, double highS, double midS, double minSB,
                        double minB, double minStatSB, double minStatB);
  int generateInvMassBinMap(TH1D* bkgd, TH1D* signal, double highS, double midS, double minSB,
                            double minB, double minStatSB, double minStatB);
  TH1D* getBinMapROOT(const TH1D* in, const char* name);

  double getSlope(double integral, double minMed, int lBins) {
    double base = integral - lBins * minMed * 1.5;
    if (base < 0)
      base = 0;
    base /= (0.5 * lBins);
    return base / lBins;
  }

private:
  ///general control booleans
  bool verb_;
  bool init_;
  bool smooth_;
  bool noviceFlag_;
  bool systOvrd_;
  bool usingBinMap_;
  bool usingSystPriors_;
  TH1D* systPriors_;

  ///histogram treatment
  double cutLowX_;
  double cutHighX_;
  double cutLowY_;
  double cutHighY_;
  int rebinX_;
  int rebinY_;
  double histMinX_;
  double histMaxX_;
  int histBinsX_;
  double histMinY_;
  double histMaxY_;
  int histBinsY_;

  ///Input Distributions
  unsigned int nBkgd_;
  unsigned int nSig_;
  map<int, map<string, CollieHistogram*>> bkgd_;
  map<int, map<string, CollieHistogram*>> sig_;
  map<int, CollieHistogram*> allBkgd_;
  map<int, CollieHistogram*> data_;

  map<int, map<string, double>> bkgdAlpha_;
  map<int, map<string, double>> sigAlpha_;

  map<int, int> parList_; //2D input variable tracker...
  vector<TH1*> systHists_;
  map<string, int> systNames_;

  ///other parameters
  int smIncr_;
  double histNorm_;
  map<int, int> binMap_;
  vector<double> binEdges_;

  //Cross section calculator
  CollieXsec xsec;

  ///IO parameters
  TFile* outFile_;
  string outName_;
  CollieChannel* channel_;
  string chanName_;
};

#endif // CollieIOFile_hh_included
