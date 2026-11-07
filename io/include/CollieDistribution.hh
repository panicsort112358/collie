#ifndef CollieDistribution_hh_included
#define CollieDistribution_hh_included

#include <TObjString.h>
#include <TNamed.h>
#include <TArrayD.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <CollieHistogram.hh>
#include "InfoForEfficiencyCalculation.hh"

typedef unsigned int uint;

class CollieMasspoint;

/** \brief Final variable distribution (1d or 2d) in the CollieIO system.

This class represents a single final variable distribtion.  The
CollieMasspoint class combines several distributions of different types
together.

To maximize the use of the limited precision of doubles, the distribution
is stored as an array which is normalized to 1.0.  The efficiency
and cross-sections are stored separately.

The distribution also contains systematic error estimates, which
are currently under development.
*/

class CollieDistribution : public TNamed {

public:
  /// special value for the model rate where the luminosity should be ignored (-1000.0f)
  static const double rate_IgnoreLuminosity;

  /// constuctor for streaming
  CollieDistribution(); // for streaming
  /// copy constructor
  CollieDistribution(const CollieDistribution& d);
  /// copy assignment operator
  CollieDistribution& operator=(const CollieDistribution& d);

  /// if the distribution is read/write (mutable), then set the efficiency for the distribution
  void setEfficiency(double efficiency);
  inline double getEfficiency() const {
    return fEfficiency;
  }

  /// if the distribution is read/write (mutable), then set the normalized contents of a given bin of the distribution
  void setNormalizedBinValue(double value, int i, int j = -1);
  /// if the distribution is read/write (mutable), then set the statistical uncertainty of a given bin of the dist
  void setBinStatErr(double value, int i, int j = -1);
  /// if the distribution is read/write (mutable), then set the cross-section for this distribution for a given model (by model id number)
  void setModelXsec(double value, int nmodel = 0);
  /// if the distribution is read/write (mutable), then set the cross-section for this distribution for a given model (by model name)
  void setModelXsec(double value, const char* name);

  /// creates a new Root Histogram (1D) and fills it with the contents of this distribution
  TH1* draw(string title = "") const;
  /// creates a new Root Histogram (2D) and fills it with the contents of this distribution
  TH2* draw2D(string title = "") const;

  /** \brief Fill the distribution using the contents of the given 1d histogram

      This routine assumes the efficiency to be the sum of the histogram bins.
      If this is not the case, it can be corrected using setEfficiency().

      If the histogram is two bins smaller than the distribution, the overflow
      contents are placed in the first and last bin.  If the histogram is the
      same size as the distribution, the contents of the overflow bins are
      added to the first and last bins of the distribution.

      \return False if there is a problem such as bin-count mismatching or an invalid histogram.
              Indicate the rebinning factor specified on the IO interface.
   */
  int fillFromHistogram(const TH1* histo, int rebin = -1);
  int fillFromHistogram(const CollieHistogram* histo, int rebin = -1);

  /** \brief Fill the distribution using the contents of the given 2d histogram
      This routine assumes the efficiency to be the sum of the histogram bins.
      If this is not the case, it can be corrected using setEfficiency().

      The overflow bins are <i>thrown away</i> in this routine.  The user must
      handle adding their content to visible bin.

      \return False if there is a problem such as bin-count mismatching or an invalid histogram.
              Indicate the rebinning factor specified on the IO interface.
   */
  int fillFromHistogram(const TH2* histo, int rebinX = -1, int rebinY = -1);
  int fillFromHistogram(const CollieHistogram2d* histo, int rebinX = -1, int rebinY = -1);

  /// add up the efficiency bins explicitly and multiply it by the recorded efficiency to give a total efficiency
  double sumEfficiency();
  // partial sum
  double sumEfficiency(int i, int j);

  /// get the number of bins in the X direction
  inline int getNXbins() const {
    return fNXbins;
  }
  /// get the number of bins in the Y direction
  inline int getNYbins() const {
    return fNYbins;
  }

  /// get the maximum edge value in the X direction in each distribution
  inline double getMaxX() const {
    return fMaxX;
  }
  /// get the minimum edge value in the X direction in each distribution
  inline double getMinX() const {
    return fMinX;
  }
  /// get the maximum edge value in the Y direction in each distribution

  inline double getMaxY() const {
    return fMaxY;
  }
  /// get the minimum edge value in the Y direction in each distribution
  inline double getMinY() const {
    return fMinY;
  }

  inline bool linearized() const {
    return fLinearized;
  }

  /// get the normalized bin value for a given bin
  //  double getNormalizedBinValue(int i, int j) const;
  // A. Haas 7/10/08
  inline double getNormalizedBinValue(int i, int j) const {
    if (i < 0 || i >= fNXbins)
      return -1;
    if (fNYbins <= 1) {
      if (fLinearized)
        return fLinBins[i];
      else
        return fBins[i];
    }

    if (j >= fNYbins)
      return -1;
    if (fLinearized)
      return fLinBins[i + j * fNXbins];
    else
      return fBins[i + j * fNXbins];
  }

  inline double getBinErrValue(int i, int j) const {
    if (i < 0 || i >= fNXbins)
      return -1;
    if (fNYbins <= 1) {
      if (fLinearized)
        return fLinBinStat[i];
      else
        return fBinStat[i];
    }

    if (j >= fNYbins)
      return -1;
    if (fLinearized)
      return fLinBinStat[i + j * fNXbins];
    else
      return fBinStat[i + j * fNXbins];
  }

  /// get the efficiency value (normalized bin value * efficiency) for a given bin
  inline double getEfficiency(int i, int j = -1) const {
    return getNormalizedBinValue(i, j) * fEfficiency;
  }

  // get the smeared efficiency for fluctuations
  inline double getSmearedEfficiency(double rS, double sP, double sN) const {
    if (rS < 0)
      return getSmearedEfficiencyN(rS, sP, sN);
    else
      return getSmearedEfficiencyP(rS, sP, sN);
  }

  // get the smeared efficiency for positive fluctuations
  inline double getSmearedEfficiencyP(double rS, double sP, double sN) const {
    return getSmearBridge(1.0 / (1.0 + 3.0 * rS), rS, sP, sP, sN);
  }

  // get the smeared efficiency for negative fluctuations
  inline double getSmearedEfficiencyN(double rS, double sP, double sN) const {
    return getSmearBridge(1.0 / (1.0 - 3.0 * rS), rS, sN, sP, sN);
  }

  // bridge asymmetric uncertainties with a rolled off quadratic matching solution
  inline double getSmearBridge(double rF, double rS, double sD, double sP, double sN) const {

    return ((1 - rF) * rS * sD + rF * (rS * (sP + sN) * 0.5 + rS * rS * (sP - sN) * 0.5));
  }

  inline double getEfficiencyVaried(int i, int j, const double* fluctList) const {
    if (fNYbins <= 1)
      return getNormalizedBinValueVaried(i, fluctList) *
             fEfficiency; //if we're 1D, just use the 1D function-- ACH
    return getNormalizedBinValueVaried(i, j, fluctList) * fEfficiency;
  }

  /// get the normalized bin value for a given bin (1D)
  double getNormalizedBinValueVaried(const int i, const double* fluctList) const;

  /// get the normalized bin value for a given bin (2D)
  double getNormalizedBinValueVaried(const int i, const int j, const double* fluctList) const;

  /// Accumulate the scaled efficiency value for each bin
  // (normalized bin value * efficiency * overall scaling factor)
  void addBinEfficiencies(const double* fluctMap, const double sigScale, double* sig,
                          const uint n_bins);

  /// Accumulate the scaled efficiency value for each bin
  // (normalized bin value * efficiency * overall scaling factor)
  // using shortcut made possible by knowing the sole changing s
  // relative to a base point.
  void addBinEfficiencies(const int perturbed_s, const double s_value, const double sigScale,
                          double* sig, const uint n_bins);

  /// Prepare data structures for fast computation using cached partial products
  void prepareExclusionSums();

  /// get the systematic value for a given bin
  double getBinSystValue(int syst, int i, int j, const double* in) const;

  /// get the statistical uncertainty given by the input histogram...if any
  double getBinStatErr(int i, int j) const;

  //  inline int getNmodels() const { return fNmodels; }
  /// get the value of the cross-section in a given model (by model id number)
  inline double getModelXsec(int i) const {
    return (i < 0 || i >= fNmodels) ? -1 : fModelXsecs[i];
  }
  /// get the value of the cross-section in a given model (by model name)
  inline double getModelXsec(const char* name) const {
    return getModelXsec(lookupModel(name));
  }
  /// get the number of model cross-sections
  inline int getNModelXsec() const {
    return fNmodels;
  }

  /// access the mass point object which contains this distribution
  const CollieMasspoint* getMasspoint() const {
    return p_MassPoint;
  }

  inline int getNsystematics() const {
    return fSystNames.GetEntriesFast();
  }

  void addSystematic(const char* systname, TH1D* pos, TH1D* neg);
  void addSystematic2D(const char* systname, TH2D* pos, TH2D* neg);

  ///Flags for floating systematics
  void setFloatFlag(string systname, bool floatIt);
  bool getFloatFlag(string systname);
  void getFloatFlagList(vector<string> inputNames, bool* floatMap);

  //Flags for using log-normal vs gaussian
  void setLogNormalFlag(string systname, bool lnIt);
  bool getLogNormalFlag(string systname);
  void getLogNormalFlagList(vector<string> inputNames, bool* LogNormalMap);

  void setPoissonFlag(double norm, double errP, double errN);

  void linearize(vector<string> inputNames);

  int getSystIndex(string syst) const;
  bool hasSystematic(string syst) const;

  inline TString getSystNameT(int i) const {
    return (i < 0 || i >= fSystNames.GetEntriesFast()) ? TString("")
                                                       : ((TObjString*)fSystNames[i])->String();
  }
  inline string getSystName(int i) const {
    return (i < 0 || i >= fSystNames.GetEntriesFast())
               ? string("")
               : string(((TObjString*)fSystNames[i])->String().Data());
  }

  inline TH1D* getPositiveSyst(int i) const {
    return (i < 0 || i >= fSystNames.GetEntriesFast()) ? NULL : (TH1D*)fSystematicsPos[i];
  }
  inline TH1D* getNegativeSyst(int i) const {
    return (i < 0 || i >= fSystNames.GetEntriesFast()) ? NULL : (TH1D*)fSystematicsNeg[i];
  }

  inline TH2D* getPositiveSyst2D(int i) const {
    return (i < 0 || i >= fSystNames.GetEntriesFast()) ? NULL : (TH2D*)fSystematicsPos[i];
  }
  inline TH2D* getNegativeSyst2D(int i) const {
    return (i < 0 || i >= fSystNames.GetEntriesFast()) ? NULL : (TH2D*)fSystematicsNeg[i];
  }

  inline void setMutable(bool mut) {
    fMutable = mut;
  }

  // print information about this distribution
  void print(map<string, int>& count, map<string, double>& posSyst,
             map<string, double>& negSyst) const;

  // methods for capturing and restoring fLinSyst sections when
  // comparing results of shortcut to full computation
  bool s_is_active(int s) {
    return (fSystIndexInner[s] >= 0);
  }
  std::vector<InfoForEfficiencyCalculation>& effInfoVector(int s) {
    return fLinSyst[fSystIndexInner[s]];
  }

  /// destructor.
  ~CollieDistribution();

private:
  friend class CollieMasspoint;

  CollieDistribution(const char* name, int nX, double mX, double MX, int nY, double mY, double MY,
                     int nModels);
  bool fMutable;                      //! whether this distribution can be changed
  const CollieMasspoint* p_MassPoint; //! pointer to the mass point (set by CollieMasspoint)

  // Name is based on CollieMasspoint name, distribution type, and more
  Double_t fMinX;       // minimum value of distribution in X
  Double_t fMaxX;       // maximum value of distribution in X
  Double_t fMinY;       // minimum value of distribution in Y (if relevant)
  Double_t fMaxY;       // maximum value of distribution in Y (if relevant)
  Double_t fEfficiency; // total efficiency

  Int_t fNXbins;       // number of bins in the X direction
  Int_t fNYbins;       // number of bins in the Y direction (<=1 implies one D)
  Int_t fTrueBinCount; // count of total bins
  TArrayD
      fBins; // the distribution itself (linearized)  Kept normalized to 1.0 (see fEfficiency for efficiency normalization)
  TArrayD
      fBinStat; // the statistical uncertainties on each bin (linearized), as obtained from the input histogram
  Int_t fNmodels;      // count of the number of models
  TArrayD fModelXsecs; // xsecs for the models

  // names of systematics, for correlations
  TObjArray fSystNames;
  // systematics values (histos) in positive direction
  TObjArray fSystematicsPos;
  // systematics values (histos) in negative direction
  TObjArray fSystematicsNeg;
  //Should this systematic be floated or use the prior info?
  TObjArray fFloatFlag;
  //Should this systematic be floated or use the prior info?
  TObjArray fLogNormalFlag;

  // Is this bkgd from a Poisson source??
  bool fPoissonFlag;
  double fPoissonNorm;
  double fPoissonErrPos;
  double fPoissonErrNeg;

  ///non-ROOT tools for fast linearized systematic calculations...
  Int_t fNsyst;
  int* fSystIndexOuter; //!
  int* fSystIndexInner; //!

  bool* fLinFloat;      //!
  bool* fLinLogN;       //!
  double** fLinSystPos; //!
  double** fLinSystNeg; //!

  // replacement for fLinSystPos and fLinSystNeg, plus remembered values
  std::vector<std::vector<InfoForEfficiencyCalculation>> fLinSyst; //!
  std::vector<double> fEfficiencySums;                             //!
  //  double fRF;
  double* fLinBins;    //!
  double* fLinBinStat; //!
  bool fLinearized;

  // Function to do a few checks on Histo stat errors, M Owen 16-3-08
  bool checkStats(const TH1* h, int verbose = 1) const;
  bool checkStats(const CollieHistogram* h, int verbose = 1) const;

  int lookupModel(const char* name) const;

  // Inner-loop computation functions for various interpolation strategies
  void binEfficienciesQbridge(double* product, const double scale, const double active_s,
                              std::vector<InfoForEfficiencyCalculation>& e);

  void binEfficienciesQbridgeLinLogN(double* product, const double scale, const double active_s,
                                     std::vector<InfoForEfficiencyCalculation>& e);

  inline double assurePositive(double x) {
    return x > 0.0 ? x : 1.0e-5;
  }

  ClassDef(CollieDistribution, 1)
};

#endif // CollieDistribution_hh_included
