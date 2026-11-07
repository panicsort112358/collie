#ifndef ProfileLH_2D_HH_INCLUDED
#define ProfileLH_2D_HH_INCLUDED

#include "SigBkgdDist.hh"
#include <memory>

namespace collie {
class ProfileLikelihoodEngine;
}

/** \brief Signal, background, and data together as distributions. */
class ProfileLH_2D {

public:
  /// create an empty ProfileLH_2D
  ProfileLH_2D();
  /// create a ProfileLH_2D associated with a SigBkgdDist
  explicit ProfileLH_2D(SigBkgdDist* asbd);
  /// destructor
  ~ProfileLH_2D();

  // FORTRAN LUN for MINUIT printouts
  void setLUN(int lun) {
    m_lun = lun;
  }

  //specify the model to be fit
  void setModel(SigBkgdDist* asbd);

  //set/get maximum number of iterations for warning messages
  inline void setMaxIterations(int iter) {
    m_maxIter = iter;
  }
  inline int getMaxIterations(void) const {
    return m_maxIter;
  }

  //Do you wish to fit signal?
  inline void fitSignal(bool value) {
    m_fitSig = value;
  }

  //Do you wish to remove bins with log(1+s/b)>X?
  inline void sigLLR(double value) {
    m_sigLLR = value;
  }

  //verbosity switch
  inline void setVerbose(bool val) {
    m_verbose = val;
  }

  //perform a fit for this model
  void fitProfile();

  //parameter values and errors, indexed over all params
  double getParamVal(uint i);
  double getParamError(uint i);

  //parameter values and errors, indexed over those fit
  inline double getFitParamErr(uint i) const {
    return m_systErrors[i];
  }
  inline double getFitParamVal(uint i) const {
    return m_systValues[i];
  }
  double getFitParamErr(string name) const;
  double getFitParamVal(string name) const;
  std::vector<double> const& getFitParams(); // USED TO BE: const double*

  //MINOS parameter errors...
  bool getMinosError(uint i, double& errP, double& errM, double& parab, double& globCC);

  //Covariance matrix
  bool getErrorMatrix(double* emat) {
    return getEMat(emat);
  }

  //Get name of fitted systematic index i
  string getFitSystName(uint i);

  //Get the index of this syst name
  int getFitSystIndex(string name);

  //How many fitted parameters are there?
  inline int getNfitSyst() const {
    return m_fitSig ? m_systNamesFit.size() : m_systNamesFitBkgd.size();
  }

  //Get the gaussian chi2 constraint term for the most recent fit
  double getSystSum2();

  //MINUIT information
  int getNiterations();
  int getStatus() {
    return m_status;
  }
  double getFuncMin() {
    return m_chi2min;
  }
  double getEstDistMin() {
    return m_edm;
  }

  //Methods used in signal floating for xsec calculations
  void floatSignal(bool floatit) {
    m_sigFloat = floatit;
  }
  double getSignalSF1() {
    return m_sigScale1;
  }
  double getSignalSF2() {
    return m_sigScale2;
  }
  void setSignalSF1(double sf) {
    m_sigScaleNom1 = sf;
  }
  void setSignalSF2(double sf) {
    m_sigScaleNom2 = sf;
  }
  double getSignalSFerr1() {
    return m_sigScaleErr1;
  }
  double getSignalSFerr2() {
    return m_sigScaleErr2;
  }

  //Indicate that a fit test is occurring
  void setFitTest(bool test) {
    m_fitTest = test;
  }

  //Print out some useful info.
  void print();

private:
  void setup();
  void fillArrays();
  void fillSyst();
  void initFitParams();
  double evaluate(const double* parameters);

  bool getEMat(double* emat);
  void getMNerrs();

  bool m_init;
  bool m_fitTest;
  bool m_fitSig;
  bool m_verbose;
  bool m_sigFloat;
  double m_sigScale1;
  double m_sigScaleNom1;
  double m_sigScaleErr1;
  double m_sigScale2;
  double m_sigScaleNom2;
  double m_sigScaleErr2;
  double m_sigLLR; // default = 0.005
  double m_chi2min;
  double m_edm;
  int m_lun;
  int m_status;
  int m_maxIter;

  std::map<std::string, int> m_systNamesFitMap;
  std::map<std::string, int> m_systNamesFitBkgdMap;
  std::vector<std::string> m_systNamesFit;
  std::vector<std::string> m_systNamesFitBkgd;
  std::vector<std::string> m_systNames;
  std::vector<double> m_systValues;
  std::vector<double> m_systErrors;
  std::vector<double> pf2d_paramErrs;
  std::vector<double> m_systErrorsUp;
  std::vector<double> m_systErrorsDn;
  std::vector<double> m_systErrorsParab;
  std::vector<double> m_systErrorsGlobCC;

  SigBkgdDist* pf2d_sbd;
  std::vector<double> pf2d_params;
  std::vector<uint> pf2d_index;
  uint pf2d_nParm;
  uint pf2d_iter;
  bool pf2d_fitSig;
  bool pf2d_sigFloat;
  double pf2d_sigLLR;
  double pf2d_scl1;
  double pf2d_scl2;
  TrialPoint pf2d_lastMinuitTrialPoint;
  std::unique_ptr<collie::ProfileLikelihoodEngine> m_engine;

  struct {
    double sig;
    double bkgd;
    double data;
    double chi2;
    double chi2_bins;
    double chi2_syst;
  } m_condI, m_condF;
};

#endif // ProfileLH_2D_HH_INCLUDED
