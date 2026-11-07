#ifndef CLpoint_hh_included
#define CLpoint_hh_included

class TTree;
class TGraph;

/** Structure containing the many values produced by a confidence level
    calculation.
*/
struct CLpoint {
private:
  int m_var1, m_var2, m_var3;

public:
  /// constructor
  CLpoint(int var1 = 0, int var2 = 0, int var3 = 0);

  /// clear the CLpoint and reset the index variables
  void reset(int var1, int var2 = 0, int var3 = 0);

  /// print out information on this CLpoint
  void print();

  /// create the appropriate branches on the tree for this CLpoint
  void branch(TTree* tree);

  /// map this CLpoint's input variables onto the tree
  void mapToTree(TTree* tree);

  /// get the first independent variable associated with this CLpoint
  inline int getVar1() const {
    return m_var1;
  }
  /// get the second independent variable associated with this CLpoint
  inline int getVar2() const {
    return m_var2;
  }
  /// get the third independent variable associated with this CLpoint
  inline int getVar3() const {
    return m_var3;
  }

  /// calculate lumifactor
  double getLumiFactor();

  /// observed confidence level for signal (from data)
  double cls_obs;
  /// background confidence level
  double clb_obs;
  /// signal + background confidence level
  double clsb_obs;

  /// background confidence level if data was equal to sig+bkgd
  double clb_sb;
  /// background confidence level if data was equal to sig+bkgd +1 sigma
  double clb_sb_p1s;
  /// background confidence level if data was equal to sig+bkgd +2 sigma
  double clb_sb_p2s;
  /// background confidence level if data was equal to sig+bkgd -1 sigma
  double clb_sb_m1s;
  /// background confidence level if data was equal to sig+bkgd -2 sigma
  double clb_sb_m2s;

  /// median expected background confidence level
  double clb_med;
  /// median expected S+B confidence level
  double clsb_med;
  /// median expected confidence level in the presence of signal
  double cls_med;
  /// +1 sigma deviation in the median expected confidence level in the presence of signal
  double cls_med_p1s;
  /// +2 sigma deviation in the median expected confidence level in the presence of signal
  double cls_med_p2s;
  /// -1 sigma deviation in the median expected confidence level in the presence of signal
  double cls_med_m1s;
  /// -2 sigma deviation in the median expected confidence level in the presence of signal
  double cls_med_m2s;
  /// +1 sigma deviation in the median expected confidence level in the presence of signal
  double clsb_med_p1s;
  /// +2 sigma deviation in the median expected confidence level in the presence of signal
  double clsb_med_p2s;
  /// -1 sigma deviation in the median expected confidence level in the presence of signal
  double clsb_med_m1s;
  /// -2 sigma deviation in the median expected confidence level in the presence of signal
  double clsb_med_m2s;

  /// observed log-likelihood ratio (from data)
  double llrobs;
  /// expected log-likelihood ratio in the absence of signal
  double llrb;
  /// +1 sigma deviation in the expected log-likelihood ratio in the absence of signal
  double llrb_p1s;
  /// +2 sigma deviation in the expected log-likelihood ratio in the absence of signal
  double llrb_p2s;
  /// -1 sigma deviation in the expected log-likelihood ratio in the absence of signal
  double llrb_m1s;
  /// -2 sigma deviation in the expected log-likelihood ratio in the absence of signal
  double llrb_m2s;

  /// expected log-likelihood ratio in the presence of signal
  double llrsb;
  /// +1 sigma deviation in the expected log-likelihood ratio in the presence of signal
  double llrsb_p1s;
  /// +2 sigma deviation in the expected log-likelihood ratio in the presence of signal
  double llrsb_p2s;
  /// -1 sigma deviation in the expected log-likelihood ratio in the presence of signal
  double llrsb_m1s;
  /// -2 sigma deviation in the expected log-likelihood ratio in the presence of signal
  double llrsb_m2s;

  /// confidence level sought for in CrossSectionLimit code
  double xsec_cl;
  /// factor which the signal must be scaled by to obtain an observed limit of [xsec_cl]
  double xsec_obsfactor;

  /// factor which the signal must be scaled by to obtain a median expected limit of [xsec_cl]
  double xsec_medfactor_p2s;
  double xsec_medfactor_p1s;
  double xsec_medfactor;
  double xsec_medfactor_m1s;
  double xsec_medfactor_m2s;

  /// scaling factor applied to signal when the application was run (default=1.0)
  double signal_scale;

  /// number of trials included to confidence level calculations
  int nTrials_exp;
  int nTrials_obs;

  /// Results from the xsec fitter!
  double fit_sigScale;
  double fit_sigScale_Err;
  double fit_sigScale_ErrP;
  double fit_sigScale_ErrM;
  double fit_sigScale_ErrStat;
  double fit_sigScale_ErrStatP;
  double fit_sigScale_ErrStatM;
  double fit_sigScale_ErrSyst;
  double fit_sigScale_ErrSystP;
  double fit_sigScale_ErrSystM;
  double fit_chi2_b;
  double fit_chi2_sb;
};

#endif
