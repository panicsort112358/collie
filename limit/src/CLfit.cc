#include <math.h>
#include "CLfit.hh"
#include <algorithm>
#include "ProfileLH.hh"

template <class T> inline const T& MAXof(const T& a, const T& b) {
  return (a > b) ? a : b;
}

CLfit::CLfit() {
  iterations_ = 0;
  sigLLR_ = 0.005;
  nmcdone_ = 0;

  //Do we want to set this randomly or dump to /dev/null.  Here choose the second option.
  //  MTwistEngine twist((int)(a.tv_usec%10000+a.tv_sec%100));
  //  lun_ = 1+(int)(twist.flat()*98);
  lun_ = 99;

  medianExpected_ = true;
  useHistoStats_ = false;
  fitSignal_ = false;

  setNoviceFlag(true);
}

bool CLfit::calculateCLs(const SigBkgdDist& sbd, CLpoint& CLs, int effort) {
  int effortIterations = 10000;
  if (effort == LEVEL_VERYVERYFAST)
    effortIterations = 5000;
  else if (effort == LEVEL_VERYFAST)
    effortIterations = 15000;
  else if (effort == LEVEL_FAST)
    effortIterations = 25000;
  else if (effort == LEVEL_STANDARD)
    effortIterations = 50000;
  else if (effort == LEVEL_FINE)
    effortIterations = 100000;
  else if (effort == LEVEL_VERYFINE)
    effortIterations = 200000;
  else if (effort == LEVEL_VERYVERYFINE)
    effortIterations = 500000;
  return doCLs(sbd, CLs, iterations_ > 0 ? iterations_ : effortIterations);
}

bool CLfit::doCLs(const SigBkgdDist& sbd, CLpoint& CLs, int its) {

  double llr_d = 0, llr_b_expected = 0, llr_sb_expected = 0, llr_sb = 0, llr_b = 0;
  double zback = 0, zsig = 0;
  int numer = 0;
  int denom = 0;
  int denom_bexp = 0;
  int nbins = sbd.nbins();
  int nmc = 0;

  //Obtain the baseline sig,bkg,data distributions
  SigBkgdDist sbBase(sbd), sblike(sbd), blike(sbd), rand(sbd);
  sblike.setBaselineModel();
  blike.setBaselineModel();
  sbBase.setBaselineModel();
  rand.setBaselineModel();
  rand.useHistoStats(useHistoStats_);

  //Set up the fitter
  ProfileLH pfLH;
  pfLH.sigLLR(sigLLR_);
  pfLH.setLUN(lun_);
  pfLH.fitSignal(fitSignal_);
  pfLH.setModel(&sbBase);

  //perform fit to data in s+b hypothesis (test hypo)
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  pfLH.fitProfile();         //fit twice on the first fit to make sure MINUIT is happy
  pfLH.fitProfile();
  llr_d = sbBase.calculateLLR();

  //perform fit to s+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = sblike.bkgd(i) + sblike.signal(i);
  pfLH.fitProfile();
  llr_sb_expected = sbBase.calculateLLR();

  //perform fit to b-only hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = blike.bkgd(i);
  pfLH.fitProfile();
  llr_b_expected = sbBase.calculateLLR();

  CLs.llrobs = llr_d;
  CLs.llrb = llr_b_expected;
  CLs.llrsb = llr_sb_expected;

  for (nmc = 0; nmc < 2 * its && (denom < its / 2.0 || nmc < its); nmc++) {
    rand.fluctuate();
    for (int j = 0; j < nbins; j++) {
      sblike.data(j) = 0;
      blike.data(j) = 0;
      if (rand.bkgd(j) == 0 && rand.signal(j) == 0)
        continue; //skip zero bins...
      if (rand.bkgd(j) > 0) {
        zback = rp_.fire(rand.bkgd(j));
      } else
        zback = 0;
      if ((rand.signal(j) + rand.bkgd(j)) > 0) {
        zsig = rp_.fire(rand.signal(j));
      } else
        zsig = 0;

      blike.data(j) = zback;
      sblike.data(j) = zsig + zback;
    }

    // first fit to the s+b pseudo-data
    sbBase.setBaselineModel();
    for (int b = 0; b < nbins; b++)
      sbBase.data(b) = sblike.data(b);
    pfLH.fitProfile();
    llr_sb = sbBase.calculateLLR();

    //  now fit to the b-only pseudo-data
    sbBase.setBaselineModel();
    for (int b = 0; b < nbins; b++)
      sbBase.data(b) = blike.data(b);
    pfLH.fitProfile();
    llr_b = sbBase.calculateLLR();

    if (llr_d <= llr_b) {
      denom++; // background "less than" data
      if (llr_d <= llr_sb)
        numer++; // s+b "less than" data
    }
    if (llr_b_expected <= llr_b)
      denom_bexp++; // background "less than" med expected b
    if (nmc < MaxKeep_) {
      btrial_[nmc] = llr_b;
      sbtrial_[nmc] = llr_sb;
    }
  }

  nmcdone_ = nmc;
  CLs.clb_obs = (double(denom) / double(nmcdone_));
  CLs.clsb_obs = (double(numer) / double(nmcdone_));
  CLs.clb_med = (double(denom_bexp) / double(nmcdone_));

  double cls_obs = 1.0;
  if (denom > 0)
    cls_obs = 1.0 - double(numer) / double(denom); // fraction of trials where data < s+b
  if (sbd.totSignal() < 0.0001)
    cls_obs = 0.0;
  CLs.cls_obs = cls_obs;

  int lim = (nmcdone_ < MaxKeep_) ? nmcdone_ : MaxKeep_;

  CLs.nTrials_obs = nmcdone_;
  CLs.nTrials_exp = lim;

  std::sort(btrial_, btrial_ + lim);
  std::sort(sbtrial_, sbtrial_ + lim);

  CLs.cls_med = -1.0;
  CLs.cls_med_m1s = -1.0;
  CLs.cls_med_m2s = -1.0;
  CLs.cls_med_p1s = -1.0;
  CLs.cls_med_p2s = -1.0;

  const int medianDenominator = std::max(1, lim / 2);
  int m1s_pos = std::min(lim - 1, int(lim * (1.0 - 0.3173 / 2.0) + 0.5) + 1);
  int p1s_pos = std::min(lim - 1, int(lim * (0.3173 / 2.0) + 0.5) + 1);
  int m2s_pos = std::min(lim - 1, int(lim * (1.0 - 0.0455 / 2.0) + 0.5) + 1);
  int p2s_pos = std::min(lim - 1, int(lim * (0.0455 / 2.0) + 0.5) + 1);

  CLs.llrsb_m2s = sbtrial_[m2s_pos];
  CLs.llrsb_m1s = sbtrial_[m1s_pos];
  CLs.llrsb_p1s = sbtrial_[p1s_pos];
  CLs.llrsb_p2s = sbtrial_[p2s_pos];

  CLs.llrb_m2s = btrial_[m2s_pos];
  CLs.llrb_m1s = btrial_[m1s_pos];
  CLs.llrb_p1s = btrial_[p1s_pos];
  CLs.llrb_p2s = btrial_[p2s_pos];

  double val = btrial_[lim / 2];
  for (int i = 0; i < lim; i++) {
    if (sbtrial_[i] >= val && CLs.cls_med < 0)
      CLs.cls_med = MAXof(0.0, double(i) / medianDenominator - 1.0);
    if (sbtrial_[i] >= btrial_[m1s_pos] && CLs.cls_med_m1s < 0)
      CLs.cls_med_m1s = 1.0 - double(lim - i) / double(lim - m1s_pos);
    if (sbtrial_[i] >= btrial_[m2s_pos] && CLs.cls_med_m2s < 0)
      CLs.cls_med_m2s = 1.0 - double(lim - i) / double(lim - m2s_pos);
    if (sbtrial_[i] >= btrial_[p1s_pos] && CLs.cls_med_p1s < 0)
      CLs.cls_med_p1s = 1.0 - double(lim - i) / double(lim - p1s_pos);
    if (sbtrial_[i] >= btrial_[p2s_pos] && CLs.cls_med_p2s < 0)
      CLs.cls_med_p2s = 1.0 - double(lim - i) / double(lim - p2s_pos);
  }
  if (CLs.cls_med == -1.0)
    CLs.cls_med = 1.0;
  if (CLs.cls_med_m1s == -1.0)
    CLs.cls_med_m1s = 1.0;
  if (CLs.cls_med_m2s == -1.0)
    CLs.cls_med_m2s = 1.0;
  if (CLs.cls_med_p1s == -1.0)
    CLs.cls_med_p1s = 1.0;
  if (CLs.cls_med_p2s == -1.0)
    CLs.cls_med_p2s = 1.0;

  if (CLs.cls_med > 1.0 && CLs.cls_med < 1.001)
    CLs.cls_med = 1.0;
  if (CLs.cls_med < 0.0)
    CLs.cls_med = 0.0;
  CLs.clsb_med = (1.0 - CLs.cls_med) / 2.0;

  if (!medianExpected_) {
    val = llr_b_expected;
    CLs.cls_med = -1.0;
    for (int i = 0; i < lim; i++) {
      if (sbtrial_[i] >= val && CLs.cls_med < 0) {
        CLs.clsb_med = 1.0 - double(i) / double(lim);
        CLs.cls_med = MAXof(0.0, 1 - CLs.clsb_med / CLs.clb_med);
      }
    }
  }

  CLs.clb_sb = -1.0;
  CLs.clb_sb_p1s = -1.0;
  CLs.clb_sb_p2s = -1.0;
  CLs.clb_sb_m1s = -1.0;
  CLs.clb_sb_m2s = -1.0;
  for (int i = 0; i < lim; i++) {
    if (btrial_[i] > CLs.llrsb && CLs.clb_sb < 0)
      CLs.clb_sb = 1.0 - double(i) / double(lim);
    if (btrial_[i] > CLs.llrsb_p1s && CLs.clb_sb_p1s < 0)
      CLs.clb_sb_p1s = 1.0 - double(i) / double(lim);
    if (btrial_[i] > CLs.llrsb_p2s && CLs.clb_sb_p2s < 0)
      CLs.clb_sb_p2s = 1.0 - double(i) / double(lim);
    if (btrial_[i] > CLs.llrsb_m1s && CLs.clb_sb_m1s < 0)
      CLs.clb_sb_m1s = 1.0 - double(i) / double(lim);
    if (btrial_[i] > CLs.llrsb_m2s && CLs.clb_sb_m2s < 0)
      CLs.clb_sb_m2s = 1.0 - double(i) / double(lim);
  }
  if (CLs.clb_sb == -1.0)
    CLs.clb_sb = 0;
  if (CLs.clb_sb_m1s == -1.0)
    CLs.clb_sb_m1s = 0;
  if (CLs.clb_sb_m2s == -1.0)
    CLs.clb_sb_m2s = 0;
  if (CLs.clb_sb_p1s == -1.0)
    CLs.clb_sb_p1s = 0;
  if (CLs.clb_sb_p2s == -1.0)
    CLs.clb_sb_p2s = 0;

  return true;
}

void CLfit::configure(const char* options) {
  parseIterations(options, iterations_);
}

TH1D* CLfit::getLLRdist_sb(const char* title, int bins, double min, double max) {

  TH1D* out = new TH1D(title, title, bins, min, max);
  int lim = (nmcdone_ < MaxKeep_) ? nmcdone_ : MaxKeep_;

  for (int i = 0; i < lim; i++)
    out->Fill(sbtrial_[i]);

  return out;
}

TH1D* CLfit::getLLRdist_b(const char* title, int bins, double min, double max) {

  TH1D* out = new TH1D(title, title, bins, min, max);
  int lim = (nmcdone_ < MaxKeep_) ? nmcdone_ : MaxKeep_;

  for (int i = 0; i < lim; i++)
    out->Fill(btrial_[i]);

  return out;
}

void CLfit::noviceWarning() {

  printf("\n*****************************************************************\n");
  printf(
      "CLfit Warning:  You are using a fitting CL calculator.  Have you\nverified your systematics and fit model using the FitTest class?\n");
  printf("*****************************************************************\n\n");

  if (useHistoStats_) {
    printf("\n*****************************************************************\n");
    printf(
        "CLfit Warning:  You're operating in novice mode and you're using\nstatistical uncertainties from your input histograms.  Are you\nsure you've set the values properly in the histograms?\n");
    printf("*****************************************************************\n\n");
  } else {
    printf("\n*****************************************************************\n");
    printf(
        "CLfit Warning:  The histogram statistical uncertainty flag is turned\noff by default.  This calculation does not include statistical uncertainties.\n");
    printf("*****************************************************************\n\n");
  }

  const char* hypo = "Signal+Background";
  if (!fitSignal_)
    hypo = "Background-Only";
  printf("\n*****************************************************************\n");
  printf("CLfit Warning:  You have chosen to fit all pseudo-experiments\nin the %s hypothesis\n",
         hypo);
  printf("*****************************************************************\n\n");

  if (!fitSignal_) {
    printf("\n*****************************************************************\n");
    printf("CLfit Warning:  You have chosen to exclude bins with log(1+s/b)>%f\n", sigLLR_);
    printf("*****************************************************************\n\n");
  }
}
