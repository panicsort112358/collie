#include <math.h>
#include "CLfit2.hh"
#include <algorithm>
#include "ProfileLH.hh"

template <class T> inline const T& MAXof(const T& a, const T& b) {
  return (a > b) ? a : b;
}

double appAWW(double inv) {
  return 0.5 * (1 + TMath::Erf(inv / sqrt(2)));
}

CLfit2::CLfit2() {
  iterations_ = 0;
  nmcdone_ = 0;

  //Do we want to set this randomly or dump to /dev/null.  Here choose the second option.
  //  MTwistEngine twist((int)(a.tv_usec%10000+a.tv_sec%100));
  //  lun_ = 1+(int)(twist.flat()*98);
  lun_ = 99;

  useHistoStats_ = false;
  medianExpected_ = true;
  doAWW_ = false;

  setNoviceFlag(true);
}

bool CLfit2::calculateCLs(const SigBkgdDist& sbd, CLpoint& CLs, int effort) {
  int effortIterations = 10000;
  if (effort == LEVEL_NONE)
    effortIterations = 1;
  else if (effort == LEVEL_TEST)
    effortIterations = 100;
  else if (effort == LEVEL_COMBO3)
    effortIterations = 250;
  else if (effort == LEVEL_COMBO2)
    effortIterations = 1500;
  else if (effort == LEVEL_COMBO)
    effortIterations = 10000;
  else if (effort == LEVEL_VERYVERYFAST)
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

bool CLfit2::doCLsAWW(const SigBkgdDist& sbd, CLpoint& CLs) {
  printf("\nPerforming CLfit2 AWW Limit Approximation\n");

  double llr_b_expected = 0, llr_sb_expected = 0, llr_d = 0;
  int nbins = sbd.nbins();

  //Obtain the baseline sig,bkg,data distributions
  SigBkgdDist sbBase(sbd);
  sbBase.setBaselineModel();

  //Set up the fitter
  ProfileLH pfLH;
  pfLH.sigLLR(1e10);
  pfLH.setLUN(lun_);
  pfLH.fitSignal(true);
  pfLH.setModel(&sbBase);

  //perform fit to data in s+b hypothesis (test hypo)
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  pfLH.fitProfile();
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3) {
    printf("CLfit2::doCLs, Fit didn't converge(d1): %d\n", pfLH.getStatus());
    sbBase.setBaselineModel(); //zeros nuisance param central values...
    pfLH.fitProfile();
  }
  llr_d = pfLH.getFuncMin();

  //perform s+b fit to s+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = sbd.bkgd(i) + sbd.signal(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfit2::doCLs, Fit didn't converge(sb1): %d\n", pfLH.getStatus());
  llr_sb_expected = pfLH.getFuncMin();

  //perform s+b fit to s+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = sbd.bkgd(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfit2::doCLs, Fit didn't converge(b1): %d\n", pfLH.getStatus());
  llr_b_expected = pfLH.getFuncMin();

  //perform fit to data in b-only hypothesis
  pfLH.fitSignal(false);
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  pfLH.fitProfile();
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfit2::doCLs, Fit didn't converge(d2): %d\n", pfLH.getStatus());
  llr_d -= pfLH.getFuncMin();

  //perform b-only fit to s+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = sbd.bkgd(i) + sbd.signal(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfit2::doCLs, Fit didn't converge(sb2): %d\n", pfLH.getStatus());
  llr_sb_expected -= pfLH.getFuncMin();

  //perform b-only fit to b-only hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = sbd.bkgd(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfit2::doCLs, Fit didn't converge(b2): %d\n", pfLH.getStatus());
  llr_b_expected -= pfLH.getFuncMin();

  if (std::isnan(llr_b_expected))
    llr_b_expected = 1e-6;
  if (std::isnan(llr_sb_expected))
    llr_sb_expected = -1e-6;

  llr_b_expected += 1e-6;
  llr_sb_expected -= 1e-6;

  CLs.llrobs = llr_d;
  CLs.llrb = llr_b_expected;
  CLs.llrsb = llr_sb_expected;

  double bsig = 2 * sqrt(llr_b_expected) + 1e-6;
  double ssig = 2 * sqrt(-1 * llr_sb_expected) + 1e-6;
  if (std::isnan(bsig))
    bsig = ssig;
  if (std::isnan(ssig))
    ssig = bsig;
  if (std::isnan(bsig)) {
    bsig = 1e-6;
    ssig = 1e-6;
  }

  CLs.llrb_p1s = llr_b_expected - bsig;
  CLs.llrb_p2s = llr_b_expected - 2 * bsig;
  CLs.llrb_m1s = llr_b_expected + bsig;
  CLs.llrb_m2s = llr_b_expected + 2 * bsig;

  CLs.llrsb_p1s = llr_sb_expected - ssig;
  CLs.llrsb_p2s = llr_sb_expected - 2 * ssig;
  CLs.llrsb_m1s = llr_sb_expected + ssig;
  CLs.llrsb_m2s = llr_sb_expected + 2 * ssig;

  //Get Obs p-values and CLs
  CLs.clsb_obs = 1 - appAWW((llr_d - llr_sb_expected) / ssig) + 1e-6;
  CLs.clb_obs = 1 - appAWW((llr_d - llr_b_expected) / bsig) + 1e-6;
  CLs.cls_obs = 1 - CLs.clsb_obs / CLs.clb_obs;
  if (CLs.cls_obs < 0)
    CLs.cls_obs = 0;
  if (CLs.cls_obs > 1)
    CLs.cls_obs = 1;
  //  printf("llrd: %f, ssig: %f, bsig: %f, llrsbexp: %f, llrbexp: %f, clsb: %f, clb: %f, cls: %f\n",llr_d,ssig,bsig,llr_sb_expected,llr_b_expected,CLs.clsb_obs,CLs.clb_obs,CLs.cls_obs);

  //Get Exp p-values and CLs
  CLs.clsb_med = 1 - appAWW((llr_b_expected - llr_sb_expected) / ssig) + 1e-6;
  CLs.clb_med = 1 - appAWW((llr_b_expected - llr_b_expected) / bsig) + 1e-6;
  CLs.cls_med = 1 - CLs.clsb_med / CLs.clb_med;

  // +2 sigma
  CLs.clsb_med_p2s = 1 - appAWW((CLs.llrb_p2s - llr_sb_expected) / ssig) + 1e-6;
  double pB = 1 - appAWW((CLs.llrb_p2s - llr_b_expected) / bsig) + 1e-6;
  CLs.cls_med_p2s = 1 - CLs.clsb_med_p2s / pB;

  // +1 sigma
  CLs.clsb_med_p1s = 1 - appAWW((CLs.llrb_p1s - llr_sb_expected) / ssig) + 1e-6;
  pB = 1 - appAWW((CLs.llrb_p1s - llr_b_expected) / bsig) + 1e-6;
  CLs.cls_med_p1s = 1 - CLs.clsb_med_p1s / pB;

  // -1 sigma
  CLs.clsb_med_m1s = 1 - appAWW((CLs.llrb_m1s - llr_sb_expected) / ssig) + 1e-6;
  pB = 1 - appAWW((CLs.llrb_m1s - llr_b_expected) / bsig) + 1e-6;
  CLs.cls_med_m1s = 1 - CLs.clsb_med_m1s / pB;

  // -2 sigma
  CLs.clsb_med_m2s = 1 - appAWW((CLs.llrb_m2s - llr_sb_expected) / ssig) + 1e-6;
  pB = 1 - appAWW((CLs.llrb_m2s - llr_b_expected) / bsig) + 1e-6;
  CLs.cls_med_m2s = 1 - CLs.clsb_med_m2s / pB;

  //  printf("CLSB: %f, %f, %f, %f\n",CLs.clsb_med_m2s,CLs.clsb_med_m1s,CLs.clsb_med_p1s,CLs.clsb_med_p2s);

  CLs.clb_sb = 1 - appAWW((CLs.llrsb - llr_b_expected) / bsig) + 1e-6;
  CLs.clb_sb_p1s = 1 - appAWW((CLs.llrsb_p1s - llr_b_expected) / bsig) + 1e-6;
  CLs.clb_sb_p2s = 1 - appAWW((CLs.llrsb_p2s - llr_b_expected) / bsig) + 1e-6;
  CLs.clb_sb_m1s = 1 - appAWW((CLs.llrsb_m1s - llr_b_expected) / bsig) + 1e-6;
  CLs.clb_sb_m2s = 1 - appAWW((CLs.llrsb_m2s - llr_b_expected) / bsig) + 1e-6;

  CLs.nTrials_obs = 1;
  CLs.nTrials_exp = 1;

  return true;
}

bool CLfit2::doCLs(const SigBkgdDist& sbd, CLpoint& CLs, int its) {

  if (doAWW_) {
    doCLsAWW(sbd, CLs);
    return true;
  }

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
  pfLH.sigLLR(1e10);
  pfLH.setLUN(lun_);
  pfLH.fitSignal(true);
  pfLH.setModel(&sbBase);

  //perform fit to data in s+b hypothesis (test hypo)
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  pfLH.fitProfile();
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3) {
    printf("CLfit2::doCLs, Fit didn't converge(d1): %d\n", pfLH.getStatus());
    sbBase.setBaselineModel(); //zeros nuisance param central values...
    pfLH.fitProfile();
  }
  llr_d = pfLH.getFuncMin();

  //perform s+b fit to s+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = sblike.bkgd(i) + sblike.signal(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfit2::doCLs, Fit didn't converge(sb1): %d\n", pfLH.getStatus());
  llr_sb_expected = pfLH.getFuncMin();

  //perform s+b fit to s+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = sblike.bkgd(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfit2::doCLs, Fit didn't converge(b1): %d\n", pfLH.getStatus());
  llr_b_expected = pfLH.getFuncMin();

  //perform fit to data in b-only hypothesis
  pfLH.fitSignal(false);
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfit2::doCLs, Fit didn't converge(d2): %d\n", pfLH.getStatus());
  llr_d -= pfLH.getFuncMin();

  //perform b-only fit to s+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = sblike.bkgd(i) + sblike.signal(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfit2::doCLs, Fit didn't converge(sb2): %d\n", pfLH.getStatus());
  llr_sb_expected -= pfLH.getFuncMin();

  //perform b-only fit to b-only hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = blike.bkgd(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfit2::doCLs, Fit didn't converge(b2): %d\n", pfLH.getStatus());
  llr_b_expected -= pfLH.getFuncMin();

  CLs.llrobs = llr_d;
  CLs.llrb = llr_b_expected;
  CLs.llrsb = llr_sb_expected;
  printf("CLfit2:  LLRobs: %f, LLRb: %f, LLRsb: %f\n", CLs.llrobs, CLs.llrb, CLs.llrsb);
  for (int ift = 0; ift < 2 * its && (denom < its / 2.0 || ift < its); ift++) {
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

    // s+b fit to the s+b pseudo-data
    pfLH.fitSignal(true);
    sbBase.setBaselineModel();
    for (int b = 0; b < nbins; b++)
      sbBase.data(b) = sblike.data(b);
    pfLH.fitProfile();
    if (pfLH.getStatus() == 0) {
      printf("CLfit2::doCLs, Fit didn't converge(l): %d\n", pfLH.getStatus());
      continue;
    }
    llr_sb = pfLH.getFuncMin();

    //  s+b fit to the b-only pseudo-data
    sbBase.setBaselineModel();
    for (int b = 0; b < nbins; b++)
      sbBase.data(b) = blike.data(b);
    pfLH.fitProfile();
    if (pfLH.getStatus() == 0) {
      printf("CLfit2::doCLs, Fit didn't converge(l): %d\n", pfLH.getStatus());
      continue;
    }
    llr_b = pfLH.getFuncMin();

    // b-only fit to the s+b pseudo-data
    pfLH.fitSignal(false);
    sbBase.setBaselineModel();
    for (int b = 0; b < nbins; b++)
      sbBase.data(b) = sblike.data(b);
    pfLH.fitProfile();
    if (pfLH.getStatus() == 0) {
      printf("CLfit2::doCLs, Fit didn't converge(l): %d\n", pfLH.getStatus());
      continue;
    }
    llr_sb -= pfLH.getFuncMin();

    // b-only fit to the b-only pseudo-data
    sbBase.setBaselineModel();
    for (int b = 0; b < nbins; b++)
      sbBase.data(b) = blike.data(b);
    pfLH.fitProfile();
    if (pfLH.getStatus() == 0) {
      printf("CLfit2::doCLs, Fit didn't converge(l): %d\n", pfLH.getStatus());
      continue;
    }
    llr_b -= pfLH.getFuncMin();

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
    nmc++;
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

void CLfit2::configure(const char* options) {
  parseIterations(options, iterations_);
}

TH1D* CLfit2::getLLRdist_sb(const char* title, int bins, double min, double max) {

  TH1D* out = new TH1D(title, title, bins, min, max);
  int lim = (nmcdone_ < MaxKeep_) ? nmcdone_ : MaxKeep_;

  for (int i = 0; i < lim; i++)
    out->Fill(sbtrial_[i]);

  return out;
}

TH1D* CLfit2::getLLRdist_b(const char* title, int bins, double min, double max) {

  TH1D* out = new TH1D(title, title, bins, min, max);
  int lim = (nmcdone_ < MaxKeep_) ? nmcdone_ : MaxKeep_;

  for (int i = 0; i < lim; i++)
    out->Fill(btrial_[i]);

  return out;
}

void CLfit2::noviceWarning() {

  printf("\n*****************************************************************\n");
  printf(
      "CLfit2 Warning:  You are using a fitting CL calculator.  Have you\nverified your systematics and fit model using the FitTest class?\n");
  printf("*****************************************************************\n\n");

  if (useHistoStats_) {
    printf("\n*****************************************************************\n");
    printf(
        "CLfit2 Warning:  You're operating in novice mode and you're using\nstatistical uncertainties from your input histograms.  Are you\nsure you've set the values properly in the histograms?\n");
    printf("*****************************************************************\n\n");
  } else {
    printf("\n*****************************************************************\n");
    printf(
        "CLfit2 Warning:  The histogram statistical uncertainty flag is turned\noff by default.  This calculation does not include statistical uncertainties.\n");
    printf("*****************************************************************\n\n");
  }
}
