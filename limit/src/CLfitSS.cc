#include <math.h>
#include "CLfitSS.hh"
#include <algorithm>
#include "ProfileLH.hh"

template <class T> inline const T& MAXof(const T& a, const T& b) {
  return (a > b) ? a : b;
}

double appAWW2(double inv) {
  return 0.5 * (1 + TMath::Erf(inv / sqrt(2)));
}

CLfitSS::CLfitSS() {
  iterations_ = 0;
  nmcdone_ = 0;

  //Do we want to set this randomly or dump to /dev/null.  Here choose the second option.
  //  MTwistEngine twist((int)(a.tv_usec%10000+a.tv_sec%100));
  //  lun_ = 1+(int)(twist.flat()*98);
  lun_ = 99;

  useHistoStats_ = false;
  medianExpected_ = true;
  doAWW_ = false;

  sigFraction_ = 1.0;
  bsmScale_ = 1.0;

  setNoviceFlag(true);
}

void CLfitSS::setSignalFraction(double frac) {
  if (frac < 0.0)
    sigFraction_ = 0.0;
  else if (frac > 1.0)
    sigFraction_ = 1.0;
  else
    sigFraction_ = frac;
}

void CLfitSS::setBSMScale(double scl) {
  if (scl < 0.0)
    bsmScale_ = 0.0;
  else if (scl > 1.0)
    bsmScale_ = 1.0;
  else
    bsmScale_ = scl;
}

bool CLfitSS::calculateCLs(const SigBkgdDist& sbd, CLpoint& CLs, int effort) {
  int effortIterations = 10000;
  if (effort == LEVEL_NONE)
    effortIterations = 1;
  else if (effort == LEVEL_TEST)
    effortIterations = 100;
  else if (effort == LEVEL_COMBO3)
    effortIterations = 350;
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

bool CLfitSS::doCLsAWW(const SigBkgdDist& sbd, CLpoint& CLs) {
  printf("\nPerforming CLfitSS AWW Limit Approximation\n");

  double nomScale = sbd.getSignalScale();
  printf("Nominal Signal Scale = %f\n", nomScale);
  double llr_b_expected = 0, llr_sb_expected = 0, llr_d = 0;
  int nbins = sbd.nbins();

  //Obtain the baseline sig,bkg,data distributions
  SigBkgdDist sbBase(sbd), s0like(sbd), s1like(sbd);

  int smsig, bsmsig;
  smsig = (sbd.getSigName(0)).find("vh") != std::string::npos ? 0 : 1;  //vvbb
  bsmsig = (sbd.getSigName(1)).find("vh") == std::string::npos ? 1 : 0; //vvbb

  // smsig  = (sbd.getSigName(0)).find("SM") != std::string::npos ? 0 : 1; //lvbb
  // bsmsig = (sbd.getSigName(1)).find("SM") == std::string::npos ? 1 : 0; //lvbb

  printf("SM  signal is index %d with name %s\n", smsig, (sbd.getSigName(smsig)).c_str());
  printf("BSM signal is index %d with name %s\n", bsmsig, (sbd.getSigName(bsmsig)).c_str());

  sbBase.setBaselineModel();
  s0like.setBaselineModel();
  s1like.setBaselineModel();

  s0like.setSignalScale("graviton", 1.0 * nomScale * bsmScale_);
  s0like.setSignalScale("pseudoscalar", 1.0 * nomScale * bsmScale_);
  s0like.setSignalScale("vh", 0.0);

  s0like.setSignalScale("GC", 1.0 * nomScale * bsmScale_);
  s0like.setSignalScale("PS", 1.0 * nomScale * bsmScale_);
  s0like.setSignalScale("SM", 0.0);

  s0like.setSignalScale("grav", 1.0 * nomScale * bsmScale_);
  s0like.setSignalScale("ps", 1.0 * nomScale * bsmScale_);
  s0like.setSignalScale("signal", 0.0);

  s0like.setBaselineModel();

  //kill signal 0...
  s1like.setSignalScale("graviton", nomScale * bsmScale_ * (1 - sigFraction_));
  s1like.setSignalScale("pseudoscalar", nomScale * bsmScale_ * (1 - sigFraction_));
  s1like.setSignalScale("vh", nomScale * sigFraction_);

  s1like.setSignalScale("GC", nomScale * bsmScale_ * (1 - sigFraction_));
  s1like.setSignalScale("PS", nomScale * bsmScale_ * (1 - sigFraction_));
  s1like.setSignalScale("SM", nomScale * sigFraction_);

  s1like.setSignalScale("grav", nomScale * bsmScale_ * (1 - sigFraction_));
  s1like.setSignalScale("ps", nomScale * bsmScale_ * (1 - sigFraction_));
  s1like.setSignalScale("signal", nomScale * sigFraction_);

  s1like.setBaselineModel();

  //Set up the fitter
  ProfileLH pfLH;
  pfLH.sigLLR(1e10);
  pfLH.setLUN(lun_);
  pfLH.fitSignal(true);
  pfLH.setModel(&sbBase);

  //perform fit to data in s0+b hypothesis (test hypo)
  sbBase.setSignalScale("graviton", nomScale * bsmScale_);
  sbBase.setSignalScale("pseudoscalar", nomScale * bsmScale_);
  sbBase.setSignalScale("vh", 0.0);

  sbBase.setSignalScale("GC", nomScale * bsmScale_);
  sbBase.setSignalScale("PS", nomScale * bsmScale_);
  sbBase.setSignalScale("SM", 0.0);

  sbBase.setSignalScale("grav", nomScale * bsmScale_);
  sbBase.setSignalScale("ps", nomScale * bsmScale_);
  sbBase.setSignalScale("signal", 0.0);

  sbBase.setBaselineModel(); //zeros nuisance param central values...
  pfLH.fitProfile();
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3) {
    printf("CLfitSS::doCLs, Fit didn't converge(d1): %d\n", pfLH.getStatus());
    sbBase.setBaselineModel(); //zeros nuisance param central values...
    pfLH.fitProfile();
    if (pfLH.getStatus() != 3)
      printf("CLfitSS::doCLs, Fit didn't converge(d1 2): %d\n", pfLH.getStatus());
  }
  llr_d = pfLH.getFuncMin();

  //perform s0+b fit to s0+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = s0like.bkgd(i) + s0like.signal(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfitSS::doCLs, Fit didn't converge(sb1): %d\n", pfLH.getStatus());
  llr_sb_expected = pfLH.getFuncMin();

  //perform s0+b fit to s1+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = s1like.bkgd(i) + s1like.signal(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfitSS::doCLs, Fit didn't converge(b1): %d\n", pfLH.getStatus());
  llr_b_expected = pfLH.getFuncMin();

  //perform fit to data in s1+b hypothesis
  sbBase.setSignalScale("graviton", nomScale * bsmScale_ * (1 - sigFraction_));
  sbBase.setSignalScale("pseudoscalar", nomScale * bsmScale_ * (1 - sigFraction_));
  sbBase.setSignalScale("vh", nomScale * sigFraction_);

  sbBase.setSignalScale("GC", nomScale * bsmScale_ * (1 - sigFraction_));
  sbBase.setSignalScale("PS", nomScale * bsmScale_ * (1 - sigFraction_));
  sbBase.setSignalScale("SM", nomScale * sigFraction_);

  sbBase.setSignalScale("grav", nomScale * bsmScale_ * (1 - sigFraction_));
  sbBase.setSignalScale("ps", nomScale * bsmScale_ * (1 - sigFraction_));
  sbBase.setSignalScale("signal", nomScale * sigFraction_);

  sbBase.setBaselineModel(); //zeros nuisance param central values...
  pfLH.fitProfile();
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3) {
    printf("CLfitSS::doCLs, Fit didn't converge(d2): %d\n", pfLH.getStatus());
    sbBase.setBaselineModel(); //zeros nuisance param central values...
    pfLH.fitProfile();
    if (pfLH.getStatus() != 3)
      printf("CLfitSS::doCLs, Fit didn't converge(d2 2): %d\n", pfLH.getStatus());
  }
  llr_d -= pfLH.getFuncMin();

  //perform s1+b fit to s0+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = s0like.bkgd(i) + s0like.signal(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfitSS::doCLs, Fit didn't converge(sb2): %d\n", pfLH.getStatus());
  llr_sb_expected -= pfLH.getFuncMin();

  //perform s1+b fit to s1+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = s1like.bkgd(i) + s1like.signal(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfitSS::doCLs, Fit didn't converge(b2): %d\n", pfLH.getStatus());
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
  CLs.clsb_obs = 1 - appAWW2((llr_d - llr_sb_expected) / ssig) + 1e-6;
  CLs.clb_obs = 1 - appAWW2((llr_d - llr_b_expected) / bsig) + 1e-6;
  CLs.cls_obs = 1 - CLs.clsb_obs / CLs.clb_obs;
  if (CLs.cls_obs < 0)
    CLs.cls_obs = 0;
  if (CLs.cls_obs > 1)
    CLs.cls_obs = 1;
  //  printf("llrd: %f, ssig: %f, bsig: %f, llrsbexp: %f, llrbexp: %f, clsb: %f, clb: %f, cls: %f\n",llr_d,ssig,bsig,llr_sb_expected,llr_b_expected,CLs.clsb_obs,CLs.clb_obs,CLs.cls_obs);

  //Get Exp p-values and CLs
  CLs.clsb_med = 1 - appAWW2((llr_b_expected - llr_sb_expected) / ssig) + 1e-6;
  CLs.clb_med = 1 - appAWW2((llr_b_expected - llr_b_expected) / bsig) + 1e-6;
  CLs.cls_med = 1 - CLs.clsb_med / CLs.clb_med;

  // +2 sigma
  CLs.clsb_med_p2s = 1 - appAWW2((CLs.llrb_p2s - llr_sb_expected) / ssig) + 1e-6;
  double pB = 1 - appAWW2((CLs.llrb_p2s - llr_b_expected) / bsig) + 1e-6;
  CLs.cls_med_p2s = 1 - CLs.clsb_med_p2s / pB;

  // +1 sigma
  CLs.clsb_med_p1s = 1 - appAWW2((CLs.llrb_p1s - llr_sb_expected) / ssig) + 1e-6;
  pB = 1 - appAWW2((CLs.llrb_p1s - llr_b_expected) / bsig) + 1e-6;
  CLs.cls_med_p1s = 1 - CLs.clsb_med_p1s / pB;

  // -1 sigma
  CLs.clsb_med_m1s = 1 - appAWW2((CLs.llrb_m1s - llr_sb_expected) / ssig) + 1e-6;
  pB = 1 - appAWW2((CLs.llrb_m1s - llr_b_expected) / bsig) + 1e-6;
  CLs.cls_med_m1s = 1 - CLs.clsb_med_m1s / pB;

  // -2 sigma
  CLs.clsb_med_m2s = 1 - appAWW2((CLs.llrb_m2s - llr_sb_expected) / ssig) + 1e-6;
  pB = 1 - appAWW2((CLs.llrb_m2s - llr_b_expected) / bsig) + 1e-6;
  CLs.cls_med_m2s = 1 - CLs.clsb_med_m2s / pB;

  CLs.clb_sb = 1 - appAWW2((CLs.llrsb - llr_b_expected) / bsig) + 1e-6;
  CLs.clb_sb_p1s = 1 - appAWW2((CLs.llrsb_p1s - llr_b_expected) / bsig) + 1e-6;
  CLs.clb_sb_p2s = 1 - appAWW2((CLs.llrsb_p2s - llr_b_expected) / bsig) + 1e-6;
  CLs.clb_sb_m1s = 1 - appAWW2((CLs.llrsb_m1s - llr_b_expected) / bsig) + 1e-6;
  CLs.clb_sb_m2s = 1 - appAWW2((CLs.llrsb_m2s - llr_b_expected) / bsig) + 1e-6;

  CLs.nTrials_obs = 1;
  CLs.nTrials_exp = 1;

  return true;
}

bool CLfitSS::doCLs(const SigBkgdDist& sbd, CLpoint& CLs, int its) {

  if (doAWW_) {
    doCLsAWW(sbd, CLs);
    return true;
  }

  double nomScale = sbd.getSignalScale();

  double llr_d = 0, llr_b_expected = 0, llr_sb_expected = 0, llr_sb = 0, llr_b = 0;
  double zsb0 = 0, zsb1 = 0;
  int numer = 0;
  int denom = 0;
  int denom_bexp = 0;
  int nbins = sbd.nbins();
  int nmc = 0;

  //Obtain the baseline sig,bkg,data distributions
  SigBkgdDist sbBase(sbd), s0like(sbd), s1like(sbd), rand0(sbd), rand1(sbd);

  int smsig, bsmsig;
  smsig = (sbd.getSigName(0)).find("vh") != std::string::npos ? 0 : 1;
  bsmsig = (sbd.getSigName(1)).find("vh") == std::string::npos ? 1 : 0;

  printf("SM signal is index %d with name %s\n", smsig, (sbd.getSigName(smsig)).c_str());
  printf("BSM signal is index %d with name %s\n", bsmsig, (sbd.getSigName(bsmsig)).c_str());

  s0like.setBaselineModel();
  s1like.setBaselineModel();
  sbBase.setBaselineModel();
  rand0.setBaselineModel();

  printf("A Sig: %f/%f/%f, Bkgd: %f/%f/%f\n", sbBase.totSignal(), s0like.totSignal(),
         s1like.totSignal(), sbBase.totBkgd(), s0like.totBkgd(), s1like.totBkgd());

  rand1.setBaselineModel();
  rand0.useHistoStats(useHistoStats_);
  rand1.useHistoStats(useHistoStats_);

  ///kill signal 1
  s0like.setSignalScale("graviton", nomScale * bsmScale_);
  s0like.setSignalScale("pseudoscalar", nomScale * bsmScale_);
  s0like.setSignalScale("vh", 0.0);

  s0like.setSignalScale("GC", nomScale * bsmScale_);
  s0like.setSignalScale("PS", nomScale * bsmScale_);
  s0like.setSignalScale("SM", 0.0);

  s0like.setSignalScale("grav", nomScale * bsmScale_);
  s0like.setSignalScale("ps", nomScale * bsmScale_);
  s0like.setSignalScale("signal", 0.0);
  s0like.setBaselineModel();

  rand0.setSignalScale("graviton", nomScale * bsmScale_);
  rand0.setSignalScale("pseudoscalar", nomScale * bsmScale_);
  rand0.setSignalScale("vh", 0.0);

  rand0.setSignalScale("GC", nomScale * bsmScale_);
  rand0.setSignalScale("PS", nomScale * bsmScale_);
  rand0.setSignalScale("SM", 0.0);

  rand0.setSignalScale("grav", nomScale * bsmScale_);
  rand0.setSignalScale("ps", nomScale * bsmScale_);
  rand0.setSignalScale("signal", 0.0);
  rand0.setBaselineModel();

  //kill signal 0...
  s1like.setSignalScale("graviton", nomScale * bsmScale_ * (1 - sigFraction_));
  s1like.setSignalScale("pseudoscalar", nomScale * bsmScale_ * (1 - sigFraction_));
  s1like.setSignalScale("vh", nomScale * sigFraction_);

  s1like.setSignalScale("GC", nomScale * bsmScale_ * (1 - sigFraction_));
  s1like.setSignalScale("PS", nomScale * bsmScale_ * (1 - sigFraction_));
  s1like.setSignalScale("SM", nomScale * sigFraction_);

  s1like.setSignalScale("grav", nomScale * bsmScale_ * (1 - sigFraction_));
  s1like.setSignalScale("ps", nomScale * bsmScale_ * (1 - sigFraction_));
  s1like.setSignalScale("signal", nomScale * sigFraction_);
  s1like.setBaselineModel();

  rand1.setSignalScale("graviton", nomScale * bsmScale_ * (1 - sigFraction_));
  rand1.setSignalScale("pseudoscalar", nomScale * bsmScale_ * (1 - sigFraction_));
  rand1.setSignalScale("vh", nomScale * sigFraction_);

  rand1.setSignalScale("GC", nomScale * bsmScale_ * (1 - sigFraction_));
  rand1.setSignalScale("PS", nomScale * bsmScale_ * (1 - sigFraction_));
  rand1.setSignalScale("SM", nomScale * sigFraction_);

  rand1.setSignalScale("grav", nomScale * bsmScale_ * (1 - sigFraction_));
  rand1.setSignalScale("ps", nomScale * bsmScale_ * (1 - sigFraction_));
  rand1.setSignalScale("signal", nomScale * sigFraction_);
  //  rand1.setSignalScale(bsmsig,nomScale*(1-sigFraction_));
  //  rand1.setSignalScale(smsig,nomScale*sigFraction_);
  rand1.setBaselineModel();

  printf("B Sig: %f/%f/%f, Bkgd: %f/%f/%f\n", sbBase.totSignal(), s0like.totSignal(),
         s1like.totSignal(), sbBase.totBkgd(), s0like.totBkgd(), s1like.totBkgd());

  //Set up the fitter
  ProfileLH pfLH;
  pfLH.sigLLR(1e10);
  pfLH.setLUN(lun_);
  pfLH.fitSignal(true);
  pfLH.setModel(&sbBase);

  //perform fit to data in s1+b hypothesis (test hypo)
  sbBase.setSignalScale("graviton", nomScale * bsmScale_);
  sbBase.setSignalScale("pseudoscalar", nomScale * bsmScale_);
  sbBase.setSignalScale("vh", 0.0);

  sbBase.setSignalScale("GC", nomScale * bsmScale_);
  sbBase.setSignalScale("PS", nomScale * bsmScale_);
  sbBase.setSignalScale("SM", 0.0);

  sbBase.setSignalScale("grav", nomScale * bsmScale_);
  sbBase.setSignalScale("ps", nomScale * bsmScale_);
  sbBase.setSignalScale("signal", 0.0);

  sbBase.setBaselineModel(); //zeros nuisance param central values...
  pfLH.fitProfile();
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3) {
    printf("CLfitSS::doCLs, Fit didn't converge(d1): %d\n", pfLH.getStatus());
    sbBase.setBaselineModel(); //zeros nuisance param central values...
    pfLH.fitProfile();
    if (pfLH.getStatus() != 3)
      printf("CLfitSS::doCLs, Fit didn't converge(d1 2): %d\n", pfLH.getStatus());
  }
  llr_d = pfLH.getFuncMin();

  //perform s0+b fit to s0+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = s0like.bkgd(i) + s0like.signal(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfitSS::doCLs, Fit didn't converge(sb1): %d\n", pfLH.getStatus());
  llr_sb_expected = pfLH.getFuncMin();

  //perform s0+b fit to s1+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = s1like.bkgd(i) + s1like.signal(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfitSS::doCLs, Fit didn't converge(b1): %d\n", pfLH.getStatus());
  llr_b_expected = pfLH.getFuncMin();

  //perform fit to data in null hypothesis
  sbBase.setSignalScale("graviton", nomScale * bsmScale_ * (1 - sigFraction_));
  sbBase.setSignalScale("pseudoscalar", nomScale * bsmScale_ * (1 - sigFraction_));
  sbBase.setSignalScale("vh", nomScale * sigFraction_);

  sbBase.setSignalScale("GC", nomScale * bsmScale_ * (1 - sigFraction_));
  sbBase.setSignalScale("PS", nomScale * bsmScale_ * (1 - sigFraction_));
  sbBase.setSignalScale("SM", nomScale * sigFraction_);

  sbBase.setSignalScale("grav", nomScale * bsmScale_ * (1 - sigFraction_));
  sbBase.setSignalScale("ps", nomScale * bsmScale_ * (1 - sigFraction_));
  sbBase.setSignalScale("signal", nomScale * sigFraction_);

  //  sbBase.setSignalScale(smsig,0.); //test

  sbBase.setBaselineModel(); //zeros nuisance param central values...

  pfLH.fitProfile();
  if (pfLH.getStatus() != 3) {
    printf("CLfitSS::doCLs, Fit didn't converge(d2): %d\n", pfLH.getStatus());
    sbBase.setBaselineModel(); //zeros nuisance param central values...
    pfLH.fitProfile();
    if (pfLH.getStatus() != 3)
      printf("CLfitSS::doCLs, Fit didn't converge(d2 2): %d\n", pfLH.getStatus());
  }
  llr_d -= pfLH.getFuncMin();

  //perform s1+b fit to s0+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = s0like.bkgd(i) + s0like.signal(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfitSS::doCLs, Fit didn't converge(sb2): %d\n", pfLH.getStatus());
  llr_sb_expected -= pfLH.getFuncMin();

  //perform s1+b fit to s1+b hypothesis
  //set reference data..
  sbBase.setBaselineModel(); //zeros nuisance param central values...
  for (int i = 0; i < nbins; i++)
    sbBase.data(i) = s1like.bkgd(i) + s1like.signal(i);
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("CLfitSS::doCLs, Fit didn't converge(b2): %d\n", pfLH.getStatus());
  llr_b_expected -= pfLH.getFuncMin();

  CLs.llrobs = llr_d;
  CLs.llrb = llr_b_expected;
  CLs.llrsb = llr_sb_expected;
  printf("CLfitSS:  LLRobs: %f, LLRb: %f, LLRsb: %f\n", CLs.llrobs, CLs.llrb, CLs.llrsb);
  double sum0p = 0;
  double sum1p = 0;
  double sum0m = 0;
  double sum1m = 0;
  double c0p = 0;
  double c1p = 0;
  double c0m = 0;
  double c1m = 0;

  for (int ift = 0; ift < 2 * its && (denom < its / 2.0 || ift < its); ift++) {

    rand0.fluctuate();
    rand1.fluctuate();

    for (int j = 0; j < nbins; j++) {
      s0like.data(j) = 0;
      s1like.data(j) = 0;
      if ((rand0.bkgd(j) + rand0.signal(j)) > 0) {
        zsb0 = rp_.fire(rand0.bkgd(j) + rand0.signal(j));
      } else
        zsb0 = 0;
      if ((rand1.bkgd(j) + rand1.signal(j)) > 0) {
        zsb1 = rp_.fire(rand1.bkgd(j) + rand1.signal(j));
      } else
        zsb1 = 0;

      s1like.data(j) = zsb1;
      s0like.data(j) = zsb0;
    }

    if (s0like.totData() > rand0.totData()) {
      sum0p += s0like.totData();
      c0p++;
    } else {
      sum0m += s0like.totData();
      c0m++;
    }

    if (s1like.totData() > rand0.totData()) {
      sum1p += s1like.totData();
      c1p++;
    } else {
      c1m++;
      sum1m += s1like.totData();
    }

    // s0+b fit to the s0+b pseudo-data
    sbBase.setSignalScale("graviton", nomScale * bsmScale_);
    sbBase.setSignalScale("pseudoscalar", nomScale * bsmScale_);
    sbBase.setSignalScale("vh", 0.0);

    sbBase.setSignalScale("GC", nomScale * bsmScale_);
    sbBase.setSignalScale("PS", nomScale * bsmScale_);
    sbBase.setSignalScale("SM", 0.0);

    sbBase.setSignalScale("grav", nomScale * bsmScale_);
    sbBase.setSignalScale("ps", nomScale * bsmScale_);
    sbBase.setSignalScale("signal", 0.0);

    sbBase.setBaselineModel(); //zeros nuisance param central values...
    sbBase.fluctuate();        //zeros nuisance param central values...
    sbBase.setBaselineModel(); //zeros nuisance param central values...

    for (int b = 0; b < nbins; b++)
      sbBase.data(b) = s0like.data(b);
    pfLH.fitProfile();
    if (pfLH.getStatus() != 3) {
      printf("CLfitSS::doCLs, Fit didn't converge(la): %d\n", pfLH.getStatus());
      continue;
    }
    llr_sb = pfLH.getFuncMin();

    //  s0+b fit to the s1+b pseudo-data
    sbBase.setBaselineModel();
    for (int b = 0; b < nbins; b++)
      sbBase.data(b) = s1like.data(b);
    pfLH.fitProfile();

    if (pfLH.getStatus() != 3) {
      printf("CLfitSS::doCLs, Fit didn't converge(lb): %d\n", pfLH.getStatus());
      continue;
    }
    llr_b = pfLH.getFuncMin();

    // s1+b fit to the s0+b pseudo-data
    sbBase.setSignalScale("graviton", nomScale * bsmScale_ * (1 - sigFraction_));
    sbBase.setSignalScale("pseudoscalar", nomScale * bsmScale_ * (1 - sigFraction_));
    sbBase.setSignalScale("vh", nomScale * sigFraction_);

    sbBase.setSignalScale("GC", nomScale * bsmScale_ * (1 - sigFraction_));
    sbBase.setSignalScale("PS", nomScale * bsmScale_ * (1 - sigFraction_));
    sbBase.setSignalScale("SM", nomScale * sigFraction_);

    sbBase.setSignalScale("grav", nomScale * bsmScale_ * (1 - sigFraction_));
    sbBase.setSignalScale("ps", nomScale * bsmScale_ * (1 - sigFraction_));
    sbBase.setSignalScale("signal", nomScale * sigFraction_);

    sbBase.setBaselineModel(); //zeros nuisance param central values...
    for (int b = 0; b < nbins; b++)
      sbBase.data(b) = s0like.data(b);
    pfLH.fitProfile();

    if (pfLH.getStatus() != 3) {
      printf("CLfitSS::doCLs, Fit didn't converge(lc): %d\n", pfLH.getStatus());
      continue;
    }
    llr_sb -= pfLH.getFuncMin();

    // s1+b fit to the s1+b pseudo-data
    sbBase.setBaselineModel();
    for (int b = 0; b < nbins; b++)
      sbBase.data(b) = s1like.data(b);

    pfLH.fitProfile();

    if (pfLH.getStatus() != 3) {
      printf("CLfitSS::doCLs, Fit didn't converge(ld): %d\n", pfLH.getStatus());
      continue;
    }
    llr_b -= pfLH.getFuncMin();

    if (llr_d <= llr_b)
      denom++; // background "less than" data

    if (llr_d <= llr_sb)
      numer++; // s+b "less than" data

    if (llr_b_expected <= llr_b)
      denom_bexp++; // background "less than" med expected b

    if (nmc < MaxKeep_) {
      btrial_[nmc] = llr_b;
      sbtrial_[nmc] = llr_sb;
    }
    nmc++;
  }

  printf("Getting Rates: %d, %f, %f/%f, %f/%f\n", nmc, rand0.totData(), sum0p, sum0m, sum1p, sum1m);
  sum0p /= (c0p * rand0.totData());
  sum0m /= (c0m * rand0.totData());
  sum1p /= (c1p * rand0.totData());
  sum1m /= (c1m * rand0.totData());

  printf("Rates: %f/%f, %f/%f\n", sum0p, sum0m, sum1p, sum1m);

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

void CLfitSS::configure(const char* options) {
  parseIterations(options, iterations_);
}

TH1D* CLfitSS::getLLRdist_sb(const char* title, int bins, double min, double max) {

  TH1D* out = new TH1D(title, title, bins, min, max);
  int lim = (nmcdone_ < MaxKeep_) ? nmcdone_ : MaxKeep_;

  for (int i = 0; i < lim; i++)
    out->Fill(sbtrial_[i]);

  return out;
}

TH1D* CLfitSS::getLLRdist_b(const char* title, int bins, double min, double max) {

  TH1D* out = new TH1D(title, title, bins, min, max);
  int lim = (nmcdone_ < MaxKeep_) ? nmcdone_ : MaxKeep_;

  for (int i = 0; i < lim; i++)
    out->Fill(btrial_[i]);

  return out;
}

void CLfitSS::noviceWarning() {

  printf("\n*****************************************************************\n");
  printf(
      "CLfitSS Warning:  You are using a fitting CL calculator.  Have you\nverified your systematics and fit model using the FitTest class?\n");
  printf("*****************************************************************\n\n");

  if (useHistoStats_) {
    printf("\n*****************************************************************\n");
    printf(
        "CLfitSS Warning:  You're operating in novice mode and you're using\nstatistical uncertainties from your input histograms.  Are you\nsure you've set the values properly in the histograms?\n");
    printf("*****************************************************************\n\n");
  } else {
    printf("\n*****************************************************************\n");
    printf(
        "CLfitSS Warning:  The histogram statistical uncertainty flag is turned\noff by default.  This calculation does not include statistical uncertainties.\n");
    printf("*****************************************************************\n\n");
  }
}
