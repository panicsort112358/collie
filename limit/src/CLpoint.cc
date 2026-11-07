#include "CLpoint.hh"
#include <TTree.h>

CLpoint::CLpoint(int var1, int var2, int var3) {
  reset(var1, var2, var3);
}

void CLpoint::reset(int var1, int var2, int var3) {
  m_var1 = var1;
  m_var2 = var2;
  m_var3 = var3;
  xsec_cl = 0;
  xsec_obsfactor = -1.0;
  xsec_medfactor = -1.0;
  xsec_medfactor_p2s = -1.0;
  xsec_medfactor_p1s = -1.0;
  xsec_medfactor_m1s = -1.0;
  xsec_medfactor_m2s = -1.0;
  signal_scale = 1.0;

  cls_obs = -1.0;
  clb_obs = -1.0;
  clsb_obs = -1.0;

  clb_sb = -1.0;
  clb_sb_p1s = -1.0;
  clb_sb_p2s = -1.0;
  clb_sb_m1s = -1.0;
  clb_sb_m2s = -1.0;

  clb_med = -1.0;
  clsb_med = -1.0;
  cls_med = -1.0;
  cls_med_p1s = -1.0;
  cls_med_p2s = -1.0;
  cls_med_m1s = -1.0;
  cls_med_m2s = -1.0;
  clsb_med_p1s = -1.0;
  clsb_med_p2s = -1.0;
  clsb_med_m1s = -1.0;
  clsb_med_m2s = -1.0;

  llrobs = -1.0;
  llrb = -1.0;
  llrb_p1s = -1.0;
  llrb_p2s = -1.0;
  llrb_m1s = -1.0;
  llrb_m2s = -1.0;

  llrsb = -1.0;
  llrsb_p1s = -1.0;
  llrsb_p2s = -1.0;
  llrsb_m1s = -1.0;
  llrsb_m2s = -1.0;

  nTrials_exp = 0;
  nTrials_obs = 0;

  fit_sigScale = -1;
  fit_sigScale_Err = -1;
  fit_sigScale_ErrP = -1;
  fit_sigScale_ErrM = -1;
  fit_sigScale_ErrStat = -1;
  fit_sigScale_ErrStatP = -1;
  fit_sigScale_ErrStatM = -1;
  fit_sigScale_ErrSyst = -1;
  fit_sigScale_ErrSystP = -1;
  fit_sigScale_ErrSystM = -1;
  fit_chi2_b = -1;
  fit_chi2_sb = -1;
}

double CLpoint::getLumiFactor() {
  if (llrb != -1 && llrsb != -1) {
    return 7.82 / (llrb - llrsb);
  }
  return -1;
}

void CLpoint::print() {
  printf("\n****************************************:\n");
  printf("CLcompute Results:\n");
  printf("      CLs_obs: %.4f,  CLs_med: %.4f\n", cls_obs, cls_med);
  printf("     CLsb_obs: %.4f, CLsb_med: %.4f\n", clsb_obs, clsb_med);
  printf("      CLb_obs: %.4f,  CLb_med: %.4f\n", clb_obs, clb_med);
  printf("     LLRobs: %.4f, LLRb: %.4f, LLRsb: %.4f, lumifactor: %.3f\n", llrobs, llrb, llrsb,
         getLumiFactor());
  printf("     LLRb_m2s: %.4f, LLRb_m1s: %.4f\n", llrb_m2s, llrb_m1s);
  printf("     LLRb_p1s: %.4f, LLRb_p2s: %.4f\n", llrb_p1s, llrb_p2s);
  printf("****************************************:\n");

  return;
}

void CLpoint::branch(TTree* tree) {
  tree->Branch("var1", &m_var1, "var1/I");
  tree->Branch("var2", &m_var2, "var2/I");
  tree->Branch("var3", &m_var3, "var3/I");

  tree->Branch("cls_obs", &cls_obs, "cls_obs/D");
  tree->Branch("clb_obs", &clb_obs, "clb_obs/D");
  tree->Branch("clsb_obs", &clsb_obs, "clsb_obs/D");
  tree->Branch("clb_med", &clb_med, "clb_med/D");

  tree->Branch("clb_sb", &clb_sb, "clb_sb/D");
  tree->Branch("clb_sb_p2s", &clb_sb_p2s, "clb_sb_p2s/D");
  tree->Branch("clb_sb_p1s", &clb_sb_p1s, "clb_sb_p1s/D");
  tree->Branch("clb_sb_m1s", &clb_sb_m1s, "clb_sb_m1s/D");
  tree->Branch("clb_sb_m2s", &clb_sb_m2s, "clb_sb_m2s/D");

  tree->Branch("cls_med", &cls_med, "cls_med/D");
  tree->Branch("cls_med_p2s", &cls_med_p2s, "cls_med_p2s/D");
  tree->Branch("cls_med_p1s", &cls_med_p1s, "cls_med_p1s/D");
  tree->Branch("cls_med_m1s", &cls_med_m1s, "cls_med_m1s/D");
  tree->Branch("cls_med_m2s", &cls_med_m2s, "cls_med_m2s/D");

  tree->Branch("clsb_med", &clsb_med, "clsb_med/D");
  tree->Branch("clsb_med_p2s", &clsb_med_p2s, "clsb_med_p2s/D");
  tree->Branch("clsb_med_p1s", &clsb_med_p1s, "clsb_med_p1s/D");
  tree->Branch("clsb_med_m1s", &clsb_med_m1s, "clsb_med_m1s/D");
  tree->Branch("clsb_med_m2s", &clsb_med_m2s, "clsb_med_m2s/D");

  tree->Branch("llrobs", &llrobs, "llrobs/D");
  tree->Branch("llrb", &llrb, "llrb/D");
  tree->Branch("llrb_p2s", &llrb_p2s, "llrb_p2s/D");
  tree->Branch("llrb_p1s", &llrb_p1s, "llrb_p1s/D");
  tree->Branch("llrb_m1s", &llrb_m1s, "llrb_m1s/D");
  tree->Branch("llrb_m2s", &llrb_m2s, "llrb_m2s/D");

  tree->Branch("llrsb", &llrsb, "llrsb/D");
  tree->Branch("llrsb_p2s", &llrsb_p2s, "llrsb_p2s/D");
  tree->Branch("llrsb_p1s", &llrsb_p1s, "llrsb_p1s/D");
  tree->Branch("llrsb_m1s", &llrsb_m1s, "llrsb_m1s/D");
  tree->Branch("llrsb_m2s", &llrsb_m2s, "llrsb_m2s/D");

  tree->Branch("xsec_cl", &xsec_cl, "xsec_cl/D");
  tree->Branch("xsec_obsfactor", &xsec_obsfactor, "xsec_obsfactor/D");
  tree->Branch("xsec_medfactor", &xsec_medfactor, "xsec_medfactor/D");
  tree->Branch("xsec_medfactor_p2s", &xsec_medfactor_p2s, "xsec_medfactor_p2s/D");
  tree->Branch("xsec_medfactor_p1s", &xsec_medfactor_p1s, "xsec_medfactor_p1s/D");
  tree->Branch("xsec_medfactor_m1s", &xsec_medfactor_m1s, "xsec_medfactor_m1s/D");
  tree->Branch("xsec_medfactor_m2s", &xsec_medfactor_m2s, "xsec_medfactor_m2s/D");
  tree->Branch("signal_scale", &signal_scale, "signal_scale/D");

  tree->Branch("nTrials_exp", &nTrials_exp, "nTrials_exp/I");
  tree->Branch("nTrials_obs", &nTrials_obs, "nTrials_obs/I");

  tree->Branch("fit_sigScale", &fit_sigScale, "fit_sigScale/D");
  tree->Branch("fit_sigScale_Err", &fit_sigScale_Err, "fit_sigScale_Err/D");
  tree->Branch("fit_sigScale_ErrP", &fit_sigScale_ErrP, "fit_sigScale_ErrP/D");
  tree->Branch("fit_sigScale_ErrM", &fit_sigScale_ErrM, "fit_sigScale_ErrM/D");
  tree->Branch("fit_sigScale_ErrStat", &fit_sigScale_ErrStat, "fit_sigScale_ErrStat/D");
  tree->Branch("fit_sigScale_ErrStatP", &fit_sigScale_ErrStatP, "fit_sigScale_ErrStatP/D");
  tree->Branch("fit_sigScale_ErrStatM", &fit_sigScale_ErrStatM, "fit_sigScale_ErrStatM/D");
  tree->Branch("fit_sigScale_ErrSyst", &fit_sigScale_ErrSyst, "fit_sigScale_ErrSyst/D");
  tree->Branch("fit_sigScale_ErrSystP", &fit_sigScale_ErrSystP, "fit_sigScale_ErrSystP/D");
  tree->Branch("fit_sigScale_ErrSystM", &fit_sigScale_ErrSystM, "fit_sigScale_ErrSystM/D");
  tree->Branch("fit_chi2_b", &fit_chi2_b, "fit_chi2_b/D");
  tree->Branch("fit_chi2_sb", &fit_chi2_sb, "fit_chi2_sb/D");
}

void CLpoint::mapToTree(TTree* tree) {

  tree->SetBranchAddress("var1", &m_var1);
  tree->SetBranchAddress("var2", &m_var2);
  tree->SetBranchAddress("var3", &m_var3);

  tree->SetBranchAddress("cls_obs", &cls_obs);
  tree->SetBranchAddress("clb_obs", &clb_obs);
  tree->SetBranchAddress("clsb_obs", &clsb_obs);
  tree->SetBranchAddress("clb_med", &clb_med);

  tree->SetBranchAddress("clb_sb", &clb_sb);
  tree->SetBranchAddress("clb_sb_p2s", &clb_sb_p2s);
  tree->SetBranchAddress("clb_sb_p1s", &clb_sb_p1s);
  tree->SetBranchAddress("clb_sb_m1s", &clb_sb_m1s);
  tree->SetBranchAddress("clb_sb_m2s", &clb_sb_m2s);

  tree->SetBranchAddress("cls_med", &cls_med);
  tree->SetBranchAddress("cls_med_p2s", &cls_med_p2s);
  tree->SetBranchAddress("cls_med_p1s", &cls_med_p1s);
  tree->SetBranchAddress("cls_med_m1s", &cls_med_m1s);
  tree->SetBranchAddress("cls_med_m2s", &cls_med_m2s);

  tree->SetBranchAddress("clsb_med", &clsb_med);
  tree->SetBranchAddress("clsb_med_p2s", &clsb_med_p2s);
  tree->SetBranchAddress("clsb_med_p1s", &clsb_med_p1s);
  tree->SetBranchAddress("clsb_med_m1s", &clsb_med_m1s);
  tree->SetBranchAddress("clsb_med_m2s", &clsb_med_m2s);

  tree->SetBranchAddress("llrobs", &llrobs);
  tree->SetBranchAddress("llrb", &llrb);
  tree->SetBranchAddress("llrb_p2s", &llrb_p2s);
  tree->SetBranchAddress("llrb_p1s", &llrb_p1s);
  tree->SetBranchAddress("llrb_m1s", &llrb_m1s);
  tree->SetBranchAddress("llrb_m2s", &llrb_m2s);

  tree->SetBranchAddress("llrsb", &llrsb);
  tree->SetBranchAddress("llrsb_p2s", &llrsb_p2s);
  tree->SetBranchAddress("llrsb_p1s", &llrsb_p1s);
  tree->SetBranchAddress("llrsb_m1s", &llrsb_m1s);
  tree->SetBranchAddress("llrsb_m2s", &llrsb_m2s);

  tree->SetBranchAddress("xsec_cl", &xsec_cl);
  tree->SetBranchAddress("xsec_obsfactor", &xsec_obsfactor);
  tree->SetBranchAddress("xsec_medfactor", &xsec_medfactor);
  tree->SetBranchAddress("xsec_medfactor_p2s", &xsec_medfactor_p2s);
  tree->SetBranchAddress("xsec_medfactor_p1s", &xsec_medfactor_p1s);
  tree->SetBranchAddress("xsec_medfactor_m1s", &xsec_medfactor_m1s);
  tree->SetBranchAddress("xsec_medfactor_m2s", &xsec_medfactor_m2s);
  tree->SetBranchAddress("signal_scale", &signal_scale);

  tree->SetBranchAddress("nTrials_exp", &nTrials_exp);
  tree->SetBranchAddress("nTrials_obs", &nTrials_obs);

  tree->SetBranchAddress("fit_sigScale", &fit_sigScale);
  tree->SetBranchAddress("fit_sigScale_Err", &fit_sigScale_Err);
  tree->SetBranchAddress("fit_sigScale_ErrP", &fit_sigScale_ErrP);
  tree->SetBranchAddress("fit_sigScale_ErrM", &fit_sigScale_ErrM);
  tree->SetBranchAddress("fit_sigScale_ErrStat", &fit_sigScale_ErrStat);
  tree->SetBranchAddress("fit_sigScale_ErrStatP", &fit_sigScale_ErrStatP);
  tree->SetBranchAddress("fit_sigScale_ErrStatM", &fit_sigScale_ErrStatM);
  tree->SetBranchAddress("fit_sigScale_ErrSyst", &fit_sigScale_ErrSyst);
  tree->SetBranchAddress("fit_sigScale_ErrSystP", &fit_sigScale_ErrSystP);
  tree->SetBranchAddress("fit_sigScale_ErrSystM", &fit_sigScale_ErrSystM);
  tree->SetBranchAddress("fit_chi2_b", &fit_chi2_b);
  tree->SetBranchAddress("fit_chi2_sb", &fit_chi2_sb);
}
