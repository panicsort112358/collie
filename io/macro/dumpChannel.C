
void dumpChannel(const char* infile) {
  gROOT->LoadMacro("style2004.C");
  style();

  TFile* f0 = new TFile(infile);
  TTree* t0 = (TTree*)f0->Get("SCAN");
  int len0 = t0->GetEntries();
  double* v0 = new double[len0];

  int a_v0;
  t0->SetBranchAddress("var1", &a_v0);

  double med0, obs0;
  t0->SetBranchAddress("xsec_medfactor", &med0);
  t0->SetBranchAddress("xsec_obsfactor", &obs0);

  double clMed, clObs, clBmed, clB, clBsb;
  t0->SetBranchAddress("cls_obs", &clObs);
  t0->SetBranchAddress("clb_obs", &clB);
  t0->SetBranchAddress("clb_sb", &clBsb);
  t0->SetBranchAddress("clb_med", &clBmed);
  t0->SetBranchAddress("cls_med", &clMed);

  double llrObs, llrB, llrSB;
  t0->SetBranchAddress("llrobs", &llrObs);
  t0->SetBranchAddress("llrb", &llrB);
  t0->SetBranchAddress("llrsb", &llrSB);

  for (int i = 0; i < len0; i++) {
    t0->GetEntry(i);
    printf("mass: %d ", a_v0);
    printf("CLobs: %f, CLmed: %f, CLBmed: %f, CLBobs: %f, CLBsb: %f\n", clObs, clMed, clBmed,
           1.0 - clB, clBsb);
    printf("LLRobs:%.3f, LLRb: %.5f, LLRsb: %.5f\n", llrObs, llrB, llrSB);
    printf("DeltaLLR: %.7f, lumifactor: %.2f\n", llrB - llrSB, 7.82 / (llrB - llrSB));
    printf("med: %.3f, obs: %.3f\n", med0, obs0);
  }
}
