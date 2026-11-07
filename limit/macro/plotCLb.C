void plotCLb(const char* filename) {
  TFile f(filename);
  TTree* t = (TTree*)f.Get("SCAN");
  if (t == NULL) {
    printf("Cannot open the 'SCAN' Tree in '%s'\n", filename);
    return;
  }

  int a_v1;
  t->SetBranchAddress("var1", &a_v1);

  double clb1;
  t->SetBranchAddress("clb_med", &clb1);

  int len = t->GetEntries();
  double* v1 = new double[len];
  double* clb = new double[len];

  double lmin = 10000, lmax = -10000;

  for (int i = 0; i < len; i++) {
    t->GetEntry(i);
    v1[i] = a_v1;
    clb[i] = 1.0 - clb1;
  }
  f.Close();

  TGraph* g_clb = new TGraph(len, v1, clb);

  TCanvas* c1 = new TCanvas("CLb", "CLb");

  TH2* h = new TH2F("h", "Confidence Level for Bkgd", 10, 100.0, 140.0, 10, 0.0005, 1.0);
  h->SetStats(kFALSE);
  h->SetYTitle("1-CL_{b}");
  h->SetXTitle("m_{H} (GeV)");
  c1->cd();
  c1->SetLogy();
  h->Draw("");

  g_clb->SetLineWidth(2);
  g_clb->Draw("L");

  c1->SetGridx();
  c1->SetGridy();

  c1->Modified();
}
