void plotCL(const char* filename) {
  TFile f(filename);
  TTree* t = (TTree*)f.Get("SCAN");
  if (t == NULL) {
    printf("Cannot open the 'SCAN' Tree in '%s'\n", filename);
    return;
  }

  int a_v1;
  t->SetBranchAddress("var1", &a_v1);
  double cl_obs;
  t->SetBranchAddress("clobs", &cl_obs);

  double clmed1, clmed_p2s, clmed_p1s, clmed_m1s, clmed_m2s;
  t->SetBranchAddress("cls_med", &clmed1);
  t->SetBranchAddress("cls_med_p2s", &clmed_p2s);
  t->SetBranchAddress("cls_med_p1s", &clmed_p1s);
  t->SetBranchAddress("cls_med_m1s", &clmed_m1s);
  t->SetBranchAddress("cls_med_m2s", &clmed_m2s);

  int len = t->GetEntries();
  double* v1 = new double[len];
  double* v1_sigma = new double[2 * len + 1];
  double* clobs = new double[len];
  double* clmed = new double[len];
  double* clmed_1s = new double[2 * len + 1];
  double* clmed_2s = new double[2 * len + 1];
  double* llr_b = new double[len];
  double* llr_b_1s = new double[2 * len + 1];
  double* llr_b_2s = new double[2 * len + 1];

  double lmin = 10000, lmax = -10000;

  for (int i = 0; i < len; i++) {
    t->GetEntry(i);
    v1[i] = a_v1;
    clobs[i] = cl_obs;
    clmed[i] = clmed1;
    v1_sigma[i] = a_v1;
    v1_sigma[2 * len - 1 - i] = a_v1;

    clmed_1s[i] = clmed_m1s;
    clmed_1s[2 * len - 1 - i] = clmed_p1s;
    clmed_2s[i] = clmed_m2s;
    clmed_2s[2 * len - 1 - i] = clmed_p2s;

    lmin = TMath::Min(
        lmin, TMath::Min(clmed_p2s, TMath::Min(clmed_p1s, TMath::Min(clmed_m1s, clmed_m2s))));
    lmax = TMath::Max(
        lmax, TMath::Max(clmed_p2s, TMath::Max(clmed_p1s, TMath::Max(clmed_m1s, clmed_m2s))));
  }
  f.Close();

  v1_sigma[2 * len] = v1_sigma[0];
  clmed_1s[2 * len] = clmed_1s[0];
  clmed_2s[2 * len] = clmed_2s[0];

  TGraph* g_clobs = new TGraph(len, v1, clobs);

  TGraph* g_clmed = new TGraph(len, v1, clmed);
  TGraph* g_clmed_1s = new TGraph(2 * len + 1, v1_sigma, clmed_1s);
  TGraph* g_clmed_2s = new TGraph(2 * len + 1, v1_sigma, clmed_2s);

  TH2* h = new TH2F("h", "Confidence Level", 10, v1[0], v1[len - 1], 10, lmin, lmax);
  h->SetStats(kFALSE);
  h->SetYTitle("1-CL_{s}");
  h->Draw("");

  g_clmed_2s->SetFillColor(5);
  g_clmed_2s->Draw("F");

  g_clmed_1s->SetFillColor(3);
  g_clmed_1s->Draw("F");

  g_clmed->SetLineWidth(2);
  g_clmed->SetLineStyle(2);
  g_clmed->Draw("L");
  g_clobs->SetLineWidth(2);
  g_clobs->Draw("L");
}
