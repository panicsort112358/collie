// ... Global attributes go here ...
void style() {

  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);

  gStyle->SetTitleColor(1); // 0
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(1);
  /*
   gStyle->SetTitleX(0.10);
   gStyle->SetTitleY(0.94);
   gStyle->SetTitleW(0.5);
   gStyle->SetTitleH(0.06);

   gStyle->SetLabelOffset(1e-04);
   gStyle->SetLabelSize(0.2);
   gStyle->SetTitleSize(0.2);
   */
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.90);
  gStyle->SetStatY(0.90);
  gStyle->SetStatW(0.30);
  gStyle->SetStatH(0.10);

  gStyle->SetErrorX(0.0); // Horizontal error bar size
  gStyle->SetOptStat(0);
  //   gStyle->SetPaperSize(10.,12.);   // Printout size
}
void h2_style(TH2* h) {

  h->SetLabelFont(62, "X");      // 42
  h->SetLabelFont(62, "Y");      // 42
  h->SetLabelOffset(0.000, "X"); // D=0.005
  h->SetLabelOffset(0.005, "Y"); // D=0.005
  h->SetLabelSize(0.055, "X");
  h->SetLabelSize(0.055, "Y");
  h->SetTitleOffset(0.85, "X");
  h->SetTitleOffset(0.85, "Y");
  h->SetTitleSize(0.06, "X");
  h->SetTitleSize(0.06, "Y");
  h->SetTitle(0);
}

void plotLLR(const char* filename) {
  style();
  TFile f(filename);

  double xmin = 115;
  double xmax = 200;
  double ymax = 1;
  int type = 99; //0=full combo
                 //1=wh combo
                 //2=zh combo
                 //3=HWW combo
                 //4=WWW combo
                 //5=whev combo
                 //6=whmv combo
                 //7=zhvv combo
                 //8=zhll combo
                 //99=no extra text

  TTree* t = (TTree*)f.Get("SCAN");
  if (t == NULL) {
    printf("Cannot open the 'SCAN' Tree in '%s'\n", filename);
    return;
  }

  int a_v1;
  t->SetBranchAddress("var1", &a_v1);
  double llrobs;
  t->SetBranchAddress("llrobs", &llrobs);

  double llrsb, llrsb_p2s, llrsb_p1s, llrsb_m1s, llrsb_m2s;
  t->SetBranchAddress("llrsb", &llrsb);
  t->SetBranchAddress("llrsb_p2s", &llrsb_p2s);
  t->SetBranchAddress("llrsb_p1s", &llrsb_p1s);
  t->SetBranchAddress("llrsb_m1s", &llrsb_m1s);
  t->SetBranchAddress("llrsb_m2s", &llrsb_m2s);
  double llrb, llrb_p2s, llrb_p1s, llrb_m1s, llrb_m2s, clb1;
  t->SetBranchAddress("llrb", &llrb);
  t->SetBranchAddress("llrb_p2s", &llrb_p2s);
  t->SetBranchAddress("llrb_p1s", &llrb_p1s);
  t->SetBranchAddress("llrb_m1s", &llrb_m1s);
  t->SetBranchAddress("llrb_m2s", &llrb_m2s);

  int len = t->GetEntries();
  int lenR = 0;
  for (int i = 0; i < len; i++) {
    t->GetEntry(i);
    if (a_v1 > xmax)
      continue;
    if (a_v1 < xmin)
      continue;
    lenR++;
  }

  double* v1 = new double[lenR];
  double* llr_obs = new double[lenR];
  double* llr_sb = new double[lenR];
  double* llr_b = new double[lenR];
  lenR = 0;
  for (int i = 0; i < len; i++) {
    t->GetEntry(i);
    if (a_v1 > xmax)
      continue;
    if (a_v1 < xmin)
      continue;
    v1[lenR] = a_v1;
    llr_obs[lenR] = llrobs;
    llr_sb[lenR] = llrsb;
    llr_b[lenR] = llrb;
    lenR++;
  }
  double* v1_sigma = new double[2 * lenR + 1];
  double* llr_sb_1s = new double[2 * lenR + 1];
  double* llr_sb_2s = new double[2 * lenR + 1];
  double* llr_b_1s = new double[2 * lenR + 1];
  double* llr_b_2s = new double[2 * lenR + 1];
  double lmin = 10000, lmax = -10000;

  int valC = 0;
  for (int i = 0; i < len; i++) {
    t->GetEntry(i);
    if (a_v1 > xmax)
      continue;
    if (a_v1 < xmin)
      continue;
    v1_sigma[valC] = a_v1;
    v1_sigma[2 * lenR - 1 - valC] = a_v1;

    llr_b_1s[valC] = llrb_m1s;
    llr_b_1s[2 * lenR - 1 - valC] = llrb_p1s;
    llr_b_2s[valC] = llrb_m2s;
    llr_b_2s[2 * lenR - 1 - valC] = llrb_p2s;

    llr_sb_1s[valC] = llrsb_m1s;
    llr_sb_1s[2 * lenR - 1 - valC] = llrsb_p1s;
    llr_sb_2s[valC] = llrsb_m2s;
    llr_sb_2s[2 * lenR - 1 - valC] = llrsb_p2s;
    valC++;
    lmin = TMath::Min(lmin,
                      TMath::Min(llrb_p2s, TMath::Min(llrb_p1s, TMath::Min(llrb_m1s, llrb_m2s))));
    lmax = TMath::Max(lmax,
                      TMath::Max(llrb_p2s, TMath::Max(llrb_p1s, TMath::Max(llrb_m1s, llrb_m2s))));
    lmin = TMath::Min(
        lmin, TMath::Min(llrsb_p2s, TMath::Min(llrsb_p1s, TMath::Min(llrsb_m1s, llrsb_m2s))));
    lmax = TMath::Max(
        lmax, TMath::Max(llrsb_p2s, TMath::Max(llrsb_p1s, TMath::Max(llrsb_m1s, llrsb_m2s))));
  }
  f.Close();

  if (lmin > (-1.0 * ymax))
    lmin = -1.0 * ymax;
  else
    lmin *= 1.1;

  if (lmax < (ymax))
    lmax = ymax;
  else
    lmax *= 1.1;

  v1_sigma[2 * lenR] = v1_sigma[0];
  llr_b_1s[2 * lenR] = llr_b_1s[0];
  llr_b_2s[2 * lenR] = llr_b_2s[0];
  llr_sb_1s[2 * lenR] = llr_sb_1s[0];
  llr_sb_2s[2 * lenR] = llr_sb_2s[0];

  TGraph* g_llrobs = new TGraph(lenR, v1, llr_obs);
  TGraph* g_llrb = new TGraph(lenR, v1, llr_b);
  TGraph* g_llrb_1s = new TGraph(2 * lenR + 1, v1_sigma, llr_b_1s);
  TGraph* g_llrb_2s = new TGraph(2 * lenR + 1, v1_sigma, llr_b_2s);
  TGraph* g_llrsb = new TGraph(lenR, v1, llr_sb);
  TCanvas* c1 = new TCanvas("ht", "ht");

  TH2* h = new TH2F("h", "LLR Distributions", 10, xmin, xmax, 10, lmin, lmax);

  h->SetTitle("");
  h->SetStats(kFALSE);
  //  h->SetYTitle("-2 ln(Q)");
  h->SetYTitle("LLR");
  h->SetXTitle("m_{H} (GeV/c^{2})");
  h->Draw("");

  g_llrb_2s->SetFillColor(5);
  g_llrb_2s->SetLineColor(5);
  g_llrb_2s->SetMarkerColor(5);
  g_llrb_2s->Draw("F");

  g_llrb_1s->SetFillColor(3);
  g_llrb_1s->SetLineColor(3);
  g_llrb_1s->SetMarkerColor(3);
  g_llrb_1s->Draw("F");

  g_llrb->SetFillColor(0);
  g_llrb->SetLineWidth(4);
  g_llrb->SetLineStyle(2);
  g_llrb->Draw("C");

  g_llrsb->SetFillColor(0);
  g_llrsb->SetLineWidth(4);
  g_llrsb->SetLineStyle(4);
  g_llrsb->SetLineColor(2);
  g_llrsb->Draw("C");

  g_llrobs->SetFillColor(0);
  g_llrobs->SetLineWidth(4);
  g_llrobs->Draw("C");

  c1->SetGridx();
  c1->SetGridy();

  TLegend* leg = new TLegend(0.73, 0.66, 0.98, 0.97, NULL, "brNDC");
  leg->AddEntry(g_llrb_2s, "LLR_{B} 2-#sigma");
  leg->AddEntry(g_llrb_1s, "LLR_{B} 1-#sigma");
  leg->AddEntry(g_llrb, "LLR_{B}");
  leg->AddEntry(g_llrsb, "LLR_{S+B}");
  leg->AddEntry(g_llrobs, "LLR_{OBS}");
  leg->SetFillColor(0);
  if (type != 99)
    leg->Draw("same");

  // Need only for outside world
  TLatex* tex0 = 0;
  if (type == 0 || type == 1 || type == 3 || type == 4 || type == 5)
    tex0 = new TLatex(125, 0, "D\349 Preliminary, L=0.4 fb^{-1}");
  if (type == 2 || type == 6 || type == 7)
    tex0 = new TLatex(125, 0, "D\349 Preliminary, L=0.4 fb^{-1}");
  if (type != 99) {
    tex0->SetTextSize(0.06);
    tex0->SetTextColor(1);
    tex0->Draw("same");
  }

  TLatex* tex3;
  if (type == 0)
    tex3 = new TLatex(105, 0, "Full D\349 Combination");
  if (type == 1)
    tex3 = new TLatex(105, 0, "WH (H#rightarrowb#bar{b}) Combination");
  if (type == 2)
    tex3 = new TLatex(105, 0, "ZH (H#rightarrowb#bar{b}) Combination");
  if (type == 3)
    tex3 = new TLatex(135, 0, "H#rightarrowW^{+}W^{-}");
  if (type == 4)
    tex3 = new TLatex(135, 0, "WH#rightarrowW W^{+}W^{-}");
  if (type == 5)
    tex3 = new TLatex(105, 0, "WH#rightarrow e#nub#bar{b}, 1+2 Tags");
  if (type == 6)
    tex3 = new TLatex(105, 0, "WH#rightarrow#mu#nub#bar{b}, 1+2 Tags");
  if (type == 7)
    tex3 = new TLatex(105, 0, "ZH#rightarrow#nu#nub#bar{b}, WH signal");
  if (type == 8)
    tex3 = new TLatex(105, 0, "ZH#rightarrow ll b#bar{b}");
  if (type != 99) {
    tex3->SetTextSize(0.06);
    tex3->SetTextColor(1);
    tex3->Draw();
  }
  c1->Modified();

  return;
}
