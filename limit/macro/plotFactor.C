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

void plotFactor(const char* filename) {
  style();
  float yMin = -1;
  float yMax = 55;
  float xMin = 115;
  float xMax = 200;
  int type = 99; //0=full combo
                 //1=wh combo
                 //2=zh combo
                 //3=HWW combo
                 //4=whev combo
                 //5=whmv combo
                 //6=zhvv combo
                 //7=zhll combo

  TFile* f6 = new TFile(filename);
  TTree* t6 = (TTree*)f6->Get("SCAN");
  int len6 = t6->GetEntries();
  int a_v6;
  t6->SetBranchAddress("var1", &a_v6);

  double med6, obs6;
  t6->SetBranchAddress("xsec_medfactor", &med6);
  t6->SetBranchAddress("xsec_obsfactor", &obs6);

  double* v6 = new double[len6];
  double* m6 = new double[len6];
  double* o6 = new double[len6];

  for (int i = 0; i < len6; i++) {
    t6->GetEntry(i);
    v6[i] = a_v6;
    m6[i] = med6;
    o6[i] = obs6;
    printf("Mass: %d, exp: %.2f, obs: %.2f\n", a_v6, med6, obs6);
  }

  TGraph* fact_obs = new TGraph(len6, v6, o6);
  TGraph* fact_med = new TGraph(len6, v6, m6);

  TCanvas* c6 = new TCanvas("factor", "factor");
  c6->SetGridx();
  c6->SetGridy();
  TH2* h6 = new TH2F("h6", "Cross-Section Limit", 20, xMin, xMax, 100, yMin, yMax);
  h6->SetTitle("");
  h6->SetStats(kFALSE);
  if (type == 0)
    h6->SetYTitle(
        "Limit / #sigma(p#bar{p}#rightarrowWH/ZH/H)#timesBR(H#rightarrowb#bar{b}/W^{+}W^{-})");
  if (type == 1 || type == 4 || type == 5)
    h6->SetYTitle("Limit / #sigma(p#bar{p}#rightarrowWH)#timesBR(H#rightarrowb#bar{b})");
  if (type == 2 || type == 6 || type == 7)
    h6->SetYTitle("Limit / #sigma(p#bar{p}#rightarrowZH)#timesBR(H#rightarrowb#bar{b})");
  if (type == 3)
    h6->SetYTitle("Limit / #sigma(p#bar{p}#rightarrowH)#times BR(H#rightarrowW^{+}W^{-})");

  h6->SetXTitle("m_{H} (GeV/c^{2})");
  h6->Draw("");

  fact_med->SetFillColor(0);
  fact_med->SetLineWidth(5);
  fact_med->SetLineStyle(2);
  fact_med->SetLineColor(2);
  fact_med->Draw("C");

  fact_obs->SetFillColor(0);
  fact_obs->SetLineWidth(5);
  fact_obs->SetLineStyle(1);
  fact_obs->SetLineColor(1);
  fact_obs->Draw("C");

  // Need only for outside world
  TText* label = new TText();
  label->SetTextFont(62);
  label->SetTextColor(1); // 4
  label->SetTextAlign(12);
  label->SetTextSize(0.06);
  //  label->DrawTextNDC(0.51,0.16,"D\349 Run II Preliminary");
  TLatex* tex0;
  if (type == 0 || type == 1 || type == 3 || type == 4 || type == 5)
    tex0 = new TLatex(125, 0, "D\349 Preliminary, L=1.0 fb^{-1}");
  if (type == 2 || type == 6 || type == 7)
    tex0 = new TLatex(125, 0, "D\349 Preliminary, L=0.9 fb^{-1}");
  tex0->SetTextSize(0.06);
  tex0->SetTextColor(1);
  tex0->Draw("same");

  // Need only for outside world
  //  TLatex *tex1 = new TLatex(105,-0.22,"363-384 pb^{-1} / 3 Analyses");
  //  tex1->SetTextSize(0.04);
  //  tex1->SetTextColor(1);
  //  tex1->Draw();             // Comment this line if result is final

  TLegend* leg4 = new TLegend(0.50, 0.75, 0.85, 0.90, NULL, "brNDC");
  leg4->AddEntry(fact_obs, "Observed Limit", "L");
  leg4->AddEntry(fact_med, "Expected Limit", "L");

  //  leg4->AddEntry(fact_obs,"P14 Exp Limit","L");
  //  leg4->AddEntry(fact_med,"P17 Exp Limit","L");

  TLatex* tex3;
  if (type == 0)
    tex3 = new TLatex(105, 10, "Full D\349 Combination");
  if (type == 1)
    tex3 = new TLatex(105, 10, "WH (H#rightarrowb#bar{b}) Combination");
  if (type == 2)
    tex3 = new TLatex(105, 10, "ZH (H#rightarrowb#bar{b}) Combination");
  if (type == 3)
    tex3 = new TLatex(105, 10, "H#rightarrowW^{+}W^{-}");
  if (type == 4)
    tex3 = new TLatex(105, 10, "WH#rightarrow e#nub#bar{b}, 3Jets 1Tag");
  if (type == 5)
    tex3 = new TLatex(105, 10, "WH#rightarrow#mu#nub#bar{b}, 3Jets");
  if (type == 6)
    tex3 = new TLatex(105, 10, "ZH#rightarrow#nu#nub#bar{b}, WH signal");
  if (type == 7)
    tex3 = new TLatex(105, 10, "ZH#rightarrow ll b#bar{b}");
  tex3->SetTextSize(0.06);
  tex3->SetTextColor(1);
  tex3->Draw();

  leg4->SetFillColor(0);
  leg4->Draw("same");

  double lenSM = 2;
  double xSM[2];
  double ySM[2];
  xSM[0] = xMin;
  xSM[1] = xMax;
  ySM[0] = 1;
  ySM[1] = 1;

  TGraph* standardModel = new TGraph(lenSM, xSM, ySM);
  standardModel->SetLineWidth(4);
  standardModel->Draw("L");
  TLatex* tex2 = new TLatex(162, 0.7, "Standard Model = 1.0");
  tex2->SetTextSize(0.05);
  tex2->SetTextColor(1);
  tex2->Draw();
}
