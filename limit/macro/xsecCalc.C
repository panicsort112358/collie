double gausSignif(double pValue, bool sm = false) {

  double n = sqrt(2.0) * TMath::ErfInverse(1.0 - 2.0 * pValue);
  if (sm == true && pValue > 0.50) {
    n = sqrt(2.0) * TMath::ErfInverse(2.0 * (pValue - 0.50));
  }

  if (!sm) {
    printf("A pValue of %f (%fe-6) corresponds\n", pValue, pValue * 1e6);
    printf("to a one-sided Gaussian significance of %0.3f sigma\n", n);
  } else {
    printf("A pValue of %f (%fe-6) corresponds\n", pValue, pValue * 1e6);
    if (pValue < 0.50)
      printf("to a one-sided Gaussian significance of %0.3f sigma above the median value\n", n);
    else
      printf("to a one-sided Gaussian significance of %0.3f sigma below the median value\n", n);
  }
  return n;
}

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
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.90);
  gStyle->SetStatY(0.90);
  gStyle->SetStatW(0.30);
  gStyle->SetStatH(0.10);
  gStyle->SetErrorX(0.0); // Horizontal error bar size
  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
}

// ... 1D histogram attributes;
void h1_style(TH1* h, int line_width = 3, int line_color = 1, int line_style = 1,
              int fill_style = 1001, int fill_color = 50, float y_min = -1111.,
              float y_max = -1111., int ndivx = 510, int ndivy = 510, int marker_style = 20,
              int marker_color = 1, float marker_size = 0.5, int optstat = 0) {
  h->SetStats(optstat);
  h->SetLabelFont(62, "X");      // 42
  h->SetLabelFont(62, "Y");      // 42
  h->SetLabelOffset(0.000, "X"); // D=0.005
  h->SetLabelOffset(0.005, "Y"); // D=0.005
  h->SetLabelSize(0.055, "X");
  h->SetLabelSize(0.055, "Y");
  h->SetTitleOffset(0.8, "X");
  h->SetTitleOffset(0.9, "Y");
  h->SetTitleSize(0.06, "X");
  h->SetTitleSize(0.06, "Y");
}

void xsecPlots(const char* filelist, double fitVal = 1.00, double sigXS = 16.1, bool isSMPE = 1,
               bool isExp = 1) {

  if (isExp)
    fitVal = 1.0;

  style();

  TFile* fin[600];
  int nf = 0;

  ifstream stream1(filelist);
  if (!stream1)
    cout << "While opening a file an error is encountered" << endl;
  else
    cout << "File is successfully opened" << endl;

  char fname[256];
  int nf = 0;

  while (stream1 >> fname) {
    cout << fname << endl;
    fin[nf++] = new TFile(fname);
  }

  TH1D* hf = (TH1D*)fin[0]->Get("Fitted Xsec");

  for (int i = 1; i < nf; i++) {
    hf->Add((TH1D*)fin[i]->Get("Fitted Xsec"));
  }
  printf("RAW mean: %f, Entries: %d\n", hf->GetMean(), hf->GetEntries());
  double obs = 0;
  for (int i = 0; i <= hf->GetNbinsX(); i++) {
    if (hf->GetBinCenter(i) <= fitVal)
      continue;
    else {
      obs = hf->Integral(i, hf->GetNbinsX() + 1);
      break;
    }
  }

  TH1D* ho = new TH1D("out", "out", hf->GetNbinsX(), 0, 8.0 * sigXS);
  h1_style(ho);

  ho->GetXaxis()->SetTitle("Fitted Cross Section (pb)");
  ho->GetXaxis()->SetRangeUser(0.0, 59.9);
  ho->GetYaxis()->SetTitle("# Pseudo-Experiments");

  for (int i = 1; i <= hf->GetNbinsX(); i++) {
    ho->SetBinContent(i, hf->GetBinContent(i));
  }

  printf("Mean value: %f, TotEvts: %.0f\n", ho->GetMean(), ho->Integral());

  TF1* gfit = new TF1("mgaus", "gaus", 0.5 * sigXS, 1.5 * sigXS);
  hoClone = (TH1D*)ho->Clone("hoClone");
  hoClone->Fit("mgaus", "r");

  char name[256];
  if (isSMPE)
    sprintf(name, "#splitline{SM Signal Pseudo-Experiments}{   Mean: %.1f pb, RMS: %.2f pb}",
            gfit->GetParameter(1), gfit->GetParameter(2));
  else
    sprintf(name, "#splitline{Zero Signal Pseudo-Experiments}{    Mean: %.1f pb, RMS: %.2f pb}",
            ho->GetMean(), ho->GetRMS());
  TLatex* tex3 = new TLatex(20, 0.9 * ho->GetMaximum(), name);
  tex3->SetTextSize(0.05);
  tex3->SetTextColor(1);

  sprintf(name, "Mean: %.1f pb, RMS: %.2f pb", ho->GetMean(), ho->GetRMS());
  TLatex* tex1 = new TLatex(25, 0.8 * ho->GetMaximum(), name);
  tex1->SetTextSize(0.05);
  tex1->SetTextColor(1);

  if (isExp)
    sprintf(name, "#splitline{%d entries above}{ SM cross section}", obs);
  else
    sprintf(name, "#splitline{   %d entries above}{observed cross section}", obs);
  TLatex* tex6 = new TLatex(35, 0.6 * ho->GetMaximum(), name);
  tex6->SetTextSize(0.05);
  tex6->SetTextColor(2);

  double pVal = 1.0 * obs / ho->Integral();
  double nsig = gausSignif(pVal);

  if (!isExp)
    sprintf(name, "#splitline{Observed p-Value: %f}{                 N-#sigma: %.2f}", pVal, nsig);
  else
    sprintf(name, "#splitline{Expected p-Value: %f}{                 N-#sigma: %.2f}", pVal, nsig);
  TLatex* tex4 = new TLatex(25, 0.4 * ho->GetMaximum(), name);
  tex4->SetTextSize(0.05);
  tex4->SetTextColor(2);

  sprintf(name, "N-#sigma: %.2f", nsig);
  TLatex* tex5 = new TLatex(40, 0.3 * ho->GetMaximum(), name);
  tex5->SetTextSize(0.05);
  tex5->SetTextColor(2);

  double xSM[2];
  double ySM[2];
  xSM[0] = fitVal * sigXS;
  xSM[1] = fitVal * sigXS;
  ySM[0] = 0;
  if (isSMPE)
    ySM[1] = 0.75 * ho->GetMaximum();
  else
    ySM[1] = 0.01 * ho->GetMaximum();

  TGraph* smLine = new TGraph(2, xSM, ySM);
  TGraph* obsLine = new TGraph(2, xSM, ySM);
  smLine->SetLineWidth(4);
  smLine->SetLineColor(2);

  ho->GetYaxis()->SetTitleOffset(1.1);
  TCanvas* c1 = new TCanvas("a", "b");
  ho->SetLineWidth(2);
  ho->SetLineColor(23);
  ho->SetFillColor(23);
  ho->Draw();
  //  tex1->Draw();
  //    tex2->Draw();
  tex3->Draw();
  tex4->Draw();
  //  tex5->Draw();
  tex6->Draw();
  //  tex7->Draw();
  smLine->Draw("L");
}
