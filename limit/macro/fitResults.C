#include <stdio.h>
#include <iostream>

/**

This ROOT macro is intended to display the results of the testFit.exe
program.  To run it, use the following procedure:

> root -l
root [0] .L fitResults.c
root [1] makePlots("yourTestFitFile.root")

The histograms that you are displayed, in the order they are
layered on the screen, are as follows:

1) Presentation: Several histograms and canvases containing the
    response of the chisquared function to the fit parameters
    specified.  The reference used is the S+B model best fit to data.
    Each systematic is varied from -10 to +10 sigma while all others
    are held at their best fit values.

    Interpretation: The response of the chisquared function is an
    indication of the nature of the constraint on the specified
    systematic.

    ==>A well-modeled systematic will return a 2nd-order
    polynomial, symmetric about its preferred minimum value.

    ==>Systematics which give a near-linear response with no clear
    minimum are underconstrained and have little impact on the fit.

    ==>Response functions which "truncate" at low fluctuations are
    zeroing the backgrounds in question.  The point of truncation
    indicates the valid response region for the specified size of the
    systematic.

    ==>By performing a fit to a 2nd-order polynomial near the minimum
    (if it exists), one can estimate the constrainable size of the
    systematic.  The variance in the region near the minimum is then
    given by 1/sqrt(C) for a fit to F(x) = A + Bx + Cx^2.  This is the
    value that corresponds to increasing the chi2 function by one
    unit. This variance is in units of the size of the systematic
    specified.  EG, if your luminosity systematic is 6.1% and
    1/sqrt(C) = 1.52, the effective constraint is 9.3%.  If the
    constraint is smaller than 1.0, your systematic may be
    overestimated.

2) Presentation: Several histograms and canvases containing the pull
    function for each systematic used in the fit.  Each distribution
    should be approximately Gaussian and have a width of 1.0.  The
    values are referenced to data and are shown for both S+B and
    B-only fits. The mean value reflects the nominal "best fit" value
    for that systematic.

    Interpretation: The pull function is a second estimate of the
    nature of the constraint on the systematic in question.  The
    previous figures test the one dimensional response, while these
    values reflect the response in N-dimensional space (N is the
    number of total systematics).

    ==>Pull functions with narrow width (<1.0) are over-constrained.
    This may be an indication that your systematic is over-estimated.

    ==>Pull functions with large width (>1.0) are under-constrained.
    This may be an indication that your systematic is under-estimated.
    Caveat: Small backgrounds generally have small contributions to
    the chisquared function and will be under-constrained by this
    measure.  This doesn't necessarily mean your systematic for this
    background is wrong.

3) Presentation: Two histograms containing the S/B values and
    LOG_10(1+S/B) values given by the nominal (non-fitted) signal and
    background distributions.

    Interpretation: If you are using the CLfit class, you must choose
    a constraint cutoff to define the background sidebands.  The
    LOG_10(1+S/B) distribution is what's used to define fitting regions.

4) Presentation: One histogram containing the S+B (red) and B-Only
    (black) central values chosen by fits to data in the S+B and
    B-Only hypotheses, respectively.

    Interpretation: These are the central values of the fitted
    systematic parameters which minimize the chisquared function.

5) Presentation: One histogram containing the S+B (red) and B-Only
    (black) chisquared distributions for fits to data in the S+B and
    B-Only hypotheses, respectively.

    Interpretation: These distributions should follow standard
    chisquared distributions for X degrees of freedom, where X is
    determined by the number of bins fit and the number of fit
    parameters.

6)  Presentation: Two canvases containing the distributions of signal
    and background systematic variations around the nominal values for
    each bin.  Each canvas contains a distribution for the non-fitted
    values (black), the values from the S+B fit (red), and the B-Only fit
    (green).

    Interpretation: These distributions indicate the total effective
    size of the signal and background systematic uncertainty before
    and after fitting.

7) Presentation: Three canvases containing the signal, background, and
    data distributions before fitting, after the S+B fit, and after the
    B-only fit.  The fourth canvas compares the background distributions
    before and after the fitting.

   Interpretation: These distributions demonstrate how the shape of
   the background changes under each "best fit" scenario.

8) Presentation: Three canvases containing the number of Gaussian
    sigma per bin based on a Data/MC comparison.  The comparison is
    performed for the baseline MC prediction, the S+B fit to data, and
    the Bkgd-Only Fit to data.  The total Chi2/ndof is also displayed
    on the figure.

9)

**/
static const int mFont = 132;

void style() {
  gStyle->SetLabelFont(mFont, "X");
  gStyle->SetLabelFont(mFont, "Y");
  gStyle->SetTitleFont(mFont, "X");
  gStyle->SetTitleFont(mFont, "Y");
  gStyle->SetTitleFont(mFont);

  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleBorderSize(0);

  gStyle->SetTitleColor(1); // 0
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.2);

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
void h1_style(TH1* h, int line_width = 3, int line_color = 1, int line_style = 1,
              int fill_style = 1001, int fill_color = 50, double y_min = -1111.,
              double y_max = -1111., int ndivx = 510, int ndivy = 510, int marker_style = 20,
              int marker_color = 1, double marker_size = 1.3, int optstat = 0) {

  h->GetXaxis()->SetTitleFont(mFont);
  h->GetYaxis()->SetTitleFont(mFont);
  h->SetStats(kFALSE);

  h->SetLabelFont(mFont, "X");   // 42
  h->SetLabelFont(mFont, "Y");   // 42
  h->SetLabelOffset(0.000, "X"); // D=0.005
  h->SetLabelOffset(0.005, "Y"); // D=0.005
  h->SetLabelSize(0.06, "X");
  h->SetLabelSize(0.06, "Y");
  h->SetTitleOffset(1.0, "X");
  h->SetTitleOffset(1.0, "Y");
  h->SetTitleSize(0.065, "X");
  h->SetTitleSize(0.065, "Y");
  //  h->SetTitle(0);
}

void axis_style(TGaxis* a) {
  a->SetTitleFont(mFont);
  a->SetLabelFont(mFont);   // 42
  a->SetLabelOffset(0.006); // D=0.005
  a->SetLabelSize(0.06);
  a->SetTitleSize(0.065);
  a->SetTitleOffset(1.2);
}

void h2_style(TH2* h) {

  h->SetLabelFont(mFont, "X");   // 42
  h->SetLabelFont(mFont, "Y");   // 42
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
void leg_style(TLegend* h) {
  h->SetTextFont(mFont);
  h->SetFillColor(0);
  h->SetFillStyle(0);
  h->SetLineColor(0);
  h->SetBorderSize(0);
}

void can_style(TCanvas* c1) {
  c1->SetLeftMargin(0.14);
  c1->SetBottomMargin(0.14);
  c1->SetTopMargin(0.06);
  c1->SetRightMargin(0.05);
}

void removeMin(TH1* hist) {
  double min = hist->GetMinimum();
  for (int i = 0; i < hist->GetNbinsX(); i++)
    hist->SetBinContent(i, hist->GetBinContent(i) - min);
  return;
}

void makePlots(const char* testFile, bool doPrint = false) {
  style();

  TFile* f = new TFile(testFile);

  TH1F* dataDist = (TH1F*)f->Get("Total Data");
  TH1F* bkgdDist = (TH1F*)f->Get("Total Background, Prefit");
  TH1F* bkgdDistSB = (TH1F*)f->Get("Total Background, S+B Fit");
  TH1F* bkgdDistB = (TH1F*)f->Get("Total Background, B-Only Fit");

  TH1F* sigDist = (TH1F*)f->Get("Total Signal, Prefit");
  TH1F* sigDistSB = (TH1F*)f->Get("Total Signal, S+B Fit");
  TH1F* sigDistB = (TH1F*)f->Get("Total Signal, B-Only Fit");
  TH1F* sbDist = (TH1F*)sigDist->Clone("sbdist");
  sbDist->Add(bkgdDist);

  TH1F* dataCDF = (TH1F*)f->Get("Data CDF");
  TH1F* bkgdCDF = (TH1F*)f->Get("Background CDF, Prefit");
  TH1F* bkgdCDFSB = (TH1F*)f->Get("Background CDF, S+B Fit");
  TH1F* bkgdCDFB = (TH1F*)f->Get("Background CDF, B-Only Fit");

  TH1F* sigCDF = (TH1F*)f->Get("Signal CDF, Prefit");
  TH1F* sigCDFSB = (TH1F*)f->Get("Signal CDF, S+B Fit");
  TH1F* sigCDFB = (TH1F*)f->Get("Signal CDF, B-Only Fit");

  TH1F* dataChi2SB = (TH1F*)f->Get("Data Chi-Square S+B");
  TH1F* dataChi2B = (TH1F*)f->Get("Data Chi-Square B-Only");
  TH1F* dataChi2SB_Ref = (TH1F*)f->Get("Data Chi-Square S+B, Ref");
  TH1F* dataChi2B_Ref = (TH1F*)f->Get("Data Chi-Square B-Only, Ref");

  TH1F* peChi2SB = (TH1F*)f->Get("PE Chi-Square S+B");
  TH1F* peChi2B = (TH1F*)f->Get("PE Chi-Square B-Only");

  TH1F* sigSyst = (TH1F*)f->Get("Signal Systematic");
  TH1F* bkgSyst = (TH1F*)f->Get("Bkgd Systematic");

  TH1F* fitsigS = (TH1F*)f->Get("Signal Systematic, S+B Fit");
  TH1F* fitbkgS = (TH1F*)f->Get("Bkgd Systematic, S+B Fit");

  TH1F* fitsigB = (TH1F*)f->Get("Signal Systematic, B-Only Fit");
  TH1F* fitbkgB = (TH1F*)f->Get("Bkgd Systematic, B-Only Fit");

  TH1F* fitsigS_PE = (TH1F*)f->Get("Signal Systematic, S+B PE Fit");
  TH1F* fitbkgS_PE = (TH1F*)f->Get("Bkgd Systematic, S+B PE Fit");

  TH1F* fitsigB_PE = (TH1F*)f->Get("Signal Systematic, B-Only PE Fit");
  TH1F* fitbkgB_PE = (TH1F*)f->Get("Bkgd Systematic, B-Only PE Fit");

  TH1F* fitParamsSB = (TH1F*)f->Get("S+B Fit Params");
  TH1F* fitParamsB = (TH1F*)f->Get("B-Only Fit Params");

  TH1F* minuitIterSB = (TH1F*)f->Get("S+B Fit Iterations");
  TH1F* minuitIterB = (TH1F*)f->Get("B-Only Fit Iterations");

  TH1F* minuitStatusSB = (TH1F*)f->Get("S+B Fit Status");
  TH1F* minuitStatusB = (TH1F*)f->Get("B-Only Fit Status");

  TH1F* minuitEmatSB = (TH1F*)f->Get("S+B Fit Error Matrix");
  TH1F* minuitEmatB = (TH1F*)f->Get("B-Only Fit Error Matrix");

  TH1F* fitErrsSB = (TH1F*)f->Get("S+B Fit Errors");
  TH1F* fitErrsB = (TH1F*)f->Get("B-Only Fit Errors");
  TH1F* fitErrsM = (TH1F*)f->Get("MINUIT Fit Errors");

  TH1F* nMinusOneChanB = (TH1F*)f->Get("B-Only N-1 Test, Channels");
  TH1F* nMinusOneSystB = (TH1F*)f->Get("B-Only N-1 Test, Systematics");
  TH1F* nMinusOneBinsB = (TH1F*)f->Get("B-Only N-1 Test, Bins");

  TH1F* nMinusOneChanSB = (TH1F*)f->Get("S+B N-1 Test, Channels");
  TH1F* nMinusOneSystSB = (TH1F*)f->Get("S+B N-1 Test, Systematics");
  TH1F* nMinusOneBinsSB = (TH1F*)f->Get("S+B N-1 Test, Bins");

  TH1F* systFinalSB = (TH1F*)f->Get("Total S+B Systematics per Bin");
  TH1F* systFinalB = (TH1F*)f->Get("Total B-Only Systematics per Bin");
  TH1F* systFinal = (TH1F*)f->Get("Total Non-Fit Systematics per Bin");

  TH1F* systDChi2_1 = (TH1F*)f->Get("Systematics DeltaChi2 Test - 1");
  TH1F* systDChi2_2 = (TH1F*)f->Get("Systematics DeltaChi2 Test - 2");
  TH1F* sourceDChi2_1 = (TH1F*)f->Get("Event Source DeltaChi2 Test - 1");
  TH1F* sourceDChi2_2 = (TH1F*)f->Get("Event Source DeltaChi2 Test - 2");

  TH1F* sOverB = new TH1F("S/B Values", "S/B Values", sigDist->GetNbinsX(),
                          sigDist->GetXaxis()->GetXmin(), sigDist->GetXaxis()->GetXmax());

  TH1F* LOGsOverB =
      new TH1F("S#timesLog(1+S/B) Values", "S#timesLog(1+S/B) Values", sigDist->GetNbinsX(),
               sigDist->GetXaxis()->GetXmin(), sigDist->GetXaxis()->GetXmax());

  TCanvas ccc("fitResults", "fitResults");
  char buffer[200];
  sprintf(buffer, "fitResults.ps[");
  if (doPrint)
    ccc.Print(buffer); //open ps file

  double sb_max = -1;
  double sb_min = 100;
  double val = 0;
  for (int b = 1; b <= sigDist->GetNbinsX(); b++) {
    if (bkgdDist->GetBinContent(b) > 0) {
      val = sigDist->GetBinContent(b) / bkgdDist->GetBinContent(b);
      sOverB->SetBinContent(b, val);
      val = sigDist->GetBinContent(b) * log(1.0 + val);
      LOGsOverB->SetBinContent(b, val);
      if (val > sb_max && sigDist->GetBinContent(b) > 1e-10)
        sb_max = val;
      if (val < sb_min && sigDist->GetBinContent(b) > 1e-10)
        sb_min = val;
    }
  }
  if (sb_max > 0)
    sb_max *= 1.01;
  else
    sb_max *= 0.99;
  sb_min *= 1.01;
  printf("SBMINMAX: %f, %f\n", sb_min, sb_max);
  int nbins = (int)(sigDist->GetNbinsX() * 0.8);
  TH1F* LOGsOverB_dat =
      new TH1F("S#timesLog(1+S/B): Data", "S#timesLog(1+S/B): Data", nbins, sb_min, sb_max);
  TH1F* LOGsOverB_bkg =
      new TH1F("S#timesLog(1+S/B): Bkgd", "S#timesLog(1+S/B): Bkgd", nbins, sb_min, sb_max);
  TH1F* LOGsOverB_fit =
      new TH1F("S#timesLog(1+S/B): Fit", "S#timesLog(1+S/B): Fit", nbins, sb_min, sb_max);

  double val = 0;
  for (int b = 1; b <= sigDist->GetNbinsX(); b++) {
    if (bkgdDist->GetBinContent(b) > 0) {
      val = sigDist->GetBinContent(b) *
            log(1 + sigDist->GetBinContent(b) / bkgdDist->GetBinContent(b));

      LOGsOverB_dat->Fill(val, dataDist->GetBinContent(b));
      LOGsOverB_bkg->Fill(val, bkgdDist->GetBinContent(b));
      LOGsOverB_fit->Fill(val, bkgdDistB->GetBinContent(b));
      //      LOGsOverB_fit->AddBinContent(val,sigDist->GetBinContent(b));
    }
  }

  LOGsOverB_bkg->SetLineWidth(4);
  LOGsOverB_fit->SetLineWidth(4);
  LOGsOverB_bkg->SetLineColor(2);
  LOGsOverB_fit->SetLineColor(4);
  LOGsOverB_dat->SetMarkerStyle(21);
  LOGsOverB_dat->SetMarkerSize(1);

  sOverB->SetLineWidth(4);
  LOGsOverB->SetLineWidth(4);

  dataDist->SetLineWidth(4);
  dataDist->SetMarkerStyle(21);
  dataDist->SetMarkerSize(1);
  dataCDF->SetLineWidth(4);
  dataCDF->SetMarkerStyle(21);
  dataCDF->SetMarkerSize(1);

  bkgdDist->SetLineWidth(4);
  bkgdDistSB->SetLineWidth(4);
  bkgdDistB->SetLineWidth(4);
  bkgdCDF->SetLineWidth(4);
  bkgdCDFSB->SetLineWidth(4);
  bkgdCDFB->SetLineWidth(4);

  bkgdDist->SetLineColor(2);
  bkgdCDF->SetLineColor(2);
  bkgdDistSB->SetLineColor(4);
  bkgdCDFSB->SetLineColor(4);
  bkgdDistB->SetLineColor(8);
  bkgdCDFB->SetLineColor(8);

  sigDist->SetLineWidth(4);
  sigDistSB->SetLineWidth(4);
  sigDistB->SetLineWidth(4);
  sigCDF->SetLineWidth(4);
  sigCDFSB->SetLineWidth(4);
  sigCDFB->SetLineWidth(4);

  sigDist->SetLineColor(14);
  sigCDF->SetLineColor(14);
  sigDistSB->SetLineColor(15);
  sigCDFSB->SetLineColor(15);
  sigDistB->SetLineColor(16);
  sigCDFB->SetLineColor(16);

  dataChi2SB->SetLineWidth(4);
  dataChi2B->SetLineWidth(4);
  dataChi2SB->SetLineColor(2);

  dataChi2SB_Ref->SetLineWidth(2);
  dataChi2B_Ref->SetLineWidth(2);
  dataChi2SB_Ref->SetLineColor(4);
  dataChi2B_Ref->SetLineColor(3);

  peChi2SB->SetLineWidth(4);
  peChi2B->SetLineWidth(4);
  peChi2SB->SetLineColor(2);

  fitParamsSB->SetLineWidth(4);
  fitParamsB->SetLineWidth(4);
  fitParamsSB->SetLineColor(2);

  if (nMinusOneChanB) {
    nMinusOneChanB->SetLineWidth(4);
    nMinusOneBinsB->SetLineWidth(4);
    nMinusOneSystB->SetLineWidth(4);

    nMinusOneChanSB->SetLineWidth(4);
    nMinusOneBinsSB->SetLineWidth(4);
    nMinusOneSystSB->SetLineWidth(4);

    nMinusOneChanSB->SetLineColor(2);
    nMinusOneBinsSB->SetLineColor(2);
    nMinusOneSystSB->SetLineColor(2);
  }
  sigSyst->SetLineWidth(4);
  bkgSyst->SetLineWidth(4);

  systFinal->SetLineWidth(4);
  systFinalSB->SetLineWidth(4);
  systFinalB->SetLineWidth(4);
  systFinalSB->SetLineColor(2);
  systFinalB->SetLineColor(3);
  systFinal->SetLineColor(4);

  fitsigS->SetLineWidth(4);
  fitbkgS->SetLineWidth(4);
  fitsigB->SetLineWidth(4);
  fitbkgB->SetLineWidth(4);
  fitsigS->SetLineColor(2);
  fitbkgS->SetLineColor(2);
  fitsigB->SetLineColor(3);
  fitbkgB->SetLineColor(3);

  fitsigS_PE->SetLineWidth(4);
  fitbkgS_PE->SetLineWidth(4);
  fitsigB_PE->SetLineWidth(4);
  fitbkgB_PE->SetLineWidth(4);
  fitsigS_PE->SetLineColor(7);
  fitbkgS_PE->SetLineColor(7);
  fitsigB_PE->SetLineColor(6);
  fitbkgB_PE->SetLineColor(6);

  minuitIterSB->SetLineWidth(4);
  minuitIterB->SetLineWidth(4);
  minuitIterB->SetLineColor(2);

  minuitStatusSB->SetLineWidth(4);
  minuitStatusB->SetLineWidth(4);
  minuitStatusB->SetLineColor(2);

  fitErrsSB->SetLineWidth(4);
  fitErrsSB->SetLineColor(2);
  fitErrsB->SetLineWidth(4);
  fitErrsB->SetLineColor(4);
  fitErrsM->SetLineWidth(4);
  fitErrsM->SetLineColor(6);

  printf("Before fit========>  Sig: %.4f, Bkgd: %.4f, Data: %.1f\n", sigDist->Integral(),
         bkgdDist->Integral(), dataDist->Integral());
  printf("After S+B fit=====>  Sig: %.4f, Bkgd: %.4f, Data: %.1f\n", sigDistSB->Integral(),
         bkgdDistSB->Integral(), dataDist->Integral());
  printf("After B-Only fit==>  Sig: %.4f, Bkgd: %.4f, Data: %.1f\n", sigDistB->Integral(),
         bkgdDistB->Integral(), dataDist->Integral());

  TH1F* nomChi2 = (TH1F*)bkgdDist->Clone("nomSig");
  nomChi2->SetTitle("# Gaussian Sigma per Bin, Pre-Fit");
  TH1F* bFitChi2 = (TH1F*)bkgdDist->Clone("bSig");
  bFitChi2->SetTitle("# Gaussian Sigma per Bin, B-Only Fit");
  TH1F* sbFitChi2 = (TH1F*)bkgdDist->Clone("sbSig");
  sbFitChi2->SetTitle("# Gaussian Sigma per Bin, S+B Fit");

  nomChi2->Scale(0);
  bFitChi2->Scale(0);
  sbFitChi2->Scale(0);
  double a1 = 0;
  double a2 = 0;
  double a3 = 0;
  TH1F* bkgdDistSB2 = (TH1F*)bkgdDistSB->Clone("tmpSB");
  bkgdDistSB2->Add(sigDistSB);
  int rBins = 0;
  for (int b = 1; b <= nomChi2->GetNbinsX(); b++) {
    double d = dataDist->GetBinContent(b);
    if (d > 0) {
      rBins++;
      a1 += fabs(bkgdDist->GetBinContent(b) - d) * fabs(bkgdDist->GetBinContent(b) - d) / d;
      a2 += fabs(bkgdDistB->GetBinContent(b) - d) * fabs(bkgdDistB->GetBinContent(b) - d) / d;
      a3 += fabs(bkgdDistSB2->GetBinContent(b) - d) * fabs(bkgdDistSB2->GetBinContent(b) - d) / d;
      nomChi2->SetBinContent(b, fabs(bkgdDist->GetBinContent(b) - d) / sqrt(d));
      bFitChi2->SetBinContent(b, fabs(bkgdDistB->GetBinContent(b) - d) / sqrt(d));
      sbFitChi2->SetBinContent(b, fabs(bkgdDistSB2->GetBinContent(b) - d) / sqrt(d));
    }
  }
  a1 /= (rBins - 1);
  a2 /= (rBins - 1);
  a3 /= (rBins - 1);

  if (nomChi2) {
    TCanvas* cA = new TCanvas("Chi2 Nominal", "Chi2 Nominal");
    can_style(cA);
    h1_style(nomChi2);
    nomChi2->GetXaxis()->SetTitle("Bin Index");
    nomChi2->GetYaxis()->SetTitle("#chi^{2}");
    nomChi2->Draw();
    char dest[256];
    sprintf(dest, "#splitline{Pre-Fit}{#chi2/ndof = %.3f}", a1);

    nomChi2->GetYaxis()->SetRangeUser(0, nomChi2->GetMaximum() * 1.25);
    double xPos = nomChi2->GetXaxis()->GetXmin() +
                  (nomChi2->GetXaxis()->GetXmax() - nomChi2->GetXaxis()->GetXmin()) * 0.5;

    TLatex* tex0 = new TLatex(xPos, nomChi2->GetMaximum() * 0.9, dest);
    tex0->SetTextFont(mFont);
    tex0->Draw("same");
    if (doPrint)
      cA->Print(buffer);
  }

  if (bFitChi2) {
    TCanvas* cB = new TCanvas("Chi2 B-Only Fit", "Chi2 B-Only Fit");
    can_style(cB);
    h1_style(bFitChi2);
    bFitChi2->GetXaxis()->SetTitle("Bin Index");
    bFitChi2->GetYaxis()->SetTitle("#chi^{2}");
    bFitChi2->Draw();
    sprintf(dest, "#splitline{Bkgd-Only Fit}{#chi2/ndof = %.3f}", a2);
    bFitChi2->GetYaxis()->SetRangeUser(0, bFitChi2->GetMaximum() * 1.25);
    TLatex* tex0 = new TLatex(xPos, bFitChi2->GetMaximum() * 0.9, dest);
    tex0->SetTextFont(mFont);
    tex0->Draw("same");
    if (doPrint)
      cB->Print(buffer);
  }

  if (sbFitChi2) {
    TCanvas* cC = new TCanvas("Chi2 S+B Fit", "Chi2 S+B Fit");
    can_style(cC);
    h1_style(sbFitChi2);
    sbFitChi2->GetXaxis()->SetTitle("Bin Index");
    sbFitChi2->GetYaxis()->SetTitle("#chi^{2}");
    sbFitChi2->Draw();
    sprintf(dest, "#splitline{S+B Fit}{#chi2/ndof = %.3f}", a3);
    sbFitChi2->GetYaxis()->SetRangeUser(0, sbFitChi2->GetMaximum() * 1.25);
    TLatex* tex0 = new TLatex(xPos, sbFitChi2->GetMaximum() * 0.9, dest);
    tex0->SetTextFont(mFont);
    tex0->Draw("same");
    if (doPrint)
      cC->Print(buffer);
  }

  double max = bkgdDist->GetMaximum();
  if (dataDist->GetMaximum() > max)
    max = dataDist->GetMaximum();
  TCanvas* c1 = new TCanvas("Pre-Fit Dists", "Pre-Fit Dists");
  can_style(c1);
  h1_style(bkgdDist);
  bkgdDist->SetMaximum(max * 1.25);
  bkgdDist->GetXaxis()->SetTitle("Final Variable");
  bkgdDist->GetYaxis()->SetTitle("Evts / Bin");
  bkgdDist->SetTitle("Data/MC Comparison, Pre-Fit");
  bkgdDist->Draw("hist");
  dataDist->Draw("PEsame");
  TH1F* sigDist2 = sigDist->Clone("tmpSig");
  sigDist2->Add(bkgdDist);
  sigDist2->Draw("samehist");
  TLegend* leg1 = new TLegend(0.75, 0.73, 0.95, 0.94);
  leg_style(leg1);
  leg1->AddEntry(dataDist, "Data");
  leg1->AddEntry(bkgdDist, "Bkgd-Only");
  leg1->AddEntry(sigDist2, "Signal+Bkgd");
  leg1->Draw("same");
  leg1->SetFillColor(0);
  if (doPrint)
    c1->Print("nomBkgdData.eps");
  if (doPrint)
    c1->Print(buffer);

  max = bkgdDistB->GetMaximum();
  if (dataDist->GetMaximum() > max)
    max = dataDist->GetMaximum();
  TCanvas* c2 = new TCanvas("B-Only Fit Dists", "B-Only Fit Dists");
  can_style(c2);
  h1_style(bkgdDistB);
  bkgdDistB->SetMaximum(max * 1.25);
  bkgdDistB->GetXaxis()->SetTitle("Final Variable");
  bkgdDistB->GetYaxis()->SetTitle("Evts / Bin");
  bkgdDistB->SetTitle("Data/MC Comparison, B-Only Fit");
  bkgdDistB->Draw("hist");
  dataDist->Draw("PEsame");
  TLegend* leg2 = new TLegend(0.75, 0.73, 0.95, 0.94);
  leg_style(leg2);
  leg2->AddEntry(dataDist, "Data");
  leg2->AddEntry(bkgdDistB, "Bkgd-Only");
  //  leg2->AddEntry(sigDistB,"Signal+Bkgd");
  leg2->SetFillColor(0);
  leg2->Draw("same");
  if (doPrint)
    c2->Print(buffer);

  sigDistSB->Add(bkgdDistSB);
  max = sigDistSB->GetMaximum();
  if (dataDist->GetMaximum() > max)
    max = dataDist->GetMaximum();
  TCanvas* c3 = new TCanvas("S+B Fit Dists", "S+B Fit Dists");
  can_style(c3);
  h1_style(sigDistSB);
  sigDistSB->SetMaximum(max * 1.25);
  sigDistSB->GetXaxis()->SetTitle("Final Variable");
  sigDistSB->GetYaxis()->SetTitle("Evts / Bin");
  sigDistSB->SetTitle("Data/MC Comparison, S+B Fit");
  sigDistSB->Draw("hist");
  dataDist->Draw("PEsame");
  TLegend* leg3 = new TLegend(0.75, 0.73, 0.95, 0.94);
  leg_style(leg3);
  leg3->AddEntry(dataDist, "Data");
  leg3->AddEntry(sigDistSB, "Signal+Bkgd");
  leg3->SetFillColor(0);
  leg3->Draw("same");
  if (doPrint)
    c3->Print(buffer);

  TCanvas* cc1 = new TCanvas("CDF Prefit", "CDF Prefit");
  can_style(cc1);
  h1_style(bkgdCDF);
  bkgdCDF->GetXaxis()->SetTitle("Final Variable");
  bkgdCDF->GetYaxis()->SetTitle("Evts / Bin");
  bkgdCDF->SetTitle("Data/MC CDF Comparison, Pre-Fit");
  bkgdCDF->Draw("hist");
  dataCDF->Draw("PEsame");
  sigCDF->Add(bkgdCDF);
  sigCDF->Draw("samehist");
  TLegend* cleg1 = new TLegend(0.75, 0.2, 0.95, 0.41);
  leg_style(cleg1);
  cleg1->AddEntry(dataCDF, "Data");
  cleg1->AddEntry(bkgdCDF, "Bkgd-Only");
  cleg1->AddEntry(sigCDF, "Signal+Bkgd");
  cleg1->Draw("same");
  cleg1->SetFillColor(0);

  TCanvas* cc2 = new TCanvas("CDF B-Only Fit", "CDF B-Only Fit");
  can_style(cc2);
  h1_style(bkgdCDFB);
  bkgdCDFB->GetXaxis()->SetTitle("Final Variable");
  bkgdCDFB->SetTitle("Data/MC CDF Comparison, B-Only Fit");
  bkgdCDFB->GetYaxis()->SetTitle("Evts / Bin");
  bkgdCDFB->Draw("hist");
  dataCDF->Draw("PEsame");
  sigCDFB->Add(bkgdCDFB);
  sigCDFB->Draw("samehist");
  TLegend* cleg2 = new TLegend(0.75, 0.2, 0.95, 0.41);
  leg_style(cleg2);
  cleg2->AddEntry(dataCDF, "Data");
  cleg2->AddEntry(bkgdCDFB, "Bkgd-Only");
  cleg2->AddEntry(sigCDFB, "Signal+Bkgd");
  cleg2->SetFillColor(0);
  cleg2->Draw("same");

  TCanvas* cc3 = new TCanvas("CDF S+B Fit", "CDF S+B Fit");
  can_style(cc3);
  h1_style(bkgdCDFSB);
  bkgdCDFSB->GetXaxis()->SetTitle("Final Variable");
  bkgdCDFSB->SetTitle("Data/MC CDF Comparison, S+B Fit");
  bkgdCDFSB->GetYaxis()->SetTitle("Evts / Bin");
  bkgdCDFSB->Draw("hist");
  dataCDF->Draw("PEsame");
  sigCDFSB->Add(bkgdCDFSB);
  sigCDFSB->Draw("samehist");
  TLegend* cleg3 = new TLegend(0.75, 0.2, 0.95, 0.41);
  leg_style(cleg3);
  cleg3->AddEntry(dataCDF, "Data");
  cleg3->AddEntry(bkgdCDFSB, "Bkgd-Only");
  cleg3->AddEntry(sigCDFSB, "Signal+Bkgd");
  cleg3->SetFillColor(0);
  cleg3->Draw("same");

  TCanvas* c4 = new TCanvas("FV Comparison", "FV Comparison");
  can_style(c4);
  h1_style(bkgdDist);
  bkgdDist->GetXaxis()->SetTitle("Final Variable");
  bkgdDist->SetTitle("Data/MC Comparison, Pre-Fit");
  bkgdDist->SetMaximum(1.25 * bkgdDist->GetMaximum());
  bkgdDist->Draw("hist");
  bkgdDistB->Draw("samehist");
  bkgdDistSB->Draw("samehist");
  TLegend* leg4 = new TLegend(0.72, 0.67, 0.94, 0.93);
  leg_style(leg4);
  leg4->AddEntry(bkgdDist, "No Fit");
  leg4->AddEntry(bkgdDistB, "Bkgd-Only Fit");
  leg4->AddEntry(bkgdDistSB, "Signal+Bkgd Fit");
  leg4->Draw("same");
  leg4->SetFillColor(0);
  if (doPrint)
    c4->Print("fitBkgds.eps");
  if (doPrint)
    c4->Print(buffer);

  if (sigSyst) {
    double RMS1 = sigSyst->GetRMS();
    double RMS2 = fitsigS->GetRMS();
    double RMS3 = fitsigB->GetRMS();
    double Mean1 = sigSyst->GetMean();
    double Mean2 = fitsigS->GetMean();
    double Mean3 = fitsigB->GetMean();

    TCanvas* c5 = new TCanvas("Signal Syst Disp", "Signal Syst Disp");
    can_style(c5);
    h1_style(sigSyst);
    sigSyst->SetTitle("Signal Per-Bin Systematic Dispersion");
    sigSyst->GetXaxis()->SetTitle("Per-Bin Dispersion");
    sigSyst->GetYaxis()->SetTitle("Frequency");

    double max = sigSyst->GetMaximum();
    if (fitsigS->GetMaximum() > max)
      max = fitsigS->GetMaximum();
    if (fitsigB->GetMaximum() > max)
      max = fitsigB->GetMaximum();

    sigSyst->Scale((max / (sigSyst->GetMaximum() + 1e-5)) / (sigSyst->Integral() + 1e-5));
    fitsigS->Scale((max / (fitsigS->GetMaximum() + 1e-5)) / (fitsigS->Integral() + 1e-5));
    fitsigB->Scale((max / (fitsigB->GetMaximum() + 1e-5)) / (fitsigB->Integral() + 1e-5));

    sigSyst->Draw("hist");
    fitsigS->Draw("samehist");
    fitsigB->Draw("samehist");

    TLegend* leg6 = new TLegend(0.73, 0.67, 0.95, 0.93);
    leg_style(leg6);
    leg6->AddEntry(sigSyst, "No Fit");
    leg6->AddEntry(fitsigB, "Bkgd-Only Fit");
    leg6->AddEntry(fitsigS, "Signal+Bkgd Fit");
    leg6->SetFillColor(0);
    leg6->Draw("same");

    char title[256];
    sprintf(
        title,
        "#splitline{Nominal RMS:   %.3f}{#splitline{S+B Fit RMS:      %.3f}{B-Only Fit RMS: %.3f}}",
        RMS1, RMS2, RMS3);
    double px = sigSyst->GetXaxis()->GetXmin() +
                0.05 * (sigSyst->GetXaxis()->GetXmax() - sigSyst->GetXaxis()->GetXmin());
    TLatex* texaa = new TLatex(px, sigSyst->GetMaximum() * 0.8, title);
    texaa->SetTextFont(mFont);
    texaa->Draw("same");

    sprintf(
        title,
        "#splitline{Nominal Mean:   %.3f}{#splitline{S+B Fit Mean:      %.3f}{B-Only Fit Mean: %.3f}}",
        Mean1, Mean2, Mean3);
    TLatex* texaa = new TLatex(px, sigSyst->GetMaximum() * 0.5, title);
    texaa->SetTextFont(mFont);
    texaa->Draw("same");
    if (doPrint)
      c5->Print("signalSyst.eps");
    if (doPrint)
      c5->Print(buffer);

    RMS1 = bkgSyst->GetRMS();
    RMS2 = fitbkgS->GetRMS();
    RMS3 = fitbkgB->GetRMS();
    Mean1 = bkgSyst->GetMean();
    Mean2 = fitbkgS->GetMean();
    Mean3 = fitbkgB->GetMean();

    TCanvas* c6 = new TCanvas("six", "six");
    can_style(c6);
    h1_style(bkgSyst);
    bkgSyst->SetTitle("Background Per-Bin Systematic Dispersion");
    bkgSyst->GetXaxis()->SetTitle("Per-Bin Dispersion");
    bkgSyst->GetYaxis()->SetTitle("Frequency");

    max = bkgSyst->GetMaximum();
    if (fitbkgS->GetMaximum() > max)
      max = fitbkgS->GetMaximum();
    if (fitbkgB->GetMaximum() > max)
      max = fitbkgB->GetMaximum();

    bkgSyst->Scale((max / (bkgSyst->GetMaximum() + 1e-5)) / (bkgSyst->Integral() + 1e-5));
    fitbkgS->Scale((max / (fitbkgS->GetMaximum() + 1e-5)) / (fitbkgS->Integral() + 1e-5));
    fitbkgB->Scale((max / (fitbkgB->GetMaximum() + 1e-5)) / (fitbkgB->Integral() + 1e-5));

    bkgSyst->Draw("hist");
    fitbkgS->Draw("samehist");
    fitbkgB->Draw("samehist");
    leg6->Draw("same");
    sprintf(
        title,
        "#splitline{Nominal RMS:   %.3f}{#splitline{S+B Fit RMS:      %.3f}{B-Only Fit RMS: %.3f}}",
        RMS1, RMS2, RMS3);
    TLatex* texaa = new TLatex(px, bkgSyst->GetMaximum() * 0.8, title);
    texaa->SetTextFont(mFont);
    texaa->Draw("same");

    sprintf(
        title,
        "#splitline{Nominal Mean:   %.3f}{#splitline{S+B Fit Mean:      %.3f}{B-Only Fit Mean: %.3f}}",
        Mean1, Mean2, Mean3);
    TLatex* texaa = new TLatex(px, bkgSyst->GetMaximum() * 0.5, title);
    texaa->SetTextFont(mFont);
    texaa->Draw("same");
    if (doPrint)
      c6->Print("bkgdSyst.eps");
    if (doPrint)
      c6->Print(buffer);
  }

  if (systFinal) {
    TCanvas* c6a = new TCanvas("sixa", "sixa");
    can_style(c6a);
    h1_style(systFinal);
    systFinal->SetTitle("Systematic Uncertainties Per Bin (events)");
    systFinal->GetXaxis()->SetTitle("Bin Index");
    systFinal->GetYaxis()->SetTitle("Events / Bin");
    systFinal->SetMinimum(0);
    systFinal->Draw("hist");
    systFinalSB->Draw("samehist");
    systFinalB->Draw("samehist");
    TLegend* leg6a = new TLegend(0.73, 0.73, 0.94, 0.94);
    leg_style(leg6a);
    leg6a->AddEntry(systFinal, "No Fit");
    leg6a->AddEntry(systFinalB, "Bkgd-Only Fit");
    leg6a->AddEntry(systFinalSB, "Signal+Bkgd Fit");
    leg6a->SetFillColor(0);
    leg6a->Draw("same");
    if (doPrint)
      c6a->Print("bkgdSyst.eps");
    if (doPrint)
      c6a->Print(buffer);

    TH1F* systFinalSB2 = (TH1F*)systFinalSB->Clone("systSB2");
    TH1F* systFinalB2 = (TH1F*)systFinalB->Clone("systB2");
    TH1F* systFinal2 = (TH1F*)systFinal->Clone("syst2");
    systFinalSB2->Divide(sigDistSB);
    systFinalB2->Divide(bkgdDistB);
    systFinal2->Divide(sbDist);

    systFinalSB2->Scale(100.0);
    systFinalB2->Scale(100.0);
    systFinal2->Scale(100.0);

    TCanvas* c6b = new TCanvas("sixb", "sixb");
    can_style(c6b);
    h1_style(systFinal2);
    systFinal2->SetTitle("Systematic Uncertainties Per Bin (%)");
    systFinal2->GetXaxis()->SetTitle("Bin Index");
    systFinal2->GetYaxis()->SetTitle("Percent / Bin");
    systFinal2->GetYaxis()->SetTitleOffset(1.2);
    systFinal2->SetMinimum(0);
    systFinal2->Draw("hist");
    systFinalSB2->Draw("samehist");
    systFinalB2->Draw("samehist");
    TLegend* leg6b = new TLegend(0.44, 0.68, 0.68, 0.92);
    leg_style(leg6b);
    leg6b->AddEntry(systFinal, "No Fit");
    leg6b->AddEntry(systFinalB2, "Bkgd-Only Fit");
    leg6b->AddEntry(systFinalSB2, "Signal+Bkgd Fit");
    leg6b->SetFillColor(0);
    leg6b->Draw("same");
    if (doPrint)
      c6b->Print("bkgdSyst.eps");
    if (doPrint)
      c6b->Print(buffer);
  }

  if (dataChi2B) {
    TCanvas* c7 = new TCanvas("seven", "seven");
    can_style(c7);
    h1_style(dataChi2B);
    dataChi2B->GetXaxis()->SetTitle("Fit #chi^{2} Value");
    dataChi2B->GetYaxis()->SetTitle("Frequency");
    dataChi2B->SetTitle("Chi2 for Fits to Data");
    dataChi2B->Draw("hist");
    dataChi2SB->Draw("samehist");
    dataChi2B_Ref->Scale(1e10);
    dataChi2B_Ref->Draw("samehist");
    dataChi2SB_Ref->Scale(1e10);
    dataChi2SB_Ref->Draw("samehist");

    TLegend* leg7 = new TLegend(0.65, 0.65, 0.85, 0.85);
    leg7->AddEntry(dataChi2B, "Background-Only");
    leg7->AddEntry(dataChi2SB, "Signal+Bkgd");
    leg7->AddEntry(dataChi2B_Ref, "Unsmeared B-Only");
    leg7->AddEntry(dataChi2SB_Ref, "Unsmeared S+B");

    leg7->SetFillColor(0);
    leg7->Draw("same");
    if (doPrint)
      c7->Print("dataChi2.eps");
    if (doPrint)
      c7->Print(buffer);
  }

  if (peChi2B) {
    TCanvas* c11 = new TCanvas("eleven", "eleven");
    can_style(c11);
    h1_style(peChi2B);
    peChi2B->GetXaxis()->SetTitle("Fit #chi^{2} Value");
    peChi2B->GetYaxis()->SetTitle("Frequency");
    peChi2B->SetTitle("Chi2 for Fits to Pseudo-Experiments");
    peChi2B->Draw("hist");
    peChi2SB->Draw("samehist");
    TLegend* leg11 = new TLegend(0.65, 0.65, 0.85, 0.85);
    leg11->AddEntry(peChi2B, "Background-Only");
    leg11->AddEntry(peChi2SB, "Signal+Bkgd");
    leg11->SetFillColor(0);
    leg11->Draw("same");
    if (doPrint)
      c11->Print("peChi2.eps");
    if (doPrint)
      c11->Print(buffer);
  }

  if (minuitEmatSB) {
    TCanvas* cme = new TCanvas("cme", "cme");
    minuitEmatSB->Draw("Lego2");
    if (doPrint)
      cme->Print(buffer);

    TCanvas* cme2 = new TCanvas("cme2", "cme2");
    can_style(cme2);
    //  minuitEmatB->GetXaxis()->LabelsOption("d");
    //  minuitEmatB->GetYaxis()->LabelsOption("d");
    //  minuitEmatB->GetXaxis()->CenterLabels(false);
    minuitEmatB->Draw("Lego2");
    if (doPrint)
      cme2->Print(buffer);
  }

  if (fitErrsSB) {
    TCanvas* cee = new TCanvas("cee", "cee");
    can_style(cee);
    cee->SetGridy();
    cee->SetGridx();
    h1_style(fitErrsSB);
    fitErrsSB->SetTitle("Post-Fit Error Size Estimate");
    fitErrsSB->GetYaxis()->SetTitle("N-Sigma");
    fitErrsSB->GetXaxis()->LabelsOption("v");
    fitErrsSB->SetMaximum(1.5);
    fitErrsSB->SetMinimum(0);
    fitErrsSB->Draw("hist");
    fitErrsM->SetMaximum(1.5);
    fitErrsM->SetMinimum(0);
    //  fitErrsM->Draw("hist");
    fitErrsM->Draw("samehist");
    fitErrsB->Draw("samehist");

    TLegend* legcee = new TLegend(0.2, 0.65, 0.45, 0.85);
    legcee->AddEntry(fitErrsB, "Collie B-Only");
    legcee->AddEntry(fitErrsSB, "Collie S+B");
    legcee->AddEntry(fitErrsM, "MINUIT");
    legcee->SetFillColor(0);
    legcee->Draw("same");
    if (doPrint)
      cee->Print(buffer);

    TCanvas* c8 = new TCanvas("eight", "eight");
    can_style(c8);
    h1_style(fitParamsB);
    fitParamsB->SetMaximum(1.0);
    fitParamsB->SetMinimum(-1.0);
    fitParamsB->GetYaxis()->SetTitle("N-Sigma");
    fitParamsB->SetTitle("Best Fit to Data Parameters");
    fitParamsB->GetXaxis()->LabelsOption("v");
    fitParamsB->Draw("hist");
    fitParamsSB->Draw("samehist");
    c8->SetGridx();
    c8->SetGridy();
    TLegend* leg8 = new TLegend(0.4, 0.65, 0.65, 0.85);
    leg8->AddEntry(fitParamsB, "Bkgd-Only Fit");
    leg8->AddEntry(fitParamsSB, "Sig+Bkgd Fit");
    leg8->SetFillColor(0);
    leg8->Draw("same");
    if (doPrint)
      c8->Print("fitParams.eps");
    if (doPrint)
      c8->Print(buffer);
  }
  if (systDChi2_1) {
    TCanvas* ccdc = new TCanvas("Syst  Delta Chi^{2}", "Syst Delta Chi^{2}");
    can_style(ccdc);
    h1_style(systDChi2_1);
    systDChi2_1->SetLineWidth(3);
    systDChi2_1->SetLineColor(2);
    ccdc->SetGridx();
    ccdc->SetGridy();
    systDChi2_1->SetTitle("Per-Systematic Impact to Expected Limits");
    systDChi2_1->GetYaxis()->SetTitle("Fractional Loss in Limit");
    systDChi2_1->Draw("hist");
    if (doPrint)
      ccdc->Print(buffer);

    TCanvas* ccdc2 = new TCanvas("Syst Delta Chi^{2}, 2D", "Syst Delta Chi^{2}, 2D");
    can_style(ccdc2);
    h1_style(systDChi2_2);
    ccdc2->SetGridx();
    ccdc2->SetGridy();
    systDChi2_2->SetTitle("2D Per-Systematic Impact to Expected Limits");
    systDChi2_2->GetZaxis()->SetTitle("Fractional Loss in Limit");
    systDChi2_2->GetXaxis()->LabelsOption("v");
    systDChi2_2->GetYaxis()->LabelsOption("v");
    systDChi2_2->Draw("lego2");
    if (doPrint)
      ccdc2->Print(buffer);

    TCanvas* csdc = new TCanvas("Source Delta Chi^{2}", "Source  Delta Chi^{2}");
    can_style(csdc);
    h1_style(sourceDChi2_1);
    sourceDChi2_1->SetLineWidth(3);
    sourceDChi2_1->SetLineColor(2);
    csdc->SetGridx();
    csdc->SetGridy();
    sourceDChi2_1->GetXaxis()->LabelsOption("v");
    sourceDChi2_1->SetTitle("Per-Source Impact to Expected Limits");
    sourceDChi2_1->GetYaxis()->SetTitle("Fractional Loss in Limit");
    sourceDChi2_1->Draw("hist");
    if (doPrint)
      csdc->Print(buffer);

    TCanvas* csdc2 = new TCanvas("Source Delta Chi^{2}, 2D", "Scource  Delta Chi^{2}, 2D");
    can_style(csdc2);
    h1_style(sourceDChi2_2);
    csdc2->SetGridx();
    csdc2->SetGridy();
    sourceDChi2_2->SetTitle("2D Per-Source Impact to Expected Limits");
    sourceDChi2_2->GetZaxis()->SetTitle("Fractional Loss in Limit");
    sourceDChi2_2->Draw("lego2");
    sourceDChi2_2->GetXaxis()->LabelsOption("v");
    sourceDChi2_2->GetYaxis()->LabelsOption("v");
    if (doPrint)
      csdc2->Print(buffer);
  }

  if (nMinusOneChanB) {
    TCanvas* cnchan = new TCanvas("n-1 chans", "n-1 chans");
    can_style(cnchan);
    h1_style(nMinusOneChanB);
    cnchan->SetGridx();
    cnchan->SetGridy();
    nMinusOneChanB->SetTitle("N-1 Test, Channels");
    nMinusOneChanB->GetYaxis()->SetTitle("#Delta#chi^{2}");
    nMinusOneChanB->GetXaxis()->LabelsOption("v");
    nMinusOneChanB->Draw("hist");
    nMinusOneChanSB->Draw("samehist");
    leg8->Draw("same");
    if (doPrint)
      cnchan->Print(buffer);

    TCanvas* cnsyst = new TCanvas("n-1 syst", "n-1 syst");
    can_style(cnsyst);
    h1_style(nMinusOneSystB);
    cnsyst->SetGridx();
    cnsyst->SetGridy();
    max = TMath::Max(nMinusOneSystB->GetMaximum(), nMinusOneSystSB->GetMaximum());
    nMinusOneSystB->SetMaximum(max * 1.05);
    nMinusOneSystB->SetTitle("N-1 Test, Systemtatics");
    nMinusOneSystB->GetYaxis()->SetTitle("#Delta#chi^{2}");
    nMinusOneSystB->GetXaxis()->LabelsOption("v");
    nMinusOneSystB->Draw("hist");
    nMinusOneSystSB->Draw("samehist");
    leg8->Draw("same");
    if (doPrint)
      cnsyst->Print(buffer);

    TCanvas* cnbins = new TCanvas("n-1 bins", "n-1 bins");
    can_style(cnbins);
    h1_style(nMinusOneBinsB);
    cnbins->SetGridx();
    cnbins->SetGridy();
    double min = TMath::Min(nMinusOneBinsB->GetMinimum(), nMinusOneBinsSB->GetMinimum());
    nMinusOneBinsB->SetMinimum(min * 1.05);
    nMinusOneBinsB->SetTitle("N-1 Test, Bins");
    nMinusOneBinsB->GetXaxis()->SetTitle("Bin Index");
    nMinusOneBinsB->GetYaxis()->SetTitle("#Delta#chi^{2}");
    nMinusOneBinsB->Draw("hist");
    nMinusOneBinsSB->Draw("samehist");
    leg8->Draw("same");
    if (doPrint)
      cnbins->Print(buffer);
  }

  if (sOverB) {
    TCanvas* c9 = new TCanvas("nine", "nine");
    can_style(c9);
    h1_style(sOverB);
    sOverB->GetXaxis()->SetTitle("Final Variable");
    sOverB->GetYaxis()->SetTitle("s/b");
    sOverB->Draw("hist");
    if (doPrint)
      c9->Print(buffer);

    TCanvas* c10 = new TCanvas("ten", "ten");
    can_style(c10);
    h1_style(LOGsOverB);
    LOGsOverB->GetXaxis()->SetTitle("Final Variable");
    LOGsOverB->GetYaxis()->SetTitle("S#timesLog(1+S/B)");
    LOGsOverB->Draw("hist");
    if (doPrint)
      c10->Print("logSoverB.eps");
    if (doPrint)
      c10->Print(buffer);

    TCanvas* c14 = new TCanvas("14", "14");
    can_style(c14);
    c14->SetLogy();
    h1_style(LOGsOverB_bkg);
    LOGsOverB_bkg->GetXaxis()->SetTitle("S#timesLog(1+S/B)");
    LOGsOverB_bkg->GetYaxis()->SetTitle("Events/Bin");
    LOGsOverB_bkg->SetTitle("Rebinned S*Log(1+S/B)");

    LOGsOverB_dat->SetMarkerStyle(22);
    LOGsOverB_dat->SetMarkerSize(1.4);
    LOGsOverB_bkg->Draw("hist");
    LOGsOverB_fit->Draw("samehist");
    LOGsOverB_dat->Draw("PEsame");
    TLegend* leg14 = new TLegend(0.71, 0.68, 0.98, 0.94);
    leg_style(leg14);
    leg14->AddEntry(LOGsOverB_bkg, "Nominal Model");
    leg14->AddEntry(LOGsOverB_fit, "Best Fit Model");
    leg14->AddEntry(LOGsOverB_dat, "Data");
    leg14->SetFillColor(0);
    leg14->Draw("same");
    if (doPrint)
      c14->Print(buffer);
  }

  if (minuitIterB) {
    TCanvas* c12 = new TCanvas("twelve", "twelve");
    can_style(c12);
    h1_style(minuitIterB);

    double maxI = 1;
    double maxII = 1;
    for (int i = 1; i <= minuitIterB->GetNbinsX(); i++) {
      if (minuitIterB->GetBinContent(i) > 0)
        maxI = i * minuitIterB->GetBinWidth(i);
      if (minuitIterSB->GetBinContent(i) > 0)
        maxII = i * minuitIterB->GetBinWidth(i);
    }

    minuitIterB->GetXaxis()->SetRangeUser(0, TMath::Max(maxI, maxII) * 1.1);
    minuitIterB->GetXaxis()->SetTitle("# Iterations");
    minuitIterB->SetTitle("MINUIT Fit Iterations");
    minuitIterB->Draw("hist");
    minuitIterSB->Draw("samehist");

    TLegend* leg12 = new TLegend(0.65, 0.65, 0.85, 0.85);
    leg12->AddEntry(minuitIterB, "Background-Only");
    leg12->AddEntry(minuitIterSB, "Signal+Bkgd");
    leg12->SetFillColor(0);
    leg12->Draw("same");
    if (doPrint)
      c12->Print("minuitIter.eps");
    if (doPrint)
      c12->Print(buffer);

    TCanvas* c13 = new TCanvas("thirteen", "thirteen");
    can_style(c13);
    h1_style(minuitStatusB);
    TLegend* leg13 = new TLegend(0.2, 0.65, 0.45, 0.85);
    leg13->AddEntry(minuitStatusB, "Background-Only");
    leg13->AddEntry(minuitStatusSB, "Signal+Bkgd");
    c13->SetLogy();
    minuitStatusB->GetXaxis()->SetTitle("Status Bit");
    minuitStatusB->SetTitle("MINUIT Fit Status");
    minuitStatusB->Draw("hist");
    minuitStatusSB->Draw("samehist");
    leg13->SetFillColor(0);
    leg13->Draw("same");
    if (doPrint)
      c13->Print("minuitStatus.eps");
    if (doPrint)
      c13->Print(buffer);
  }

  int nsb = fitParamsB->GetNbinsX();
  const int maxSyst = nsb;
  TH1F* sigPull_D[maxSyst];
  TH1F* bkgPull_D[maxSyst];
  TH1F* sigPull_PE[maxSyst];
  TH1F* bkgPull_PE[maxSyst];
  TH1F* chi2_sb[maxSyst];
  TH1F* chi2_b[maxSyst];

  TH1F* p1S_bkgd[maxSyst];
  TH1F* m1S_bkgd[maxSyst];
  TH1F* p1S_sig[maxSyst];
  TH1F* m1S_sig[maxSyst];

  TH1F* dataParSB[maxSyst];
  TH1F* dataParB[maxSyst];

  char dest[256];
  for (int i = 0; i < maxSyst; i++) {

    sprintf(dest, "S+B Fit to Data: pull %d", i);
    sigPull_D[i] = (TH1F*)f->Get(dest);
    if (!sigPull_D[i])
      continue;
    sprintf(dest, "B-Only Fit to Data: pull %d", i);
    bkgPull_D[i] = (TH1F*)f->Get(dest);

    sprintf(dest, "S+B Fit to PEs: pull %d", i);
    sigPull_PE[i] = (TH1F*)f->Get(dest);
    sprintf(dest, "B-Only Fit to PEs: pull %d", i);
    bkgPull_PE[i] = (TH1F*)f->Get(dest);

    sprintf(dest, "S+B Fit chi2: %d", i);
    chi2_sb[i] = (TH1F*)f->Get(dest);

    sprintf(dest, "B-Only Fit chi2: %d", i);
    chi2_b[i] = (TH1F*)f->Get(dest);

    sprintf(dest, "Plus 1 Sigma Bkgd: %d", i);
    p1S_bkgd[i] = (TH1F*)f->Get(dest);

    sprintf(dest, "Minus 1 Sigma Bkgd: %d", i);
    m1S_bkgd[i] = (TH1F*)f->Get(dest);

    sprintf(dest, "Plus 1 Sigma Signal: %d", i);
    p1S_sig[i] = (TH1F*)f->Get(dest);

    sprintf(dest, "Minus 1 Sigma Signal: %d", i);
    m1S_sig[i] = (TH1F*)f->Get(dest);

    p1S_bkgd[i]->SetLineWidth(4);
    m1S_bkgd[i]->SetLineWidth(4);

    p1S_bkgd[i]->SetLineColor(2);
    m1S_bkgd[i]->SetLineColor(4);

    p1S_sig[i]->SetLineWidth(4);
    m1S_sig[i]->SetLineWidth(4);

    p1S_sig[i]->SetLineColor(2);
    m1S_sig[i]->SetLineColor(4);

    sigPull_D[i]->SetLineColor(2);
    sigPull_D[i]->SetLineWidth(4);
    bkgPull_D[i]->SetLineWidth(4);

    sigPull_PE[i]->SetLineColor(2);
    sigPull_PE[i]->SetLineWidth(4);
    bkgPull_PE[i]->SetLineWidth(4);
    chi2_sb[i]->SetLineWidth(4);
    chi2_b[i]->SetLineWidth(4);
    chi2_sb[i]->SetLineColor(2);

    sprintf(dest, "CLdpsb%d", i);
    dataParSB[i] = (TH1F*)sigPull_PE[i]->Clone(dest);
    dataParSB[i]->Scale(0);
    dataParSB[i]->Fill(fitParamsSB->GetBinContent(i + 1), 1e9);

    sprintf(dest, "CLdpb%d", i);
    dataParB[i] = (TH1F*)sigPull_PE[i]->Clone(dest);
    dataParB[i]->Scale(0);
    dataParB[i]->Fill(fitParamsB->GetBinContent(i + 1), 1e9);

    dataParSB[i]->SetLineColor(3);
    dataParB[i]->SetLineColor(4);
    dataParB[i]->SetFillColor(0);
    dataParSB[i]->SetFillColor(0);
    dataParSB[i]->SetLineWidth(2);
    dataParB[i]->SetLineWidth(2);
  }

  double variance_B[maxSyst][2];
  double minVal_B[maxSyst];
  double gWidthD_B[maxSyst];
  double gWidthPE_B[maxSyst];
  double gMeanD_B[maxSyst];
  double gMeanPE_B[maxSyst];

  double variance_SB[maxSyst][2];
  double minVal_SB[maxSyst];
  double gWidthD_SB[maxSyst];
  double gWidthPE_SB[maxSyst];
  double gMeanD_SB[maxSyst];
  double gMeanPE_SB[maxSyst];

  TCanvas* cT = new TCanvas("fit", "fit");

  for (int i = 0; i < maxSyst; i++) {
    if (!sigPull_D[i])
      continue;
    int minBinSB = chi2_sb[i]->GetMinimumBin();
    minVal_SB[i] = chi2_sb[i]->GetBinCenter(minBinSB);

    TF1* f2 = new TF1("f2", "gaus", -2.5, 2.5);
    sigPull_PE[i]->Fit("f2", "RQ1");
    gWidthPE_SB[i] = f2->GetParameter(2);
    gMeanPE_SB[i] = f2->GetParameter(1);

    sigPull_D[i]->Fit("f2", "RQ2");
    gWidthD_SB[i] = f2->GetParameter(2);
    gMeanD_SB[i] = f2->GetParameter(1);

    int minBinB = chi2_b[i]->GetMinimumBin();
    minVal_B[i] = chi2_b[i]->GetBinCenter(minBinB);

    bkgPull_PE[i]->Fit("f2", "RQ1");
    gWidthPE_B[i] = f2->GetParameter(2);
    gMeanPE_B[i] = f2->GetParameter(1);

    bkgPull_D[i]->Fit("f2", "RQ2");
    gWidthD_B[i] = f2->GetParameter(2);
    gMeanD_B[i] = f2->GetParameter(1);

    double p1s = 1000;
    double dchiP = 10000;
    double m1s = 1000;
    double dchiM = 10000;
    double mval = chi2_sb[i]->GetBinContent(minBinSB);
    for (int b = 1; b < chi2_sb[i]->GetNbinsX(); b++) {
      if (b < minBinSB) {
        if (fabs((chi2_sb[i]->GetBinContent(b) - mval) - 1) < dchiM) {
          dchiM = fabs((chi2_sb[i]->GetBinContent(b) - mval) - 1);
          m1s = chi2_sb[i]->GetBinCenter(b);
        }
      } else {
        if (fabs((chi2_sb[i]->GetBinContent(b) - mval) - 1) < dchiP) {
          dchiP = fabs((chi2_sb[i]->GetBinContent(b) - mval) - 1);
          p1s = chi2_sb[i]->GetBinCenter(b);
        }
      }
    }
    variance_SB[i][0] = minVal_SB[i] - m1s;
    variance_SB[i][1] = p1s - minVal_SB[i];

    p1s = 1000;
    dchiP = 10000;
    m1s = 1000;
    dchiM = 10000;
    mval = chi2_b[i]->GetBinContent(minBinB);
    for (int b = 1; b < chi2_b[i]->GetNbinsX(); b++) {
      if (b < minBinB) {
        if (fabs((chi2_b[i]->GetBinContent(b) - mval) - 1) < dchiM) {
          dchiM = fabs((chi2_b[i]->GetBinContent(b) - mval) - 1);
          m1s = chi2_b[i]->GetBinCenter(b);
        }
      } else {
        if (fabs((chi2_b[i]->GetBinContent(b) - mval) - 1) < dchiP) {
          dchiP = fabs((chi2_b[i]->GetBinContent(b) - mval) - 1);
          p1s = chi2_b[i]->GetBinCenter(b);
        }
      }
    }
    variance_B[i][0] = minVal_B[i] - m1s;
    variance_B[i][1] = p1s - minVal_B[i];
  }

  for (int i = 0; i < maxSyst; i++) {
    sprintf(dest, "Syst Canvas %d", i);
    if (!sigPull_D[i])
      continue;
    TCanvas* c1 = new TCanvas(dest, dest);
    can_style(c1);
    c1->Divide(2, 2);
    c1->cd(1);
    TVirtualPad* p1 = c1->GetPad(1);
    p1->SetLeftMargin(0.13);
    p1->SetBottomMargin(0.14);
    p1->SetTopMargin(0.075);
    p1->SetRightMargin(0.01);
    h1_style(sigPull_D[i]);
    sigPull_D[i]->GetXaxis()->SetTitle("N #sigma");
    sigPull_D[i]->GetYaxis()->SetTitle("Frequency");
    max = TMath::Max(sigPull_D[i]->GetMaximum(), bkgPull_D[i]->GetMaximum());
    sigPull_D[i]->Scale((max / (sigPull_D[i]->GetMaximum() + 1e-5)) /
                        (sigPull_D[i]->Integral() + 1));
    bkgPull_D[i]->Scale((max / (bkgPull_D[i]->GetMaximum() + 1e-5)) /
                        (bkgPull_D[i]->Integral() + 1));
    sigPull_D[i]->Draw("hist");
    bkgPull_D[i]->Draw("samehist");
    TLegend* leg = new TLegend(0.65, 0.65, 0.85, 0.85);
    leg_style(leg);
    leg->AddEntry(sigPull_D[i], "S+B Fit");
    leg->AddEntry(bkgPull_D[i], "B-Only Fit");
    leg->SetFillColor(0);
    leg->Draw("same");
    TLatex* tex0 = new TLatex(-4, sigPull_D[i]->GetMaximum() * 0.9, "(A)");
    tex0->SetTextSize(0.07);
    tex0->Draw("same");

    c1->cd(2);
    p1 = c1->GetPad(2);
    p1->SetLeftMargin(0.13);
    p1->SetBottomMargin(0.14);
    p1->SetTopMargin(0.075);
    p1->SetRightMargin(0.01);
    h1_style(sigPull_PE[i]);
    sigPull_PE[i]->GetXaxis()->SetTitle("N-Sigma");
    sigPull_PE[i]->GetYaxis()->SetTitle("Frequency");
    max = TMath::Max(sigPull_PE[i]->GetMaximum(), bkgPull_PE[i]->GetMaximum());
    sigPull_PE[i]->Scale((max / (sigPull_PE[i]->GetMaximum() + 1e-5)) /
                         (sigPull_PE[i]->Integral() + 1));
    bkgPull_PE[i]->Scale((max / (bkgPull_PE[i]->GetMaximum() + 1e-5)) /
                         (bkgPull_PE[i]->Integral() + 1));
    sigPull_PE[i]->Draw("hist");
    bkgPull_PE[i]->Draw("samehist");
    dataParSB[i]->Draw("samehist");
    dataParB[i]->Draw("samehist");
    TLegend* leg2 = new TLegend(0.65, 0.65, 0.85, 0.85);
    leg_style(leg2);
    leg2->AddEntry(sigPull_PE[i], "S+B Fit");
    leg2->AddEntry(bkgPull_PE[i], "B-Only Fit");
    leg2->AddEntry(dataParB[i], "B-Only Fit to Data");
    leg2->AddEntry(dataParSB[i], "S+B Fit to Data");
    leg2->SetFillColor(0);
    leg2->Draw("same");
    TLatex* tex0 = new TLatex(-2, sigPull_PE[i]->GetMaximum() * 0.9, "(B)");
    tex0->SetTextSize(0.07);
    tex0->Draw("same");

    c1->cd(3);
    p1 = c1->GetPad(3);
    p1->SetLeftMargin(0.13);
    p1->SetBottomMargin(0.14);
    p1->SetTopMargin(0.075);
    p1->SetRightMargin(0.01);

    h1_style(chi2_sb[i]);
    removeMin(chi2_sb[i]);
    removeMin(chi2_b[i]);
    chi2_sb[i]->GetXaxis()->SetRangeUser(-5.9, 5.9);
    double max = TMath::Max(chi2_sb[i]->GetMaximum(), chi2_b[i]->GetMaximum());
    chi2_sb[i]->SetMaximum(max * 1.5);
    chi2_sb[i]->GetXaxis()->SetTitle("N-Sigma");
    chi2_sb[i]->GetYaxis()->SetTitle("#Delta #chi^{2}");
    chi2_sb[i]->Draw("hist");
    chi2_b[i]->Draw("samehist");
    leg->Draw("same");

    double val1 =
        chi2_sb[i]->GetMinimum() + (chi2_sb[i]->GetMaximum() - chi2_sb[i]->GetMinimum()) * 0.9;
    TLatex* tex0 = new TLatex(0.0, val1, "(C)");
    tex0->SetTextSize(0.07);
    tex0->Draw("same");

    c1->cd(4);
    char dest[256];

    TLatex* tex0 = new TLatex(0.1, 0.9, "[A] Smeared B-Only (S+B) Fit to Data");
    tex0->SetTextSize(0.07);
    tex0->Draw();

    sprintf(dest, "        Mean: %.3f (%.3f)", gMeanD_B[i], gMeanD_SB[i]);
    TLatex* tex0 = new TLatex(0.1, 0.8, dest);
    tex0->SetTextSize(0.07);
    tex0->Draw();
    sprintf(dest, "        Width: %.3f (%.3f)", gWidthD_B[i], gWidthD_SB[i]);
    TLatex* tex1 = new TLatex(0.1, 0.7, dest);
    tex1->SetTextSize(0.07);
    tex1->Draw();

    TLatex* tex0 = new TLatex(0.1, 0.6, "[B] B-Only (S+B) Fit to PseudoData");
    tex0->SetTextSize(0.07);
    tex0->Draw();

    sprintf(dest, "        Mean: %.3f (%.3f)", gMeanPE_B[i], gMeanPE_SB[i]);
    TLatex* tex0 = new TLatex(0.1, 0.5, dest);
    tex0->SetTextSize(0.07);
    tex0->Draw();
    sprintf(dest, "        Width: %.3f (%.3f)", gWidthPE_B[i], gWidthPE_SB[i]);
    TLatex* tex1 = new TLatex(0.1, 0.4, dest);
    tex1->SetTextSize(0.07);
    tex1->Draw();

    TLatex* tex0 = new TLatex(0.1, 0.3, "[C] B-Only (S+B) #chi^{2} Response Function");
    tex0->SetTextSize(0.07);
    tex0->Draw();

    sprintf(dest, "        Minimum: %.2f (%.2f)", minVal_B[i], minVal_SB[i]);
    TLatex* tex0 = new TLatex(0.1, 0.2, dest);
    tex0->SetTextSize(0.07);
    tex0->Draw();
    sprintf(dest, "        +1#sigma: %.3f (%.3f)", variance_B[i][0], variance_SB[i][0]);
    TLatex* tex1 = new TLatex(0.1, 0.1, dest);
    tex1->SetTextSize(0.07);
    tex1->Draw();
    sprintf(dest, "        +1#sigma: %.3f (%.3f)", variance_B[i][1], variance_SB[i][1]);
    TLatex* tex2 = new TLatex(0.1, 0.0, dest);
    tex2->SetTextSize(0.07);
    tex2->Draw();

    sprintf(dest, "systCanvas%d.eps", i);
    if (doPrint)
      c1->Print(dest);
    if (doPrint)
      c1->Print(buffer);
  }

  for (int i = 0; i < maxSyst; i++) {
    sprintf(dest, "Bkgd Fluct Canvas %d", i);

    TString pch = p1S_bkgd[i]->GetTitle();
    pch.ReplaceAll("Plus 1 Sigma", "#pm1 Sigma");
    p1S_bkgd[i]->SetTitle(pch.Data());
    p1S_bkgd[i]->Add(bkgdDist, -1);
    m1S_bkgd[i]->Add(bkgdDist, -1);

    p1S_bkgd[i]->Divide(bkgdDist);
    m1S_bkgd[i]->Divide(bkgdDist);

    max = fabs(TMath::Max(p1S_bkgd[i]->GetMaximum(), m1S_bkgd[i]->GetMaximum()));
    min = fabs(TMath::Min(p1S_bkgd[i]->GetMinimum(), m1S_bkgd[i]->GetMinimum()));
    if (min > max)
      max = min;

    const char* ttitle = p1S_bkgd[i]->GetTitle();
    h1_style(p1S_bkgd[i]);
    p1S_bkgd[i]->GetYaxis()->SetRangeUser(-1.5 * max, 1.5 * max);
    p1S_bkgd[i]->GetYaxis()->SetTitle("Fractional Uncertainty");
    p1S_bkgd[i]->GetXaxis()->SetTitle("Final Variable");

    TCanvas* c1 = new TCanvas(dest, dest);
    c1->SetLeftMargin(0.14);
    c1->SetBottomMargin(0.14);
    c1->SetTopMargin(0.06);
    c1->SetRightMargin(0.05);
    c1->SetGridx();
    c1->SetGridy();
    p1S_bkgd[i]->Draw("hist");
    m1S_bkgd[i]->Draw("samehist");

    TLegend* leg = new TLegend(0.82, 0.73, 0.998, 0.945, NULL, "brNDC");
    leg_style(leg);
    leg->AddEntry(p1S_bkgd[i], "+1#sigma");
    leg->AddEntry(m1S_bkgd[i], "-1#sigma");
    leg->SetFillColor(0);
    leg->Draw("same");
    if (doPrint)
      c1->Print(buffer);

    sprintf(dest, "Signal Fluct Canvas %d", i);
    pch = p1S_sig[i]->GetTitle();
    pch.ReplaceAll("Plus 1 Sigma", "#pm1 Sigma");
    p1S_sig[i]->SetTitle(pch.Data());
    p1S_sig[i]->Add(sigDist, -1);
    m1S_sig[i]->Add(sigDist, -1);

    p1S_sig[i]->Divide(sigDist);
    m1S_sig[i]->Divide(sigDist);

    max = fabs(TMath::Max(p1S_sig[i]->GetMaximum(), m1S_sig[i]->GetMaximum()));
    min = fabs(TMath::Min(p1S_sig[i]->GetMinimum(), m1S_sig[i]->GetMinimum()));
    if (min > max)
      max = min;

    ttitle = p1S_sig[i]->GetTitle();
    h1_style(p1S_sig[i]);
    p1S_sig[i]->GetYaxis()->SetRangeUser(-1.5 * max, 1.5 * max);
    p1S_sig[i]->GetYaxis()->SetTitle("Fractional Uncertainty");
    p1S_sig[i]->GetXaxis()->SetTitle("Final Variable");

    TCanvas* c2 = new TCanvas(dest, dest);
    c2->SetLeftMargin(0.14);
    c2->SetBottomMargin(0.14);
    c2->SetTopMargin(0.06);
    c2->SetRightMargin(0.05);
    c2->SetGridx();
    c2->SetGridy();
    p1S_sig[i]->Draw("hist");
    m1S_sig[i]->Draw("samehist");

    TLegend* leg2 = new TLegend(0.82, 0.73, 0.998, 0.945, NULL, "brNDC");
    leg_style(leg2);
    leg2->AddEntry(p1S_sig[i], "+1#sigma");
    leg2->AddEntry(m1S_sig[i], "-1#sigma");
    leg2->SetFillColor(0);
    leg2->Draw("same");
    if (doPrint)
      c2->Print(buffer);
  }

  sprintf(buffer, "fitResults.ps]");
  if (doPrint)
    ccc.Print(buffer); //close ps file
}

void fitResults(const char* testFile = "myFitResults.root") {

  makePlots(testFile, 1);
}
