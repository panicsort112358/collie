static const int mFont = 132;

// ... Global attributes go here ...
void style() {
  gStyle->SetTextFont(mFont);
  gStyle->SetLabelFont(mFont, "X");
  gStyle->SetLabelFont(mFont, "Y");
  gStyle->SetTitleFont(mFont, "X");
  gStyle->SetTitleFont(mFont, "Y");

  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);

  gStyle->SetPadLeftMargin(0.145);
  gStyle->SetPadRightMargin(0.025);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetPadTopMargin(0.025);

  gStyle->SetTitleColor(1); // 0
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
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
  //  gStyle->SetErrorX(0.0);   // Horizontal error bar size
  gStyle->SetOptStat(0);
  //   gStyle->SetPaperSize(10.,12.);   // Printout size
}

// ... 1D histogram attributes; (2D to come) ...
void h1_style(TH1* h, int line_width = 3, int line_color = 1, int line_style = 1,
              int fill_style = 1001, int fill_color = 50, float y_min = -1111.,
              float y_max = -1111., int ndivx = 510, int ndivy = 510, int marker_style = 20,
              int marker_color = 1, float marker_size = 1.3, int optstat = 0) {
  /*
   h->SetLineWidth(line_width);
   h->SetLineColor(line_color);
   h->SetLineStyle(line_style);
   h->SetFillColor(fill_color);
   h->SetFillStyle(fill_style);
   h->SetMaximum(y_max);
   h->SetMinimum(y_min);
   h->GetXaxis()->SetNdivisions(ndivx);
   h->GetYaxis()->SetNdivisions(ndivy);

   h->SetMarkerStyle(marker_style);
   h->SetMarkerColor(marker_color);
   h->SetMarkerSize(marker_size);
   h->SetStats(optstat);
  */
  h->GetXaxis()->SetTitleFont(mFont);
  h->GetYaxis()->SetTitleFont(mFont);
  h->SetStats(kFALSE);

  h->SetLabelFont(mFont, "X");   // 42
  h->SetLabelFont(mFont, "Y");   // 42
  h->SetLabelOffset(0.000, "X"); // D=0.005
  h->SetLabelOffset(0.005, "Y"); // D=0.005
  h->SetLabelSize(0.06, "X");
  h->SetLabelSize(0.06, "Y");
  h->SetTitleOffset(0.95, "X");
  h->SetTitleOffset(0.95, "Y");
  h->SetTitleSize(0.065, "X");
  h->SetTitleSize(0.065, "Y");
  h->SetTitle(0);
}

void axis_style(TGaxis* a) {
  a->SetTitleFont(mFont);
  a->SetLabelFont(mFont);   // 42
  a->SetLabelOffset(0.006); // D=0.005
  a->SetLabelSize(0.06);
  a->SetTitleSize(0.065);
  a->SetTitleOffset(1.2);
}

void hstack_style(THStack* h) {
  if (h == NULL) {
    printf("This stack is null!\n");
    return;
  }
  h->GetXaxis()->SetTitleFont(mFont);
  h->GetYaxis()->SetTitleFont(mFont);
  h->GetXaxis()->SetLabelFont(mFont); // 42
  h->GetYaxis()->SetLabelFont(mFont); // 42

  h->GetXaxis()->SetLabelOffset(0.006); // D=0.005
  h->GetYaxis()->SetLabelOffset(0.006); // D=0.005

  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->SetTitleSize(0.065);
  h->GetYaxis()->SetTitleSize(0.065);
  h->SetTitle(0);
}

void h2_style(TH2* h) {
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
  h->SetTitleOffset(1.025, "Y");
  h->SetTitleSize(0.07, "X");
  h->SetTitleSize(0.07, "Y");
  h->SetTitle(0);
}

void leg_style(TLegend* h) {
  h->SetTextFont(mFont);
  h->SetTextSize(0.065);
  h->SetBorderSize(0);
  h->SetFillColor(0);
  h->SetFillStyle(0);
  h->SetLineColor(0);
}

void plotFitHistos(const char* infile, bool nullFit = false, int rebin = 1, int ntags = -1,
                   const char* varType = "My Variable") {

  style();

  TFile* fin = new TFile(infile);

  TString type("NULL Fit,");
  if (!nullFit)
    type = TString("TEST Fit,");

  TObjArray sigNames;
  TObjArray dataNames;
  TObjArray bkgdNames;
  TObjArray posNames;
  TObjArray negNames;
  TFile* fin = new TFile(infile, "READ");

  TList* myKeys = fin->GetListOfKeys();

  TString tagString1 = "";
  if (ntags == 0)
    tagString1 = "0tag";
  if (ntags == 1)
    tagString1 = "1tag";
  if (ntags == 2)
    tagString1 = "2tag";

  TString tagString2 = "";
  if (ntags == 0)
    tagString2 = "0Tag";
  if (ntags == 1)
    tagString2 = "1Tag";
  if (ntags == 2)
    tagString2 = "2Tag";

  for (int i = 0; i < myKeys->GetEntries(); i++) {
    TObjString* ostr = myKeys->At(i);
    TString* astr = (TString*)ostr->GetString();
    TString* str = new TString(astr->Data());
    if (!str->Contains(tagString1) && !str->Contains(tagString2))
      continue;

    if (str->Contains(type)) {
      if (str->Contains("Sig")) {
        if (sigNames.FindObject(str->Data()) == NULL)
          sigNames.Add(new TObjString(str->Data()));
        printf("Adding %s\n", str->Data());
      } else if (str->Contains("Bkgd") && !str->Contains("Sum")) {
        if (bkgdNames.FindObject(str->Data()) == NULL)
          bkgdNames.Add(new TObjString(str->Data()));
        printf("Adding %s\n", str->Data());
      } else if (str->Contains("Data")) {
        if (dataNames.FindObject(str->Data()) == NULL)
          dataNames.Add(new TObjString(str->Data()));
        printf("Adding %s\n", str->Data());
      }
    }
    if (str->Contains("Total")) {
      if (str->Contains("Positive Error Matrix Systematics")) {
        if (posNames.FindObject(str->Data()) == NULL)
          posNames.Add(new TObjString(str->Data()));
        printf("Adding %s\n", str->Data());
      } else if (str->Contains("Negative Error Matrix Systematics")) {
        if (negNames.FindObject(str->Data()) == NULL)
          negNames.Add(new TObjString(str->Data()));
        printf("Adding %s\n", str->Data());
      }
    } //contains mass
  } //key entries

  TObjArray bkgdHists;
  TObjArray dataHists;
  TObjArray sigHists;
  TObjArray p1sHists;
  TObjArray m1sHists;

  TH1D* totData = 0;
  TH1D* totBkgd = 0;
  TH1D* totSig = 0;
  TH1D* totSig0 = 0;
  TH1D* totSig1 = 0;
  TH1D* totP1S = 0;
  TH1D* totM1S = 0;
  //  TH1D* sigScales = (TH1D*)fin->Get("Signal Scale Factor");

  for (int b = 0; b < dataNames.GetEntries(); b++) {
    const char* tName = ((TString*)(((TObjString*)(dataNames.At(b)))->GetString()))->Data();
    TH1D* d = (TH1D*)fin->Get(tName);
    d->Rebin(rebin);
    dataHists.Add(d);
    //    printf("%s: %d\n",tName,d->GetNbinsX());
    if (totData == NULL)
      totData = d;
    else
      totData->Add(d);
  }

  for (int b = 0; b < bkgdNames.GetEntries(); b++) {
    const char* tName = ((TString*)(((TObjString*)(bkgdNames.At(b)))->GetString()))->Data();
    TH1D* h = (TH1D*)fin->Get(tName);
    h->Rebin(rebin);
    bkgdHists.Add(h);

    if (totBkgd == NULL)
      totBkgd = h;
    else
      totBkgd->Add(h);
  }

  //  TAxis* as = sigScales->GetXaxis();
  for (int b = 0; b < sigNames.GetEntries(); b++) {
    const char* tName = ((TString*)(((TObjString*)(sigNames.At(b)))->GetString()))->Data();
    TH1D* h = (TH1D*)fin->Get(tName);
    h->Rebin(rebin);
    sigHists.Add(h);
    TString sname(tName);
    if (sname.Contains("WW")) {
      if (totSig0 == NULL)
        totSig0 = (TH1D*)h->Clone("WW Sig");
      else
        totSig0->Add(h);
    }
    if (sname.Contains("VZ")) {
      if (totSig1 == NULL)
        totSig1 = (TH1D*)h->Clone("VZ Sig");
      else
        totSig1->Add(h);
    }

    if (totSig == NULL)
      totSig = h;
    else
      totSig->Add(h);
  }

  for (int b = 0; b < negNames.GetEntries(); b++) {
    const char* tName = ((TString*)(((TObjString*)(negNames.At(b)))->GetString()))->Data();
    TH1D* h = (TH1D*)fin->Get(tName);
    h->Rebin(rebin);
    m1sHists.Add(h);
    if (totM1S == NULL)
      totM1S = h;
    else
      totM1S->Add(h);
  }

  for (int b = 0; b < posNames.GetEntries(); b++) {
    const char* tName = ((TString*)(((TObjString*)(posNames.At(b)))->GetString()))->Data();

    TH1D* h = (TH1D*)fin->Get(tName);
    h->Rebin(rebin);
    p1sHists.Add(h);
    if (totP1S == NULL)
      totP1S = h;
    else
      totP1S->Add(h);
  }

  //  totP1S->Add(totBkgd,-1);
  //  totM1S->Add(totBkgd,-1);
  totM1S->Scale(-1);

  TH1D* dataSub = (TH1D*)totData->Clone("Data-Bkgd");
  dataSub->Add(totBkgd, -1);

  totP1S->SetLineColor(4);
  totM1S->SetLineColor(4);
  totP1S->SetLineWidth(3);
  totM1S->SetLineWidth(3);

  totSig->SetLineColor(2);
  totSig->SetFillColor(2);
  totSig0->SetLineColor(2);
  totSig0->SetFillColor(2);
  totSig1->SetLineColor(7);
  totSig1->SetFillColor(7);

  totBkgd->SetLineColor(5);
  totBkgd->SetFillColor(5);

  char buffer[256];
  sprintf(buffer, "Final Variable");
  THStack* hstack = new THStack(buffer, buffer);
  sprintf(buffer, "Final Variable2");
  THStack* hstack2 = new THStack(buffer, buffer);
  hstack->Add(totBkgd);
  //  hstack->Add(totSig);
  hstack->Add(totSig0);
  hstack->Add(totSig1);
  hstack->Draw();
  hstack2->Add(totSig0);
  hstack2->Add(totSig1);
  hstack2->Draw();

  h1_style(totData);
  totData->SetMarkerSize(1.1);
  totData->SetLineWidth(2.0);
  totData->SetMarkerStyle(20);

  h1_style(dataSub);
  dataSub->SetMarkerSize(1.1);
  dataSub->SetLineWidth(2.0);
  dataSub->SetMarkerStyle(20);

  double max = totBkgd->GetMaximum();
  if (totData->GetMaximum() > max)
    max = totData->GetMaximum();
  max *= 1.75;
  hstack->SetMaximum(max);

  TLegend* leg = new TLegend(0.65, 0.65, 0.9, 0.9);
  leg_style(leg);
  leg->AddEntry(totData, "Data", "PE");
  leg->AddEntry(totBkgd, "Bkgd", "F");
  leg->AddEntry(totSig0, "WW");
  leg->AddEntry(totSig1, "WZ/ZZ");
  leg->SetFillColor(0);
  leg->SetBorderSize(1);

  sprintf(buffer, "Final Variable Canvas");
  TCanvas* c2 = new TCanvas(buffer, buffer);
  hstack->Draw("hist");
  totData->Draw("E1P0same");
  //  totSig->Draw("samehist");
  leg->Draw("same");

  TLegend* leg2 = new TLegend(0.65, 0.65, 0.9, 0.9);
  leg_style(leg2);
  leg2->AddEntry(totData, "Data-Bkgd", "PE");
  leg2->AddEntry(totSig0, "WW");
  leg2->AddEntry(totSig1, "WZ/ZZ");
  leg2->AddEntry(totM1S, "Bkgd #pm 1#sigma");
  leg2->SetFillColor(0);
  leg2->SetBorderSize(1);

  max = dataSub->GetMaximum();
  if (totM1S->GetMaximum() > max)
    max = totM1S->GetMaximum();
  if (totP1S->GetMaximum() > max)
    max = totP1S->GetMaximum();
  max *= 1.75;

  double min = dataSub->GetMinimum();
  if (totM1S->GetMinimum() < min)
    min = totM1S->GetMinimum();
  if (totP1S->GetMinimum() < min)
    min = totP1S->GetMinimum();
  min *= 1.75;
  hstack2->SetMinimum(min);
  hstack2->SetMaximum(max);

  sprintf(buffer, "Subtracted Canvas");
  hstack2->SetTitle("");
  TCanvas* c3 = new TCanvas(buffer, buffer);
  //  dataSub->Draw("E1P0");
  hstack2->Draw("hist");
  dataSub->Draw("E1P0same");

  totM1S->Draw("samehist");
  totP1S->Draw("samehist");
  leg2->Draw("same");

  printf("Sig: %f (%f/%f), Bkgd: %f, Data: %f (%f)\n", totSig->Integral(), totSig0->Integral(),
         totSig1->Integral(), totBkgd->Integral(), totData->Integral(), dataSub->Integral());
}
