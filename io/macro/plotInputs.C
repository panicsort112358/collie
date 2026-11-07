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

void plotInputs(const char* infile, int massPoint = -1, const char* channelName = "My Channel",
                const char* varType = "MyVariable") {
  style();

  if (massPoint == -1) {
    printf("You must specify a masspoint!\n");
    return;
  }

  TObjArray sigNames;
  TObjArray bkgdNames;

  TFile* fin = new TFile(infile, "READ");

  TList* myKeys = fin->GetListOfKeys();
  char mass[256];
  sprintf(mass, "%d", massPoint);
  for (int i = 0; i < myKeys->GetEntries(); i++) {
    TObjString* ostr = myKeys->At(i);
    TString* astr = (TString*)ostr->GetString();
    TString* str = new TString(astr->Data());

    if (str->Contains(mass)) {
      if (str->Contains("Signal")) {
        if (str->Contains("Final Variable")) {
          str->Replace(0, 7, "");
          str->ReplaceAll("Final Variable - ", "");
          str->ReplaceAll(mass, "");
          str->ReplaceAll(" ", "");
          if (sigNames.FindObject(str->Data()) == NULL)
            sigNames.Add(new TObjString(str->Data()));
        }
      } else {
        if (str->Contains("Final Variable") && !str->Contains("All") && !str->Contains("Data")) {

          str->ReplaceAll("Final Variable - ", "");
          str->ReplaceAll(mass, "");
          str->ReplaceAll(" ", "");
          if (bkgdNames.FindObject(str->Data()) == NULL)
            bkgdNames.Add(new TObjString(str->Data()));
        }
      }

    } //contains mass
  } //key entries

  char buffer[256];
  sprintf(buffer, "Data Final Variable - %d", massPoint);
  TH1D* data = (TH1D*)fin->Get(buffer);

  sprintf(buffer, "All Bkgd Final Variable - %d", massPoint);
  TH1D* bkgd = (TH1D*)fin->Get(buffer);

  TH1D* signal = NULL;
  for (int b = 0; b < sigNames.GetEntries(); b++) {
    const char* tName = ((TString*)(((TObjString*)(sigNames.At(b)))->GetString()))->Data();
    sprintf(buffer, "Signal %s Final Variable - %d", tName, massPoint);
    printf("Getting: %s\n", buffer);
    if (signal == NULL)
      signal = (TH1D*)fin->Get(buffer);
    else
      signal->Add((TH1D*)fin->Get(buffer));
  }
  signal->Scale(10.0);

  TObjArray bkgdHists;
  for (int b = 0; b < bkgdNames.GetEntries(); b++) {
    const char* tName = ((TString*)(((TObjString*)(bkgdNames.At(b)))->GetString()))->Data();
    sprintf(buffer, "%s Final Variable - %d", tName, massPoint);
    printf("Getting: %s\n", buffer);
    TH1D* h = (TH1D*)fin->Get(buffer);
    if (h == NULL)
      printf("NULL histo!\n");
    bkgdHists.Add((TH1D*)fin->Get(buffer));
  }

  sprintf(buffer, "Final Variable %d", massPoint);
  THStack* hstack = new THStack(buffer, buffer);

  signal->SetLineColor(2);
  signal->SetLineWidth(3);

  int colors[20];
  colors[0] = 3;
  colors[1] = 4;
  colors[2] = 5;
  colors[3] = 6;
  colors[4] = 18;
  colors[5] = 38;
  colors[6] = 5;

  for (int b = 0; b < bkgdHists.GetEntries(); b++) {
    TH1D* h = (TH1D*)bkgdHists.At(b);
    printf("Adding: %s\n", h->GetName());
    h->SetFillColor(colors[b]);
    hstack->Add(h);
  }
  hstack->Draw();

  hstack_style(hstack);
  sprintf(buffer, "");
  hstack->SetTitle(buffer);
  hstack->GetXaxis()->SetTitle("Final Variable");
  hstack->GetYaxis()->SetTitle("Events");
  hstack->GetYaxis()->SetTitleOffset(1.00);

  h1_style(data);
  data->SetMarkerSize(1.1);
  data->SetLineWidth(2.0);
  data->SetMarkerStyle(20);

  double max = bkgd->GetMaximum();
  if (data->GetMaximum() > max)
    max = data->GetMaximum();
  max *= 1.75;
  hstack->SetMaximum(max);

  TLegend* leg = new TLegend(0.65, 0.65, 0.9, 0.9);
  leg_style(leg);
  leg->AddEntry(data, "Data", "PE");
  leg->AddEntry(signal, "Signal #times10");
  for (int b = 0; b < bkgdHists.GetEntries(); b++) {
    TH1D* h = (TH1D*)bkgdHists.At(b);
    const char* tName = ((TString*)(((TObjString*)(bkgdNames.At(b)))->GetString()))->Data();
    leg->AddEntry(h, tName, "f");
  }
  leg->SetFillColor(0);
  leg->SetBorderSize(1);

  sprintf(buffer, "Final Variable Canvas");
  TCanvas* c2 = new TCanvas(buffer, buffer);
  hstack->Draw("hist");
  data->Draw("E1P0same");
  signal->Draw("samehist");
  leg->Draw("same");
}
