static const int mFont = 132;

void style() {
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
  gStyle->SetPadTopMargin(0.075);

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
  //   gStyle->SetPaperSize(10.,12.);   // Printout size
}

// ... 1D histogram attributes; (2D to come) ...
void h1_style(TH1* h, int line_width = 3, int line_color = 1, int line_style = 1,
              int fill_style = 1001, int fill_color = 50, float y_min = -1111.,
              float y_max = -1111., int ndivx = 510, int ndivy = 510, int marker_style = 20,
              int marker_color = 1, float marker_size = 1.3, int optstat = 0) {
  h->GetXaxis()->SetTitleFont(mFont);
  h->GetYaxis()->SetTitleFont(mFont);
  h->SetStats(kFALSE);

  h->SetLabelFont(mFont, "X");   // 42
  h->SetLabelFont(mFont, "Y");   // 42
  h->SetLabelOffset(0.000, "X"); // D=0.005
  h->SetLabelOffset(0.005, "Y"); // D=0.005
  h->SetLabelSize(0.05, "X");
  h->SetLabelSize(0.05, "Y");
  h->SetTitleOffset(0.8, "X");
  h->SetTitleOffset(1.2, "Y");
  h->SetTitleSize(0.06, "X");
  h->SetTitleSize(0.06, "Y");
  h->SetTitle(0);
}

void h2_style(TH2* h) {

  h->SetLabelFont(mFont, "X");   // 42
  h->SetLabelFont(mFont, "Y");   // 42
  h->SetLabelOffset(0.000, "X"); // D=0.005
  h->SetLabelOffset(0.02, "Y");  // D=0.005
  h->SetLabelSize(0.055, "X");
  h->SetLabelSize(0.055, "Y");
  h->SetTitleOffset(0.85, "X");
  h->SetTitleOffset(1.1, "Y");
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

/*********************************************************
Systematics Plotting Macro:

Input variables:
const char* infile = input file name. Must be fv_xyz.root file from
collieIO production step.

int massPoint = mass point for which plots should be generated

const char* channelName = name specifier for plots

const char* varType = variable name for plots (ie, x-axis label)

All systematics for all signal and background sources will be plotted.

********************************************************/

void plotSystematic(const char* infile, int massPoint = -1, const char* channelName = "My Channel",
                    const char* varType = "MyVariable") {
  style();

  if (massPoint == -1) {
    printf("You must specify a masspoint!\n");
    return;
  }

  TObjArray sigSysts;
  TObjArray sigNames;
  TObjArray bkgdSysts;
  TObjArray bkgdNames;

  TFile* fin = new TFile(infile, "READ");

  TList* myKeys = fin->GetListOfKeys();
  char mass[256];
  sprintf(mass, "%d", massPoint);
  for (int i = 0; i < myKeys->GetEntries(); i++) {
    TObjString* ostr = myKeys->At(i);
    TString* astr = (TString*)ostr->GetString();
    TString* str = new TString(astr->Data());
    //    printf("Key: %s\n",str->Data());

    if (str->Contains(mass)) {

      if (str->Contains("Signal")) {

        if (str->Contains("systematic") && str->Contains("Pos")) {

          str->Replace(0, str->Index(":") + 2, "");
          str->ReplaceAll(", Pos", "");
          str->ReplaceAll(" ", "");
          if (sigSysts.FindObject(str->Data()) == NULL)
            sigSysts.Add(new TObjString(str->Data()));

        } else if (str->Contains("Final Variable")) {

          str->Replace(0, 7, "");
          str->ReplaceAll("Final Variable - ", "");
          str->ReplaceAll(mass, "");
          str->ReplaceAll(" ", "");
          if (sigNames.FindObject(str->Data()) == NULL)
            sigNames.Add(new TObjString(str->Data()));
        }
      } else {

        if (str->Contains("systematic") && str->Contains("Pos")) {
          str->Replace(0, str->Index(":") + 2, "");
          str->ReplaceAll(", Pos", "");
          str->ReplaceAll(" ", "");
          if (bkgdSysts.FindObject(str->Data()) == NULL)
            bkgdSysts.Add(new TObjString(str->Data()));

        } else if (str->Contains("Final Variable") && !str->Contains("All") &&
                   !str->Contains("Data")) {

          str->ReplaceAll("Final Variable - ", "");
          str->ReplaceAll(mass, "");
          str->ReplaceAll(" ", "");
          if (bkgdNames.FindObject(str->Data()) == NULL)
            bkgdNames.Add(new TObjString(str->Data()));
        }
      }

    } //contains mass
  } //key entries

  char title[256];
  for (int s = 0; s < sigSysts.GetEntries(); s++) {
    cout << "Signal: " << s + 1 << " of " << sigSysts.GetEntries() << ": " << endl;

    const char* sName = ((TString*)(((TObjString*)(sigSysts.At(s)))->GetString()))->Data();

    for (int b = 0; b < sigNames.GetEntries(); b++) {

      bool isShape = true;

      const char* tName = ((TString*)(((TObjString*)(sigNames.At(b)))->GetString()))->Data();
      sprintf(title, "Shape Signal %s systematic %d: %s, Pos", tName, massPoint, sName);
      TH1D* pos = (TH1D*)fin->Get(title);
      sprintf(title, "Shape Signal %s systematic %d: %s, Neg", tName, massPoint, sName);
      TH1D* neg = (TH1D*)fin->Get(title);

      if (!neg && !pos) {
        sprintf(title, "Signal %s systematic %d: %s, Pos", tName, massPoint, sName);
        pos = (TH1D*)fin->Get(title);
        sprintf(title, "Signal %s systematic %d: %s, Neg", tName, massPoint, sName);
        neg = (TH1D*)fin->Get(title);
        isShape = false;
      }

      if (neg && pos) {
        h1_style(pos);
        h1_style(neg);
        pos->GetXaxis()->SetTitle(varType);
        pos->GetYaxis()->SetTitle("Fractional Uncertainty");
        pos->GetYaxis()->SetTitleOffset(1.1);
        if (isShape)
          sprintf(title, "Signal %s Shape systematic: %s", tName, sName);
        else
          sprintf(title, "Signal %s systematic: %s", tName, sName);
        pos->SetTitle(title);
        neg->Scale(-1.0);

        pos->SetLineWidth(5);
        pos->SetLineColor(4);
        neg->SetLineWidth(5);
        neg->SetLineColor(2);
        sprintf(title, "Sig Canvas: %s %d", sName, b);

        TCanvas* c1 = new TCanvas(title, title);
        can_style(c1);
        c1->SetGridx();
        c1->SetGridy();
        double maxp, minp, maxn, minn, max;
        maxp = pos->GetMaximum();
        minp = fabs(pos->GetMinimum());
        maxn = neg->GetMaximum();
        minn = fabs(neg->GetMinimum());
        max = 0.05;
        for (; maxp > max || minp > max || maxn > max || minn > max; max += 0.05) {
        }
        max *= 1.75;
        pos->SetMaximum(max);
        pos->SetMinimum(-max);
        if (isShape) {
          pos->Draw("hist");
          neg->Draw("samehist");
        } else {
          pos->Draw("hist");
          neg->Draw("samehist");
        }

        TLegend* leg = new TLegend(0.8, 0.44, 1.0, 0.64, NULL, "brNDC");
        leg_style(leg);
        leg->AddEntry(pos, "+1 #sigma");
        leg->AddEntry(neg, "-1 #sigma");
        leg->Draw("same");

        sprintf(title, "%s_%s_%s_%s.eps", channelName, varType, tName, sName);
        c1->Print(title);
      } else {
        printf("Could not find %s systematic %d: %s\n", tName, massPoint, sName);
      }
    }
  }

  for (int s = 0; s < bkgdSysts.GetEntries(); s++) {
    cout << "Background: " << s + 1 << " of " << bkgdSysts.GetEntries() << ": " << endl;

    const char* sName = ((TString*)(((TObjString*)(bkgdSysts.At(s)))->GetString()))->Data();

    for (int b = 0; b < bkgdNames.GetEntries(); b++) {

      bool isShape = true;
      const char* tName = ((TString*)(((TObjString*)(bkgdNames.At(b)))->GetString()))->Data();

      sprintf(title, "Shape %s systematic %d: %s, Pos", tName, massPoint, sName);
      TH1D* pos = (TH1D*)fin->Get(title);
      sprintf(title, "Shape %s systematic %d: %s, Neg", tName, massPoint, sName);
      TH1D* neg = (TH1D*)fin->Get(title);

      if (!neg && !pos) {
        sprintf(title, "%s systematic %d: %s, Pos", tName, massPoint, sName);
        pos = (TH1D*)fin->Get(title);
        sprintf(title, "%s systematic %d: %s, Neg", tName, massPoint, sName);
        neg = (TH1D*)fin->Get(title);
        isShape = false;
      }

      if (neg && pos) {
        h1_style(pos);
        h1_style(neg);
        pos->GetXaxis()->SetTitle(varType);
        pos->GetYaxis()->SetTitle("Fractional Uncertainty");
        pos->GetYaxis()->SetTitleOffset(1.1);
        if (isShape)
          sprintf(title, "Bkgd %s Shape systematic: %s", tName, sName);
        else
          sprintf(title, "Bkgd %s systematic: %s", tName, sName);

        pos->SetTitle(title);
        neg->Scale(-1.0);

        pos->SetLineWidth(5);
        pos->SetLineColor(4);
        neg->SetLineWidth(5);
        neg->SetLineColor(2);
        sprintf(title, "Bkgd Canvas: %s %d", sName, b);

        TCanvas* c1 = new TCanvas(title, title);
        can_style(c1);
        c1->SetGridx();
        c1->SetGridy();
        double maxp, minp, maxn, minn, max;
        maxp = pos->GetMaximum();
        minp = fabs(pos->GetMinimum());
        maxn = neg->GetMaximum();
        minn = fabs(neg->GetMinimum());
        max = 0.05;
        for (; maxp > max || minp > max || maxn > max || minn > max; max += 0.05) {
        }
        max *= 1.75;
        pos->SetMaximum(max);
        pos->SetMinimum(-max);
        if (isShape) {
          pos->Draw("hist");
          neg->Draw("samehist");
        } else {
          pos->Draw("hist");
          neg->Draw("samehist");
        }

        TLegend* leg = new TLegend(0.8, 0.44, 1.0, 0.64, NULL, "brNDC");
        leg_style(leg);
        leg->AddEntry(pos, "+1 #sigma");
        leg->AddEntry(neg, "-1 #sigma");
        leg->Draw("same");

        sprintf(title, "%s_%s_%s_%s.eps", channelName, varType, tName, sName);
        TString outTitle(title);
        outTitle.ReplaceAll(" ", "_");
        c1->Print(outTitle);
      } else {
        printf("Could not find %s systematic %d: %s\n", tName, massPoint, sName);
      }
    }
  }
}
