#ifndef COLLIE_PROFILE_LIKELIHOOD_PLOT_STYLE_HH
#define COLLIE_PROFILE_LIKELIHOOD_PLOT_STYLE_HH

#include <algorithm>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TColor.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

namespace CollieProfileLikelihoodPlot {

inline constexpr int kCollieFont = 42;
inline constexpr int kSignalColor = 2;
inline constexpr int kObservedColor = 1;
inline constexpr int kExpectedBandColor = 5;
inline constexpr int kExpectedOneSigmaColor = 3;
inline constexpr int kLlrExpectedLineStyle = 1;

inline int llr_background_color() {
  return TColor::GetColor("#0000ff");
}

inline int llr_signal_color() {
  return TColor::GetColor("#ff0000");
}

inline int llr_observed_color() {
  return TColor::GetColor("#000000");
}

inline void apply_collie_style() {
  const int grid_color = TColor::GetColor("#e5e5e5");

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(kCollieFont);
  gStyle->SetEndErrorSize(0.0);
  gStyle->SetPadLeftMargin(0.095);
  gStyle->SetPadRightMargin(0.095);
  gStyle->SetPadTopMargin(0.125);
  gStyle->SetPadBottomMargin(0.115);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(grid_color);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetLineStyleString(11, "20 12");
  gStyle->SetLineStyleString(12, "84 28");
  gStyle->SetLineStyleString(13, "10 22");
  gStyle->SetLineStyleString(14, "54 24 12 24");
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(255);
  gStyle->SetLabelSize(0.035, "XYZ");
  gStyle->SetTitleSize(0.040, "XYZ");
  gStyle->SetTitleOffset(1.05, "X");
  gStyle->SetTitleOffset(1.10, "Y");
  gStyle->SetTitleOffset(1.15, "Z");
  gStyle->SetNdivisions(510, "X");
  gStyle->SetNdivisions(506, "Y");

  gStyle->SetTextFont(kCollieFont);
  gStyle->SetLabelFont(kCollieFont, "XYZ");
  gStyle->SetTitleFont(kCollieFont, "XYZ");
  gStyle->SetTitleFont(kCollieFont);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetErrorX(0.0);
  gStyle->SetPaintTextFormat(".2f");
}

inline void style_collie_canvas(TCanvas& canvas, bool right_axis = false) {
  canvas.SetCanvasSize(2800, 1575);
  canvas.SetFillColor(0);
  canvas.SetFrameFillColor(0);
  canvas.SetLeftMargin(0.095);
  canvas.SetRightMargin(right_axis ? 0.135 : 0.095);
  canvas.SetBottomMargin(0.115);
  canvas.SetTopMargin(0.125);
  canvas.SetTicks(1, right_axis ? 0 : 1);
  canvas.SetGrid(0, 0);
}

inline void ensure_output_directory(const std::string& output_stem) {
  const auto slash = output_stem.find_last_of('/');
  if (slash != std::string::npos) {
    gSystem->mkdir(output_stem.substr(0, slash).c_str(), kTRUE);
  }
}

inline void save_canvas(TCanvas& canvas, const std::string& output_stem) {
  ensure_output_directory(output_stem);
  canvas.SaveAs((output_stem + ".pdf").c_str());
  canvas.SaveAs((output_stem + ".png").c_str());
}

inline void style_axis_histogram(TH1D& hist) {
  hist.SetStats(kFALSE);
  hist.SetTitle("");
  hist.GetXaxis()->SetTitleFont(kCollieFont);
  hist.GetYaxis()->SetTitleFont(kCollieFont);
  hist.GetXaxis()->SetLabelFont(kCollieFont);
  hist.GetYaxis()->SetLabelFont(kCollieFont);
  hist.GetXaxis()->SetTitleSize(0.040);
  hist.GetYaxis()->SetTitleSize(0.040);
  hist.GetXaxis()->SetLabelSize(0.035);
  hist.GetYaxis()->SetLabelSize(0.035);
  hist.GetXaxis()->SetLabelOffset(0.000);
  hist.GetYaxis()->SetLabelOffset(0.005);
  hist.GetXaxis()->SetTitleOffset(1.05);
  hist.GetYaxis()->SetTitleOffset(1.10);
  hist.GetXaxis()->SetNdivisions(510);
  hist.GetYaxis()->SetNdivisions(506);
}

inline void style_density(TH1D& hist, int color, int line_style = 1) {
  hist.SetStats(kFALSE);
  hist.SetTitle("");
  hist.SetFillStyle(0);
  hist.SetLineColor(color);
  hist.SetLineStyle(line_style);
  hist.SetLineWidth(2);
  hist.GetXaxis()->SetTitle("NLLR");
  hist.GetYaxis()->SetTitle("Probability Density");
  style_axis_histogram(hist);
}

inline void style_legend(TLegend& legend, double text_size = 0.028) {
  legend.SetTextFont(kCollieFont);
  legend.SetTextSize(text_size);
  legend.SetFillColor(0);
  legend.SetFillStyle(0);
  legend.SetLineColor(0);
  legend.SetBorderSize(0);
}

inline void normalise_density(TH1D& hist) {
  const double integral = hist.Integral("width");
  if (integral > 0.0) {
    hist.Scale(1.0 / integral);
  }
}

inline double max_bin_content(const TH1D& hist) {
  double value = 0.0;
  for (int bin = 1; bin <= hist.GetNbinsX(); ++bin) {
    value = std::max(value, hist.GetBinContent(bin));
  }
  return value;
}

inline TH1D make_tail_histogram(const TH1D& source, const char* name, double threshold,
                                bool keep_right_tail, int fill_color) {
  TH1D tail(source);
  tail.SetName(name);
  tail.SetFillColor(fill_color);
  tail.SetFillStyle(1001);
  tail.SetLineColor(fill_color);
  tail.SetLineWidth(0);
  for (int bin = 1; bin <= tail.GetNbinsX(); ++bin) {
    const double center = tail.GetBinCenter(bin);
    const bool keep = keep_right_tail ? center >= threshold : center <= threshold;
    if (!keep) {
      tail.SetBinContent(bin, 0.0);
    }
  }
  return tail;
}

inline void make_band(const std::vector<double>& x, const std::vector<double>& low,
                      const std::vector<double>& high, std::vector<double>& band_x,
                      std::vector<double>& band_y) {
  const int n = static_cast<int>(x.size());
  band_x.assign(2 * n + 1, 0.0);
  band_y.assign(2 * n + 1, 0.0);
  for (int i = 0; i < n; ++i) {
    band_x[i] = x[i];
    band_y[i] = low[i];
    band_x[2 * n - 1 - i] = x[i];
    band_y[2 * n - 1 - i] = high[i];
  }
  band_x[2 * n] = band_x[0];
  band_y[2 * n] = band_y[0];
}

} // namespace CollieProfileLikelihoodPlot

#endif
