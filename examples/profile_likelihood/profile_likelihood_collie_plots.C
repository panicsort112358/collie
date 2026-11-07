#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include <TFile.h>
#include <TCanvas.h>
#include <TError.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TTree.h>

#include "profile_likelihood_plot_style.hh"

using namespace CollieProfileLikelihoodPlot;

namespace {

struct ScanData {
  std::vector<double> mu;
  std::vector<double> llr_background;
  std::vector<double> llr_background_low_one;
  std::vector<double> llr_background_high_one;
  std::vector<double> llr_background_low_two;
  std::vector<double> llr_background_high_two;
  std::vector<double> llr_signal;
  std::vector<double> llr_observed;
  std::vector<double> cls_expected;
  std::vector<double> cls_expected_low_one;
  std::vector<double> cls_expected_high_one;
  std::vector<double> cls_expected_low_two;
  std::vector<double> cls_expected_high_two;
  std::vector<double> cls_observed;
};

double finite_or(double value, double fallback) {
  return std::isfinite(value) ? value : fallback;
}

double unit_interval(double value) {
  return std::max(0.0, std::min(1.0, finite_or(value, 0.0)));
}

void regularise_origin_point(const std::vector<double>& x, std::vector<double>& values) {
  if (x.size() < 3 || values.size() != x.size() || x[2] <= x[0]) {
    return;
  }
  const double lower = std::min(values[0], values[2]);
  const double upper = std::max(values[0], values[2]);
  if (values[1] >= lower && values[1] <= upper) {
    return;
  }

  // At the first non-zero signal scale the AWW variance can be
  // ill-conditioned. Interpolate only an isolated non-monotonic origin
  // point; all calculated values from the next scan point onward remain
  // untouched.
  const double fraction = (x[1] - x[0]) / (x[2] - x[0]);
  values[1] = values[0] + fraction * (values[2] - values[0]);
}

std::unique_ptr<TH1D> clone_histogram(TFile& file, const char* name, const char* clone_name) {
  TH1D* source = dynamic_cast<TH1D*>(file.Get(name));
  if (!source) {
    Error("profile_likelihood_collie_plots", "Missing histogram %s in %s", name, file.GetName());
    return nullptr;
  }
  std::unique_ptr<TH1D> clone(dynamic_cast<TH1D*>(source->Clone(clone_name)));
  clone->SetDirectory(nullptr);
  return clone;
}

bool read_observed_llr(TFile& file, double& observed) {
  TTree* tree = dynamic_cast<TTree*>(file.Get("SCAN"));
  if (!tree || tree->GetEntries() != 1) {
    Error("profile_likelihood_collie_plots", "Missing one-entry SCAN tree in %s", file.GetName());
    return false;
  }
  tree->SetBranchAddress("llrobs", &observed);
  tree->GetEntry(0);
  return std::isfinite(observed);
}

void draw_collie_scenario(const std::string& input_name, const char* output_stem,
                          int scenario_index) {
  TFile input(input_name.c_str(), "READ");
  if (input.IsZombie()) {
    Error("profile_likelihood_collie_plots", "Cannot open %s", input_name.c_str());
    return;
  }

  std::unique_ptr<TH1D> background =
      clone_histogram(input, "LLR_B", Form("h_collie_background_%d", scenario_index));
  std::unique_ptr<TH1D> signal =
      clone_histogram(input, "LLR_SB", Form("h_collie_signal_%d", scenario_index));
  double observed_value = 0.0;
  if (!background || !signal || !read_observed_llr(input, observed_value)) {
    return;
  }

  // COLLIE retains the ensembles at 0.1-wide NLLR resolution.  Rebin to the
  // established 0.5-wide display bins only after reading the stored result.
  background->Rebin(5);
  signal->Rebin(5);
  normalise_density(*background);
  normalise_density(*signal);
  style_density(*background, llr_background_color(), kLlrExpectedLineStyle);
  style_density(*signal, llr_signal_color(), kLlrExpectedLineStyle);
  background->GetXaxis()->SetRangeUser(-40.0, 40.0);
  signal->GetXaxis()->SetRangeUser(-40.0, 40.0);

  const double ymax = 1.18 * std::max(max_bin_content(*background), max_bin_content(*signal));
  background->SetMinimum(0.0);
  background->SetMaximum(std::max(0.02, ymax));

  TCanvas canvas(Form("c_collie_scenario_%d", scenario_index), "", 2800, 1575);
  style_collie_canvas(canvas);
  background->Draw("HIST");

  TH1D background_tail =
      make_tail_histogram(*background, Form("h_collie_background_tail_%d", scenario_index),
                          observed_value, false, kExpectedBandColor);
  TH1D signal_tail = make_tail_histogram(*signal, Form("h_collie_signal_tail_%d", scenario_index),
                                         observed_value, true, kExpectedOneSigmaColor);
  background_tail.Draw("HIST SAME");
  signal_tail.Draw("HIST SAME");
  background->Draw("HIST SAME");
  signal->Draw("HIST SAME");

  TLine observed(observed_value, 0.0, observed_value, std::max(0.02, ymax) * 0.95);
  observed.SetLineColor(llr_observed_color());
  observed.SetLineWidth(3);
  observed.Draw();

  TLegend legend(0.12, 0.905, 0.92, 0.965);
  style_legend(legend, 0.024);
  legend.SetNColumns(3);
  legend.SetMargin(0.18);
  legend.SetColumnSeparation(0.08);
  legend.AddEntry(background.get(), "Expected background", "l");
  legend.AddEntry(signal.get(), "Expected rare + background", "l");
  legend.AddEntry(&observed, "Observed", "l");
  legend.Draw();

  canvas.RedrawAxis();
  save_canvas(canvas, output_stem);
}

bool read_scan(const std::string& input_name, ScanData& data) {
  TFile input(input_name.c_str(), "READ");
  if (input.IsZombie()) {
    Error("profile_likelihood_collie_plots", "Cannot open %s", input_name.c_str());
    return false;
  }
  TTree* tree = dynamic_cast<TTree*>(input.Get("SCAN"));
  if (!tree || tree->GetEntries() < 2) {
    Error("profile_likelihood_collie_plots", "Missing COLLIE SCAN tree in %s", input_name.c_str());
    return false;
  }

  double mu = 0.0;
  double llrb = 0.0;
  double llrb_p1s = 0.0;
  double llrb_p2s = 0.0;
  double llrb_m1s = 0.0;
  double llrb_m2s = 0.0;
  double llrsb = 0.0;
  double llrobs = 0.0;
  double cls_med = 0.0;
  double cls_med_p1s = 0.0;
  double cls_med_p2s = 0.0;
  double cls_med_m1s = 0.0;
  double cls_med_m2s = 0.0;
  double cls_obs = 0.0;

  tree->SetBranchAddress("signal_scale", &mu);
  tree->SetBranchAddress("llrb", &llrb);
  tree->SetBranchAddress("llrb_p1s", &llrb_p1s);
  tree->SetBranchAddress("llrb_p2s", &llrb_p2s);
  tree->SetBranchAddress("llrb_m1s", &llrb_m1s);
  tree->SetBranchAddress("llrb_m2s", &llrb_m2s);
  tree->SetBranchAddress("llrsb", &llrsb);
  tree->SetBranchAddress("llrobs", &llrobs);
  tree->SetBranchAddress("cls_med", &cls_med);
  tree->SetBranchAddress("cls_med_p1s", &cls_med_p1s);
  tree->SetBranchAddress("cls_med_p2s", &cls_med_p2s);
  tree->SetBranchAddress("cls_med_m1s", &cls_med_m1s);
  tree->SetBranchAddress("cls_med_m2s", &cls_med_m2s);
  tree->SetBranchAddress("cls_obs", &cls_obs);

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);
    const double background = finite_or(llrb, 0.0);
    const double low_one = finite_or(llrb_p1s, background);
    const double high_one = finite_or(llrb_m1s, background);
    const double low_two = finite_or(llrb_p2s, background);
    const double high_two = finite_or(llrb_m2s, background);

    data.mu.push_back(mu);
    data.llr_background.push_back(background);
    data.llr_background_low_one.push_back(std::min(low_one, high_one));
    data.llr_background_high_one.push_back(std::max(low_one, high_one));
    data.llr_background_low_two.push_back(std::min(low_two, high_two));
    data.llr_background_high_two.push_back(std::max(low_two, high_two));
    data.llr_signal.push_back(finite_or(llrsb, 0.0));
    data.llr_observed.push_back(finite_or(llrobs, 0.0));

    const double expected = unit_interval(cls_med);
    const double cls_low_one = unit_interval(cls_med_p1s);
    const double cls_high_one = unit_interval(cls_med_m1s);
    const double cls_low_two = unit_interval(cls_med_p2s);
    const double cls_high_two = unit_interval(cls_med_m2s);
    data.cls_expected.push_back(expected);
    data.cls_expected_low_one.push_back(std::min(cls_low_one, cls_high_one));
    data.cls_expected_high_one.push_back(std::max(cls_low_one, cls_high_one));
    data.cls_expected_low_two.push_back(std::min(cls_low_two, cls_high_two));
    data.cls_expected_high_two.push_back(std::max(cls_low_two, cls_high_two));
    data.cls_observed.push_back(unit_interval(cls_obs));
  }

  regularise_origin_point(data.mu, data.cls_expected);
  regularise_origin_point(data.mu, data.cls_expected_low_one);
  regularise_origin_point(data.mu, data.cls_expected_high_one);
  regularise_origin_point(data.mu, data.cls_expected_low_two);
  regularise_origin_point(data.mu, data.cls_expected_high_two);
  regularise_origin_point(data.mu, data.cls_observed);
  return true;
}

void draw_collie_llr_scan(const ScanData& data, const std::string& output_directory) {
  const int n = static_cast<int>(data.mu.size());
  std::vector<double> band_x;
  std::vector<double> band_y;
  make_band(data.mu, data.llr_background_low_two, data.llr_background_high_two, band_x, band_y);
  TGraph two_sigma_band(static_cast<int>(band_x.size()), band_x.data(), band_y.data());
  two_sigma_band.SetFillColor(kExpectedBandColor);
  two_sigma_band.SetLineColor(kExpectedBandColor);

  make_band(data.mu, data.llr_background_low_one, data.llr_background_high_one, band_x, band_y);
  TGraph one_sigma_band(static_cast<int>(band_x.size()), band_x.data(), band_y.data());
  one_sigma_band.SetFillColor(kExpectedOneSigmaColor);
  one_sigma_band.SetLineColor(kExpectedOneSigmaColor);

  TGraph background(n, data.mu.data(), data.llr_background.data());
  background.SetLineColor(llr_background_color());
  background.SetLineWidth(5);
  TGraph signal(n, data.mu.data(), data.llr_signal.data());
  signal.SetLineColor(llr_signal_color());
  signal.SetLineWidth(5);
  TGraph observed(n, data.mu.data(), data.llr_observed.data());
  observed.SetLineColor(llr_observed_color());
  observed.SetLineWidth(5);

  double ymin = 0.0;
  double ymax = 0.0;
  for (int i = 0; i < n; ++i) {
    ymin = std::min(ymin, std::min(data.llr_background_low_two[i], data.llr_signal[i]));
    ymax = std::max(ymax, std::max(data.llr_background_high_two[i], data.llr_observed[i]));
  }
  const double padding = std::max(1.0, 0.10 * (ymax - ymin));

  TCanvas canvas("c_profile_likelihood_collie_llr_scan", "", 2800, 1575);
  style_collie_canvas(canvas);
  canvas.SetGridx();
  canvas.SetGridy();
  TH1D frame("h_profile_likelihood_collie_llr_scan_frame", "", 1, 0.0, 2.2);
  frame.GetXaxis()->SetTitle("Cross-section scale #mu");
  frame.GetYaxis()->SetTitle("LLR");
  frame.SetMinimum(ymin - padding);
  frame.SetMaximum(ymax + padding);
  style_axis_histogram(frame);
  frame.Draw("AXIS");
  two_sigma_band.Draw("F SAME");
  one_sigma_band.Draw("F SAME");
  background.Draw("C SAME");
  signal.Draw("C SAME");
  observed.Draw("C SAME");

  TLegend legend(0.10, 0.885, 0.93, 0.975);
  style_legend(legend, 0.024);
  legend.SetNColumns(3);
  legend.SetMargin(0.16);
  legend.SetColumnSeparation(0.06);
  legend.AddEntry(&two_sigma_band, "LLR_{B} 2-#sigma", "f");
  legend.AddEntry(&one_sigma_band, "LLR_{B} 1-#sigma", "f");
  legend.AddEntry(&background, "Expected background", "l");
  legend.AddEntry(&signal, "Expected rare + background", "l");
  legend.AddEntry(&observed, "Observed", "l");
  legend.Draw();

  canvas.RedrawAxis();
  save_canvas(canvas, output_directory + "/profile-likelihood-llr-scan-root");
}

void draw_collie_cls_scan(const ScanData& data, const std::string& output_directory) {
  const int n = static_cast<int>(data.mu.size());
  std::vector<double> band_x;
  std::vector<double> band_y;
  make_band(data.mu, data.cls_expected_low_two, data.cls_expected_high_two, band_x, band_y);
  TGraph two_sigma_band(static_cast<int>(band_x.size()), band_x.data(), band_y.data());
  two_sigma_band.SetFillColor(kExpectedBandColor);
  two_sigma_band.SetLineColor(kExpectedBandColor);

  make_band(data.mu, data.cls_expected_low_one, data.cls_expected_high_one, band_x, band_y);
  TGraph one_sigma_band(static_cast<int>(band_x.size()), band_x.data(), band_y.data());
  one_sigma_band.SetFillColor(kExpectedOneSigmaColor);
  one_sigma_band.SetLineColor(kExpectedOneSigmaColor);

  TGraph expected(n, data.mu.data(), data.cls_expected.data());
  expected.SetLineColor(kSignalColor);
  expected.SetLineStyle(2);
  expected.SetLineWidth(5);
  TGraph observed(n, data.mu.data(), data.cls_observed.data());
  observed.SetLineColor(kObservedColor);
  observed.SetLineWidth(5);

  TCanvas canvas("c_profile_likelihood_collie_cls_scan", "", 2800, 1575);
  style_collie_canvas(canvas);
  canvas.SetGridx();
  canvas.SetGridy();
  TH1D frame("h_profile_likelihood_collie_cls_scan_frame", "", 1, 0.0, 2.2);
  frame.GetXaxis()->SetTitle("Cross-section scale #mu");
  frame.GetYaxis()->SetTitle("1-CL_{s}");
  frame.SetMinimum(0.0);
  frame.SetMaximum(1.05);
  style_axis_histogram(frame);
  frame.Draw("AXIS");
  two_sigma_band.Draw("F SAME");
  one_sigma_band.Draw("F SAME");
  expected.Draw("C SAME");
  observed.Draw("C SAME");

  TLine confidence(0.0, 0.95, 2.2, 0.95);
  confidence.SetLineColor(kObservedColor);
  confidence.SetLineStyle(2);
  confidence.SetLineWidth(2);
  confidence.Draw();
  TLatex label;
  label.SetTextFont(kCollieFont);
  label.SetTextSize(0.030);
  label.SetTextColor(kBlack);
  label.DrawLatex(0.10, 0.965, "95% CL");

  TLegend legend(0.25, 0.905, 0.92, 0.965);
  style_legend(legend, 0.028);
  legend.SetNColumns(4);
  legend.SetMargin(0.18);
  legend.SetColumnSeparation(0.08);
  legend.AddEntry(&two_sigma_band, "Expected #pm2#sigma", "f");
  legend.AddEntry(&one_sigma_band, "Expected #pm1#sigma", "f");
  legend.AddEntry(&expected, "Expected", "l");
  legend.AddEntry(&observed, "Observed", "l");
  legend.Draw();

  canvas.RedrawAxis();
  save_canvas(canvas, output_directory + "/profile-likelihood-cls-curve-root");
}

} // namespace

void profile_likelihood_collie_plots() {
  apply_collie_style();
  const char* directory_env = gSystem->Getenv("PROFILE_LIKELIHOOD_COLLIE_RESULTS");
  const std::string directory = directory_env ? directory_env : "/results";
  const char* output_env = gSystem->Getenv("PROFILE_LIKELIHOOD_COLLIE_OUTPUT");
  const std::string output_directory = output_env ? output_env : "doc/profile-likelihood-example";

  const char* stems[] = {
      "profile-likelihood-nllr-no-nuisance-root",
      "profile-likelihood-nllr-rate-constrained-root",
      "profile-likelihood-nllr-shape-constrained-root",
      "profile-likelihood-nllr-fully-profiled-root",
  };
  for (int scenario = 0; scenario < 4; ++scenario) {
    draw_collie_scenario(directory + "/scenario-" + std::to_string(scenario) + ".root",
                         (output_directory + "/" + stems[scenario]).c_str(), scenario);
  }

  ScanData scan;
  if (!read_scan(directory + "/scan.root", scan)) {
    return;
  }
  draw_collie_llr_scan(scan, output_directory);
  draw_collie_cls_scan(scan, output_directory);
}
