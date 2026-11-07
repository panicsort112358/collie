#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TList.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>

#include <CLfast.hh>
#include <CLfit2.hh>
#include <CLpoint.hh>
#include <CollieIOFile.hh>
#include <CollieLoader.hh>
#include <RunContext.hh>
#include <SigBkgdDist.hh>

#include "profile_likelihood_toy_model.hh"

namespace {

constexpr int kMassPoint = 1000;
constexpr int kLlrBins = 1200;
constexpr double kLlrMin = -60.0;
constexpr double kLlrMax = 60.0;

struct Scenario {
  const char* name;
  bool rate;
  bool shape;
};

const Scenario kScenarios[] = {
    {"no-nuisance", false, false},
    {"rate-constrained", true, false},
    {"shape-constrained", false, true},
    {"fully-profiled", true, true},
};

void usage(const char* executable) {
  std::cerr << "Usage:\n"
            << "  " << executable << " create <collie-input.root>\n"
            << "  " << executable
            << " scenario <collie-input.root> <result.root> <0..3> <effort> <seed>\n"
            << "  " << executable << " scan <collie-input.root> <result.root>\n"
            << "Effort: test, veryveryfast, veryfast, fast, standard, fine, veryfine\n";
}

void set_bin(TH1D& hist, int bin, double value) {
  hist.SetBinContent(bin, value);
  hist.SetBinError(bin, 0.0);
}

int create_input(const char* output_name) {
  const ProfileLikelihoodToy::ToyModel model;
  const auto& templates = model.templates();

  TH1D data("collie_toy_data", "", ProfileLikelihoodToy::kNBins, ProfileLikelihoodToy::kXMin,
            ProfileLikelihoodToy::kXMax);
  TH1D signal("collie_toy_signal", "", ProfileLikelihoodToy::kNBins, ProfileLikelihoodToy::kXMin,
              ProfileLikelihoodToy::kXMax);
  TH1D background1("collie_toy_detector_background", "", ProfileLikelihoodToy::kNBins,
                   ProfileLikelihoodToy::kXMin, ProfileLikelihoodToy::kXMax);
  TH1D background2("collie_toy_conventional_neutrino", "", ProfileLikelihoodToy::kNBins,
                   ProfileLikelihoodToy::kXMin, ProfileLikelihoodToy::kXMax);
  TH1D signal_shape_positive("collie_toy_signal_shape_positive", "", ProfileLikelihoodToy::kNBins,
                             ProfileLikelihoodToy::kXMin, ProfileLikelihoodToy::kXMax);
  TH1D signal_shape_negative("collie_toy_signal_shape_negative", "", ProfileLikelihoodToy::kNBins,
                             ProfileLikelihoodToy::kXMin, ProfileLikelihoodToy::kXMax);
  TH1D background_shape_positive("collie_toy_background_shape_positive", "",
                                 ProfileLikelihoodToy::kNBins, ProfileLikelihoodToy::kXMin,
                                 ProfileLikelihoodToy::kXMax);
  TH1D background_shape_negative("collie_toy_background_shape_negative", "",
                                 ProfileLikelihoodToy::kNBins, ProfileLikelihoodToy::kXMin,
                                 ProfileLikelihoodToy::kXMax);

  for (int bin = 1; bin <= ProfileLikelihoodToy::kNBins; ++bin) {
    const int index = bin - 1;
    set_bin(data, bin, model.data()[index]);
    set_bin(signal, bin, templates.signal[index]);
    set_bin(background1, bin, templates.background1[index]);
    set_bin(background2, bin, templates.background2[index]);

    // COLLIE expects signed fractional excursions. Supplying the same
    // signed shape for its positive and negative widths gives the
    // symmetric +theta/-theta morph used by the thesis toy.
    signal_shape_positive.SetBinContent(bin, templates.signal_shape[index]);
    signal_shape_negative.SetBinContent(bin, templates.signal_shape[index]);
    background_shape_positive.SetBinContent(bin, templates.background2_shape[index]);
    background_shape_negative.SetBinContent(bin, templates.background2_shape[index]);
  }

  CollieIOFile output;
  output.initFile(output_name, "ProfileLikelihoodToy");
  output.setNoviceFlag(false);
  output.setSystematicsOverride(true);
  output.setInputHist(ProfileLikelihoodToy::kXMin, ProfileLikelihoodToy::kXMax,
                      ProfileLikelihoodToy::kNBins);
  output.setRebin(1);
  output.setSmooth(false);

  std::vector<std::string> background_names = {
      "DetectorBackground",
      "ConventionalNeutrino",
  };
  output.createChannel(background_names);

  std::vector<TH1D*> backgrounds = {&background1, &background2};
  std::vector<double> smoothing = {-1.0, -1.0};
  output.createMassPoint(kMassPoint, &data, &signal, -1.0, backgrounds, smoothing);

  // Correlated rate nuisances. Reusing a name correlates the corresponding
  // fractional response across signal and background components.
  output.createFlatSigSystematic("FluxNorm", 0.060, 0.060, kMassPoint);
  output.createFlatBkgdSystematic(0, "FluxNorm", 0.060, 0.060, kMassPoint);
  output.createFlatBkgdSystematic(1, "FluxNorm", 0.060, 0.060, kMassPoint);

  output.createFlatSigSystematic("Selection", 0.080, 0.080, kMassPoint);
  output.createFlatBkgdSystematic(0, "Selection", 0.045, 0.045, kMassPoint);
  output.createFlatBkgdSystematic(1, "Selection", 0.045, 0.045, kMassPoint);

  output.createFlatBkgdSystematic(0, "DetectorBackgroundRate", 0.140, 0.140, kMassPoint);
  output.createFlatBkgdSystematic(1, "ConventionalNeutrinoRate", 0.125, 0.125, kMassPoint);

  output.createSigSystematic("SignalShape", &signal_shape_positive, &signal_shape_negative,
                             kMassPoint);
  output.createBkgdSystematic(1, "ConventionalNeutrinoShape", &background_shape_positive,
                              &background_shape_negative, kMassPoint);

  output.storeFile();
  std::cout << "Created COLLIE input " << output_name << std::endl;
  return 0;
}

std::unique_ptr<SigBkgdDist> load_model(const char* input_name) {
  TFile input(input_name, "READ");
  if (input.IsZombie() || !input.GetListOfKeys() || input.GetListOfKeys()->GetSize() == 0) {
    std::cerr << "Cannot read COLLIE input " << input_name << std::endl;
    return nullptr;
  }
  const std::string channel_name = input.GetListOfKeys()->At(0)->GetName();
  input.Close();

  CollieLoader loader;
  const std::string options = "name='" + channel_name + "'";
  if (!loader.open(input_name, options.c_str())) {
    std::cerr << "COLLIE failed to open " << input_name << std::endl;
    return nullptr;
  }

  SigBkgdDist* model = loader.get(kMassPoint, -1, -1);
  if (!model) {
    std::cerr << "COLLIE mass point " << kMassPoint << " is missing" << std::endl;
    return nullptr;
  }
  return std::unique_ptr<SigBkgdDist>(model);
}

bool is_rate_systematic(const std::string& name) {
  return name == "FluxNorm" || name == "Selection" || name == "DetectorBackgroundRate" ||
         name == "ConventionalNeutrinoRate";
}

bool is_shape_systematic(const std::string& name) {
  return name == "SignalShape" || name == "ConventionalNeutrinoShape";
}

void apply_scenario(SigBkgdDist& model, int scenario_index) {
  const Scenario& scenario = kScenarios[scenario_index];
  for (int i = 0; i < model.getNsyst(); ++i) {
    const std::string name = model.getSystName(i);
    const bool enabled = (scenario.rate && is_rate_systematic(name)) ||
                         (scenario.shape && is_shape_systematic(name));
    model.setSystFitFlag(i, enabled);
  }
  model.setSignalScale(1.0);
  model.useHistoStats(false);
}

int effort_code(const std::string& effort) {
  if (effort == "test")
    return CLcompute::LEVEL_TEST;
  if (effort == "veryveryfast")
    return CLcompute::LEVEL_VERYVERYFAST;
  if (effort == "veryfast")
    return CLcompute::LEVEL_VERYFAST;
  if (effort == "fast")
    return CLcompute::LEVEL_FAST;
  if (effort == "standard")
    return CLcompute::LEVEL_STANDARD;
  if (effort == "fine")
    return CLcompute::LEVEL_FINE;
  if (effort == "veryfine")
    return CLcompute::LEVEL_VERYFINE;
  return 999999;
}

int requested_toys(const std::string& effort) {
  if (effort == "test")
    return 100;
  if (effort == "veryveryfast")
    return 5000;
  if (effort == "veryfast")
    return 15000;
  if (effort == "fast")
    return 25000;
  if (effort == "standard")
    return 50000;
  if (effort == "fine")
    return 100000;
  if (effort == "veryfine")
    return 200000;
  return 0;
}

int run_scenario(const char* input_name, const char* output_name, int scenario_index,
                 const std::string& effort, unsigned int seed) {
  if (scenario_index < 0 || scenario_index >= 4 || effort_code(effort) == 999999) {
    std::cerr << "Invalid scenario or effort" << std::endl;
    return 1;
  }

  // The modern run context replaces COLLIE's historical wall-clock seed.
  // One worker matches the serial thesis reproduction exactly and keeps
  // every generated ensemble deterministic.
  collie::configureGlobalRun(seed, 1);
  std::unique_ptr<SigBkgdDist> model = load_model(input_name);
  if (!model)
    return 1;
  apply_scenario(*model, scenario_index);

  TFile output(output_name, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Cannot create " << output_name << std::endl;
    return 1;
  }
  TTree scan("SCAN", "COLLIE profile-likelihood scenario");
  CLpoint result;
  result.branch(&scan);
  result.reset(scenario_index);

  bool ok = false;
  TH1D* background = nullptr;
  TH1D* signal = nullptr;
  const int level = effort_code(effort);
  if (scenario_index == 0) {
    std::unique_ptr<CLfast> calculator(new CLfast());
    calculator->setNoviceFlag(false);
    calculator->useHistoStats(false);
    ok = calculator->calculateCLs(*model, result, level);
    background = calculator->getLLRdist_b("LLR_B", kLlrBins, kLlrMin, kLlrMax);
    signal = calculator->getLLRdist_sb("LLR_SB", kLlrBins, kLlrMin, kLlrMax);
  } else {
    std::unique_ptr<CLfit2> calculator(new CLfit2());
    calculator->setNoviceFlag(false);
    calculator->useHistoStats(false);
    ok = calculator->calculateCLs(*model, result, level);
    background = calculator->getLLRdist_b("LLR_B", kLlrBins, kLlrMin, kLlrMax);
    signal = calculator->getLLRdist_sb("LLR_SB", kLlrBins, kLlrMin, kLlrMax);
  }

  if (!ok || !background || !signal) {
    std::cerr << "COLLIE failed for scenario " << scenario_index << std::endl;
    return 1;
  }
  const int requested_trials = requested_toys(effort);
  if (result.nTrials_exp < requested_trials) {
    std::cerr << "COLLIE retained only " << result.nTrials_exp << " of " << requested_trials
              << " requested trials for scenario " << scenario_index << std::endl;
    return 1;
  }

  scan.Fill();
  scan.Write();
  background->Write();
  signal->Write();
  TNamed engine("collie_engine", scenario_index == 0 ? "CLfast" : "CLfit2");
  engine.Write();
  TNamed scenario_name("scenario", kScenarios[scenario_index].name);
  scenario_name.Write();
  TNamed effort_name("effort", effort.c_str());
  effort_name.Write();
  TParameter<int> requested("requested_toys_per_hypothesis", requested_trials);
  requested.Write();
  TParameter<int> seed_parameter("seed", static_cast<int>(seed));
  seed_parameter.Write();
  output.Write();
  output.Close();

  std::cout << "Finished " << kScenarios[scenario_index].name << " with " << result.nTrials_exp
            << " stored COLLIE trials" << std::endl;
  return 0;
}

void set_zero_signal_result(CLpoint& result) {
  result.cls_obs = 0.0;
  result.clb_obs = 0.5;
  result.clsb_obs = 0.5;
  result.clb_med = 0.5;
  result.clsb_med = 0.5;
  result.cls_med = 0.0;
  result.cls_med_p2s = 0.0;
  result.cls_med_p1s = 0.0;
  result.cls_med_m1s = 0.0;
  result.cls_med_m2s = 0.0;
  result.llrobs = 0.0;
  result.llrb = 0.0;
  result.llrb_p2s = 0.0;
  result.llrb_p1s = 0.0;
  result.llrb_m1s = 0.0;
  result.llrb_m2s = 0.0;
  result.llrsb = 0.0;
  result.llrsb_p2s = 0.0;
  result.llrsb_p1s = 0.0;
  result.llrsb_m1s = 0.0;
  result.llrsb_m2s = 0.0;
  result.nTrials_exp = 1;
  result.nTrials_obs = 1;
}

int run_scan(const char* input_name, const char* output_name) {
  collie::configureGlobalRun(314159, 1);
  std::unique_ptr<SigBkgdDist> model = load_model(input_name);
  if (!model)
    return 1;
  apply_scenario(*model, 3);

  TFile output(output_name, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Cannot create " << output_name << std::endl;
    return 1;
  }
  TTree scan("SCAN", "COLLIE Asimov-Wilks-Wald signal-scale scan");
  CLpoint result;
  result.branch(&scan);

  std::unique_ptr<CLfit2> calculator(new CLfit2());
  calculator->setNoviceFlag(false);
  calculator->useHistoStats(false);
  calculator->useAWWapproximation(true);

  constexpr int n_steps = 110;
  constexpr double step = 0.02;
  for (int i = 0; i <= n_steps; ++i) {
    const double mu = step * i;
    result.reset(i);
    result.signal_scale = mu;
    if (i == 0) {
      set_zero_signal_result(result);
    } else {
      model->setSignalScale(mu);
      if (!calculator->calculateCLs(*model, result, CLcompute::LEVEL_NONE)) {
        std::cerr << "COLLIE scan failed at mu=" << mu << std::endl;
        return 1;
      }
      result.signal_scale = mu;
    }
    scan.Fill();
  }

  scan.Write();
  TNamed engine("collie_engine", "CLfit2 AWW approximation");
  engine.Write();
  output.Write();
  output.Close();
  std::cout << "Finished COLLIE scan with " << (n_steps + 1) << " signal-scale points" << std::endl;
  return 0;
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    usage(argv[0]);
    return 1;
  }

  const std::string mode = argv[1];
  if (mode == "create" && argc == 3) {
    return create_input(argv[2]);
  }
  if (mode == "scenario" && argc == 7) {
    return run_scenario(argv[2], argv[3], std::atoi(argv[4]), argv[5],
                        static_cast<unsigned int>(std::strtoul(argv[6], nullptr, 10)));
  }
  if (mode == "scan" && argc == 4) {
    return run_scan(argv[2], argv[3]);
  }

  usage(argv[0]);
  return 1;
}
