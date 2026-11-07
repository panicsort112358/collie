#include "CollieIOFile.hh"

#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TObject.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
  if (argc < 2 || argc > 3)
    return 2;
  const bool nonfinite = argc == 3 && std::string(argv[2]) == "--nonfinite";

  TH1::AddDirectory(false);
  CollieIOFile output;
  output.initFile(argv[1], "test_channel");
  output.setInputHist(0.0, 2.0, 2);
  output.setRebin(1);
  output.createChannel(std::vector<std::string>{"background"});

  for (const int point : {125, 130}) {
    TH1D data(("data_" + std::to_string(point)).c_str(), "", 2, 0.0, 2.0);
    TH1D signal(("signal_" + std::to_string(point)).c_str(), "", 2, 0.0, 2.0);
    TH1D background(("background_" + std::to_string(point)).c_str(), "", 2, 0.0, 2.0);
    data.Sumw2();
    signal.Sumw2();
    background.Sumw2();

    data.SetBinContent(1, 11.0);
    data.SetBinContent(2, 8.0);
    signal.SetBinContent(1, point == 125 ? 2.0 : 1.6);
    signal.SetBinContent(2, point == 125 ? 1.0 : 1.3);
    background.SetBinContent(1, 9.0);
    background.SetBinContent(2, 7.0);
    for (int bin = 1; bin <= 2; ++bin) {
      data.SetBinError(bin, std::sqrt(data.GetBinContent(bin)));
      signal.SetBinError(bin, 0.1);
      background.SetBinError(bin, 0.2);
    }

    std::vector<TH1D*> backgrounds{&background};
    std::vector<double> alphas{-1.0};
    output.createMassPoint(point, &data, &signal, -1.0, backgrounds, alphas);

    TH1D up(("up_" + std::to_string(point)).c_str(), "", 2, 0.0, 2.0);
    TH1D down(("down_" + std::to_string(point)).c_str(), "", 2, 0.0, 2.0);
    up.SetBinContent(1, 0.10);
    up.SetBinContent(2, 0.05);
    down.SetBinContent(1, 0.08);
    down.SetBinContent(2, 0.04);
    output.createSigSystematic("calibration", &up, &down, point);
    output.createFlatBkgdSystematic(0, "calibration", 0.03, 0.03, point);
  }

  output.storeFile();

  // CollieIOFile correctly rejects non-finite inputs. Corrupt a persisted
  // systematic deliberately so the validator's read-side defence is tested.
  if (nonfinite) {
    TFile file(argv[1], "UPDATE");
    CollieChannel* channel = CollieChannel::loadChannel("test_channel", &file);
    CollieMasspoint* point = channel == nullptr ? nullptr : channel->getMasspoint(125);
    const CollieDistribution* distribution = point == nullptr ? nullptr : point->getSignalDist(0);
    TH1D* positive = distribution == nullptr ? nullptr : distribution->getPositiveSyst(0);
    TDirectory* directory = file.GetDirectory("test_channel");
    if (file.IsZombie() || point == nullptr || positive == nullptr || directory == nullptr) {
      return 1;
    }
    positive->SetBinContent(1, std::numeric_limits<double>::quiet_NaN());
    directory->WriteObject(point, point->GetName(), "Overwrite");
    file.Close();
  }
  return 0;
}
