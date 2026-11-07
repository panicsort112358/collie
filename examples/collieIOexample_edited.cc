#include "CollieIOFile.hh"
#include "TRandom.h"

int main(int argc, char* argv[]) {

  /////////////////////////////////////////
  ///Create IO file with input parameters
  /////////////////////////////////////////
  CollieIOFile* cfile = new CollieIOFile();
  // Specify outputfile and channel name
  cfile->initFile("exampleCollieIOfile.root", "My Channel");

  // Define your input histograms

  //Otherwise, get your input histograms from an external file

  TFile infile("../data/sbnd_0.310GeV_0.1_smearing.root");
  //TH1D* data = (TH1D*)infile.Get("data");
  TH1D* signal = (TH1D*)infile.Get("invmass_signal");
  TH1D* bkgd1 = (TH1D*)infile.Get("invmass_bkg");

  TRandom r(1234);
  double niter = 5e5;

  //We'll make three mass points...
  for (int m = 100; m <= 150; m += 25) {

    bkgd1->Scale(0);
    bkgd2->Scale(0);
    sig->Scale(0);
    data->Scale(0);
    bkgd2SystP->Scale(0);
    bkgd2SystN->Scale(0);
    sigSystP->Scale(0);
    sigSystN->Scale(0);

    // For this example, randomly fill histograms
    for (int i = 0; i < niter; i++) {
      //Flat background
      bkgd1->Fill(r.Rndm(), 150.0 / niter);

      //Fill data with 5% higher bkgd1 rate
      data->Fill(r.Rndm(), 157.5 / niter);

      //Falling background
      bkgd2->Fill(r.Gaus(0, 0.30), 600.0 / niter);

      //Shape-dependent bkgd systematic
      bkgd2SystP->Fill(r.Gaus(0, 0.305), 605.0 / niter);
      bkgd2SystN->Fill(r.Gaus(0, 0.295), 595.0 / niter);

      //Fill data with 2% lower rate and 1/2 sigma morphed shape
      data->Fill(r.Gaus(0, 0.2975), 588.0 / niter);

      //Peaked signal
      sig->Fill(r.Gaus(1, 0.15), 15 / niter);

      //Shape-dependent signal systematic
      sigSystP->Fill(r.Gaus(1, 0.145), 15 / niter);
      sigSystN->Fill(r.Gaus(1, 0.155), 15 / niter);
    }

    //Backgrounds are passed in via vector
    vector<TH1D*> vbkgd;
    vbkgd.push_back(bkgd1);
    vbkgd.push_back(bkgd2);

    //Alpha parameters only matter when smoothing is utilized
    //  Input values don't matter if you're not smoothing.
    //  Don't smooth unless you know what you're doing.
    vector<double> valpha;
    valpha.push_back(-1);
    valpha.push_back(-1);

    ///Use this tool to allow collie to generate a low-stats safe
    // binning for your histograms.  Histogram binning should be
    // larger than your desired final number of bins.  For more
    // details see CollieIOfile.h
    // Should be done for each mass point in the file, but with
    // the same number of output bins (specified above).  IE, if
    // you have a mass point loop, put this inside the loop.
    //
    //         TH1D* btotal = (TH1D*)bkgd1->Clone("btotal");
    //         btotal->Add(bkgd2);
    //         cfile->generateBinMap(btotal,sig,"MVA");
    //         delete btotal; btotal = NULL;

    //Each parameter point has a signal histo, data histo, and an array of backgrounds...
    //  Smoothing parameters are also passed in.
    cfile->createMassPoint(m, data, sig, -1, vbkgd, valpha);

    // If you have more than one mass point, you may choose to interpolate on some parameter grid
    //cfile->interpolateMassGrid(5,100,110);

    // Add systematics...either flat or by shape (ie, function of final variable)
    //   if by shape, must supply a histogram of the values in percent(%) fluctuations...
    //   Signal requires no index, but backgrounds must be specifically indexed (0->N bkgds)
    //   Read the instructions in collie/io/include/CollieIOFile.hh if you're in doubt
    cfile->createFlatSigSystematic("Lumi", 0.06, 0.06, m);
    cfile->createShapeSigSystematic("SigShape", sigSystP, sigSystN, m);

    cfile->createFlatBkgdSystematic(0, "Lumi", 0.06, 0.06, m);
    cfile->createFlatBkgdSystematic(1, "Lumi", 0.06, 0.06, m);

    cfile->createFlatSigSystematic("Eff", 0.10, 0.10, m);
    cfile->createFlatBkgdSystematic(0, "Eff", 0.10, 0.10, m);
    cfile->createFlatBkgdSystematic(1, "Eff", 0.10, 0.10, m);

    cfile->createFlatBkgdSystematic(0, "Xsec1", 0.15, 0.15, m);
    cfile->createFlatBkgdSystematic(1, "Xsec2", 0.15, 0.15, m);
    cfile->createShapeBkgdSystematic(1, "Bkgd2Shape", bkgd2SystP, bkgd2SystN, m);

    // Example of systematics input as histograms, can be flat or function of final variable
    //==>Use this method if you're inputing fractional shape systematics
    //    TH1D* systP = (TH1D*)infile.Get("signal_Systematic_positive");
    //    TH1D* systN = (TH1D*)infile.Get("signal_Systematic_negative");
    //    cfile->createSigSystematic("ShapeSyst",systP,systN,m);

    //==>Use this method if you're inputing a different shape template
    //    systP = (TH1D*)infile.Get("bkgd_BkgdShape_positive");
    //    systN = (TH1D*)infile.Get("bkgd_BkgdShape_negative");
    //    cfile->createShapeBkgdSystematic(0,"BkgdShape",systP,systN,m);
    //    cfile->createShapeBkgdSystematic(1,"BkgdShape",systP,systN,m);

    //  ==>Option to remove prior constraint on systematic uncertainty PDF.
    //     Floating makes a parameter a free parameter in the fit.
    //
    //  cfile->setBkgdFloatFlag(0,"Eff",true,m);
    //  cfile->setBkgdFloatFlag(1,"Eff",true,m);
    //  cfile->setSigFloatFlag("Eff",true,m);

    //  ==>For large uncertainties (eg, >30%) use a log-normal PDF to avoid
    //     problems in Gaussian PDF modeling.
    //
    //  cfile->setLogNormalFlag("Eff",true,m);
  }
  ///store and output channel information
  cfile->storeFile();
}
