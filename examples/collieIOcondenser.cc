#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TArrayI.h>
#include <TCanvas.h>

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sys/time.h>
#include <vector>

#include "CollieLoader.hh"
#include "CollieIOFile.hh"
#include "SigBkgdDist.hh"

/*
  CollieIOCondenser takes four user arguments and three "expert" arguments.  Users should feel
  free to play around with modifying the expert arguments, but it should not be necessary in general.

  User Arguments:
  ---------------

  char* outFile:  The name of the output file you wish to create.

  char* outName:  The Collie channel name of the new output file you're going make.

  char* inList:   A text file containing a list of the channels you wish to combine.  The default limit is 200, but this
                  can be changed in the code if one needs to do so.

  char* m:        The mass (or model) point you wish to combine.  The program will currently only condense one point
                  at a time.  The first input channel in the list must contain this point, but the others don't need to
      (they will just be skipped).

  Expert Arguments:
  -----------------

  char* binfact:  The ratio of input to output bins.  This is more of a guideline for the code rather than an explicit
                  rebinning factor, but will result in a ratio (in bins)/(out bins) ~ binfact.  Default = 2.

  char* cf:       The cutoff in log(1+s/b) that marks the definition of low and high S/B regions.  This cutoff is a fraction
                  of the total sum of log(1+s/b).  Bins below this fraction are grouped together, bins above it are not.
      The default = 0.10 (10%).

  char* lG:       The grouping factor for the low S/B regions.  In the region below the log(1+s/b) cutoff, bins are merged
                  until they satisfy the condition that they are all within the grouping factor of the maximum bin.  The
      default is 0.75 (all low S/B bins within 75% of the maximum bin).
*/

void collieCondenser(char* outFile, char* outName, char* inList, char* m, char* binfact, char* cf,
                     char* lG) {

  int mass = atoi(m);
  double binFactor = atof(binfact);
  double cutOff = atof(cf);
  double lowGroup = atof(lG);

  CollieLoader loaders[200];
  string chanNames[200];
  int nld = 0;

  ifstream streamIn(inList);
  if (!streamIn) {
    cout << "Error: Could not open the input list:" << inList << endl;
    return;
  }

  bool ok = true;
  char fname[1024];
  char options[1024];
  while (streamIn >> fname) {

    TFile* ftest = new TFile(fname);
    TList* aList = ftest->GetListOfKeys();
    if (aList->GetEntries() != 1) {
      cout << "Incorrect key length for file: " << fname << endl;
      return;
    }
    chanNames[nld] = aList->At(0)->GetName();

    cout << endl << "Reading file: " << fname << " with channel name: " << chanNames[nld] << endl;

    aList->Delete();
    ftest->Close();
    ftest->Delete();

    sprintf(options, "name='%s'", chanNames[nld].c_str());
    if (!loaders[nld].open(fname, options)) {
      std::cout << "Failed to open " << fname << " using " << options << "!\n";
      ok = false;
    }

    nld++;
  }

  cout << endl << "************************************************" << endl;
  cout << "Collie IO Condenser" << endl;
  cout << nld << " channels to be combined" << endl;
  cout << "************************************************" << endl;

  if (!ok)
    return;

  int len = loaders[0].getNMasspoints();

  if (len <= 0) {
    std::cout << "Cannot handle loader with " << len << "masspoints...\n";
    return;
  }
  std::vector<int> v1(len);
  std::vector<int> v2(len);
  std::vector<int> v3(len);
  loaders[0].getMasspointList(len, v1.data(), v2.data(), v3.data());

  TFile f("channelComboInspection.root", "RECREATE");

  for (int i = 0; i < len; i++) {
    if (v1[i] == mass || mass == -1) {
      cout << endl << "Mass: " << v1[i] << endl;

      SigBkgdDist* sbd = loaders[0].get(v1[i], v2[i], v3[i]);
      cout << "==>Adding " << chanNames[0].c_str() << endl;
      for (int l = 1; l < nld; l++) {
        SigBkgdDist* sbd2 = loaders[l].get(v1[i], v2[i], v3[i]);
        if (sbd2) {
          cout << "==>Adding " << chanNames[l].c_str() << endl;
          sbd->append(*sbd2);
        } else
          cout << "==>Not adding " << chanNames[l].c_str() << endl;
      }

      cout << endl
           << "Original Sig: " << sbd->totSignal() << ", Bkgd: " << sbd->totBkgd()
           << ", Data: " << sbd->totData() << endl;
      sbd->fluctuate();
      sbd->setBaselineModel();
      cout << "N Systematics: " << sbd->getNsyst() << endl;
      cout << "N Original Bins: " << sbd->nbins() << endl;
      int noutbins = sbd->nbins() * 100;

      char title[256];
      sprintf(title, "Original Data %d", v1[i]);
      TH1D* dataDist = new TH1D(title, title, sbd->nbins(), 0, 1);
      sprintf(title, "Original Bkgd%d", v1[i]);
      TH1D* bkgdDist = new TH1D(title, title, sbd->nbins(), 0, 1);
      sprintf(title, "Original Signal %d", v1[i]);
      TH1D* sigDist = new TH1D(title, title, sbd->nbins(), 0, 1);

      for (int b = 0; b < sbd->nbins(); b++) {
        dataDist->SetBinContent(b + 1, sbd->data(b));
        bkgdDist->SetBinContent(b + 1, sbd->bkgd(b));
        sigDist->SetBinContent(b + 1, sbd->signal(b));

        bkgdDist->SetBinError(b + 1, sbd->getBkgdStatErr(b));
        sigDist->SetBinError(b + 1, sbd->getSignalStatErr(b));
      }

      bkgdDist->SetLineColor(2);
      bkgdDist->SetLineWidth(2);
      sprintf(title, "Original Distributions %d", v1[i]);
      TCanvas* c0 = new TCanvas(title, title);
      dataDist->Draw("E0P1");
      bkgdDist->Draw("samehist");
      sigDist->Draw("samehist");
      c0->Draw();
      c0->Write();

      map<int, int> mapD_DLLR;
      map<int, int> mapDLLR_RB;

      ///Find DLLR ranges and maximum
      double dLLR_max = -1000;
      double dLLR_min = 1000;
      double dLLR_tot = 0;
      for (int b = 1; b <= sigDist->GetNbinsX(); b++) {
        sigDist->AddBinContent(b, 1.0e-6);
        bkgdDist->AddBinContent(b, 1.0e-6);
        if ((bkgdDist->GetBinContent(b) == 1.0e-6) && (sigDist->GetBinContent(b) == 1.0e-6))
          continue;

        double sbt = log(1.0 + sigDist->GetBinContent(b) / bkgdDist->GetBinContent(b));

        if (bkgdDist->GetBinContent(b) > 1e-4 && sbt > dLLR_max)
          dLLR_max = sbt;
        if (sbt < dLLR_min)
          dLLR_min = sbt;
        dLLR_tot += sigDist->GetBinContent(b) * sbt;
      }
      printf("Max/Min S/B: %f/%f\nTot dLLR: %f\n", dLLR_max, dLLR_min, dLLR_tot);

      //Reformat into DLLR histograms
      sprintf(title, "Data S/B Dist %d", v1[i]);
      TH1D* dataDLLR = new TH1D(title, title, noutbins, 0, dLLR_max * 1.01);
      sprintf(title, "Bkgd S/B Dist %d", v1[i]);
      TH1D* bkgdDLLR = new TH1D(title, title, noutbins, 0, dLLR_max * 1.01);
      sprintf(title, "Signal S/B Dist %d", v1[i]);
      TH1D* sigDLLR = new TH1D(title, title, noutbins, 0, dLLR_max * 1.01);

      double dLLRThis = 0;
      for (int b = 1; b <= sigDist->GetNbinsX(); b++) {
        if (bkgdDist->GetBinContent(b) <= 1.1e-6 && sigDist->GetBinContent(b) <= 1.1e-6) {
          if (dataDist->GetBinContent(b) > 0)
            assert(false);
          continue;
        }
        dLLRThis = log(1.0 + sigDist->GetBinContent(b) / bkgdDist->GetBinContent(b));

        if (dLLRThis > dLLR_max)
          dLLRThis = dLLR_max * 0.999;

        double sigIn = sigDist->GetBinContent(b) - 1.0e-6;
        double sigErr = sigDist->GetBinError(b);
        if (sigIn < 0)
          sigIn = 0;

        double bkgdIn = bkgdDist->GetBinContent(b) - 1.0e-6;
        double bkgdErr = bkgdDist->GetBinError(b);
        if (bkgdIn < 0)
          bkgdIn = 0;

        if (sigIn <= 0 && bkgdIn <= 0) {
          if (dataDist->GetBinContent(b) > 0)
            assert(false);
          continue;
        }

        int ibin = sigDLLR->FindBin(dLLRThis);

        double errIS = sigDLLR->GetBinError(ibin);
        double errIB = bkgdDLLR->GetBinError(ibin);

        errIS = errIS * errIS + sigErr * sigErr;
        errIB = errIB * errIB + bkgdErr * bkgdErr;

        sigDLLR->AddBinContent(ibin, sigIn);
        bkgdDLLR->AddBinContent(ibin, bkgdIn);
        dataDLLR->AddBinContent(ibin, dataDist->GetBinContent(b));

        sigDLLR->SetBinError(ibin, sqrt(errIS));
        bkgdDLLR->SetBinError(ibin, sqrt(errIB));

        mapD_DLLR[b] = ibin;
      }

      double dLLRTest = 0;
      int cutBin = -1;
      for (int b = 1; b < noutbins; b++) {
        if (bkgdDLLR->GetBinContent(b) <= 0)
          continue;
        double sbt = log(1.0 + sigDLLR->GetBinContent(b) / bkgdDLLR->GetBinContent(b));
        dLLRTest += sigDLLR->GetBinContent(b) * sbt;
        if (dLLRTest > cutOff * dLLR_tot) {
          printf("\nS/B cutoff bin: %d\n", b);
          cutBin = b;
          break;
        }
      }

      bkgdDLLR->SetLineColor(2);
      bkgdDLLR->SetLineWidth(2);
      sprintf(title, "S/B Distributions %d", v1[i]);
      TCanvas* c2 = new TCanvas(title, title);
      dataDLLR->Draw("E0P1");
      bkgdDLLR->Draw("samehist");
      sigDLLR->Draw("samehist");
      c2->Draw();
      c2->Write();

      //How many filled bins do we have?
      int nfBins = 0;
      int hiB = 0;
      for (int b = 1; b <= sigDLLR->GetNbinsX(); b++) {
        if (bkgdDLLR->GetBinContent(b) > 0 || dataDLLR->GetBinContent(b) > 0) {
          nfBins++;
          if (b > cutBin)
            hiB++;
        }
      }
      printf("Non-empty bins: %d\n", nfBins);

      printf("Bkgd Integral:  Low Region = %f, High Region = %f\n", bkgdDLLR->Integral(0, cutBin),
             bkgdDLLR->Integral(cutBin + 1, noutbins));

      int loB = int(dataDist->GetNbinsX() / binFactor) - hiB;
      double lowAVG = lowGroup * bkgdDLLR->Integral(0, cutBin) / (1.0 * loB);
      printf("N Bins: Low Region = %d, High Region = %d, Avg bkgd / Bin Low Region = %f\n", loB,
             hiB, lowAVG);

      int nFinalBins = 0;
      for (int b = 1; b < bkgdDLLR->GetNbinsX(); b++) {
        if (bkgdDLLR->GetBinContent(b) > 0 || dataDLLR->GetBinContent(b) > 0) {

          double inb = bkgdDLLR->GetBinContent(b);
          while (inb < lowAVG && b < cutBin) {
            b++;
            inb += bkgdDLLR->GetBinContent(b);
          }
          nFinalBins++;
        }
      }

      sprintf(title, "Data Final Rebin Dist %d", v1[i]);
      TH1D* dataRB = new TH1D(title, title, nFinalBins, 0, 1);
      sprintf(title, "Bkgd Final Rebin Dist %d", v1[i]);
      TH1D* bkgdRB = new TH1D(title, title, nFinalBins, 0, 1);
      sprintf(title, "signal Final Rebin Dist %d", v1[i]);
      TH1D* sigRB = new TH1D(title, title, nFinalBins, 0, 1);

      dataRB->Sumw2();
      bkgdRB->Sumw2();
      sigRB->Sumw2();

      int inbin = 1;
      for (int b = 1; b < bkgdDLLR->GetNbinsX(); b++) {
        if (bkgdDLLR->GetBinContent(b) > 0 || dataDLLR->GetBinContent(b) > 0) {

          double ins = sigDLLR->GetBinContent(b);
          double inb = bkgdDLLR->GetBinContent(b);
          double ind = dataDLLR->GetBinContent(b);

          double errs = sigDLLR->GetBinError(b);
          double errb = bkgdDLLR->GetBinError(b);
          errs = errs * errs;
          errb = errb * errb;

          int sum = b;
          int istart = b;
          while (inb < lowAVG && sum < cutBin) {
            sum++;
            b++;
            ins += sigDLLR->GetBinContent(sum);
            inb += bkgdDLLR->GetBinContent(sum);
            ind += dataDLLR->GetBinContent(sum);

            errs += sigDLLR->GetBinError(sum) * sigDLLR->GetBinError(sum);
            errb += bkgdDLLR->GetBinError(sum) * bkgdDLLR->GetBinError(sum);
          }

          dataRB->SetBinContent(inbin, ind);
          sigRB->SetBinContent(inbin, ins);
          bkgdRB->SetBinContent(inbin, inb);

          sigRB->SetBinError(inbin, sqrt(errs));
          bkgdRB->SetBinError(inbin, sqrt(errb));

          for (int bb = istart; bb <= sum; bb++)
            mapDLLR_RB[bb] = inbin;
          inbin++;
        }
      }

      printf("N Bins Final: %d\n", inbin - 1);

      bkgdRB->SetLineColor(2);
      bkgdRB->SetLineWidth(2);
      sigRB->SetLineWidth(2);
      sprintf(title, "Final Rebinned Dists %d", v1[i]);
      TCanvas* c3 = new TCanvas(title, title);
      dataRB->Draw("E0P1");
      bkgdRB->Draw("samehist");
      sigRB->Draw("samehist");
      c3->Draw();
      c3->Write();

      double dLLR_final = 0;
      for (int b = 1; b <= bkgdRB->GetNbinsX(); b++) {
        if (bkgdRB->GetBinContent(b) > 0) {
          double sbt = log(1.0 + sigRB->GetBinContent(b) / bkgdRB->GetBinContent(b));
          dLLR_final += sigRB->GetBinContent(b) * sbt;
        }
      }

      printf("\nFinal dLLR total: %f\n", dLLR_final);
      printf("        Orig    Mid RB   Final RB\n");
      printf("Data: %.3f, %.3f, %.3f\n", dataDist->Integral(), dataDLLR->Integral(),
             dataRB->Integral());
      printf("Bkgd: %.3f, %.3f, %.3f\n", bkgdDist->Integral(), bkgdDLLR->Integral(),
             bkgdRB->Integral());
      printf(" Sig:  %.3f,  %.3f,  %.3f\n\n", sigDist->Integral(), sigDLLR->Integral(),
             sigRB->Integral());

      //Now make the systematics histograms...
      const int nSyst = sbd->getNsyst();
      TH1D* bkgdSystP[nSyst];
      TH1D* sigSystP[nSyst];
      TH1D* bkgdSystN[nSyst];
      TH1D* sigSystN[nSyst];
      vector<double> systParams;
      systParams.resize(nSyst);

      for (int s = 0; s < nSyst; s++) {
        for (int a = 0; a < nSyst; a++)
          systParams[a] = 0;

        systParams[s] = 1;
        sbd->setBaselineModel();
        sbd->fluctuate(systParams);
        sprintf(title, "Bkgd Syst Plus %d-%d", s, v1[i]);
        bkgdSystP[s] = new TH1D(title, title, bkgdRB->GetNbinsX(), 0, 1);
        sprintf(title, "Signal Syst Plus %d-%d", s, v1[i]);
        sigSystP[s] = new TH1D(title, title, bkgdRB->GetNbinsX(), 0, 1);

        for (int b = 1; b <= bkgdDist->GetNbinsX(); b++) {
          int dLLRBin = mapD_DLLR.find(b)->second;
          int rbBin = mapDLLR_RB.find(dLLRBin)->second;
          bkgdSystP[s]->AddBinContent(rbBin, sbd->bkgd(b - 1));
          sigSystP[s]->AddBinContent(rbBin, sbd->signal(b - 1));
        }

        systParams[s] = -1;
        sbd->setBaselineModel();
        sbd->fluctuate(systParams);
        sprintf(title, "Bkgd Syst Neg %d-%d", s, v1[i]);
        bkgdSystN[s] = new TH1D(title, title, bkgdRB->GetNbinsX(), 0, 1);
        sprintf(title, "Signal Syst Neg %d-%d", s, v1[i]);
        sigSystN[s] = new TH1D(title, title, bkgdRB->GetNbinsX(), 0, 1);

        for (int b = 1; b <= bkgdDist->GetNbinsX(); b++) {
          int dLLRBin = mapD_DLLR.find(b)->second;
          int rbBin = mapDLLR_RB.find(dLLRBin)->second;
          bkgdSystN[s]->AddBinContent(rbBin, sbd->bkgd(b - 1));
          sigSystN[s]->AddBinContent(rbBin, sbd->signal(b - 1));
        }
      }

      /////////////////////////////////////////
      ///Create IO file with input parameters
      /////////////////////////////////////////
      CollieIOFile* cfile = new CollieIOFile();
      // Specify outputfile and channel name
      cfile->initFile(outFile, outName);

      //Do not try to modify any systematics!
      cfile->setNoviceFlag(false);
      cfile->setSystematicsOverride(true);

      // Define your input histograms
      cfile->setInputHist(0.0, 1.0, bkgdRB->GetNbinsX());
      cfile->setRebin(1);

      //Define backgrounds
      vector<string> bkgdNames;
      bkgdNames.push_back("Background");
      cfile->createChannel(bkgdNames);

      //Backgrounds are passed in via vector
      vector<TH1D*> vbkgd;
      vbkgd.push_back(bkgdRB);

      //Alpha parameters only matter when smoothing is utilized
      //  Input values don't matter if you're not smoothing.
      //  Don't smooth unless you know what you're doing.
      vector<double> valpha;
      valpha.push_back(-1);

      //Each parameter point has a signal histo, data histo, and an array of backgrounds...
      //  Smoothing parameters are also passed in.
      cfile->createMassPoint(mass, dataRB, sigRB, -1, vbkgd, valpha);

      //Add systematics...
      for (int s = 0; s < nSyst; s++) {
        cfile->createShapeSigSystematic(sbd->getSystName(s).c_str(), sigSystP[s], sigSystN[s], mass,
                                        -1, -1, 1.0, 0, 0, 0);
        cfile->createShapeBkgdSystematic(0, sbd->getSystName(s).c_str(), bkgdSystP[s], bkgdSystN[s],
                                         mass, -1, -1, 1.0, 0, 0, 0);
      }

      //store and output channel information
      cfile->storeFile();
    }
  }
  f.Write();
  f.Close();
}

void usage() {
  printf("**************** CollieIOCondenser Usage ****************\n\n");
  printf("./collieIOcondenser.exe <output file> <output channel name> <input file list>\n\n");
  printf("*********************************************************\n");
}

int main(int argc, char* argv[]) {

  timeval a, b;
  gettimeofday(&a, NULL);

  if (argc < 4) {
    usage();
    return 0;
  }

  collieCondenser(argv[1], argv[2], argv[3], argv[4], "1.0", "0.05", "0.50");

  gettimeofday(&b, NULL);
  double deltat = a.tv_sec * 1000.0 + a.tv_usec / 1000.0;
  deltat = b.tv_sec * 1000.0 + b.tv_usec / 1000.0 - deltat;
  printf(" %f sec run time\n", deltat / 1000);
  printf("\n");
  return 0;
}
