#include <TFile.h>
#include <TTree.h>
#include <CrossSectionLimit.hh>
#include <CollieLoader.hh>
#include <FitTest.hh>
#include <CLfast.hh>
#include <CLsyst.hh>
#include <CLfit.hh>
#include <CLfit2.hh>
#include <sys/time.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

void calcLimit(const string& outFile, const string& inList, const vector<string>& inFiles, int mass,
               int prec, int type, int sigma) {

  /*
  int mass = -1;
  int prec = atoi(p);
  int sigma = atoi(s);
  int type = atoi(tt);
  */

  /*
  printf("\n\ncollieLimitCalc.exe\n");
  printf("Input List: %s\n",inList);
  printf("Output File: %s\n",outFile);
  printf("Mass Point: %d\n",mass);
  */

  CollieLoader loaders[300];
  string chanNames[300];
  string fileNames[300];
  int nld = 0;
  bool ok = true;
  if (!inList.empty()) {
    ifstream streamIn(inList.c_str());
    if (!streamIn) {
      cout << "Error: Could not open " << inList << endl;
      return;
    }

    char fname[1024];
    char options[1024];
    while (streamIn >> fname) {
      cout << "Reading: " << fname << endl;

      TFile* ftest = new TFile(fname);
      TList* aList = ftest->GetListOfKeys();
      if (aList->GetEntries() < 1) {
        printf("Incorrect key length for %s\n", fname);
        return;
      }
      chanNames[nld] = aList->At(0)->GetName();
      TString name(fname);
      fileNames[nld] = name.Data();
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
  }
  for (vector<string>::const_iterator f = inFiles.begin(); f != inFiles.end(); ++f) {
    const char* fname = f->c_str();
    char options[1024];

    cout << "Reading: " << fname << endl;

    TFile* ftest = new TFile(fname);
    TList* aList = ftest->GetListOfKeys();
    if (aList->GetEntries() < 1) {
      printf("Incorrect key length for %s\n", fname);
      return;
    }
    chanNames[nld] = aList->At(0)->GetName();
    TString name(fname);
    fileNames[nld] = name.Data();
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

  /*
  printf("\n************************************************\n");
  printf("Collie Example Limit Calculation\n");
  printf("%d channel(s) available\n",nld);
  printf("************************************************\n");
  */

  if (!ok)
    return;

  ///Create an output container for the results you're about to calculate
  TFile f(outFile.c_str(), "RECREATE");
  TTree t("SCAN", "SCAN");
  CLpoint clresults;
  clresults.branch(&t);

  int sigmas[] = {0, -2, -1, 0, 1, 2};

  double results[6];

  if (sigma >= -1 && sigma < 6) {
    int lsig = sigma, hsig = sigma;
    if (lsig == -1)
      lsig = 0;
    if (hsig == -1)
      hsig = 5;

    for (int sigma_i = lsig; sigma_i <= hsig; ++sigma_i) {

      int nsigma = sigmas[sigma_i];
      bool isobs = (sigma_i == 0);
      bool isexp = (sigma_i > 0);

      // Choose a systematics treatment...
      // The CLfast computation uses no systematics.  This class should only be used for testing purposes.
      //CLfast clcompute;
      //CLfast clfast;
      CLcompute* clcompute = 0;

      if (type == 1) {
        clcompute = new CLsyst;
      } else if (type == 2) {
        clcompute = new CLfit2;
      } else if (type == 3) {
        clcompute = new CLfit;
        ((CLfit*)clcompute)->fitSignal(false);
        ((CLfit*)clcompute)->logSigExclusion(0.005);
      } else if (type == 4) {
        clcompute = new CLfast;
      } else {
        std::cerr << "Unknown type of compute: " << type << std::endl;
        return;
      }

      // The CLsyst computation applies all systematics via Gaussian distribution
      //  CLsyst clsyst;

      // Use CLfit2 for profileLH fitting of systematics-smeared distributions using two fits per pseudoexperiement
      //  CLfit2 clfit2;

      //Use CLfit for profileLH fitting of systematics-smeared distributions using just one fit per pseudoexperiment
      //  CLfit clfit;
      /**
      // If you choose the CLfit option (faster but less powerful than CLfit2), you must
      // specify whether the fit will include signal contributions.  If
      // not, you must specify at which level to exclude signal bins.
      // The cutoff is calculated in terms of log(1+s/b) and the default
      // value is 0.005 (ie, remove bins if log(1+s/b)>0.005.
      clcompute.fitSignal(false);
      clcompute.logSigExclusion(0.005);
       **/

      //  clcompute->setNoviceFlag(false);  // deactivate novice flag if you want to use stat uncertainties
      //  clcompute.useHistoStats(true);  // statistics is turned off by default, only has meaning for CLsyst, CLfit, CLfit2

      /// This is the class for computing cross section limits
      CrossSectionLimit csLim;
      csLim.setup(clcompute);
      csLim.setVerbose(false);

      //95% CL is the default value
      csLim.setCLlevel(0.9);

      //The range of CL values that will satisfy the algorithm: -0.001 < (CL-0.95) < 0.001
      csLim.setAccuracy(0.001);

      //Toggle the number of pseudo-experiments used to find the limit 0 is lowest(fastest), 4 is highest(slowest)
      csLim.setPrecision(prec);

      //Toggle expected/observed to speed things up if you wish
      csLim.calculateExpected(isexp);
      csLim.calculateObserved(isobs);

      //Calculate the expected limit in the case of -2,-1,0,1, or 2-sigma variations of the data relative to bkgd
      csLim.setNSigma(nsigma);

      //Start the cross section limit search at a cross section of 1.0 times the nominal input value
      //  Use this to shorten your calculation if you know roughly where the limit will be.
      //csLim.setSearchSeed(seed);

      // This class is used to test the fit used by the CLfit and CLfit2 classes
      //  Use this to determine the quality of your fit model.
      FitTest fitTest;
      // Set the number of pseudo-experiments to fit
      fitTest.setIterations(2000);
      // Determine if you want fitted pseudo-experiments in the tests
      fitTest.testPE(true);

      //extract the total number of masspoints in the file
      int len = loaders[0].getNMasspoints();
      if (len <= 0) {
        std::cout << "Cannot handle loader with " << len << "masspoints" << std::endl;
        return;
      }

      //create list of mass point indices
      vector<int> v1(len);
      vector<int> v2(len);
      vector<int> v3(len);
      loaders[0].getMasspointList(len, v1.data(), v2.data(), v3.data());

      //loop over all masspoints and perform calculations
      for (int i = 0; i < len; i++) {
        if (v1[i] == mass || mass == -1) {
          //tell the container what point you're working on
          clresults.reset(v1[i], v2[i], v3[i]);

          //printf("Calculating for parameters: %d/%d/%d\n",v1[i],v2[i],v3[i]);

          //Extract the signal & background distributions associated with this point
          SigBkgdDist* sbd = loaders[0].get(v1[i], v2[i], v3[i]);
          //printf("==>Adding channel %s\n",chanNames[0].c_str());
          for (int l = 1; l < nld; l++) {
            SigBkgdDist* sbd2 = loaders[l].get(v1[i], v2[i], v3[i]);
            if (sbd2) {
              //printf("==>Adding channel %s\n",chanNames[l].c_str());
              sbd->append(*sbd2);
            }
          }

          //printf("Sig: %f, Bkgd: %f, Data: %f\n",sbd->totSignal(),sbd->totBkgd(),sbd->totData());

          //If you wish to run the fit test, uncomment the next two lines
          //      fitTest.runTest(sbd,1e6);
          //      continue;

          //calculate CLs
          clcompute->calculateCLs(*sbd, clresults,
                                  CLcompute::LEVEL_VERYVERYFINE /*/LEVEL_VERYFAST*/);

          //report your results for interested observers
          clresults.print();

          //Calculate a cross section limit...
          //These results are reported in the factor by which you must
          //multiply your nominal signal cross section to obtain a 95% CL
          //upper limit for this model... IE, multiply this factor by
          //your model xsec to get your limit in barns

          //csLim.calculate(*sbd,clresults);
          //report your results for interested observers
          //csLim.print();
          //
          if (isobs)
            results[sigma_i] = clresults.xsec_obsfactor;
          if (isexp)
            results[sigma_i] = clresults.xsec_medfactor;

          t.Fill();
          delete sbd;

          if (sigma == -1 || sigma == 0)
            std::cout << "Mass " << v1[i] << " Observed scale factor: " << results[0]
                      << " mass lim: " << 500. * sqrt(results[0]) << " meV" << std::endl;
          if (sigma == -1 || sigma == 3)
            std::cout << "Mass " << v1[i] << " Expected scale factor: " << results[3]
                      << " mass lim: " << 500. * sqrt(results[3]) << " meV" << std::endl;
          if (sigma == -1 || sigma == 2 || sigma == 4)
            std::cout << "Mass " << v1[i] << " 1 sigma range: " << results[2] << " -- "
                      << results[4] << std::endl;
          if (sigma == -1 || sigma == 1 || sigma == 5)
            std::cout << "Mass " << v1[i] << " 2 sigma range: " << results[1] << " -- "
                      << results[5] << std::endl;
        }
      }
      delete clcompute;
    }
  }

  f.Write();
}
void Usage() {
  printf("Using Limits.exe:\n");
  printf(
      " Limits.exe [ Output ROOT File ] [ List of input files ] [ precision ] [ limit type ] [ seed ]\n");
  printf(" The test point input is the integer test variable you wish to look at.\n");
}

int main(int argc, char* argv[]) {

  timeval a, b;
  gettimeofday(&a, NULL);

  string outfile = "/dev/null";
  string infile = "";
  vector<string> infiles;
  int mass = -1;
  int sigma = -1;
  int type = 1;
  int prec = 1;
  char c;
  while ((c = getopt(argc, argv, "o:i:t:s:m:p:")) != -1) {
    string T = optarg;
    switch (c) {
      case 'o':
        outfile = optarg;
        break;
      case 'i':
        infile = optarg;
        break;
      case 't':
        if (T.find("syst") != string::npos)
          type = 1;
        else if (T.find("fit2") != string::npos)
          type = 2;
        else if (T.find("fit") != string::npos)
          type = 3;
        else if (T.find("stat") != string::npos)
          type = 4;
        else
          istringstream(optarg) >> type;
        break;
      case 's':
        istringstream(optarg) >> sigma;
        break;
      case 'm':
        istringstream(optarg) >> mass;
        break;
      case 'p':
        istringstream(optarg) >> prec;
      default:
        break;
    }
  }
  for (int i = optind; i < argc; ++i) {
    infiles.push_back(argv[i]);
  }

  if (infile.empty() && infiles.empty()) {
    cerr << "Need to provide input files!" << endl;
    return -1;
  }
  calcLimit(outfile, infile, infiles, mass, prec, type, sigma);

  gettimeofday(&b, NULL);
  double deltat = a.tv_sec * 1000.0 + a.tv_usec / 1000.0;
  deltat = b.tv_sec * 1000.0 + b.tv_usec / 1000.0 - deltat;
  printf(" %f sec run time\n", deltat / 1000);
  printf("\n");
  return 0;
}

/************* Lines to include for LLR histograms **********
 ************* Include after line 136 ***********************
   int bins = 500;
    double min = -15;
    double max = 15;
    TH1D* sigLLR = clcompute.getLLRdist_sb("LLR_SB",bins,min,max);
    TH1D* bkgLLR = clcompute.getLLRdist_b("LLR_B",bins,min,max);
    TH1D* LLRd = new TH1D("LLR_D","LLR_D",bins,min,max);
    TH1D* LLRsigma1 = new TH1D("LLR_B_1sigmas","LLR_B_1sigmas",bins,min,max);
    TH1D* LLRsigma2 = new TH1D("LLR_B_2sigmas","LLR_B_2sigmas",bins,min,max);

    LLRd->Fill(clresults.llrobs);
    LLRsigma2->Fill(clresults.llrb_m2s);
    LLRsigma1->Fill(clresults.llrb_m1s);
    LLRsigma1->Fill(clresults.llrb_p1s);
    LLRsigma2->Fill(clresults.llrb_p2s);
***********************************************************/
