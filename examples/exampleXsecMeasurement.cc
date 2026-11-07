#include <TFile.h>
#include <TTree.h>
#include <CrossSectionCalc.hh>
#include <CollieLoader.hh>
#include <FitTest.hh>
#include <CLfast.hh>
#include <CLsyst.hh>
#include <CLfit.hh>
#include <CLfit2.hh>
#include <sys/time.h>
#include <fstream>
#include <vector>

void calcXsec(char* outFile, char* inList, char* m) {

  int mass = atoi(m);

  printf("\n\ncollieXsecCalc.exe\n");
  printf("Input List: %s\n", inList);
  printf("Output File: %s\n", outFile);
  printf("Mass Point: %d\n", mass);

  CollieLoader loaders[300];
  string chanNames[300];
  string fileNames[300];
  int nld = 0;
  ifstream streamIn(inList);
  if (!streamIn) {
    cout << "Error: Could not open " << inList << endl;
    return;
  }

  bool ok = true;
  char fname[1024];
  char options[1024];
  while (streamIn >> fname) {
    cout << "Reading: " << fname << endl;

    TFile* ftest = new TFile(fname);
    TList* aList = ftest->GetListOfKeys();
    if (aList->GetEntries() != 1) {
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

  printf("\n************************************************\n");
  printf("Collie Example Xsec Calculation\n");
  printf("%d channel(s) available\n", nld);
  printf("************************************************\n");

  if (!ok)
    return;

  ///Create output container for the results you're about to calculate
  TFile f(outFile, "RECREATE");
  TTree t("SCAN", "SCAN");
  CLpoint clresults;
  clresults.branch(&t);

  // Choose a systematics treatment...
  // The CLfast computation uses no systematics.  This class should only be used for testing purposes.
  CLfast clcompute;

  /// This is the class for computing cross section limits
  CrossSectionCalc csCalc;
  csCalc.setup(&clcompute);
  csCalc.setVerbose(false);
  csCalc.setPrecision(0); //0 is lowest(fastest), 4 is highest(slowest)

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
  std::vector<int> v1(len);
  std::vector<int> v2(len);
  std::vector<int> v3(len);
  loaders[0].getMasspointList(len, v1.data(), v2.data(), v3.data());

  //loop over all masspoints and perform calculations
  for (int i = 0; i < len; i++) {
    if (v1[i] == mass || mass == -1) {
      //tell the container what point you're working on
      clresults.reset(v1[i], v2[i], v3[i]);

      printf("Calculating for parameters: %d/%d/%d\n", v1[i], v2[i], v3[i]);

      //Extract the signal & background distributions associated with this point
      SigBkgdDist* sbd = loaders[0].get(v1[i], v2[i], v3[i]);
      printf("==>Adding channel %s\n", chanNames[0].c_str());
      for (int l = 1; l < nld; l++) {
        SigBkgdDist* sbd2 = loaders[l].get(v1[i], v2[i], v3[i]);
        if (sbd2) {
          printf("==>Adding channel %s\n", chanNames[l].c_str());
          sbd->append(*sbd2);
        }
      }

      printf("Sig: %f, Bkgd: %f, Data: %f\n", sbd->totSignal(), sbd->totBkgd(), sbd->totData());

      //If you wish to run the fit test, uncomment the next two lines
      //      fitTest.runTest(sbd,1e6);
      //      continue;

      //Calculate a cross section
      // The signal rate is floated as a free parameter
      // The resulting fit gives you the fitted xsec in units of
      // the input cross section.
      csCalc.calculate(*sbd, clresults);
      //report your results for interested observers
      csCalc.print();
      //print out the covariance matrix if you want
      csCalc.printEMAT();

      //Perform a significance test for your signal cross section
      // You do not HAVE to do this to get a cross section calculation!
      // The second parameter determines the signal cross section size
      // to be used in the generation of pseudo-experiments (0.0 for NULL hypothesis, 1.0 for TEST hypothesis)
      // The third parameter determines how many pseudo-experiments will
      // be generated for the test.  The results are stored in the output file.
      csCalc.testFitPE(*sbd, 0.0, 15000);

      t.Fill();
      delete sbd;
    }
  }

  f.Write();
}
void Usage() {
  printf("Using collieXsecCalc.exe:\n");
  printf(" collieXsecCalc.exe [ Output ROOT File ] [ List of input files ] [ Test Point ]\n");
  printf(" The test point input is the integer test variable you wish to look at.\n");
  printf(" If you leave off the test point, the code will loop over all\n");
  printf(" available points in succession.\n");
}

int main(int argc, char* argv[]) {

  timeval a, b;
  gettimeofday(&a, NULL);
  if (argc == 1) {
    Usage();
    return 1;
  }

  if (argc < 4)
    calcXsec(argv[1], argv[2], "-1");
  else
    calcXsec(argv[1], argv[2], argv[3]);

  gettimeofday(&b, NULL);
  double deltat = a.tv_sec * 1000.0 + a.tv_usec / 1000.0;
  deltat = b.tv_sec * 1000.0 + b.tv_usec / 1000.0 - deltat;
  printf(" %f sec run time\n", deltat / 1000);
  printf("\n");
  return 0;
}
