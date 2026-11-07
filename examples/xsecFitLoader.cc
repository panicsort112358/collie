#include <TFile.h>
#include <TTree.h>
#include <CrossSectionLimit.hh>
#include <ExclusionLimit.hh>
#include <ThreeSigmaEvidence.hh>
#include <CollieLoader.hh>
#include <CLfast.hh>
#include <CLfit.hh>
#include <CLfit2.hh>
#include <CLsyst.hh>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <FitTest.hh>
#include <vector>

void newcslLoader(char* outFile, char* inList, char* fitFile, char* m) {

  int mass = atoi(m);

  printf("fitLoader.exe\n");
  printf("InputList: %s\n", inList);
  printf("FitFile: %s\n", fitFile);
  printf("OutFile: %s\n", outFile);
  printf("Mass: %d\n", mass);

  TFile* fitf = new TFile(fitFile);
  if (fitf == NULL) {
    printf("fitLoader ==> Error loading input historam file!\n");
    return;
  }
  if (fitf->IsZombie())
    return;
  TH2D* errMat = (TH2D*)fitf->Get("Error Matrix");
  TH1D* nullPar = (TH1D*)fitf->Get("Null Fit Params");
  TH1D* testPar = (TH1D*)fitf->Get("Test Fit Params");
  TH1D* sigScale = (TH1D*)fitf->Get("Signal Scale Factor");

  if (errMat == NULL || nullPar == NULL || testPar == NULL || sigScale == NULL) {
    printf("fitLoader ==> NULL input fit histograms!\n");
    return;
  }

  CollieLoader loaders[300];
  string fileNames[300];
  string chanNames[300];
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

    fileNames[nld] = fname;

    TFile* ftest = new TFile(fname);
    TList* aList = ftest->GetListOfKeys();
    if (aList->GetEntries() != 1) {
      printf("Incorrect key length for %s\n", fname);
      return;
    }
    chanNames[nld] = aList->At(0)->GetName();

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
  printf("Fit Loader\n");
  printf("%d channels available\n", nld);
  printf("************************************************\n");

  if (!ok)
    return;

  TFile f(outFile, "RECREATE");

  int len = loaders[0].getNMasspoints();

  if (len <= 0) {
    std::cout << "Cannot handle loader with " << len << "masspoints...\n";
    return;
  }
  std::vector<int> v1(len);
  std::vector<int> v2(len);
  std::vector<int> v3(len);
  loaders[0].getMasspointList(len, v1.data(), v2.data(), v3.data());

  for (int i = 0; i < len; i++) {
    if (v1[i] == mass || mass == -1) {

      SigBkgdDist* sbd = loaders[0].get(v1[i], v2[i], v3[i]);
      printf("==>Adding %s\n", chanNames[0].c_str());
      for (int l = 1; l < nld; l++) {
        SigBkgdDist* sbd2 = loaders[l].get(v1[i], v2[i], v3[i]);
        if (sbd2) {
          printf("==>Adding %s\n", chanNames[l].c_str());
          sbd->append(*sbd2);
        }
      }

      printf("Mass: %d\n", v1[i]);
      printf("Sig: %f, Bkgd: %f, Data: %f\n", sbd->totSignal(), sbd->totBkgd(), sbd->totData());

      printf("Loading Fit Results:\n");
      printf("=>Generating TEST fit histograms...\n");
      sbd->generateFitHistos(sigScale, testPar, errMat, "TEST");
      printf("===>Generating NULL fit histograms...\n");
      sbd->generateFitHistos(sigScale, nullPar, errMat, "NULL");
      printf("=====>Done!\n");
    }
  }

  f.Write();
}

int main(int argc, char* argv[]) {

  timeval a, b;
  gettimeofday(&a, NULL);

  if (argc < 5)
    newcslLoader(argv[1], argv[2], argv[3], "-1");
  else
    newcslLoader(argv[1], argv[2], argv[3], argv[4]);

  gettimeofday(&b, NULL);
  double deltat = a.tv_sec * 1000.0 + a.tv_usec / 1000.0;
  deltat = b.tv_sec * 1000.0 + b.tv_usec / 1000.0 - deltat;
  printf(" %f sec run time\n", deltat / 1000);
  printf("\n");
  return 0;
}
