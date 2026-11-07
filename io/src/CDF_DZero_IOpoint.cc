#include "CDF_DZero_IOpoint.hh"

ClassImp(CDF_DZero_IOpoint);

//-------------------------------------------
CDF_DZero_IOpoint::CDF_DZero_IOpoint() {
  fDataDistribution = NULL;
  fIndepVar1 = -1;
  fIndepVar2 = -1;
  fIndepVar3 = -1;
}

//-------------------------------------------
CDF_DZero_IOpoint::CDF_DZero_IOpoint(int var1, int var2, int var3) {
  fDataDistribution = NULL;
  fIndepVar1 = var1;
  fIndepVar2 = var2;
  fIndepVar3 = var3;
}

//-------------------------------------------
CDF_DZero_IOpoint::~CDF_DZero_IOpoint() {
  //  if(fDataDistribution != NULL) delete fDataDistribution;
  //  fDataDistribution = NULL;
}

//-------------------------------------------
void CDF_DZero_IOpoint::Book(const char* name, int nx, double minX, double maxX) {
  Book(name, nx, minX, maxX, 0, -1, -1);
}

//-------------------------------------------
void CDF_DZero_IOpoint::Book(const char* name, int nx, double minX, double maxX, int ny,
                             double minY, double maxY) {
  if (fIndepVar1 == -1 && fIndepVar2 == -1 && fIndepVar3 == -1) {
    printf("CDF_DZero_IOpoint::Book ==> This class is not properly initialized!\n");
    return;
  }

  int nvars = 1;
  if (fIndepVar2 > -1) {
    nvars++;
    if (fIndepVar3 > -1)
      nvars++;
  }

  const std::string pointName = buildName(name, fIndepVar1, fIndepVar2, fIndepVar3, nvars);
  fName = pointName.c_str();
  fNXbins = nx;
  fNYbins = (ny < 1) ? 1 : ny;
  fMaxX = maxX;
  fMinX = minX;
  fMaxY = maxY;
  fMinY = minY;
  return;
}

bool CDF_DZero_IOpoint::checkHist(TH1* hist) {

  if (hist == NULL)
    return false;
  if (hist->GetNbinsX() != fNXbins)
    return false;
  //floating point protection...
  if (fabs(hist->GetXaxis()->GetXmax() - fMaxX) > 1e-6 * (fMaxX - fMinX))
    return false;
  if (fabs(hist->GetXaxis()->GetXmin() - fMinX) > 1e-6 * (fMaxX - fMinX))
    return false;
  return true;
}

bool CDF_DZero_IOpoint::checkHist(TH2* hist) {
  if (hist == NULL)
    return false;
  if (hist->GetNbinsX() != fNXbins)
    return false;
  if (fabs(hist->GetXaxis()->GetXmax() - fMaxX) > 1e-6 * (fMaxX - fMinX))
    return false;
  if (fabs(hist->GetXaxis()->GetXmin() - fMinX) > 1e-6 * (fMaxX - fMinX))
    return false;
  if (hist->GetNbinsY() != fNYbins)
    return false;
  if (fabs(hist->GetYaxis()->GetXmax() - fMaxY) > 1e-6 * (fMaxY - fMinY))
    return false;
  if (fabs(hist->GetYaxis()->GetXmin() - fMinY) > 1e-6 * (fMaxY - fMinY))
    return false;
  return true;
}

bool CDF_DZero_IOpoint::checkDist(CDF_DZero_Distribution* dist) {
  if (dist == NULL)
    return false;
  if (dist->getNXbins() != fNXbins)
    return false;
  if (fabs(dist->getMaxX() - fMaxX) > 1e-6 * (fMaxX - fMinX))
    return false;
  if (fabs(dist->getMinX() - fMinX) > 1e-6 * (fMaxX - fMinX))
    return false;
  if (dist->getNYbins() > 1) {
    if (dist->getNYbins() != fNYbins)
      return false;
    if (fabs(dist->getMaxY() - fMaxY) > 1e-6 * (fMaxY - fMinY))
      return false;
    if (fabs(dist->getMinY() - fMinY) > 1e-6 * (fMaxY - fMinY))
      return false;
  }

  return true;
}

//-------------------------------------------
std::string CDF_DZero_IOpoint::buildName(const char* channelname, int var1, int var2, int var3,
                                         int nvar) {
  char tool[1024];
  std::string a(channelname);
  if (nvar >= 3)
    sprintf(tool, "(%d,%d,%d)", var1, var2, var3);
  else if (nvar == 2)
    sprintf(tool, "(%d,%d)", var1, var2);
  else if (nvar <= 1)
    sprintf(tool, "(%d)", var1);
  a += tool;

  return a;
}

//-------------------------------------------
void CDF_DZero_IOpoint::addSignalDist(CDF_DZero_Distribution* dist) {
  if (dist == NULL) {
    printf("CDF_DZero_IOpoint::addSignalDist ==> Cannot add a NULL distribution\n");
    return;
  }
  if (!checkDist(dist)) {
    printf("CDF_DZero_IOpoint::addSignalDist ==> Failed format check!\n");
    return;
  }

  TString nname(dist->GetName());
  if (nname.Length() <= 1) {
    char tool[1024];
    sprintf(tool, "%s-SIGNAL %d", fName.Data(), fSignals.GetEntriesFast());
    dist->setDistName(tool);
  }
  fSignals.AddAtAndExpand(dist, fSignals.GetEntriesFast());

  return;
}

//-------------------------------------------
CDF_DZero_Distribution* CDF_DZero_IOpoint::getSignalDist(int i) {
  if (i < 0 || i >= fSignals.GetEntriesFast())
    return NULL;
  return (CDF_DZero_Distribution*)fSignals[i];
}

//-------------------------------------------
void CDF_DZero_IOpoint::addBkgdDist(CDF_DZero_Distribution* dist) {
  if (dist == NULL) {
    printf("CDF_DZero_IOpoint::addBkgdDist ==> Cannot add a NULL distribution\n");
    return;
  }
  if (!checkDist(dist)) {
    printf("CDF_DZero_IOpoint::addBkgdDist ==> Failed format check!\n");
    return;
  }

  TString nname(dist->GetName());
  if (nname.Length() <= 1) {
    char tool[1024];
    if (string(dist->GetName()) == "") {
      sprintf(tool, "%s-BKGD %d", fName.Data(), fBackgrounds.GetEntriesFast());
      dist->setDistName(tool);
    } else {
      sprintf(tool, "%s-%s", fName.Data(), dist->GetName());
      dist->setDistName(tool);
    }
  }
  fBackgrounds.AddAtAndExpand(dist, fBackgrounds.GetEntriesFast());

  return;
}

//-------------------------------------------
CDF_DZero_Distribution* CDF_DZero_IOpoint::getBkgdDist(int i) {
  if (i < 0 || i >= fBackgrounds.GetEntriesFast())
    return NULL;
  return (CDF_DZero_Distribution*)fBackgrounds[i];
}

//-------------------------------------------
void CDF_DZero_IOpoint::addDataDist(TH1* data) {
  if (data == NULL) {
    printf("CDF_DZero_IOpoint::addDataDist ==> Cannot add a NULL distribution\n");
    return;
  }
  if (!checkHist(data)) {
    printf("CDF_DZero_IOpoint::addDataDist ==> Failed format check!\n");
    return;
  }

  std::string a(fName);
  a += "-DATA";

  if (fDataDistribution == NULL) {
    fDataDistribution = (TH1*)data;
    fDataDistribution->SetName(a.c_str());
    fDataDistribution->SetTitle(a.c_str());
  } else
    printf(
        "CDF_DZero_IOpoint::addDataDist ==> You already have a distribution for this instance!\n");

  return;
}

//-------------------------------------------
void CDF_DZero_IOpoint::addDataDist2D(TH2* data) {
  if (data == NULL) {
    printf("CDF_DZero_IOpoint::addDataDist ==> Cannot add a NULL distribution\n");
    return;
  }
  if (!checkHist(data)) {
    printf("CDF_DZero_IOpoint::addDataDist ==> Failed format check!\n");
    return;
  }
  std::string a(fName);
  a += "-DATA";

  if (fDataDistribution == NULL) {
    fDataDistribution = (TH2D*)data;
    fDataDistribution->SetName(a.c_str());
    fDataDistribution->SetTitle(a.c_str());
  } else
    printf(
        "CDF_DZero_IOpoint::addDataDist2D ==> You already have a distribution for this instance!\n");

  return;
}

void CDF_DZero_IOpoint::checkPoint() {

  //Zero any data bins with odd data/bkgd or data/sig ratios
  if (fNYbins <= 1) {
    for (int bx = 1; bx <= fNXbins; bx++) {
      double sumSig = 0;
      double sumBkg = 0;
      for (int i = 0; i < getNSignalDists(); i++) {
        sumSig += getSignalDist(i)->inBin(bx);
      }
      for (int i = 0; i < getNBkgdDists(); i++) {
        sumBkg += getBkgdDist(i)->inBin(bx);
      }
      if (sumSig <= 0 && sumBkg <= 0)
        getDataDist()->SetBinContent(bx, 0);
      if (sumSig > 0 && sumBkg <= 0)
        getDataDist()->SetBinContent(bx, 0);
    }
  } else {
    for (int bx = 1; bx <= fNXbins; bx++) {
      for (int by = 1; by <= fNYbins; by++) {
        double sumSig = 0;
        double sumBkg = 0;
        for (int i = 0; i < getNSignalDists(); i++) {
          sumSig += getSignalDist(i)->inBin(bx, by);
        }
        for (int i = 0; i < getNBkgdDists(); i++) {
          sumBkg += getBkgdDist(i)->inBin(bx, by);
        }
        if (sumSig <= 0 && sumBkg <= 0)
          getDataDist()->SetBinContent(bx, by, 0);
        if (sumSig > 0 && sumBkg <= 0)
          getDataDist()->SetBinContent(bx, by, 0);
      }
    }
  }
}

void CDF_DZero_IOpoint::print(map<string, int>& nSyst, map<string, double>& posSyst,
                              map<string, double>& negSyst) {

  printf("\nCDF_DZero_IOpoint: %s\n", fName.Data());
  printf("Vars (%d, %d, %d)\n", fIndepVar1, fIndepVar2, fIndepVar3);
  float sum = 0;

  for (int i = 0; i < getNSignalDists(); i++) {
    CDF_DZero_Distribution* sig = getSignalDist(i);
    if (sig == NULL) {
      printf("CDF_DZero_IOpoint::print, NULL signal!\n");
      continue;
    }
    printf("\nSignal %d: %s\n", i, sig->GetName());
    sig->print(nSyst, posSyst, negSyst);
  }
  for (int i = 0; i < getNBkgdDists(); i++) {
    CDF_DZero_Distribution* bkg = getBkgdDist(i);
    if (bkg == NULL) {
      printf("CDF_DZero_IOpoint::print, NULL background!\n");
      continue;
    }
    printf("\nBkgd %d: %s\n", i, bkg->GetName());
    bkg->print(nSyst, posSyst, negSyst);
    sum += bkg->getDistribution()->Integral();
  }

  printf("\nTotal Bkgd: %.2f\n", sum);
  printf("Total Data: %.1f\n", fDataDistribution->Integral());
}
