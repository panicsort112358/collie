#include "CDF_DZero_Distribution.hh"
#include <cmath>

ClassImp(CDF_DZero_Distribution);

//-------------------------------------------
CDF_DZero_Distribution::CDF_DZero_Distribution() {
  fDistribution = NULL;
  fName = "";
  f2Ddistribution = false;
}

//-------------------------------------------
CDF_DZero_Distribution::~CDF_DZero_Distribution() {
  //  if(fDistribution != NULL) delete fDistribution;
  //  fDistribution = NULL;
}

//-------------------------------------------
void CDF_DZero_Distribution::addDistribution(TH1* dist) {
  if (dist == NULL) {
    printf("CDF_DZero_Distribution::addDistribution ==> Cannot add a NULL distribution\n");
    return;
  }

  if (fDistribution == NULL)
    fDistribution = (TH1*)dist;
  else
    printf(
        "CDF_DZero_Distribution::addDistribution ==> You already have a distribution for this instance!\n");

  return;
}

//-------------------------------------------
void CDF_DZero_Distribution::addDistribution2D(TH2* dist) {
  if (dist == NULL) {
    printf("CDF_DZero_Distribution::addDistribution2D ==> Cannot add a NULL distribution\n");
    return;
  }

  if (fDistribution == NULL)
    fDistribution = (TH2*)dist;
  else
    printf(
        "CDF_DZero_Distribution::addDistribution2D ==> You already have a distribution for this instance!\n");
  f2Ddistribution = true;
  return;
}

//-------------------------------------------
void CDF_DZero_Distribution::addSystematic(const char* name, TH1* pos, TH1* neg) {

  if (pos == NULL || neg == NULL) {
    printf("CDF_DZero_Distribution::addSystematic ==> Cannot add a NULL systematic\n");
    return;
  }

  std::string p(fDistribution->GetName());
  p += "-PosSyst";

  std::string n(fDistribution->GetName());
  n += "-NegSyst";

  fSystematicsPos.AddAtAndExpand((TH1*)pos->Clone(p.c_str()), fSystematicsPos.GetEntriesFast());
  fSystematicsNeg.AddAtAndExpand((TH1*)neg->Clone(n.c_str()), fSystematicsNeg.GetEntriesFast());
  fSystematicsNames.AddAtAndExpand(new TObjString(name), fSystematicsNames.GetEntriesFast());

  return;
}

//-------------------------------------------
void CDF_DZero_Distribution::addSystematic2D(const char* name, TH2* pos, TH2* neg) {

  if (pos == NULL || neg == NULL) {
    printf("CDF_DZero_Distribution::addSystematic ==> Cannot add a NULL systematic\n");
    return;
  }

  std::string p(fDistribution->GetName());
  p += "-PosSyst";
  std::string n(fDistribution->GetName());
  n += "-NegSyst";

  fSystematicsPos.AddAtAndExpand((TH2*)pos->Clone(p.c_str()), fSystematicsPos.GetEntriesFast());
  fSystematicsNeg.AddAtAndExpand((TH2*)neg->Clone(n.c_str()), fSystematicsNeg.GetEntriesFast());
  fSystematicsNames.AddAtAndExpand(new TObjString(name), fSystematicsNames.GetEntriesFast());

  return;
}

//-------------------------------------------
TH1* CDF_DZero_Distribution::getPositiveSystematic(int i) {
  if (i < 0 || i >= fSystematicsPos.GetEntriesFast())
    return NULL;
  return (TH1*)fSystematicsPos[i];
}

//-------------------------------------------
TH1* CDF_DZero_Distribution::getNegativeSystematic(int i) {
  if (i < 0 || i >= fSystematicsNeg.GetEntriesFast())
    return NULL;
  return (TH1*)fSystematicsNeg[i];
}

//-------------------------------------------
TH2* CDF_DZero_Distribution::getPositiveSystematic2D(int i) {
  if (i < 0 || i >= fSystematicsPos.GetEntriesFast())
    return NULL;
  return (TH2*)fSystematicsPos[i];
}

//-------------------------------------------
TH2* CDF_DZero_Distribution::getNegativeSystematic2D(int i) {
  if (i < 0 || i >= fSystematicsNeg.GetEntriesFast())
    return NULL;
  return (TH2*)fSystematicsNeg[i];
}

//-------------------------------------------
const char* CDF_DZero_Distribution::getSystematicName(int i) {
  if (i < 0 || i >= fSystematicsNames.GetEntriesFast())
    return NULL;
  TObjString* name = (TObjString*)fSystematicsNames[i];
  return name->GetName();
}

//-------------------------------------------
void CDF_DZero_Distribution::setDistName(const char* name) {
  fName = name;
  if (fDistribution != NULL) {
    fDistribution->SetTitle(fName);
    fDistribution->SetName(fName);
  }
  return;
}

//-------------------------------------------
const char* CDF_DZero_Distribution::getDistName() {
  return fName;
}

//-------------------------------------------
int CDF_DZero_Distribution::getNXbins() const {
  if (fDistribution == NULL) {
    printf("CDF_DZero_Distribution, This distribution is NULL!\n");
    return -1;
  } else
    return fDistribution->GetNbinsX();
}

//-------------------------------------------
int CDF_DZero_Distribution::getNYbins() const {
  if (fDistribution == NULL) {
    printf("CDF_DZero_Distribution, This distribution is NULL!\n");
    return -1;
  } else
    return fDistribution->GetNbinsY();
}

//-------------------------------------------
double CDF_DZero_Distribution::getMaxX() const {
  if (fDistribution == NULL) {
    printf("CDF_DZero_Distribution, This distribution is NULL!\n");
    return -1;
  } else
    return fDistribution->GetXaxis()->GetXmax();
}

//-------------------------------------------
double CDF_DZero_Distribution::getMinX() const {
  if (fDistribution == NULL) {
    printf("CDF_DZero_Distribution, This distribution is NULL!\n");
    return -1;
  } else
    return fDistribution->GetXaxis()->GetXmin();
}

//-------------------------------------------
double CDF_DZero_Distribution::getMaxY() const {
  if (fDistribution == NULL) {
    printf("CDF_DZero_Distribution, This distribution is NULL!\n");
    return -1;
  } else
    return fDistribution->GetYaxis()->GetXmax();
}

//-------------------------------------------
double CDF_DZero_Distribution::getMinY() const {
  if (fDistribution == NULL) {
    printf("CDF_DZero_Distribution, This distribution is NULL!\n");
    return -1;
  } else
    return fDistribution->GetYaxis()->GetXmin();
}

//-------------------------------------------
double CDF_DZero_Distribution::inBin(int i, int j) const {
  if (fDistribution == NULL) {
    printf("CDF_DZero_Distribution, This distribution is NULL!\n");
    return -1;
  } else {
    if (getNYbins() <= 1)
      return fDistribution->GetBinContent(i);
    else
      return fDistribution->GetBinContent(i, j);
  }
}

void CDF_DZero_Distribution::print(map<string, int>& nSyst, map<string, double>& posSyst,
                                   map<string, double>& negSyst) {

  printf("Total evts: %.2f\n", fDistribution->Integral());
  printf("%d systematics:\n", getNSystematics());
  for (int s = 0; s < getNSystematics(); s++) {
    double meanP = 0, meanN = 0;
    double maxP = 0, maxN = 0;

    if (fDistribution->GetNbinsY() > 1) {
      for (int x = 1; x <= getPositiveSystematic(s)->GetNbinsX(); x++) {
        for (int y = 1; y <= getPositiveSystematic(s)->GetNbinsY(); y++) {
          double ref = fDistribution->GetBinContent(x, y);
          meanP += getPositiveSystematic(s)->GetBinContent(x, y);
          meanN += getNegativeSystematic(s)->GetBinContent(x, y);
          if (ref == 0)
            continue;

          double vp = (getPositiveSystematic2D(s)->GetBinContent(x, y) - ref) / ref;
          double vn = (getNegativeSystematic2D(s)->GetBinContent(x, y) - ref) / ref;

          if (fabs(vp) > maxP)
            maxP = fabs(vp);
          if (fabs(vn) > maxN)
            maxN = fabs(vn);
        }
      }
    } else {
      for (int x = 1; x <= getPositiveSystematic(s)->GetNbinsX(); x++) {
        double ref = fDistribution->GetBinContent(x);
        meanP += getPositiveSystematic(s)->GetBinContent(x);
        meanN += getNegativeSystematic(s)->GetBinContent(x);

        if (ref == 0)
          continue;
        double vp = (getPositiveSystematic(s)->GetBinContent(x) - ref) / ref;
        double vn = (getNegativeSystematic(s)->GetBinContent(x) - ref) / ref;
        if (fabs(vp) > maxP)
          maxP = fabs(vp);
        if (fabs(vn) > maxN)
          maxN = fabs(vn);
      }
    }

    maxN *= -1.0;

    double ref = fDistribution->Integral();
    //    printf("Ref: %f, MeanP: %f, MeanN: %f\n",ref,meanP,meanN);
    if (ref > 0) {
      meanP -= ref;
      meanN -= ref;
      meanP /= ref;
      meanN /= ref;
    }

    if (fabs(maxP - meanP) > fabs(1e-2 * meanP) || fabs(maxN - meanN) > fabs(1e-2 * meanN))
      printf("Shape systematic %d, %s: MeanP: %.4f, MaxP: %.4f, MeanN: %.4f, MaxN: %.4f\n", s,
             getSystematicName(s), meanP, maxP, meanN, maxN);
    else
      printf("Flat systematic %d, %s: MeanP: %.4f, MeanN: %.4f\n", s, getSystematicName(s), meanP,
             meanN);
    if (ref > 1e-5) {
      if (posSyst.find(getSystematicName(s)) == posSyst.end()) {
        posSyst[getSystematicName(s)] = meanP;
        negSyst[getSystematicName(s)] = meanN;
        nSyst[getSystematicName(s)] = 1;
      } else {
        //      if(posSyst[getSystematicName(s)]<meanP) posSyst[getSystematicName(s)] = meanP;
        //      if(negSyst[getSystematicName(s)]>meanN) negSyst[getSystematicName(s)] = meanN;
        posSyst[getSystematicName(s)] += meanP;
        negSyst[getSystematicName(s)] += meanN;
        nSyst[getSystematicName(s)] += 1;
      }
    }
  }
}
