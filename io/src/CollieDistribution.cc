#include <CollieDistribution.hh>
#include <CollieMasspoint.hh>
#include <CollieChannel.hh>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>

using namespace std;

ClassImp(CollieDistribution);

const double CollieDistribution::rate_IgnoreLuminosity = 1000.0f;

CollieDistribution::CollieDistribution()
    : TNamed(),
      fMutable(false),
      p_MassPoint(),
      fMinX(-1),
      fMaxX(-1),
      fMinY(-1),
      fMaxY(-1),
      fEfficiency(-1),
      fNXbins(1),
      fNYbins(1),
      fTrueBinCount(0),
      fBins(),
      fBinStat(),
      fNmodels(0),
      fModelXsecs(),
      fSystNames(),
      fSystematicsPos(),
      fSystematicsNeg(),
      fFloatFlag(),
      fLogNormalFlag(),
      fPoissonFlag(false),
      fPoissonNorm(0),
      fPoissonErrPos(0),
      fPoissonErrNeg(0),
      fNsyst(0),
      fSystIndexOuter(),
      fSystIndexInner(),
      fLinFloat(),
      fLinLogN(),
      fLinSystPos(),
      fLinSystNeg(),
      fLinSyst(),
      fEfficiencySums(),
      fLinBins(),
      fLinBinStat(),
      fLinearized(false) {
  fName = "Default";
  //std::cout << "CollieDistribution()\n";
}

CollieDistribution::CollieDistribution(const CollieDistribution& d) {
  fMutable = false;
  fMinX = -1;
  fMaxX = -1;
  fMinY = -1;
  fMaxY = -1;
  fNXbins = 1;
  fNYbins = 1;
  fEfficiency = -1;
  fTrueBinCount = 0;
  fNmodels = 0;
  fLinearized = false;
  fPoissonFlag = false;
  fPoissonNorm = 0;
  fPoissonErrPos = 0;
  fPoissonErrNeg = 0;
  (*this) = d;
  //std::cout << "CollieDistribution(d)\n";
}

CollieDistribution::CollieDistribution(const char* name, int nX, double mX, double MX, int nY,
                                       double mY, double MY, int nModels)
    : TNamed(),
      fMutable(true),
      p_MassPoint(),
      fMinX(mX),
      fMaxX(MX),
      fMinY(mY),
      fMaxY(MY),
      fEfficiency(),
      fNXbins(nX),
      fNYbins(nY <= 1 ? 1 : nY),
      fTrueBinCount(fNXbins * fNYbins),
      fBins(),
      fBinStat(),
      fNmodels(nModels),
      fModelXsecs(),
      fSystNames(),
      fSystematicsPos(),
      fSystematicsNeg(),
      fFloatFlag(),
      fLogNormalFlag(),
      fPoissonFlag(false),
      fPoissonNorm(0),
      fPoissonErrPos(0),
      fPoissonErrNeg(0),
      fNsyst(0),
      fSystIndexOuter(),
      fSystIndexInner(),
      fLinFloat(),
      fLinLogN(),
      fLinSystPos(),
      fLinSystNeg(),
      fLinSyst(),
      fEfficiencySums(),
      fLinBins(),
      fLinBinStat(),
      fLinearized(false) {
  fName = name;
  fBins.Set(fTrueBinCount);
  fBins.Reset();

  fBinStat.Set(fTrueBinCount);
  fBinStat.Reset();

  fModelXsecs.Set(fNmodels);
  fModelXsecs.Reset(-1);
  //std::cout << "CollieDistribution(name, ...)\n";
}

CollieDistribution& CollieDistribution::operator=(const CollieDistribution& d) {
  fEfficiency = d.fEfficiency;
  fNXbins = d.fNXbins;
  fNYbins = d.fNYbins;
  fMinX = d.fMinX;
  fMaxX = d.fMaxX;
  fMinY = d.fMinY;
  fMaxY = d.fMaxY;

  fName = d.fName;

  fTrueBinCount = d.fTrueBinCount;
  fBins = d.fBins;
  fBinStat = d.fBinStat;

  fNmodels = d.fNmodels;
  fModelXsecs = d.fModelXsecs;

  fMutable = d.fMutable;
  fLinearized = false;
  fPoissonFlag = d.fPoissonFlag;
  fPoissonNorm = d.fPoissonNorm;
  fPoissonErrPos = d.fPoissonErrPos;
  fPoissonErrNeg = d.fPoissonErrNeg;

  fSystNames.Delete();
  fFloatFlag.Delete();
  fLogNormalFlag.Delete();
  fSystematicsPos.Delete();
  fSystematicsPos.Delete();
  for (int i = 0; i < d.fSystNames.GetEntriesFast(); ++i) {
    fSystNames.Add(d.fSystNames.At(i));
    fFloatFlag.Add(d.fFloatFlag.At(i));
    fLogNormalFlag.Add(d.fLogNormalFlag.At(i));
    fSystematicsPos.Add(d.fSystematicsPos.At(i));
    fSystematicsNeg.Add(d.fSystematicsNeg.At(i));
  }

  fLinSyst = d.fLinSyst;
  fEfficiencySums = d.fEfficiencySums;

  return (*this);
}

CollieDistribution::~CollieDistribution() {}

double CollieDistribution::getBinStatErr(int i, int j) const {
  if (i < 0 || i >= fNXbins)
    return -1;
  if (!fLinearized) {
    if (fNYbins <= 1)
      return fBinStat[i];
    else {
      if (j >= fNYbins)
        return -1;
      return fBinStat[i + j * fNXbins];
    }
  } else {
    if (fNYbins <= 1)
      return fLinBinStat[i];
    else {
      if (j >= fNYbins)
        return -1;
      return fLinBinStat[i + j * fNXbins];
    }
  }
}

int CollieDistribution::getSystIndex(string syst) const {

  for (int i = 0; i < fSystNames.GetEntriesFast(); i++) {
    string ts = string(((TObjString*)fSystNames[i])->String().Data());
    if (syst.compare(ts) == 0)
      return i;
  }

  return -1;
}

bool CollieDistribution::hasSystematic(string syst) const {
  int outI = getSystIndex(syst);
  if (outI < 0)
    return false;
  return true;
}

double CollieDistribution::getBinSystValue(int syst, int i, int j, const double* fluctList) const {

  if (j >= fNYbins)
    return -1;
  if (i < 0 || i >= fNXbins)
    return -1;
  if (!fLinearized) {
    printf("CollieDistribution::getBinSystValue, this dist is not linearized...: %s\n",
           fName.Data());
    return -1;
  }

  //get fluctuated values
  double fluct = 0;
  if (fSystIndexInner[syst] == -1)
    return fluct;
  double rand = fluctList[syst];

  InfoForEfficiencyCalculation e;
  e = (fNYbins > 1) ? fLinSyst[fSystIndexInner[syst]][i + j * fNXbins]
                    : fLinSyst[fSystIndexInner[syst]][i];
  fluct = (rand < 0) ? e.sigmaN : e.sigmaP;

  return fluct;
}

void CollieDistribution::setEfficiency(double efficiency) {
  if (fMutable)
    fEfficiency = efficiency;
}

void CollieDistribution::setNormalizedBinValue(double value, int i, int j) {
  if (!fMutable || i < 0 || i >= fNXbins)
    return;
  if (fNYbins > 1) {
    if (j < 0 || j >= fNYbins)
      return;
    fBins[i + j * fNXbins] = value;
  } else
    fBins[i] = value;
}

void CollieDistribution::setBinStatErr(double value, int i, int j) {
  if (!fMutable || i < 0 || i >= fNXbins)
    return;
  if (fNYbins > 1) {
    if (j < 0 || j >= fNYbins)
      return;
    fBinStat[i + j * fNXbins] = value;
  } else
    fBinStat[i] = value;
}

void CollieDistribution::setModelXsec(double value, int nmodel) {
  if (!fMutable || nmodel < 0 || nmodel >= fNmodels)
    return;
  fModelXsecs[nmodel] = value;
}

void CollieDistribution::setModelXsec(double value, const char* name) {
  if (!fMutable)
    return;
  int nmodel = lookupModel(name);
  if (nmodel >= 0)
    fModelXsecs[nmodel] = value;
}

TH1* CollieDistribution::draw(string title) const {
  if (fNYbins > 1)
    return NULL;
  std::string a(GetName());

  if (title.size() > 0) {
    a += "-HIST-";
    a += title;
  } else {
    a += "-HIST";
  }

  TH1* retval = new TH1D(a.c_str(), a.c_str(), fNXbins, fMinX, fMaxX);
  retval->Sumw2();
  for (int i = 0; i < fNXbins; ++i) {
    retval->SetBinContent(i + 1, getEfficiency(i));
    retval->SetBinError(i + 1, getBinStatErr(i, -1));
  }
  return retval;
}

TH2* CollieDistribution::draw2D(string title) const {
  std::string a(GetName());
  a += "-HIST";
  if (title.size() > 0)
    a = title;

  TH2* retval = new TH2D(a.c_str(), a.c_str(), fNXbins, fMinX, fMaxX, fNYbins, fMinY, fMaxY);
  retval->Sumw2();
  for (int i = 0; i < fNXbins; ++i) {
    for (int j = 0; j < fNYbins; ++j) {
      retval->SetBinContent(i + 1, j + 1, getEfficiency(i, j));
      retval->SetBinError(i + 1, j + 1, getBinStatErr(i, j));
    }
  }
  return retval;
}

int CollieDistribution::lookupModel(const char* name) const {
  if (p_MassPoint == NULL || p_MassPoint->getChannel() == NULL)
    return -1;
  for (int i = 0; i < fNmodels; ++i)
    if (!strcmp(name, p_MassPoint->getChannel()->getModelName(i)))
      return i;
  return -1;
}

int CollieDistribution::fillFromHistogram(const TH1* histo, int rebin) {
  if (histo == NULL) {
    printf("CollieDistribution::fillFromHistogram, NULL histo!\n");
    return false;
  }
  TH1* htmp = (TH1*)histo->Clone("temp1d");

  if (htmp->GetNbinsX() != fNXbins) {
    if (rebin > 0) {
      if ((htmp->GetNbinsX() / rebin) == fNXbins)
        htmp->Rebin(rebin);
      else {
        printf("CollieDistribution::fillFromHistogram, Unrecoverable binning mismatch!\n");
        htmp->Delete();
        return false;
      }
    } else {
      printf("CollieDistribution::fillFromHistogram, Unrecoverable binning mismatch!\n");
      htmp->Delete();
      return false;
    }
  }

  //wf 5/27/08, ignore under/overflow
  //  double sum=htmp->GetSumOfWeights()+htmp->GetBinContent(0)+htmp->GetBinContent(fNXbins+1);
  double sum = htmp->GetSumOfWeights();
  fEfficiency = sum;
  if (sum > 0)
    sum = 1.0 / sum;

  if (htmp->GetNbinsX() == fNXbins) {

    ////Check stat errors (M. Owen 16-3-08)
    bool check = checkStats(histo, 1);
    if (!check) {
      cerr << "CollieDistribution::fillFromHistogram problem in bin stat errors" << endl;
      return false;
    }

    for (int i = 0; i < fNXbins; ++i) {
      if (htmp->GetBinContent(i + 1) * sum < 0) {
        printf(
            "CollieDistribution::fillFromHistogram, Negative input histogram bin value.  Bin=%d\n",
            i);
        return false;
      }
      setNormalizedBinValue(htmp->GetBinContent(i + 1) * sum, i, -1);
      setBinStatErr(htmp->GetBinError(i + 1), i, -1);
    }

    //wf 5/27/08, ignore under/overflow
    /*
    //catch overflow, underflow
    double uf = htmp->GetBinContent(0)*sum;
    double of = htmp->GetBinContent(fNXbins+1)*sum;

    double suf = htmp->GetBinError(0);
    double sof = htmp->GetBinError(fNXbins+1);

    fBinStat[0] = fBinStat[0]*fBinStat[0]*fBins[0] + uf*suf*suf;
    if((fBins[0] + uf)>0) fBinStat[0] /= (fBins[0] + suf);
    assert(fBinStat[0]>=0);
    fBinStat[0] = sqrt(fBinStat[0]);

    fBinStat[fNXbins-1] = fBinStat[fNXbins-1]*fBinStat[fNXbins-1]*fBins[fNXbins-1] + of*sof*sof;
    if((fBins[fNXbins-1] + of)>0) fBinStat[fNXbins-1] /= (fBins[fNXbins-1] + sof);
    assert(fBinStat[fNXbins-1]>=0);
    fBinStat[fNXbins-1] = sqrt(fBinStat[fNXbins-1]);

    fBins[0]+= uf*sum;
    fBins[fNXbins-1]+=of*sum;
    */
  } else {
    htmp->Delete();
    printf("CollieDistribution::fillFromHistogram, Unrecoverable binning mismatch(2)!\n");
    return false;
  }

  htmp->Delete();
  return true;
}

int CollieDistribution::fillFromHistogram(const CollieHistogram* histo, int rebin) {
  if (histo == NULL) {
    printf("CollieDistribution::fillFromHistogram, NULL histo!\n");
    return false;
  }

  if (histo->nBins() != fNXbins) {
    printf("CollieDistribution::fillFromHistogram, Unrecoverable binning mismatch!\n");
    return false;
  }

  //ignore under/overflow
  double sum = histo->Integral();
  fEfficiency = sum;
  if (sum > 0)
    sum = 1.0 / sum;

  ////Check stat errors (M. Owen 16-3-08)
  bool check = checkStats(histo, 1);
  if (!check) {
    cerr << "CollieDistribution::fillFromHistogram problem in bin stat errors" << endl;
    return false;
  }

  for (int i = 0; i < fNXbins; ++i) {
    if (histo->InBin(i) * sum < 0) {
      printf("CollieDistribution::fillFromHistogram, Negative input histogram bin value.  Bin=%d\n",
             i);
      return false;
    }
    setNormalizedBinValue(histo->InBin(i) * sum, i, -1);
    setBinStatErr(histo->BinErr(i), i, -1);
  }

  return true;
}

int CollieDistribution::fillFromHistogram(const TH2* histo, int rebinX, int rebinY) {
  if (histo == NULL)
    return false;
  TH2* htmp = (TH2*)histo->Clone("temp2d");

  if (histo->GetNbinsX() != fNXbins) {
    if (rebinX > 0) {
      if ((histo->GetNbinsX() / rebinX) == fNXbins)
        htmp->RebinX(rebinX);
      else {
        htmp->Delete();
        return false;
      }
    } else {
      htmp->Delete();
      return false;
    }
  }
  if (histo->GetNbinsY() != fNYbins) {
    if (rebinY > 0) {
      if ((histo->GetNbinsY() / rebinY) == fNYbins)
        htmp->RebinY(rebinY);
      else {
        htmp->Delete();
        return false;
      }
    } else {
      htmp->Delete();
      return false;
    }
  }

  double sum = htmp->GetSumOfWeights();
  fEfficiency = sum;
  if (sum > 0)
    sum = 1.0 / sum;

  if (htmp->GetNbinsX() == fNXbins) {

    ////Check stat errors (M. Owen 16-3-08)
    bool check = checkStats(histo, 1);
    if (!check) {
      cerr << "CollieDistribution::fillFromHistogram problem in bin stat errors" << endl;
      return false;
    }

    for (int i = 0; i < fNXbins; ++i)
      for (int j = 0; j < fNYbins; ++j) {
        if (htmp->GetBinContent(i + 1, j + 1) * sum < 0) {
          printf(
              "CollieDistribution::fillFromHistogram, Negative input histogram bin value.  Bin=%d,%d\n",
              i, j);
          return false;
        }
        setNormalizedBinValue(htmp->GetBinContent(i + 1, j + 1) * sum, i, j);
        setBinStatErr(htmp->GetBinError(i + 1, j + 1), i, j);
      }
  } else {
    htmp->Delete();
    return false;
  }

  htmp->Delete();
  return true;
}

int CollieDistribution::fillFromHistogram(const CollieHistogram2d* histo, int rebinX, int rebinY) {
  if (histo == NULL)
    return false;

  if (histo->nxBins() != fNXbins) {
    printf("CollieDistribution::fillFromHistogram, Unrecoverable binning mismatch!\n");
    return false;
  }
  if (histo->nyBins() != fNYbins) {
    printf("CollieDistribution::fillFromHistogram, Unrecoverable binning mismatch!\n");
    return false;
  }

  double sum = histo->Integral();
  fEfficiency = sum;
  if (sum > 0)
    sum = 1.0 / sum;

  ////Check stat errors (M. Owen 16-3-08)
  bool check = checkStats(histo, 1);
  if (!check) {
    cerr << "CollieDistribution::fillFromHistogram problem in bin stat errors" << endl;
    return false;
  }

  for (int i = 0; i < fNXbins; ++i)
    for (int j = 0; j < fNYbins; ++j) {
      if (histo->InBin(i, j) * sum < 0) {
        printf(
            "CollieDistribution::fillFromHistogram, Negative input histogram bin value.  Bin=%d,%d\n",
            i, j);
        return false;
      }
      setNormalizedBinValue(histo->InBin(i, j) * sum, i, j);
      setBinStatErr(histo->BinErr(i, j), i, j);
    }

  return true;
}

double CollieDistribution::sumEfficiency() {
  double sum = 0;
  for (int i = 0; i < fTrueBinCount; ++i)
    sum += fBins[i];
  return sum * fEfficiency;
}

double CollieDistribution::sumEfficiency(int i, int j) {
  double sum = 0;
  for (int k = 0; k < fTrueBinCount; ++k) {
    if (k >= i && k < j)
      sum += fBins[k];
  }
  return sum * fEfficiency;
}

void CollieDistribution::addSystematic(const char* name, TH1D* pos, TH1D* neg) {
  if (pos == NULL || neg == NULL) {
    printf("CollieDistribution::addSystematic, NULL histos!\n");
    return;
  }
  string aname(name);
  for (int i = 0; i < fSystNames.GetEntriesFast(); ++i)
    if (aname == (string)((TObjString*)fSystNames[i])->GetName()) {
      printf("We've already got this systematic!!\n");
      return;
    }

  char title[250];
  sprintf(title, "Pos syst %s:%s", pos->GetName(), fName.Data());
  TH1D* p = (TH1D*)pos->Clone(title);
  fSystematicsPos.Add(p);
  sprintf(title, "Neg syst %s:%s", pos->GetName(), fName.Data());
  TH1D* n = (TH1D*)neg->Clone(title);
  fSystematicsNeg.Add(n);

  fSystNames.Add(new TObjString(name));
  fFloatFlag.Add(new TObject());
  fLogNormalFlag.Add(new TObject());
}

void CollieDistribution::addSystematic2D(const char* name, TH2D* pos, TH2D* neg) {
  if (pos == NULL || neg == NULL) {
    printf("CollieDistribution::addSystematic2D, NULL histos!\n");
    return;
  }

  string aname(name);

  for (int i = 0; i < fSystNames.GetEntriesFast(); ++i)
    if (aname == (string)((TObjString*)fSystNames[i])->GetName()) {
      printf("We've already got this systematic!!\n");
      return;
    }

  char title[250];
  sprintf(title, "Pos syst %s:%s", pos->GetName(), fName.Data());
  TH2D* p = (TH2D*)pos->Clone(title);
  fSystematicsPos.Add(p);
  sprintf(title, "Neg syst %s:%s", pos->GetName(), fName.Data());
  TH2D* n = (TH2D*)neg->Clone(title);
  fSystematicsNeg.Add(n);

  fSystNames.Add(new TObjString(name));
  fFloatFlag.Add(new TObject());
  fLogNormalFlag.Add(new TObject());
}

void CollieDistribution::setPoissonFlag(double norm, double errPos, double errNeg) {
  fPoissonFlag = true;
  fPoissonNorm = norm;
  fPoissonErrPos = errPos;
  fPoissonErrNeg = errNeg;

  fSystematicsPos.Clear();
  fSystematicsNeg.Clear();
  fSystNames.Clear();
  fFloatFlag.Clear();
  fLogNormalFlag.Clear();

  return;
}

void CollieDistribution::linearize(vector<string> inputNames) {

  ///Clear it first
  fNsyst = fSystNames.GetEntries();

  if (!fLinearized) {
    //index from outside to inside
    fSystIndexOuter = new int[fNsyst];
    memset(fSystIndexOuter, 0, sizeof(int) * fNsyst);

    // THIS CHUNK WILL BE REPLACED:
    // (but for now, they are still used on non-critical paths)
    fLinSystPos = new double*[fNsyst];
    fLinSystNeg = new double*[fNsyst];
    for (int i = 0; i < fNsyst; ++i) {
      fLinSystPos[i] = new double[fTrueBinCount];
      fLinSystNeg[i] = new double[fTrueBinCount];
      memset(fLinSystPos[i], 0, sizeof(double) * fTrueBinCount);
      memset(fLinSystNeg[i], 0, sizeof(double) * fTrueBinCount);
    }

    // Eventually REPLACES ESTABLISHMENT OF fLinSystPos and fLinSystNeg:
    fLinSyst.clear();
    fLinSyst.resize(fNsyst);
    {
      std::vector<InfoForEfficiencyCalculation> tmp(fTrueBinCount);
      fill(fLinSyst.begin(), fLinSyst.end(), tmp);
      fEfficiencySums.clear();
      fEfficiencySums.resize(fTrueBinCount);
    }

    fLinBins = new double[fTrueBinCount];
    fLinBinStat = new double[fTrueBinCount];
    memset(fLinBins, 0, sizeof(double) * fTrueBinCount);
    memset(fLinBinStat, 0, sizeof(double) * fTrueBinCount);

    fLinFloat = new bool[fNsyst];
    fLinLogN = new bool[fNsyst];
  }

  for (int x = 0; x < fNXbins; ++x) {
    for (int y = 0; y < fNYbins; ++y) {
      fLinBins[x + y * fNXbins] = fBins[x + y * fNXbins];
      fLinBinStat[x + y * fNXbins] = fBinStat[x + y * fNXbins];
    }
  }

  //index from inside to outside
  fSystIndexInner = new int[inputNames.size()];
  memset(fSystIndexInner, 0, sizeof(int) * inputNames.size());
  for (uint s = 0; s < inputNames.size(); ++s)
    fSystIndexInner[s] = -1;

  for (int i = 0; i < fNsyst; ++i) {
    for (uint s = 0; s < inputNames.size(); ++s) {
      bool looking = true;
      if (looking && inputNames[s] == string(((TObjString*)fSystNames[i])->GetName())) {
        fSystIndexOuter[i] = s;
        fSystIndexInner[s] = i;
        looking = false;
      }
    }
  }
  if (fLinearized)
    return;

  for (int i = 0; i < fNsyst; ++i) {
    fLinFloat[i] = getFloatFlag(fSystNames[i]->GetName());
    fLinLogN[i] = getLogNormalFlag(fSystNames[i]->GetName());

    //Force large systematics to use logNormal instead of Gaussian
    //    if(!fLinLogN[i]) fLinLogN[i] = overLogNormThresh(i);
  }

  ///form "flat and linear" systematics dists
  if (fNYbins > 1) {
    for (int i = 0; i < fNsyst; ++i) {
      for (int bx = 0; bx < fNXbins; ++bx) {
        for (int by = 0; by < fNYbins; ++by) {
          // THIS CHUNK WILL BE REPLACED:
          // But for now, the old ones are still used in non-critical paths
          fLinSystPos[i][bx + by * fNXbins] =
              ((TH1D*)fSystematicsPos[i])->GetBinContent(bx + 1, by + 1);
          fLinSystNeg[i][bx + by * fNXbins] =
              ((TH1D*)fSystematicsNeg[i])->GetBinContent(bx + 1, by + 1);
          // REPLACEMENT:
          InfoForEfficiencyCalculation& sb = fLinSyst[i][bx + by * fNXbins];
          sb.sigmaP = ((TH1D*)fSystematicsPos[i])->GetBinContent(bx + 1, by + 1);
          sb.sigmaN = ((TH1D*)fSystematicsNeg[i])->GetBinContent(bx + 1, by + 1);

          sb.assym = true;
          if (sb.sigmaP == sb.sigmaN)
            sb.assym = false;
        }
      }
    }
  } else {
    for (int i = 0; i < fNsyst; ++i) {
      for (int b = 0; b < fNXbins; ++b) {
        // THIS CHUNK WILL BE REPLACED:
        fLinSystPos[i][b] = ((TH1D*)fSystematicsPos[i])->GetBinContent(b + 1);
        fLinSystNeg[i][b] = ((TH1D*)fSystematicsNeg[i])->GetBinContent(b + 1);
        // REPLACEMENT:
        InfoForEfficiencyCalculation& sb = fLinSyst[i][b];
        sb.sigmaP = ((TH1D*)fSystematicsPos[i])->GetBinContent(b + 1);
        sb.sigmaN = ((TH1D*)fSystematicsNeg[i])->GetBinContent(b + 1);

        sb.assym = true;
        if (sb.sigmaP == sb.sigmaN)
          sb.assym = false;
      }
    }
  }

  fLinearized = true;
  return;
}

void CollieDistribution::getFloatFlagList(vector<string> inputNames, bool* floatMap) {
  for (uint s = 0; s < inputNames.size(); ++s) {
    if (getFloatFlag(inputNames[s]))
      floatMap[s] = true;
  }
  return;
}

bool CollieDistribution::getFloatFlag(string systname) {
  for (int i = 0; i < fSystNames.GetEntriesFast(); ++i)
    if (systname == (string)((TObjString*)fSystNames[i])->GetName()) {
      return fFloatFlag[i]->TestBit(1);
    }

  return false;
}

void CollieDistribution::setFloatFlag(string systname, bool floatIt) {
  for (int i = 0; i < fSystNames.GetEntriesFast(); ++i) {
    if (systname == (string)((TObjString*)fSystNames[i])->GetName()) {
      if (floatIt)
        fFloatFlag[i]->SetBit(1, true);
      else
        fFloatFlag[i]->SetBit(1, false);
      return;
    }
  }
  return;
}

void CollieDistribution::getLogNormalFlagList(vector<string> inputNames, bool* lnMap) {
  for (uint s = 0; s < inputNames.size(); ++s) {
    if (getLogNormalFlag(inputNames[s]))
      lnMap[s] = true;
  }
  return;
}

bool CollieDistribution::getLogNormalFlag(string systname) {

  for (int i = 0; i < fSystNames.GetEntriesFast(); ++i)
    if (systname == (string)((TObjString*)fSystNames[i])->GetName()) {
      return fLogNormalFlag[i]->TestBit(1);
    }

  return false;
}

void CollieDistribution::setLogNormalFlag(string systname, bool floatIt) {
  for (int i = 0; i < fSystNames.GetEntriesFast(); ++i) {
    if (systname == (string)((TObjString*)fSystNames[i])->GetName()) {
      if (floatIt)
        fLogNormalFlag[i]->SetBit(1, true);
      else
        fLogNormalFlag[i]->SetBit(1, false);
      return;
    }
  }
  return;
}

/////// Code to perform a few checks on the stat errors
/////// of the incoming histogram.
/////// 'verbose' controls how much info is printed.
/////// Added 16-3-08 by M. Owen (markowen@fnal.gov)
bool CollieDistribution::checkStats(const TH1* h, int verbose) const {

  if (!h) {
    cerr << "CollieDistribution::checkStats, NULL histogram!" << endl;
    return false;
  }

  float hsumw2n = true;
  if (!h->GetSumw2N()) {
    cerr << "CollieDistribution::checkStats, WARNING: Input hist " << h->GetName();
    cerr << " does not have sum of weight structure - bin stat errors";
    cerr << " may be incorrect" << endl;
    hsumw2n = false;
  }

  double totstaterr_sq = 0.0;
  int nbins = h->GetNbinsX() + 2;
  if (h->GetDimension() > 1)
    nbins += h->GetNbinsY() + 2;
  if (h->GetDimension() > 2)
    nbins += h->GetNbinsZ() + 2;
  if (verbose > 2) {
    cout << "Checking stat errors for " << h->GetName() << endl;
  }
  for (int i = 0; i < nbins; ++i) {
    float err = h->GetBinError(i);
    float cont = h->GetBinContent(i);
    totstaterr_sq += err * err;
    if (err < 0) {
      cerr << "ERROR: Bin " << i << " of " << h->GetName();
      cerr << " has error < 0" << endl;
      cerr << cont << " +- " << err << endl;
      if (hsumw2n)
        return false;
    }

    if (verbose > 2)
      cout << "Bin " << i << " = " << cont << " +- " << err << endl;
    if ((err - cont) > 1e-5 && cont > 1.0) {
      cerr << "ERROR: Bin " << i << " of " << h->GetName();
      cerr << " has larger error than content" << endl;
      cerr << cont << " +- " << err << endl;
      if (hsumw2n)
        return false;
    }

  } //loop on bins

  if (verbose > 1) {
    if (totstaterr_sq == 0) {
      cout << "INFO: " << h->GetName() << " has a stat error of 0" << endl;
    } else
      cout << "INFO: " << h->GetName() << " total = " << h->Integral() << " +- "
           << sqrt(totstaterr_sq) << endl;
  }

  return true;

} //checkStats

bool CollieDistribution::checkStats(const CollieHistogram* h, int verbose) const {

  if (!h) {
    cerr << "CollieDistribution::checkStats, NULL histogram!" << endl;
    return false;
  }

  double totstaterr_sq = 0.0;
  int nbins = h->nBins();
  if (verbose > 2) {
    cout << "Checking stat errors for " << h->GetName() << endl;
  }
  for (int i = 0; i < nbins; ++i) {
    float err = h->BinErr(i);
    float cont = h->InBin(i);
    totstaterr_sq += err * err;
    if (err < 0) {
      cerr << "ERROR: Bin " << i << " of " << h->GetName();
      cerr << " has error < 0" << endl;
      cerr << cont << " +- " << err << endl;
    }

    if (verbose > 2)
      cout << "Bin " << i << " = " << cont << " +- " << err << endl;
    if ((err - cont) > 1e-5 && cont > 1.0) {
      cerr << "ERROR: Bin " << i << " of " << h->GetName();
      cerr << " has larger error than content" << endl;
      cerr << cont << " +- " << err << endl;
    }
  } //loop on bins

  if (verbose > 1) {
    if (totstaterr_sq == 0) {
      cout << "INFO: " << h->GetName() << " has a stat error of 0" << endl;
    } else
      cout << "INFO: " << h->GetName() << " total = " << h->Integral() << " +- "
           << sqrt(totstaterr_sq) << endl;
  }

  return true;
} //checkStats

double CollieDistribution::getNormalizedBinValueVaried(const int i, const int j,
                                                       const double* fluctList) const {
  if (j >= fNYbins)
    return 0;
  if (i < 0 || i >= fNXbins)
    return 0;
  if (!fLinearized) {
    printf("CollieDistribution::getNormalizedBinValueVaried, this dist is not linearized...: %s\n",
           fName.Data());
    return 0;
  }

  //get fluctuated values
  double fval = 1.0;
  double delta = 0;

  //common vars to be used in the loop
  const uint tbin = (fNYbins > 1) ? i + j * fNXbins : i;
  const double inbin = fLinBins[tbin];
  const uint ns = fNsyst;

  for (uint s = 0; s < ns; ++s) {
    //For asymmetric errors, bridge PDF discontinuity with quadratic match
    fval = getSmearedEfficiency(fluctList[fSystIndexOuter[s]], fLinSyst[s][tbin].sigmaP,
                                fLinSyst[s][tbin].sigmaN);

    //Log-Normal PDF estimation
    if (fLinLogN[s])
      fval = exp(fval) - 1.0; //Match Mean

    // sum the changes in rate
    delta += fval;
  }

  if (std::isinf(delta) || std::isnan(delta))
    delta = 0;

  if ((inbin + delta * inbin) < 0)
    return 0.0;
  else
    return inbin + delta * inbin;
}

double CollieDistribution::getNormalizedBinValueVaried(int i, const double* fluctList) const {
  if (i < 0 || i >= fNXbins)
    return 0;

  if (!fLinearized) {
    printf("CollieDistribution::getNormalizedBinValueVaried, this dist is not linearized...: %s\n",
           fName.Data());
    return 0;
  }

  //get fluctuated values
  double fval = 1.0;
  double delta = 0;

  //common vars to be used in the loop
  const uint tbin = i;
  const double inbin = fLinBins[tbin];
  const uint ns = fNsyst;

  for (uint s = 0; s < ns; ++s) {
    //For asymmetric errors, bridge PDF discontinuity with quadratic match
    fval = getSmearedEfficiency(fluctList[fSystIndexOuter[s]], fLinSyst[s][tbin].sigmaP,
                                fLinSyst[s][tbin].sigmaN);

    //Log-Normal PDF estimation
    if (fLinLogN[s])
      fval = exp(fval) - 1.0; //Match Mean

    // sum the changes in rate
    delta += fval;
  }

  if (std::isinf(delta) || std::isnan(delta))
    delta = 0;
  if ((inbin + delta * inbin) < 0)
    return 0.0;
  else
    return inbin + delta * inbin;

} // getNormalizedBinValueVaried(i, fluctList)

void CollieDistribution::addBinEfficiencies(const double* fluctMap, const double sigScale,
                                            double* sig, const uint n_bins) {

  uint s = 0;
  const double scale = sigScale * fEfficiency;
  for (double *p = fLinBins + 0, *end = fLinBins + fTrueBinCount, *dest = &fEfficiencySums[0];
       p != end; ++p, ++dest) {
    (*dest) = (*p) * scale;
  }

  // s loops over the active s-parameters for this distribution:
  for (s = 0; s < (uint)fNsyst; ++s) {
    // To optimize speed per bin, do as much computation as we can outside
    // the bin loop, and select a specific function to do the bin loop
    // with as few conditionals as possible.
    //
    if (fLinLogN[s]) {
      binEfficienciesQbridgeLinLogN(&fEfficiencySums[0], scale, fluctMap[fSystIndexOuter[s]],
                                    fLinSyst[s]);
    } else {
      binEfficienciesQbridge(&fEfficiencySums[0], scale, fluctMap[fSystIndexOuter[s]], fLinSyst[s]);
    }
    // The above call could be pulled out into a set of funtions controlled by
    // a switch, if we decide to allow each distribution (or even each dist, s)
    // to have its own policy for bridging the asymmetric sensitivities.  Here
    // we hardwire just one choice, so we don't need that switch.  But we do
    // need the conditional on fLinLogN.
  }

  for (s = 0; s != (uint)fTrueBinCount; ++s) {
    if (fEfficiencySums[s] > 0)
      sig[s + n_bins] += fEfficiencySums[s];
  }

} // addBinEfficiencies

void CollieDistribution::addBinEfficiencies(const int perturbed_s, const double s_value,
                                            const double sigScale, double* sig, const uint n_bins) {

  const int s = fSystIndexInner[perturbed_s];
  std::vector<double> sums(fTrueBinCount);
  double* es = fEfficiencySums.data(); // If s is not active for this dist,
                                       // just contribute base sums
  if (s >= 0) {
    typedef std::vector<InfoForEfficiencyCalculation>::iterator Iter;
    es = sums.data();
    int idx = 0;
    for (Iter be = fLinSyst[s].begin(), en = fLinSyst[s].end(); be != en; ++be, ++es) {
      (*es) = be->exclusionSum;
      ++idx;
    }

    if (fLinLogN[s]) {
      binEfficienciesQbridgeLinLogN(sums.data(), fEfficiency * sigScale, s_value, fLinSyst[s]);
    } else {
      binEfficienciesQbridge(sums.data(), fEfficiency * sigScale, s_value, fLinSyst[s]);
    }

    es = sums.data();
  }

  for (int b = 0; b != fTrueBinCount; ++b) {
    if (es[b] > 0)
      sig[b + n_bins] += es[b];
  }
} // addBinEfficiencies (shortcut case)

void CollieDistribution::binEfficienciesQbridge(double* sum, const double scale,
                                                const double active_s,
                                                std::vector<InfoForEfficiencyCalculation>& e) {
  double partialSum;
  typedef std::vector<InfoForEfficiencyCalculation>::iterator InfoIt;
  const InfoIt infoEnd = e.end();
  double* nom = &fLinBins[0];

  if (active_s < 0) {
    for (InfoIt info = e.begin(); info != infoEnd; ++info, ++sum, ++nom) {
      //pos/neg bridge functionality
      if (info->assym)
        partialSum = (*nom) * scale * getSmearedEfficiencyN(active_s, info->sigmaP, info->sigmaN);
      else
        partialSum = (*nom) * scale * active_s * info->sigmaN;

      info->baseDeltaEfficiency = partialSum;
      (*sum) += partialSum;
    }
  } else {
    for (InfoIt info = e.begin(); info != infoEnd; ++info, ++sum, ++nom) {
      //pos/neg bridge functionality
      if (info->assym)
        partialSum = (*nom) * scale * getSmearedEfficiencyP(active_s, info->sigmaP, info->sigmaN);
      else
        partialSum = (*nom) * scale * active_s * info->sigmaP;

      info->baseDeltaEfficiency = partialSum;
      (*sum) += partialSum;
    }
  }

} // binEfficienciesQbridge()

void CollieDistribution::binEfficienciesQbridgeLinLogN(
    double* sum, const double scale, const double active_s,
    std::vector<InfoForEfficiencyCalculation>& e) {

  double partialSum;
  typedef std::vector<InfoForEfficiencyCalculation>::iterator InfoIt;
  InfoIt infoEnd = e.end();
  double* nom = &fLinBins[0];

  if (active_s < 0.0) {
    for (InfoIt info = e.begin(); info != infoEnd; ++info, ++sum) {
      if (info->assym)
        partialSum = (*nom) * scale *
                     (std::exp(getSmearedEfficiencyN(active_s, info->sigmaP, info->sigmaN)) - 1.0);
      else
        partialSum = (*nom) * scale * (std::exp(active_s * info->sigmaN) - 1.0);

      info->baseDeltaEfficiency = partialSum;
      (*sum) += partialSum;
    }
  } else { // active_s >= 0
    for (InfoIt info = e.begin(); info != infoEnd; ++info, ++sum) {
      if (info->assym) {
        partialSum = (*nom) * scale *
                     (std::exp(getSmearedEfficiencyP(active_s, info->sigmaP, info->sigmaN)) - 1.0);
      } else
        partialSum = (*nom) * scale * (std::exp(active_s * info->sigmaP) - 1.0);

      info->baseDeltaEfficiency = partialSum;
      (*sum) += partialSum;
    }
  }
} // binEfficienciesQbridgeLinLogN()

void CollieDistribution::prepareExclusionSums() {

  typedef std::vector<InfoForEfficiencyCalculation>::iterator InfoIt;
  InfoIt info, infoEnd;

  for (int s = 0; s < fNsyst; ++s) {
    std::vector<InfoForEfficiencyCalculation>& e = fLinSyst[s];
    infoEnd = e.end();
    std::vector<double>::const_iterator fes = fEfficiencySums.begin();

    for (info = e.begin(); info != infoEnd; ++info, ++fes) {
      info->exclusionSum = (*fes) - info->baseDeltaEfficiency;
    }
  }
} // prepareExclusionSums

void CollieDistribution::print(map<string, int>& count, map<string, double>& posSyst,
                               map<string, double>& negSyst) const {

  printf("Total evts: %.2f\n", getEfficiency());
  printf("%d systematics:\n", getNsystematics());
  for (int s = 0; s < getNsystematics(); ++s) {
    double meanP = 0, meanN = 0;
    double maxP = 0, maxN = 0;
    int nBins = 0;
    if (fNYbins > 1) {
      for (int x = 1; x <= fNXbins; ++x) {
        for (int y = 1; y <= fNYbins; ++y) {
          if (getEfficiency(x, y) < 1e-6)
            continue;
          ++nBins;
          double vp = getPositiveSyst2D(s)->GetBinContent(x, y);
          double vn = getNegativeSyst2D(s)->GetBinContent(x, y);
          meanP += vp;
          meanN += vn;
          if (fabs(vp) > fabs(maxP))
            maxP = vp;
          if (fabs(vn) > fabs(maxN))
            maxN = vn;
        }
      }
    } else {
      for (int x = 1; x < fNXbins; ++x) {
        if (getEfficiency(x) < 1e-6)
          continue;
        ++nBins;
        double vp = getPositiveSyst(s)->GetBinContent(x);
        double vn = getNegativeSyst(s)->GetBinContent(x);
        meanP += vp;
        meanN += vn;
        if (fabs(vp) > fabs(maxP))
          maxP = vp;
        if (fabs(vn) > fabs(maxN))
          maxN = vn;
      }
    }
    if (nBins > 0) {
      meanP /= nBins;
      meanN /= nBins;
    }

    if (fabs(maxP - meanP) > fabs(1e-2 * meanP) || fabs(maxN - meanN) > fabs(1e-2 * meanN))
      printf("Shape systematic %d, %s: (Mean/Max) Pos : %.4f/%.4f, Neg: %.4f/%.4f\n", s,
             getSystName(s).c_str(), meanP, maxP, meanN, maxN);

    else
      printf("Flat systematic %d, %s: Pos: %.4f, Neg: %.4f\n", s, getSystName(s).c_str(), meanP,
             meanN);

    if (getEfficiency() > 1e-5) {
      if (posSyst.find(getSystName(s)) == posSyst.end()) {
        posSyst[getSystName(s)] = meanP;
        negSyst[getSystName(s)] = meanN;
        count[getSystName(s)] = 1;
      } else {
        count[getSystName(s)]++;
        int tc = count[getSystName(s)];
        posSyst[getSystName(s)] = posSyst[getSystName(s)] * (tc - 1) + meanP;
        negSyst[getSystName(s)] = negSyst[getSystName(s)] * (tc - 1) + meanN;
        posSyst[getSystName(s)] /= tc;
        negSyst[getSystName(s)] /= tc;
      }
    }
  }
  return;
}
