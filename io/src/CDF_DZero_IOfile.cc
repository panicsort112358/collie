#include <TSystem.h>
#include "CDF_DZero_IOfile.hh"

ClassImp(CDF_DZero_IOfile);

//-------------------------------------------
CDF_DZero_IOfile::CDF_DZero_IOfile() {
  fMutable = false;
  fDirectory = NULL;
}

CDF_DZero_IOfile::CDF_DZero_IOfile(const char* name) {
  fDirectory = NULL;
  fChannelName = name;
  fMutable = true;
  fNindepVars = -1;
  std::string a(fChannelName);
  a += "-CHANNEL";
  fName = a.c_str();
  fTitle = fName;
  fCreationComputer = gSystem->HostName();
  fCreationUser = gSystem->Getenv("LOGNAME");
  fCreationTime.Set();
}

CDF_DZero_IOfile* CDF_DZero_IOfile::createChannel(const char* name, TFile* rootFile,
                                                  const char* options) {
  std::string a("Data for ");
  a += name;
  if (rootFile == NULL)
    return NULL;
  rootFile->cd("/");
  TDirectory* td = rootFile->mkdir(name, a.c_str());
  if (td == 0)
    return NULL;

  CDF_DZero_IOfile* c = new CDF_DZero_IOfile(name);
  c->fDirectory = td;
  a = "/";
  a += name;
  rootFile->cd(a.c_str());

  return c;
}

CDF_DZero_IOfile* CDF_DZero_IOfile::loadChannel(const char* channelName, TFile* rootFile) {
  std::string a("/");
  a += channelName;
  a += '/';
  a += channelName;
  a += "-CHANNEL";
  CDF_DZero_IOfile* c = (CDF_DZero_IOfile*)rootFile->Get(a.c_str());
  if (c == NULL)
    return NULL;
  c->associateFile(rootFile);
  return c;
}

void CDF_DZero_IOfile::associateFile(TFile* tf) {
  std::string a(fChannelName);
  fDirectory = (TDirectory*)tf->Get(a.c_str());
  return;
}

void CDF_DZero_IOfile::Store() {
  if (!fMutable)
    return;
  doCD();
  this->storePoints();
  Write();
}

void CDF_DZero_IOfile::setNIndepVariables(int nVariables) {
  if (!fMutable || fNindepVars >= 1)
    return;
  fNindepVars = nVariables;
  return;
}

void CDF_DZero_IOfile::storePoints() {

  for (int i = 0; i < fPoints_var1.GetSize(); i++) {
    if (fPoints_var2.GetSize() > 0) {
      for (int j = 0; j < fPoints_var2.GetSize(); j++) {
        if (fPoints_var3.GetSize() > 0) {
          for (int k = 0; k < fPoints_var3.GetSize(); i++) {
            this->storePoint(this->getPoint(fPoints_var1[i], fPoints_var2[j], fPoints_var3[k]));
          }
        } else
          this->storePoint(this->getPoint(fPoints_var1[i], fPoints_var2[j]));
      }
    } else {
      this->storePoint(this->getPoint(fPoints_var1[i]));
    }
  }

  return;
}

void CDF_DZero_IOfile::storePoint(CDF_DZero_IOpoint* iop) {
  if (iop == NULL) {
    printf("\n\nCDF_DZero_IOfile::storePoint ==> This point is NULL!\n\n\n\n");
    return;
  }
  fDirectory->cd();
  iop->Write(); // store it...
  // store this...
  //  fDirectory->Write();
  return;
}

static void insertAt(TArrayI& vec, int val, int pos) {
  if (vec.GetSize() <= pos)
    vec.Set(pos + 1);
  else {
    vec.Set(vec.GetSize() + 1);
    // copy the integers...
    memmove(&(vec.GetArray()[pos + 1]), &(vec.GetArray()[pos]),
            (vec.GetSize() - pos - 1) * sizeof(int));
  }
  vec[pos] = val;
  return;
}

void CDF_DZero_IOfile::logPoint(CDF_DZero_IOpoint* iop) {

  if (iop == NULL) {
    printf("\n\nCDF_DZero_IOfile::logPoint ==> This point is NULL!\n\n\n\n");
    return;
  }

  fDirectory->Append(iop);

  // find the place in the lists...
  int i;
  for (i = 0; i < fPoints_var1.GetSize() && fPoints_var1[i] < iop->getVar1(); i++)
    ;
  if (i < fPoints_var1.GetSize() && fPoints_var1[i] == iop->getVar1()) {
    if (fNindepVars == 1)
      i = -1;
    else {
      for (; i < fPoints_var2.GetSize() && fPoints_var2[i] < iop->getVar2(); i++)
        ;
      if (i < fPoints_var2.GetSize() && fPoints_var2[i] == iop->getVar2()) {
        if (fNindepVars == 2)
          i = -1;
        else {
          for (; i < fPoints_var3.GetSize() && fPoints_var3[i] < iop->getVar3(); i++)
            ;
          if (i < fPoints_var3.GetSize() && fPoints_var3[i] == iop->getVar3())
            i = -1;
        }
      }
    }
  }

  if (i == -1) {
    printf("\n\nCDF_DZero_IOfile::logPoint ==> Could not find a spot to log!\n\n\n\n");
    return;
  }

  insertAt(fPoints_var1, iop->getVar1(), i);
  if (fNindepVars > 1)
    fPoints_var2.AddAt(iop->getVar2(), i);
  if (fNindepVars > 2)
    fPoints_var3.AddAt(iop->getVar3(), i);

  return;
}

CDF_DZero_IOpoint* CDF_DZero_IOfile::createPoint(int var1, int var2, int var3) {
  return new CDF_DZero_IOpoint(var1, var2, var3);
}

CDF_DZero_IOpoint* CDF_DZero_IOfile::getPoint(int var1, int var2, int var3) {
  if (fDirectory == NULL) {
    printf("CDF_DZero_IOfile::getPoint, NULL directory\n");
    return NULL;
  }

  std::string a(CDF_DZero_IOpoint::buildName(fChannelName, var1, var2, var3, fNindepVars));

  CDF_DZero_IOpoint* mp = (CDF_DZero_IOpoint*)fDirectory->Get(a.c_str());
  return mp;
}

void CDF_DZero_IOfile::getIndependentVariables(int idx, int& var1, int& var2, int& var3) {
  var1 = -1;
  var2 = -1;
  var3 = -1;

  if (idx >= 0 && idx < fPoints_var1.GetSize()) {
    var1 = fPoints_var1[idx];
    if (fPoints_var2.GetSize() > idx)
      var2 = fPoints_var2[idx];
    if (fPoints_var3.GetSize() > idx)
      var3 = fPoints_var3[idx];
  } else
    printf("CDF_DZero_IOfile::getIndependentVariables, This point doesn't exist: %d (%d)\n", idx,
           fPoints_var1.GetSize());
  return;
}

void CDF_DZero_IOfile::print() {

  printf("\nCDF/DZero Common IO File\n");
  printf("Channel name: %s\n", fChannelName.Data());
  printf("Created %s by  %s using %s\n", fCreationTime.AsString(), fCreationUser.Data(),
         fCreationComputer.Data());

  printf("\n%d Independent Variable(s), total for each possible: (%d, %d, %d)\n", fNindepVars,
         fPoints_var1.GetSize(), fPoints_var2.GetSize(), fPoints_var3.GetSize());

  map<string, int> nSyst;
  map<string, double> posSyst;
  map<string, double> negSyst;

  for (int i = 0; i < fPoints_var1.GetSize(); i++) {
    if (fPoints_var2.GetSize() > 0) {
      for (int j = 0; j < fPoints_var2.GetSize(); j++) {
        if (fPoints_var3.GetSize() > 0) {
          for (int k = 0; k < fPoints_var3.GetSize(); i++) {
            this->getPoint(fPoints_var1[i], fPoints_var2[j], fPoints_var3[k])
                ->print(nSyst, posSyst, negSyst);
          }
        } else
          this->getPoint(fPoints_var1[i], fPoints_var2[j])->print(nSyst, posSyst, negSyst);
      }
    } else
      this->getPoint(fPoints_var1[i])->print(nSyst, posSyst, negSyst);
  }

  map<string, double>::iterator iterP = posSyst.begin();
  map<string, double>::iterator iterN = negSyst.begin();
  map<string, int>::iterator iterT = nSyst.begin();

  printf("Systematics summary for channel %s\n", fChannelName.Data());
  printf("Reporting average mean values over all points and distributions:\n");
  for (iterP = posSyst.begin(); iterP != posSyst.end(); iterP++) {
    printf("Systematic %s: %.3f, %.3f\n", iterP->first.c_str(), iterP->second / iterT->second,
           iterN->second / iterT->second);
    iterN++;
    iterT++;
  }

  return;
}
