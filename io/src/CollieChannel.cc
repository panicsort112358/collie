#include <TSystem.h>
#include <CollieIterator.hh>
#include <CollieChannel.hh>

ClassImp(CollieChannel);

CollieChannel::CollieChannel() {
  fMutable = false;
  fDirectory = NULL;
}

CollieChannel::CollieChannel(const char* name) {

  fChannelName = name;
  fMutable = true;
  fNindepVars = -1;
  fNsignals = -1;
  fNbackgrounds = -1;
  fNmodels = -1;
  std::string a(fChannelName);
  a += "-CHANNEL";
  fName = a.c_str();
  fTitle = fName;
  fCollieVersion = "V00-04-13";
  fCreationComputer = gSystem->HostName();
  fCreationUser = gSystem->Getenv("LOGNAME");
  fCreationTime.Set();
}

CollieChannel* CollieChannel::createChannel(const char* name, TFile* rootFile,
                                            const char* options) {
  std::string a("Collie data for ");
  a += name;
  if (rootFile == NULL)
    return NULL;
  rootFile->cd("/");
  TDirectory* td = rootFile->mkdir(name, a.c_str());
  if (td == 0)
    return NULL;

  CollieChannel* c = new CollieChannel(name);
  c->fDirectory = td;
  a = "/";
  a += name;
  rootFile->cd(a.c_str());

  return c;
}

CollieChannel* CollieChannel::loadChannel(const char* channelName, TFile* rootFile) {
  std::string a("/");
  a += channelName;
  a += '/';
  a += channelName;
  a += "-CHANNEL";
  CollieChannel* c = (CollieChannel*)rootFile->Get(a.c_str());
  if (c == NULL)
    return NULL;
  c->associateFile(rootFile);
  return c;
}

void CollieChannel::associateFile(TFile* tf) {
  std::string a(fChannelName);
  fDirectory = (TDirectory*)tf->Get(a.c_str());
  return;
}

/// copy constructor
CollieChannel* CollieChannel::CopyChannel(TFile* rootFile) {
  std::string a("Collie data for ");
  a += this->fChannelName;
  if (rootFile == NULL)
    return NULL;
  rootFile->cd("/");
  TDirectory* td = rootFile->mkdir(this->fChannelName, a.c_str());
  if (td == 0)
    return NULL;

  CollieChannel* c = new CollieChannel(this->fChannelName);
  c->fDirectory = td;
  doCD();

  c->fName = this->fName;
  c->fTitle = this->fTitle;
  c->fLuminosity = this->fLuminosity;

  c->fNindepVars = this->fNindepVars;
  c->fIndepVarName[0] = this->fIndepVarName[0];
  c->fIndepVarName[1] = this->fIndepVarName[1];
  c->fIndepVarName[2] = this->fIndepVarName[2];

  c->fNsignals = this->fNsignals;
  this->fSignalNames = TObjArray(c->fSignalNames);

  c->fNbackgrounds = this->fNbackgrounds;
  this->fBkgdNames = TObjArray(c->fBkgdNames);

  c->fNmodels = this->fNmodels;
  this->fModelNames = TObjArray(c->fModelNames);

  for (int i = 0; i < fMasspoints_var1.GetSize(); i++) {

    CollieMasspoint* mp = 0;
    if (fNindepVars == 1)
      mp = (CollieMasspoint*)this->getMasspoint(fMasspoints_var1[i])->Clone();
    else if (fNindepVars == 2)
      mp = (CollieMasspoint*)this->getMasspoint(fMasspoints_var1[i], fMasspoints_var2[i])->Clone();
    else if (fNindepVars == 3)
      mp = (CollieMasspoint*)this
               ->getMasspoint(fMasspoints_var1[i], fMasspoints_var2[i], fMasspoints_var3[i])
               ->Clone();

    if (mp != NULL) {
      mp->setMutable(true);
      mp->pChannel = c;
      c->logPoint(mp);
    } else
      printf("CollieChannel::CopyChannel, NULL masspoint!\n");
  }

  return c;
}

CollieChannel::~CollieChannel() {

  //  if(fDirectory!=NULL) fDirectory->DeleteAll();// fDirectory=NULL;

  fSignalNames.Delete();
  fBkgdNames.Delete();
  fModelNames.Delete();
}

CollieIterator* CollieChannel::getIterator() {
  if (fNindepVars == 1)
    return new CollieIterator(this, fMasspoints_var1.GetSize(), fMasspoints_var1.GetArray(), NULL,
                              NULL);
  else if (fNindepVars == 2)
    return new CollieIterator(this, fMasspoints_var1.GetSize(), fMasspoints_var1.GetArray(),
                              fMasspoints_var2.GetArray(), NULL);
  else if (fNindepVars == 3)
    return new CollieIterator(this, fMasspoints_var1.GetSize(), fMasspoints_var1.GetArray(),
                              fMasspoints_var2.GetArray(), fMasspoints_var3.GetArray());
  return NULL;
}

void CollieChannel::setNIndepVariables(int nVariables) {
  if (!fMutable || fNindepVars >= 1)
    return;
  fNindepVars = nVariables;
  return;
}

void CollieChannel::setNSignals(int nSignals) {
  if (!fMutable || fNsignals >= 1)
    return;
  fNsignals = nSignals;
  for (int i = 0; i < nSignals; i++)
    fSignalNames.AddAtAndExpand(new TObjString(), i);
  fSignalNames.Compress();
  return;
}

void CollieChannel::setNBackgrounds(int n) {
  if (!fMutable || fNbackgrounds >= 1)
    return;
  fNbackgrounds = n;
  for (int i = 0; i < n; i++)
    fBkgdNames.AddAtAndExpand(new TObjString(), i);
  fBkgdNames.Compress();
  return;
}

void CollieChannel::setNmodels(int nModels) {
  if (!fMutable || fNmodels >= 0)
    return;
  fNmodels = nModels;
  for (int i = 0; i < nModels; i++)
    fModelNames.AddAtAndExpand(new TObjString(), i);
  fModelNames.Compress();
  return;
}

CollieMasspoint* CollieChannel::createMasspoint(int var1, int var2, int var3) {
  fDirectory->cd();
  return new CollieMasspoint(this, var1, var2, var3);
}

CollieMasspoint* CollieChannel::getMasspoint(int var1, int var2, int var3) {
  //  fDirectory->cd();
  std::string a(CollieMasspoint::buildName(fChannelName, var1, var2, var3, fNindepVars));

  if (fDirectory == NULL) {
    printf("CollieChannel::getMasspoint, NULL directory (%d, %d, %d)!\n", var1, var2, var3);
    return NULL;
  }

  CollieMasspoint* mp = (CollieMasspoint*)fDirectory->Get(a.c_str());
  if (mp != NULL)
    mp->pChannel = this;
  else
    printf("CollieChannel::getMasspoint, NULL mass point (%d, %d, %d)!\n", var1, var2, var3);
  return mp;
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

void CollieChannel::Store() {
  if (!fMutable)
    return;
  doCD();
  this->storePoints();
  Write();
}

void CollieChannel::storePoints() {

  int pts = 0;
  for (int i = 0; i < fMasspoints_var1.GetSize(); i++) {
    if (fNindepVars == 1)
      this->storePoint(this->getMasspoint(fMasspoints_var1[i]));
    else if (fNindepVars == 2)
      this->storePoint(this->getMasspoint(fMasspoints_var1[i], fMasspoints_var2[i]));
    else if (fNindepVars == 3)
      this->storePoint(
          this->getMasspoint(fMasspoints_var1[i], fMasspoints_var2[i], fMasspoints_var3[i]));
    pts++;
  }
  printf("Stored %d mass point(s)\n", pts);

  return;
}

void CollieChannel::storePoint(CollieMasspoint* cmp) {
  if (cmp == NULL) {
    printf("\n\nCollieChannel::storePoint, This mass point is NULL!\n\n\n\n");
    return;
  }
  fDirectory->cd();
  //  cmp->Write(); // store it...
  cmp->Store(); // store it...
  // store this...
  //  fDirectory->Write();
  return;
}

void CollieChannel::logPoint(CollieMasspoint* cmp) {

  fDirectory->cd();
  fDirectory->Add(cmp);

  // find the place in the lists...
  int i;
  for (i = 0; i < fMasspoints_var1.GetSize() && fMasspoints_var1[i] < cmp->getVar1(); i++)
    ;
  if (i < fMasspoints_var1.GetSize() && fMasspoints_var1[i] == cmp->getVar1()) {
    if (fNindepVars == 1)
      i = -1;
    else {
      for (; i < fMasspoints_var2.GetSize() && fMasspoints_var2[i] < cmp->getVar2(); i++)
        ;
      if (i < fMasspoints_var2.GetSize() && fMasspoints_var2[i] == cmp->getVar2() &&
          fMasspoints_var1[i] == cmp->getVar1()) {
        if (fNindepVars == 2)
          i = -1;
        else {
          for (; i < fMasspoints_var3.GetSize() && fMasspoints_var3[i] < cmp->getVar3(); i++)
            ;
          if (i < fMasspoints_var3.GetSize() && fMasspoints_var3[i] == cmp->getVar3() &&
              fMasspoints_var2[i] == cmp->getVar2() && fMasspoints_var1[i] == cmp->getVar1())
            i = -1;
        }
      }
    }
  }

  if (i == -1) {
    if (fNindepVars == 1)
      printf("CollieChannel::logPoint, this point already exists: %d\n", cmp->getVar1());
    else if (fNindepVars == 2)
      printf("CollieChannel::logPoint, this point already exists: %d, %d!\n", cmp->getVar1(),
             cmp->getVar2());
    else if (fNindepVars == 3)
      printf("CollieChannel::logPoint, this point already exists: %d, %d, %d!\n", cmp->getVar1(),
             cmp->getVar2(), cmp->getVar3());
    return;
  }

  insertAt(fMasspoints_var1, cmp->getVar1(), i);
  if (fNindepVars > 1)
    insertAt(fMasspoints_var2, cmp->getVar2(), i);
  if (fNindepVars > 2)
    insertAt(fMasspoints_var3, cmp->getVar3(), i);

  return;
}

void CollieChannel::scaleSignal(double sf) {

  for (int i = 0; i < fMasspoints_var1.GetSize(); i++) {
    if (fNindepVars == 1)
      this->getMasspoint(fMasspoints_var1[i])->scaleSignal(sf);
    else if (fNindepVars == 2)
      this->getMasspoint(fMasspoints_var1[i], fMasspoints_var2[i])->scaleSignal(sf);
    else if (fNindepVars == 3)
      this->getMasspoint(fMasspoints_var1[i], fMasspoints_var2[i], fMasspoints_var3[i])
          ->scaleSignal(sf);
  }
  return;
}

void CollieChannel::scaleBackground(int idx, double sf) {

  for (int i = 0; i < fMasspoints_var1.GetSize(); i++) {
    if (fNindepVars == 1)
      this->getMasspoint(fMasspoints_var1[i])->scaleBackground(idx, sf);
    else if (fNindepVars == 2)
      this->getMasspoint(fMasspoints_var1[i], fMasspoints_var2[i])->scaleBackground(idx, sf);
    else if (fNindepVars == 3)
      this->getMasspoint(fMasspoints_var1[i], fMasspoints_var2[i], fMasspoints_var3[i])
          ->scaleBackground(idx, sf);
  }
  return;
}

void CollieChannel::scaleData(double sf) {

  for (int i = 0; i < fMasspoints_var1.GetSize(); i++) {
    if (fNindepVars == 1)
      this->getMasspoint(fMasspoints_var1[i])->scaleData(sf);
    else if (fNindepVars == 2)
      this->getMasspoint(fMasspoints_var1[i], fMasspoints_var2[i])->scaleData(sf);
    else if (fNindepVars == 3)
      this->getMasspoint(fMasspoints_var1[i], fMasspoints_var2[i], fMasspoints_var3[i])
          ->scaleData(sf);
  }
  return;
}

void CollieChannel::print() {

  printf("\nCollie IO File\n");
  printf("Channel name: %s\n", fChannelName.Data());
  printf("Comments: %s\n", fComments.Data());
  printf("Created %s by  %s using %s\n", fCreationTime.AsString(), fCreationUser.Data(),
         fCreationComputer.Data());
  printf("Collie Version: %s\n", fCollieVersion.Data());

  printf("\n%d Independent Variable(s), total for each possible: (%d, %d, %d)\n", fNindepVars,
         fMasspoints_var1.GetSize(), fMasspoints_var2.GetSize(), fMasspoints_var3.GetSize());

  map<string, double> posSyst;
  map<string, double> negSyst;
  map<string, int> count;
  //  map<string,double> posSystM;
  //  map<string,double> negSystM;

  for (int i = 0; i < fMasspoints_var1.GetSize(); i++) {
    if (fMasspoints_var2.GetSize() > 0) {
      for (int j = 0; j < fMasspoints_var2.GetSize(); j++) {
        if (fMasspoints_var3.GetSize() > 0) {
          for (int k = 0; k < fMasspoints_var3.GetSize(); i++) {
            this->getMasspoint(fMasspoints_var1[i], fMasspoints_var2[j], fMasspoints_var3[k])
                ->print(count, posSyst, negSyst);
          }
        } else
          this->getMasspoint(fMasspoints_var1[i], fMasspoints_var2[j])
              ->print(count, posSyst, negSyst);
      }
    } else
      this->getMasspoint(fMasspoints_var1[i])->print(count, posSyst, negSyst);
  }

  map<string, double>::iterator iterP = posSyst.begin();
  map<string, double>::iterator iterN = negSyst.begin();

  printf("Systematics summary for channel %s\n", fChannelName.Data());
  printf("      Reporting average mean deviations over all points and distributions\n");
  for (iterP = posSyst.begin(); iterP != posSyst.end(); iterP++) {
    printf("Systematic %s: %.3f, %.3f\n", iterP->first.c_str(), iterP->second, iterN->second);
    iterN++;
  }
  iterN = negSyst.begin();
  for (iterP = posSyst.begin(); iterP != posSyst.end(); iterP++) {
    if (iterP->second == iterN->second)
      printf("%s & $\\pm %.1f$\n", iterP->first.c_str(), iterP->second * 100);
    else
      printf("%s & $^{+%.1f}_{-%.1f}$\n", iterP->first.c_str(), iterP->second * 100,
             fabs(iterN->second) * 100);
    iterN++;
  }
  return;
}
