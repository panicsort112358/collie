#include <CollieMasspoint.hh>
#include <CollieChannel.hh>

ClassImp(CollieMasspoint);

CollieMasspoint::~CollieMasspoint() {
  if (fDataDistribution != NULL)
    delete fDataDistribution;
  if (fDataEventList != NULL)
    delete fDataEventList;

  fSignals.Delete();
  fBackgrounds.Delete();
}

CollieMasspoint::CollieMasspoint() {
  pChannel = NULL;
  fDataDistribution = NULL;
  fDataEventList = NULL;
  fMutable = false;
}

CollieMasspoint::CollieMasspoint(CollieChannel* theChannel, int var1, int var2, int var3) {

  fDataDistribution = NULL;
  fDataEventList = NULL;

  string name =
      buildName(theChannel->getChannelName(), var1, var2, var3, theChannel->getNIndepVariables());
  fName = name.c_str();
  fTitle = fName;

  fIndepVar1 = var1;
  fIndepVar2 = var2;
  fIndepVar3 = var3;
  pChannel = theChannel;
  fMutable = true;
}

void CollieMasspoint::ConfigEventList(int ndoubleVars, int nintVars) {
  if (!fMutable)
    return;
  int dims = (fNYbins > 1) ? 2 : 1;
  std::string name = GetName();
  name += ":";
  name += "DATA";
  fDataEventList = new CollieEventList(name.c_str(), dims, ndoubleVars, nintVars);
  fDataEventList->p_MassPoint = this;
}

void CollieMasspoint::BookData() {
  if (!fMutable)
    return;
  if (fDataDistribution != NULL)
    return;
  fDataDistribution = BookDataIntl();
}

CollieDistribution* CollieMasspoint::BookDataIntl() const {
  std::string name = GetName();
  name += ":";
  name += "DATA";
  CollieDistribution* d = new CollieDistribution(name.c_str(), fNXbins, fMinX, fMaxX, fNYbins,
                                                 fMinY, fMaxY, pChannel->getNModels());
  d->p_MassPoint = this;
  return d;
}

void CollieMasspoint::Book(int nx, double minX, double maxX) {
  Book(nx, minX, maxX, 0, -1, -1);
}

void CollieMasspoint::Book(int nx, double minX, double maxX, int ny, double minY, double maxY) {
  if (!fMutable)
    return;
  fNXbins = nx;
  fNYbins = (ny < 1) ? 1 : ny;
  fMaxX = maxX;
  fMinX = minX;
  fMaxY = maxY;
  fMinY = minY;

  if (pChannel == NULL)
    return;

  /// create the distributions
  for (int i = 0; i < pChannel->getNSignals(); i++) {
    std::string name = pChannel->getSignalName(i);
    CollieDistribution* d = new CollieDistribution(name.c_str(), nx, minX, maxX, ny, minY, maxY,
                                                   pChannel->getNModels());
    fSignals.AddAtAndExpand(d, i);
  }
  for (int i = 0; i < pChannel->getNBackgrounds(); i++) {
    std::string name = pChannel->getBackgroundName(i);
    CollieDistribution* d = new CollieDistribution(name.c_str(), nx, minX, maxX, ny, minY, maxY,
                                                   pChannel->getNModels());
    fBackgrounds.AddAtAndExpand(d, i);
  }
}

std::string CollieMasspoint::buildName(const char* channelname, int var1, int var2, int var3,
                                       int nvar) {
  char tool[1024];

  if (nvar >= 3)
    sprintf(tool, "%s(%d,%d,%d)", channelname, var1, var2, var3);
  else if (nvar == 2)
    sprintf(tool, "%s(%d,%d)", channelname, var1, var2);
  else if (nvar <= 1)
    sprintf(tool, "%s(%d)", channelname, var1);

  return tool;
}

void CollieMasspoint::Store() {
  Write();
}

void CollieMasspoint::logPoint() {
  if (pChannel != NULL)
    pChannel->logPoint(this);
  else
    printf("CollieMasspoint::logPoint, NULL channel!\n");
}

const CollieDistribution* CollieMasspoint::getSignalDist(int i) const {
  if (i < 0 || i >= getNSignalDists())
    return NULL;
  CollieDistribution* dist = (CollieDistribution*)fSignals[i];
  dist->p_MassPoint = this;
  return dist;
}
const CollieDistribution* CollieMasspoint::getBkgdDist(int i) const {
  if (i < 0 || i >= getNBkgdDists())
    return NULL;
  CollieDistribution* dist = (CollieDistribution*)fBackgrounds[i];
  dist->p_MassPoint = this;
  return dist;
}
CollieDistribution* CollieMasspoint::getSignalDistMutable(int i) {
  if (!fMutable)
    return NULL;
  if (i < 0 || i >= getNSignalDists())
    return NULL;
  CollieDistribution* dist = (CollieDistribution*)fSignals[i];
  dist->p_MassPoint = this;
  return dist;
}

CollieDistribution* CollieMasspoint::getBkgdDistMutable(int i) {
  if (!fMutable)
    return NULL;
  if (i < 0 || i >= getNBkgdDists())
    return NULL;
  CollieDistribution* dist = (CollieDistribution*)fBackgrounds[i];
  dist->p_MassPoint = this;
  return dist;
}

CollieDistribution* CollieMasspoint::getDataDistMutable() {
  if (!fMutable || fDataDistribution == NULL)
    return NULL;
  fDataDistribution->p_MassPoint = this;
  return fDataDistribution;
}

const CollieDistribution* CollieMasspoint::getDataDist() const {
  if (fDataDistribution != NULL)
    fDataDistribution->p_MassPoint = this;
  else if (fDataEventList != NULL) {
    CollieDistribution* d = BookDataIntl(); // setup a distribution
    d->p_MassPoint = this;
    fDataEventList->fillDistribution(d);
    ((CollieMasspoint*)this)->fDataDistribution = d; // "violate" the const
    return d;
  }
  return fDataDistribution;
}

CollieEventList* CollieMasspoint::getDataEventListMutable() {
  if (!fMutable || fDataEventList == NULL)
    return NULL;
  fDataEventList->p_MassPoint = this;
  return fDataEventList;
}

const CollieEventList* CollieMasspoint::getDataEventList() const {
  if (fDataEventList != NULL)
    fDataEventList->p_MassPoint = this;
  return fDataEventList;
}

void CollieMasspoint::scaleSignal(double sf) {
  CollieDistribution* cd;
  for (int i = 0; i < this->getNSignalDists(); i++) {
    cd = this->getSignalDistMutable(i);
    if (cd != NULL) {
      cd->setMutable(true);
      cd->setEfficiency(sf * cd->getEfficiency());
    }
  }
  return;
}

void CollieMasspoint::scaleBackground(int idx, double sf) {
  CollieDistribution* cd;
  if (idx < 0) {
    for (int i = 0; i < this->getNBkgdDists(); i++) {
      cd = this->getBkgdDistMutable(i);
      if (cd != NULL) {
        cd->setMutable(true);
        cd->setEfficiency(sf * cd->getEfficiency());
      }
    }
  } else {
    cd = this->getBkgdDistMutable(idx);
    if (cd != NULL) {
      cd->setMutable(true);
      cd->setEfficiency(sf * cd->getEfficiency());
    }
  }
  return;
}

void CollieMasspoint::scaleData(double sf) {
  CollieDistribution* cd = this->getDataDistMutable();
  if (cd != NULL) {
    cd->setMutable(true);
    cd->setEfficiency(sf * cd->getEfficiency());
  }
  return;
}

int CollieMasspoint::lookupBkgd(const char* name) const {
  if (pChannel == NULL)
    return -1;
  for (int i = 0; i < pChannel->getNBackgrounds(); i++)
    if (!strcmp(name, pChannel->getBackgroundName(i)))
      ;
  return -1;
}

int CollieMasspoint::lookupSignal(const char* name) const {
  if (pChannel == NULL)
    return -1;
  for (int i = 0; i < pChannel->getNSignals(); i++)
    if (!strcmp(name, pChannel->getSignalName(i)))
      ;
  return -1;
}

void CollieMasspoint::print(map<string, int>& count, map<string, double>& posSyst,
                            map<string, double>& negSyst) {

  printf("\nCollie Masspoint: %s\n", fName.Data());
  printf("Vars (%d, %d, %d)\n", fIndepVar1, fIndepVar2, fIndepVar3);
  float sum = 0;

  for (int i = 0; i < getNSignalDists(); i++) {
    const CollieDistribution* sig = getSignalDist(i);
    if (sig == NULL) {
      printf("CDF_DZero_IOpoint::print, NULL signal!\n");
      continue;
    }
    printf("\nSignal %d\n", i);

    sig->print(count, posSyst, negSyst);
  }
  for (int i = 0; i < getNBkgdDists(); i++) {
    const CollieDistribution* bkg = getBkgdDist(i);
    if (bkg == NULL) {
      printf("CDF_DZero_IOpoint::print, NULL background!\n");
      continue;
    }
    printf("\nBkgd %d: %s\n", i, bkg->GetName());
    bkg->print(count, posSyst, negSyst);
    sum += bkg->getEfficiency();
  }

  printf("\nTotal Bkgd: %.2f\n", sum);
  printf("Total Data: %.1f\n", fDataDistribution->getEfficiency());
}
