#include "CollieEventList.hh"

CollieEventList::CollieEventList() {
  fMutable = false;
  p_MassPoint = NULL;
  fNevents = 0;
  fDimensionality = 0;
  fNdoubles = 0;
  fNints = 0;
}

CollieEventList::~CollieEventList() {}

CollieEventList::CollieEventList(const char* name, int dims, int nF, int nI) {
  fMutable = true;
  fName = name;
  fDimensionality = dims;
  fNdoubles = nF + dims;
  fNints = nI;
  fNevents = 0;
  p_MassPoint = NULL;

  int i = 0;
  char temp[32];
  fValueNames.AddAtAndExpand(new TObjString("X"), i++);
  if (fDimensionality > 1)
    fValueNames.AddAtAndExpand(new TObjString("Y"), i++);
  for (int j = 0; j < nF; j++) {
    sprintf(temp, "Double%d", j);
    fValueNames.AddAtAndExpand(new TObjString(temp), i++);
  }
  for (int j = 0; j < nI; j++) {
    sprintf(temp, "Integer%d", j);
    fValueNames.AddAtAndExpand(new TObjString(temp), i++);
  }
}

void CollieEventList::addEvent(double x, double* dvars, int* ivars) {
  if (fDimensionality != 1 || !fMutable)
    return;
  // expand the arrays
  fDoubles.Set(fDoubles.GetSize() + fNdoubles);
  if (fNints > 0)
    fIntegers.Set(fIntegers.GetSize() + fNints);

  if (fNevents > 0) {
    // move/copy the required blocks around...
    for (int i = fNdoubles - 1; i > 0; i--)
      memmove(fDoubles.GetArray() + fNevents * i + i, fDoubles.GetArray() + fNevents * i,
              fNevents * sizeof(double));
    for (int i = fNints - 1; i > 0; i--)
      memmove(fIntegers.GetArray() + fNevents * i + i, fIntegers.GetArray() + fNevents * i,
              fNevents * sizeof(int));
  }

  // store the values...
  fDoubles[fNevents] = x;
  for (int i = 1; i < fNdoubles; i++)
    fDoubles[fNevents + (fNevents + 1) * i] = dvars[i - 1];
  for (int i = 0; i < fNints; i++)
    fIntegers[fNevents + (fNevents + 1) * i] = ivars[i];

  fNevents++;
}

void CollieEventList::addEvent(double x, double y, double* dvars, int* ivars) {
  if (fDimensionality != 2 || !fMutable)
    return;
  // expand the arrays
  fDoubles.Set(fDoubles.GetSize() + fNdoubles);
  if (fNints > 0)
    fIntegers.Set(fIntegers.GetSize() + fNints);

  if (fNevents > 0) {
    // move/copy the required blocks around...
    for (int i = fNdoubles - 1; i > 0; i--)
      memmove(fDoubles.GetArray() + fNevents * i + i, fDoubles.GetArray() + fNevents * i,
              fNevents * sizeof(double));
    for (int i = fNints - 1; i > 0; i--)
      memmove(fIntegers.GetArray() + fNevents * i + i, fIntegers.GetArray() + fNevents * i,
              fNevents * sizeof(int));
  }

  // store the values...
  fDoubles[fNevents] = x;
  fDoubles[fNevents + fNevents + 1] = y;
  for (int i = 2; i < fNdoubles; i++)
    fDoubles[fNevents + (fNevents + 1) * i] = dvars[i - 2];
  for (int i = 0; i < fNints; i++)
    fIntegers[fNevents + (fNevents + 1) * i] = ivars[i];

  fNevents++;
}

int CollieEventList::lookupIVarName(const char* name) const {
  for (int i = 0; i < fNints; i++)
    if (((TObjString*)(fValueNames[i + fNdoubles]))->String() == name)
      return i;
  return -1;
}
int CollieEventList::lookupDVarName(const char* name) const {
  for (int i = 0; i < fNdoubles; i++)
    if (((TObjString*)(fValueNames[i]))->String() == name)
      return i;
  return -1;
}
void CollieEventList::setDVarName(int nvar, const char* name) {
  if (!fMutable || nvar < 0 || nvar >= fNdoubles)
    return;
  ((TObjString*)(fValueNames[nvar]))->SetString((char*)name);
}
void CollieEventList::setIVarName(int nvar, const char* name) {
  if (!fMutable || nvar < 0 || nvar >= fNints)
    return;
  ((TObjString*)(fValueNames[fNdoubles + nvar]))->SetString((char*)name);
}

void CollieEventList::fillDistribution(CollieDistribution* cd, int overflow) const {
  double minX = cd->getMinX();
  int nX = cd->getNXbins();
  double deltaX = double(nX) / (cd->getMaxX() - minX);
  int nY = cd->getNYbins();
  double minY = cd->getMinY();

  if (nY < 2) {
    for (int n = 0; n < fNevents; n++) {
      int i = int((getX(n) - minX) * deltaX);
      if (i < 0 && overflow)
        i = 0;
      if (i >= nX && overflow)
        i = nX - 1;
      if (i >= 0 && i < nX)
        cd->setNormalizedBinValue(cd->getNormalizedBinValue(i, -1) + 1, i);
    }
  } else {
    double deltaY = double(nY) / (cd->getMaxY() - minY);
    for (int n = 0; n < fNevents; n++) {
      int i = int((getX(n) - minX) * deltaX);
      if (i < 0 && overflow)
        i = 0;
      if (i >= nX && overflow)
        i = nX - 1;
      int j = int((getY(n) - minY) * deltaY);
      if (j < 0 && overflow)
        j = 0;
      if (j >= nY && overflow)
        j = nY - 1;
      if (i >= 0 && j >= 0 && i < nX && j < nY)
        cd->setNormalizedBinValue(cd->getNormalizedBinValue(i, j) + 1, i, j);
    }
  }
  for (int i = 0; i < cd->getNModelXsec(); i++)
    cd->setModelXsec(CollieDistribution::rate_IgnoreLuminosity, i);
  cd->setEfficiency(1.0);
}

TH1D* CollieEventList::draw(const char* title, int nbins, double xlo, double xhi) {
  if (fDimensionality != 1)
    return NULL;

  TH1D* out = new TH1D(title, title, nbins, xlo, xhi);
  out->Sumw2();
  for (int i = 0; i < getNEvents(); i++)
    out->Fill(getX(i));

  return out;
}

TH2D* CollieEventList::draw2D(const char* title, int nbinsx, double xlo, double xhi, int nbinsy,
                              double ylo, double yhi) {
  if (fDimensionality != 2)
    return NULL;

  TH2D* out = new TH2D(title, title, nbinsx, xlo, xhi, nbinsy, ylo, yhi);
  out->Sumw2();
  for (int i = 0; i < getNEvents(); i++)
    out->Fill(getX(i), getY(i));

  return out;
}
