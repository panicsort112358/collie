#include "CollieIOFile.hh"

////////////////////////////
// Constructors/Destructor
////////////////////////////
CollieIOFile::CollieIOFile() {
  cutLowX_ = UNDEF;
  cutHighX_ = UNDEF;
  histMinX_ = UNDEF;
  histMaxX_ = UNDEF;
  histBinsX_ = 0;
  rebinX_ = 1;
  cutLowY_ = UNDEF;
  cutHighY_ = UNDEF;
  histMinY_ = UNDEF;
  histMaxY_ = UNDEF;
  histBinsY_ = 0;
  rebinY_ = 1;

  histNorm_ = HMAX;
  smIncr_ = 0;
  nBkgd_ = 0;
  nSig_ = 0;
  outName_ = "";
  chanName_ = "";
  channel_ = NULL;

  smooth_ = false;
  usingBinMap_ = false;
  init_ = false;
  verb_ = false;
  noviceFlag_ = true;
  systOvrd_ = false;

  data_.clear();
  sig_.clear();
  bkgd_.clear();
  allBkgd_.clear();
  sigAlpha_.clear();
  bkgdAlpha_.clear();
  systHists_.clear();
  systNames_.clear();
  parList_.clear();

  usingSystPriors_ = false;
  systPriors_ = NULL;

  outFile_ = NULL;
}

CollieIOFile::CollieIOFile(string outname, string channame) {

  if (outname.size() < 3) {
    printf("CollieIOFile::CollieIOFile, Error: Output file name is too short: %s\n",
           outname.c_str());
    return;
  }
  if (channame.size() < 3) {
    printf("CollieIOFile::CollieIOFile, Error: Channel name is too short: %s\n", channame.c_str());
    return;
  }
  outFile_ = NULL;

  usingSystPriors_ = false;
  systPriors_ = NULL;

  cutLowX_ = UNDEF;
  cutHighX_ = UNDEF;
  histMinX_ = UNDEF;
  histMaxX_ = UNDEF;
  histBinsX_ = 0;
  rebinX_ = 1;
  cutLowY_ = UNDEF;
  cutHighY_ = UNDEF;
  histMinY_ = UNDEF;
  histMaxY_ = UNDEF;
  histBinsY_ = 0;
  rebinY_ = 1;

  histNorm_ = HMAX;
  verb_ = false;
  smooth_ = false;
  noviceFlag_ = true;
  systOvrd_ = false;
  usingBinMap_ = false;

  smIncr_ = 0;
  nBkgd_ = 0;
  nSig_ = 0;

  data_.clear();
  sig_.clear();
  bkgd_.clear();
  sigAlpha_.clear();
  bkgdAlpha_.clear();
  allBkgd_.clear();
  parList_.clear();
  systHists_.clear();
  systNames_.clear();
  outName_ = outname;
  outFile_ = new TFile(outName_.c_str(), "RECREATE");

  usingSystPriors_ = false;
  systPriors_ = NULL;

  chanName_ = channame;
  channel_ = CollieChannel::createChannel(chanName_.c_str(), outFile_);
}

CollieIOFile::~CollieIOFile() {

  for (map<int, CollieHistogram*>::iterator it = allBkgd_.begin(); it != allBkgd_.end(); ++it)
    delete it->second;
  for (map<int, CollieHistogram*>::iterator it = data_.begin(); it != data_.end(); ++it)
    delete it->second;

  for (map<int, map<string, CollieHistogram*>>::iterator it = bkgd_.begin(); it != bkgd_.end();
       ++it) {
    for (map<string, CollieHistogram*>::iterator it2 = (it->second).begin();
         it2 != (it->second).end(); ++it2)
      delete it2->second;
  }

  for (map<int, map<string, CollieHistogram*>>::iterator it = sig_.begin(); it != sig_.end();
       ++it) {
    for (map<string, CollieHistogram*>::iterator it2 = (it->second).begin();
         it2 != (it->second).end(); ++it2)
      delete it2->second;
  }

  data_.clear();
  sig_.clear();
  bkgd_.clear();
  allBkgd_.clear();

  delete channel_;
  channel_ = NULL;

  if (outFile_) {
    outFile_->Close("R");
    outFile_->Flush();
    outFile_->Delete("T*;*");
  }
}

void CollieIOFile::initFile(string outname, string channame) {

  if (outname.size() < 3) {
    printf("CollieIOFile::initFile, Error: Output file name is too short: %s\n", outname.c_str());
    return;
  }
  if (channame.size() < 3) {
    printf("CollieIOFile::initFile, Error: Channel name is too short: %s\n", channame.c_str());
    return;
  }

  if (outName_ == outname && chanName_ == channame) {
    printf("CollieIOFile::initFile, Error: This file has already been initialized!\n");
    return;
  }
  if (outFile_) {
    printf("CollieIOFile::initFile, Error: This file has already been initialized!\n");
    return;
  }
  if (channel_) {
    printf("CollieIOFile::initFile, Error: This file has already been initialized!\n");
    return;
  }

  systHists_.clear();
  systNames_.clear();
  outName_ = outname;
  outFile_ = new TFile(outName_.c_str(), "RECREATE");

  usingSystPriors_ = false;
  systPriors_ = NULL;

  CollieHistogram* t = new CollieHistogram();
  t->ClearHistograms();

  chanName_ = channame;
  channel_ = CollieChannel::createChannel(chanName_.c_str(), outFile_);

  return;
}

void CollieIOFile::storeFile() {
  if (!check()) {
    printf("CollieIOFile::storeFile, Error: The file has not been properly initialized!\n");
    return;
  }

  //smooth if requested
  if (smooth_)
    generateSmoothedHistos();

  report();

  ///M.Owen 9/11/2007
  //  string outHists_ = "fv_"+outName_;
  string outHists_ = outName_;
  string::size_type pos = outHists_.find_last_of('/');
  if (pos != string::npos)
    outHists_.insert(pos + 1, "fv_");
  else
    outHists_ = "fv_" + outName_;

  printf("==>Saving inspection histos to %s (%zu systs)\n", outHists_.c_str(), systHists_.size());
  fflush(stdout);
  if (systHists_.size() > 0)
    CollieHistogram::StoreHistograms(outHists_.c_str(), "ROOT", systHists_);
  else
    CollieHistogram::StoreHistograms(outHists_.c_str());

  //store channel
  printf("==>Saving channel data to %s....\n", outName_.c_str());
  fflush(stdout);
  channel_->Store();
  outFile_->Close();

  return;
}

void CollieIOFile::generateSmoothedHistos() {
  if (!smooth_) {
    printf("CollieIOFile::generateSmoothedHistos, Error: Smoothing parameters not set!!\n");
    return;
  }
  if (!check()) {
    printf("CollieIOFile::generateSmoothedHistos, Error: File parameters not set!!\n");
    return;
  }

  if (data_.size() == 0)
    printf("CollieIOFile::generateSmoothedHistos, No data points registered!\n");

  //  list<int> massList;
  map<int, CollieHistogram*>::iterator iter, iter2;
  map<int, map<string, CollieHistogram*>>::iterator iter3;
  map<int, map<string, double>>::iterator iterA;
  map<string, double>::iterator iterA2;
  char title[256];

  CollieMasspoint* mp = 0;
  CollieDistribution* dist = 0;
  for (iter = data_.begin(); iter != data_.end(); iter++) {
    int massX = -1;
    int massY = -1;
    int massZ = -1;
    mLookUp(iter->first, massX, massY, massZ);

    iter3 = sig_.find(iter->first);
    if (iter3 != sig_.end()) {

      for (map<string, CollieHistogram*>::iterator inner = iter3->second.begin();
           inner != iter3->second.end(); inner++) {
        iterA = sigAlpha_.find(iter->first);
        iterA2 = iterA->second.find(inner->first);
        if (iterA2->second > 0 || iterA2->second < -1) {

          if (iterA2->second > 0)
            printf("==>Smoothing Signal %s for mass point %d/%d\n", inner->first.c_str(), massX,
                   massY);
          else
            printf("==>Grouping bins in Signal %s for mass point %d/%d\n", inner->first.c_str(),
                   massX, massY);

          CollieHistogram* smHist = new CollieHistogram();

          const char* style = "Smoothed";
          if (iterA2->second < -1)
            style = "Grouped Bin";
          if (massY > 0)
            sprintf(title, "%s Signal %s Final Variable - %d,%d", style, inner->first.c_str(),
                    massX, massY);
          else
            sprintf(title, "%s Signal %s Final Variable - %d", style, inner->first.c_str(), massX);
          smHist->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);
          //Create new smoothed signal
          if (iterA2->second > 0) {

            if (!getSmoothedCollie(iterA2->second, inner->second, smHist)) {
              printf("CollieIOFile::generateSmoothedHistos, signal %s smoothing failure!\n",
                     inner->first.c_str());
              return;
            }
          } else {
            if (!getGroupedBinCollie(iterA2->second, inner->second, smHist)) {
              printf("CollieIOFile::generateSmoothedHistos, signal %s grouping failure!\n",
                     inner->first.c_str());
              return;
            }
          }
          //Reset reference signal
          applyCuts(smHist);
          inner->second = smHist;
          mp = channel_->getMasspoint(massX, massY);
          dist = mp->getSignalDistMutable(getSigIndex(inner->first));
          if (!dist->fillFromHistogram(smHist, rebinX_)) {
            printf("CollieIOFile::generateSmoothedHistos, Error extracting signal dist %s!\n",
                   inner->first.c_str());
            return;
          }
        }
      }
    }

    iter3 = bkgd_.find(iter->first);
    if (iter3 != bkgd_.end()) {

      //  reset all bkgd histo
      iter2 = allBkgd_.find(iter->first);
      (*(iter2->second)) *= 0;

      for (map<string, CollieHistogram*>::iterator inner = iter3->second.begin();
           inner != iter3->second.end(); inner++) {
        iterA = bkgdAlpha_.find(iter->first);
        iterA2 = iterA->second.find(inner->first);
        if (iterA2->second > 0 || iterA2->second < -1) {
          if (iterA2->second > 0)
            printf("==>Smoothing %s for mass point %d/%d\n", inner->first.c_str(), massX, massY);
          else
            printf("==>Grouping bins in %s for mass point %d/%d\n", inner->first.c_str(), massX,
                   massY);
          CollieHistogram* smHist = new CollieHistogram();

          const char* style = (iterA2->second > 0) ? "Smoothed" : "Grouped Bin";
          if (massY > 0)
            sprintf(title, "%s %s Final Variable - %d,%d", style, inner->first.c_str(), massX,
                    massY);
          else
            sprintf(title, "%s %s Final Variable - %d", style, inner->first.c_str(), massX);
          smHist->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

          //Create new smoothed background
          if (iterA2->second > 0) {
            if (!getSmoothedCollie(iterA2->second, inner->second, smHist)) {
              printf("CollieIOFile::generateSmoothedHistos, background %s smoothing failure!\n",
                     inner->first.c_str());
              return;
            }
          } else {
            if (!getGroupedBinCollie(iterA2->second, inner->second, smHist)) {
              printf("CollieIOFile::generateSmoothedHistos, background %s grouping failure!\n",
                     inner->first.c_str());
              return;
            }
          }

          //Reset reference background
          applyCuts(smHist);
          inner->second = smHist;
          mp = channel_->getMasspoint(massX, massY);
          dist = mp->getBkgdDistMutable(getBkgdIndex(inner->first));
          if (!dist->fillFromHistogram(smHist, rebinX_)) {
            printf("CollieIOFile::generateSmoothedHistos, Error extracting bkgd dist %s!\n",
                   inner->first.c_str());
            return;
          }
        }
        iter2->second->Add(*(inner->second));
      }
    }
  }

  return;
}

void CollieIOFile::report() {

  printf("\n//*****************************************//\n");
  printf("CollieIOFile::report\n");

  if (!check()) {
    printf("CollieIOFile::report, Error: File parameters not set!!\n");
    return;
  }

  if (data_.size() == 0)
    printf("CollieIOFile::report, No data points registered!\n");

  //  channel_->print();

  printf("\n");
  map<int, CollieHistogram*>::iterator iter, iter2;
  map<int, map<string, CollieHistogram*>>::iterator iter3;
  for (iter = data_.begin(); iter != data_.end(); iter++) {
    int mx = -1;
    int my = -1;
    int mz = -1;
    mLookUp(iter->first, mx, my, mz);
    if (mz > 0 && my > 0)
      printf("Mass: %d, %d, %d\n", mx, my, mz);
    else if (my > 0)
      printf("Mass: %d, %d\n", mx, my);
    else
      printf("Mass: %d\n", mx);
    printf("     Data: %d\n", (int)(iter->second->Integral()));

    iter3 = sig_.find(iter->first);
    if (iter3 != sig_.end()) {
      for (map<string, CollieHistogram*>::iterator inner = iter3->second.begin();
           inner != iter3->second.end(); inner++)
        printf("   Signal: %s, %.3f\n", inner->first.c_str(), inner->second->Integral());
    }

    iter3 = bkgd_.find(iter->first);
    if (iter3 != bkgd_.end()) {
      for (map<string, CollieHistogram*>::iterator inner = iter3->second.begin();
           inner != iter3->second.end(); inner++)
        printf("     Bkgd: %s, %.3f\n", inner->first.c_str(), inner->second->Integral());

      iter2 = allBkgd_.find(iter->first);
      printf("     Allbkgd: %.2f\n", iter2->second->Integral());
    }
  }
  printf("\n//*****************************************//\n");

  return;
}

void CollieIOFile::createChannel(vector<string> bkgdNames, vector<string> signalNames) {
  vector<string> varNames;
  varNames.push_back("Final Variable");
  createChannel(bkgdNames, signalNames, varNames);
}

void CollieIOFile::createChannel(vector<string> bkgdNames) {
  vector<string> varNames;
  vector<string> sigNames;
  varNames.push_back("Final Variable");
  sigNames.push_back("Signal");

  createChannel(bkgdNames, sigNames, varNames);
}

void CollieIOFile::createChannel(vector<string> bkgdNames, vector<string> sigNames,
                                 vector<string> varNames) {

  if (!channel_) {
    printf("CollieIOFile::createChannel, Error: The file is not initialized!!\n");
    return;
  }
  nBkgd_ = bkgdNames.size();
  nSig_ = sigNames.size();

  if (varNames.size() > 1) {
    channel_->setNIndepVariables(varNames.size());
    for (uint n = 0; n < varNames.size(); n++)
      channel_->setIndepVarName(n, varNames[n].c_str());
  } else {
    channel_->setNIndepVariables(1);
    channel_->setIndepVarName(0, "Final Variable");
  }

  channel_->setNSignals(nSig_);
  channel_->setNBackgrounds(nBkgd_);
  for (uint i = 0; i < nSig_; i++)
    channel_->setSignalName(i, sigNames[i].c_str());
  for (uint i = 0; i < nBkgd_; i++)
    channel_->setBackgroundName(i, bkgdNames[i].c_str());

  channel_->setNmodels(1);
  channel_->setModelName(0, "SM");
  channel_->setLuminosity(1.0);

  init_ = true;

  return;
}

bool CollieIOFile::logMassPoint(int pX, int pY, int pZ) {
  map<int, int>::iterator iter = parList_.find(mIndex(pX, pY, pZ));
  if (iter == parList_.end()) {
    parList_[mIndex(pX, pY, pZ)] = 1;
  } else {
    printf("collieIOFile::logMassPoint, Warning: Already have a mass point for %d,%d,%d\n", pX, pY,
           pZ);
    return false;
  }

  return true;
}

const CollieEventList* CollieIOFile::getDataEventList(int massX, int massY) const {

  CollieMasspoint* point;
  if (massY > 0)
    point = channel_->getMasspoint(massX, massY);
  else
    point = channel_->getMasspoint(massX);

  if (point)
    return point->getDataEventList();

  printf("CollieIOFile::getDataEventList, This mass point doesn't exist: %d,%d\n", massX, massY);
  return NULL;
}

CollieEventList* CollieIOFile::getDataEventListMutable(int massX, int massY) {
  CollieMasspoint* point;
  if (massY > 0)
    point = channel_->getMasspoint(massX, massY);
  else
    point = channel_->getMasspoint(massX);

  if (point)
    return point->getDataEventListMutable();

  printf("CollieIOFile::getDataEventListMutable, This mass point doesn't exist: %d,%d\n", massX,
         massY);
  return NULL;
}

/// if possible, store data event list
void CollieIOFile::finalizeDataEventList(int massX, int massY) {

  channel_->doCD();

  CollieMasspoint* point;
  if (massY > 0)
    point = channel_->getMasspoint(massX, massY);
  else
    point = channel_->getMasspoint(massX);

  if (point) {
    point->getDataDist();
    CollieEventList* evtList = point->getDataEventListMutable();

    if (histBinsY_ <= 0) {
      CollieHistogram* mdata = new CollieHistogram();
      char title[256];
      if (massY > 0)
        sprintf(title, "Data Final Variable - %d,%d", massX, massY);
      else
        sprintf(title, "Data Final Variable - %d", massX);
      mdata->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

      TH1D* data = evtList->draw(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);
      if (!getCollie(data, mdata)) {
        printf("CollieIOFile::finalizeDataEventList, data histo error!\n");
        return;
      }

      CollieDistribution* dist = point->getDataDistMutable();
      if (!dist->fillFromHistogram(data, rebinX_)) {
        printf("CollieIOFile::finalizeDataEventList, Error extracting data dist!\n");
        return;
      }

      //log data dist
      data_[mIndex(massX, massY)] = mdata;
    } else {
      CollieHistogram2d* mdata = new CollieHistogram2d();
      char title[256];
      if (massY > 0)
        sprintf(title, "Data Final Variable - %d,%d", massX, massY);
      else
        sprintf(title, "Data Final Variable - %d", massX);
      mdata->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), (int)(histBinsY_ * 1.0 / rebinY_),
                  histMinX_, histMaxX_, histMinY_, histMaxY_);

      TH2D* data = evtList->draw2D(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_,
                                   (int)(histBinsY_ * 1.0 / rebinY_), histMinY_, histMaxY_);
      if (!getCollie2D(data, mdata)) {
        printf("CollieIOFile::finalizeDataEventList, data histo error!\n");
        return;
      }

      CollieDistribution* dist = point->getDataDistMutable();
      if (!dist->fillFromHistogram(data, rebinX_, rebinY_)) {
        printf("CollieIOFile::finalizeDataEventList, Error extracting data dist!\n");
        return;
      }

      //log data dist
      data_[mIndex(massX, massY)] = mdata;
    }

    return;
  }

  printf("CollieIOFile::finalizeDataEventList, This mass point doesn't exist: %d,%d\n", massX,
         massY);
  return;
}

int CollieIOFile::getBkgdIndex(const string name) {
  for (uint i = 0; i < nBkgd_; i++) {
    if (string(channel_->getBackgroundName(i)) == name)
      return i;
  }
  return -1;
}

int CollieIOFile::getSigIndex(const string name) {
  for (uint i = 0; i < nSig_; i++) {
    if (string(channel_->getSignalName(i)) == name)
      return i;
  }
  return -1;
}

void CollieIOFile::setDataDist(TH1D* data, int parX, int parY, int parZ) {

  ///fill data dist
  char title[250];
  CollieHistogram* mdata = new CollieHistogram();
  if (parZ > 0)
    sprintf(title, "Data Final Variable - %d,%d,%d", parX, parY, parZ);
  else if (parY > 0)
    sprintf(title, "Data Final Variable - %d,%d", parX, parY);
  else
    sprintf(title, "Data Final Variable - %d", parX);
  mdata->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

  getCollie(data, mdata);
  data_[mIndex(parX, parY, parZ)] = mdata;
  CollieMasspoint* point;
  if (parZ > 0)
    point = channel_->getMasspoint(parX, parY, parZ);
  else if (parY > 0)
    point = channel_->getMasspoint(parX, parY);
  else
    point = channel_->getMasspoint(parX);
  CollieDistribution* dist = point->getDataDistMutable();
  //  fillDist(mdata,dist);//WF 1/17/08
  if (!dist->fillFromHistogram(data, rebinX_)) {
    printf("CollieIOFile::setDataDist, Error extracting data dist!\n");
    return;
  }
  return;
}

void CollieIOFile::setDataDist2D(TH2D* data, int parX, int parY, int parZ) {
  ///fill data data
  char title[250];
  CollieHistogram2d* mdata = new CollieHistogram2d();
  if (parZ > 0)
    sprintf(title, "Data Final Variable - %d,%d,%d", parX, parY, parZ);
  else if (parY > 0)
    sprintf(title, "Data Final Variable - %d,%d", parX, parY);
  else
    sprintf(title, "Data Final Variable - %d", parX);
  mdata->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), (int)(histBinsY_ * 1.0 / rebinY_),
              histMinX_, histMaxX_, histMinY_, histMaxY_);

  getCollie2D(data, mdata);
  data_[mIndex(parX, parY, parZ)] = mdata;
  CollieMasspoint* point;
  if (parZ > 0)
    point = channel_->getMasspoint(parX, parY, parZ);
  else if (parY > 0)
    point = channel_->getMasspoint(parX, parY);
  else
    point = channel_->getMasspoint(parX);
  CollieDistribution* dist = point->getDataDistMutable();
  //  fillDist(mdata,dist);//WF 1/17/08
  if (!dist->fillFromHistogram(data, rebinX_, rebinY_)) {
    printf("CollieIOFile::setDataDist2D, Error extracting data dist!\n");
    return;
  }
  return;
}

//////////////////////////////
//MassPoint Creation Methods
//////////////////////////////
void CollieIOFile::createMassPoint(int par, TH1D* data, TH1D* sig, double alphaS,
                                   vector<TH1D*>& bkg, vector<double>& alphaB) {
  vector<TH1D*> sigV;
  sigV.push_back(sig);
  vector<double> sigA;
  sigA.push_back(alphaS);
  return createMassPoint(par, -1, -1, data, sigV, sigA, bkg, alphaB);
}

void CollieIOFile::createMassPoint(int parX, int parY, TH1D* data, TH1D* sig, double alphaS,
                                   vector<TH1D*>& bkg, vector<double>& alphaB) {
  vector<TH1D*> sigV;
  sigV.push_back(sig);
  vector<double> sigA;
  sigA.push_back(alphaS);
  return createMassPoint(parX, parY, -1, data, sigV, sigA, bkg, alphaB);
}

void CollieIOFile::createMassPoint(int parX, int parY, int parZ, TH1D* data, TH1D* sig,
                                   double alphaS, vector<TH1D*>& bkg, vector<double>& alphaB) {
  vector<TH1D*> sigV;
  sigV.push_back(sig);
  vector<double> sigA;
  sigA.push_back(alphaS);
  return createMassPoint(parX, parY, parZ, data, sigV, sigA, bkg, alphaB);
}

void CollieIOFile::createMassPoint(int par, TH1D* data, vector<TH1D*>& sig, vector<double>& alphaS,
                                   vector<TH1D*>& bkg, vector<double>& alphaB) {
  return createMassPoint(par, -1, -1, data, sig, alphaS, bkg, alphaB);
}

void CollieIOFile::createMassPoint(int parX, int parY, TH1D* data, vector<TH1D*>& sig,
                                   vector<double>& alphaS, vector<TH1D*>& bkg,
                                   vector<double>& alphaB) {
  return createMassPoint(parX, parY, -1, data, sig, alphaS, bkg, alphaB);
}

void CollieIOFile::createMassPoint2D_EvtList(int mass, int nVarD, int nVarI, TH2D* sig,
                                             double alphaS, vector<TH2D*>& bkgd,
                                             vector<double>& alphaB) {
  vector<TH2D*> sigV;
  sigV.push_back(sig);
  vector<double> sigA;
  sigA.push_back(alphaS);
  createMassPoint2D_EvtList(mass, -1, -1, nVarD, nVarI, sigV, sigA, bkgd, alphaB);
}

void CollieIOFile::createMassPoint2D_EvtList(int massX, int massY, int nVarD, int nVarI, TH2D* sig,
                                             double alphaS, vector<TH2D*>& bkgd,
                                             vector<double>& alphaB) {
  vector<TH2D*> sigV;
  sigV.push_back(sig);
  vector<double> sigA;
  sigA.push_back(alphaS);
  createMassPoint2D_EvtList(massX, massY, -1, nVarD, nVarI, sigV, sigA, bkgd, alphaB);
}

void CollieIOFile::createMassPoint2D_EvtList(int massX, int massY, int massZ, int nVarD, int nVarI,
                                             TH2D* sig, double alphaS, vector<TH2D*>& bkgd,
                                             vector<double>& alphaB) {
  vector<TH2D*> sigV;
  sigV.push_back(sig);
  vector<double> sigA;
  sigA.push_back(alphaS);
  createMassPoint2D_EvtList(massX, massY, massZ, nVarD, nVarI, sigV, sigA, bkgd, alphaB);
}

void CollieIOFile::createMassPoint2D_EvtList(int mass, int nVarD, int nVarI, vector<TH2D*>& sig,
                                             vector<double>& alphaS, vector<TH2D*>& bkgd,
                                             vector<double>& alphaB) {
  printf("CollieIOFile::createMassPoint2D_EvtList Needs Implementation!\n");
}

void CollieIOFile::createMassPoint2D_EvtList(int massX, int massY, int nVarD, int nVarI,
                                             vector<TH2D*>& sig, vector<double>& alphaS,
                                             vector<TH2D*>& bkgd, vector<double>& alphaB) {
  printf("CollieIOFile::createMassPoint2D_EvtList Needs Implementation!\n");
}

void CollieIOFile::createMassPoint2D_EvtList(int massX, int massY, int massZ, int nVarD, int nVarI,
                                             vector<TH2D*>& sig, vector<double>& alphaS,
                                             vector<TH2D*>& bkgd, vector<double>& alphaB) {
  printf("CollieIOFile::createMassPoint2D_EvtList Needs Implementation!\n");
}

void CollieIOFile::createMassPoint_EvtList(int mass, int nVarD, int nVarI, TH1D* sig, double alphaS,
                                           vector<TH1D*>& bkgd, vector<double>& alphaB) {
  vector<TH1D*> sigV;
  sigV.push_back(sig);
  vector<double> sigA;
  sigA.push_back(alphaS);
  createMassPoint_EvtList(mass, -1, -1, nVarD, nVarI, sigV, sigA, bkgd, alphaB);
}

void CollieIOFile::createMassPoint_EvtList(int massX, int massY, int nVarD, int nVarI, TH1D* sig,
                                           double alphaS, vector<TH1D*>& bkgd,
                                           vector<double>& alphaB) {
  vector<TH1D*> sigV;
  sigV.push_back(sig);
  vector<double> sigA;
  sigA.push_back(alphaS);
  createMassPoint_EvtList(massX, massY, -1, nVarD, nVarI, sigV, sigA, bkgd, alphaB);
}

void CollieIOFile::createMassPoint_EvtList(int massX, int massY, int massZ, int nVarD, int nVarI,
                                           TH1D* sig, double alphaS, vector<TH1D*>& bkgd,
                                           vector<double>& alphaB) {
  vector<TH1D*> sigV;
  sigV.push_back(sig);
  vector<double> sigA;
  sigA.push_back(alphaS);
  createMassPoint_EvtList(massX, massY, massZ, nVarD, nVarI, sigV, sigA, bkgd, alphaB);
}

void CollieIOFile::createMassPoint_EvtList(int mass, int nVarD, int nVarI, vector<TH1D*>& sigV,
                                           vector<double>& sigA, vector<TH1D*>& bkgd,
                                           vector<double>& alphaB) {
  createMassPoint_EvtList(mass, -1, -1, nVarD, nVarI, sigV, sigA, bkgd, alphaB);
}

void CollieIOFile::createMassPoint_EvtList(int massX, int massY, int nVarD, int nVarI,
                                           vector<TH1D*>& sigV, vector<double>& sigA,
                                           vector<TH1D*>& bkgd, vector<double>& alphaB) {
  createMassPoint_EvtList(massX, massY, -1, nVarD, nVarI, sigV, sigA, bkgd, alphaB);
}

void CollieIOFile::createMassPoint_EvtList(int parX, int parY, int parZ, int nVarD, int nVarI,
                                           vector<TH1D*>& sig, vector<double>& alphaS,
                                           vector<TH1D*>& bkg, vector<double>& alphaB) {
  if (!check()) {
    printf("CollieIOFile::createMassPoint_EvtList, Error: The file is not properly initialized!\n");
    return;
  }
  if (checkMassPoint(parX, parY, parZ)) {
    printf(
        "CollieIOFile::createMassPoint_EvtList, Error: This mass point has already been created: %d, %d, %d\n",
        parX, parY, parZ);
    return;
  }
  if (bkg.size() != nBkgd_) {
    printf("CollieIOFile::createMassPoint_EvtList, Error: Mismatched bkgd sizes: %zu, %u\n",
           bkg.size(), nBkgd_);
    return;
  }
  if (nVarD < 0 || nVarI < 0) {
    printf("CollieIOFile::createMassPoint_EvtList, Error: Negative variable sizes: %d, %d\n", nVarD,
           nVarI);
    return;
  }

  char title[250];
  CollieMasspoint* point = channel_->createMasspoint(parX, parY, parZ);
  point->Book((int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);
  point->ConfigEventList(nVarD, nVarI);

  ///fill signal dists
  map<string, CollieHistogram*> inputSig;
  map<string, double> inputAlphaS;
  for (uint i = 0; i < sig.size(); i++) {
    if (sig[i] == NULL)
      continue;
    const char* tmpName = channel_->getSignalName(i);

    CollieDistribution* dist = point->getSignalDistMutable(i);
    CollieHistogram* sigHist = new CollieHistogram();
    if (parZ > 0 && parY > 0)
      sprintf(title, "Signal %s Final Variable - %d,%d,%d", tmpName, parX, parY, parZ);
    else if (parY > 0)
      sprintf(title, "Signal %s Final Variable - %d,%d", tmpName, parX, parY);
    else
      sprintf(title, "Signal %s Final Variable - %d", tmpName, parX);
    sigHist->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

    if (!getCollie(sig[i], sigHist)) {
      printf("CollieIOFile::createMassPoint, signal histo error!\n");
      return;
    }
    inputSig[tmpName] = sigHist;

    if (alphaS.size() == sig.size())
      inputAlphaS[tmpName] = alphaS[i];
    else
      inputAlphaS[tmpName] = -1;

    if (!dist->fillFromHistogram(sig[i], rebinX_)) {
      printf("CollieIOFile::createMassPoint, Error extracting signal dist!\n");
      return;
    }
    dist->setModelXsec(1.0);
  }
  //  nSig_ = inputSig.size();
  sigAlpha_[mIndex(parX, parY, parZ)] = inputAlphaS;
  sig_[mIndex(parX, parY, parZ)] = inputSig;

  ///fill bkgd dists
  CollieHistogram* abkgd = new CollieHistogram();
  if (parZ > 0 && parY > 0)
    sprintf(title, "All Bkgd Final Variable - %d,%d,%d", parX, parY, parZ);
  else if (parY > 0)
    sprintf(title, "All Bkgd Final Variable - %d,%d", parX, parY);
  else
    sprintf(title, "All Bkgd Final Variable - %d", parX);
  abkgd->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

  map<string, CollieHistogram*> inputB;
  map<string, double> inputAlphaB;
  for (uint i = 0; i < bkg.size(); i++) {
    if (bkg[i] == NULL)
      continue;
    const char* tmpName = channel_->getBackgroundName(i);
    CollieDistribution* dist = point->getBkgdDistMutable(i);

    CollieHistogram* hist = new CollieHistogram();
    if (parZ > 0 && parY > 0)
      sprintf(title, "%s Final Variable - %d,%d,%d", tmpName, parX, parY, parZ);
    else if (parY > 0)
      sprintf(title, "%s Final Variable - %d,%d", tmpName, parX, parY);
    else
      sprintf(title, "%s Final Variable - %d", tmpName, parX);
    hist->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

    if (!getCollie(bkg[i], hist)) {
      printf("CollieIOFile::createMassPoint, bkgd histo error!\n");
      return;
    }
    inputB[tmpName] = hist;

    if (alphaB.size() == bkg.size())
      inputAlphaB[tmpName] = alphaB[i];
    else
      inputAlphaB[tmpName] = -1;

    if (!dist->fillFromHistogram(bkg[i], rebinX_)) {
      printf("CollieIOFile::createMassPoint, Error extracting bkgd dist %d!\n", i);
      return;
    }
    dist->setModelXsec(1.0);
    abkgd->Add(*hist);
  }
  bkgd_[mIndex(parX, parY, parZ)] = inputB;
  bkgdAlpha_[mIndex(parX, parY, parZ)] = inputAlphaB;
  allBkgd_[mIndex(parX, parY, parZ)] = abkgd;
  nBkgd_ = inputB.size();

  //log mass point
  if (!logMassPoint(parX, parY, parZ)) {
    printf("CollieIOFile::createMassPoint, Error logging mass point!\n");
    return;
  }
  if (parZ > 0 && parY > 0)
    printf("==>Created mass point %d, %d, %d\n", parX, parY, parZ);
  else if (parY > 0)
    printf("==>Created mass point %d, %d\n", parX, parY);
  else
    printf("==>Created mass point %d\n", parX);

  return;
}

void CollieIOFile::createMassPoint(int parX, int parY, int parZ, TH1D* data, vector<TH1D*>& sig,
                                   vector<double>& alphaS, vector<TH1D*>& bkg,
                                   vector<double>& alphaB) {

  if (!check()) {
    printf("CollieIOFile::createMassPoint, Error: The file is not properly initialized!\n");
    return;
  }
  if (checkMassPoint(parX, parY, parZ)) {
    printf(
        "CollieIOFile::createMassPoint, Error: This mass point has already been created: %d, %d, %d\n",
        parX, parY, parZ);
    return;
  }
  if (bkg.size() != nBkgd_) {
    printf("CollieIOFile::createMassPoint, Error: Mismatched bkgd sizes: %zu, %u\n", bkg.size(),
           nBkgd_);
    return;
  }
  if (sig.size() != nSig_) {
    printf("CollieIOFile::createMassPoint, Error: Mismatched signal sizes: %zu, %u\n", sig.size(),
           nSig_);
    return;
  }
  checkBins(data, sig, bkg);

  //  printf("Internal rate1: %f, %f, %f\n",bkg[0]->Integral(), data->Integral(),sig[0]->Integral());
  //  printf("information0: %d, %d, %d, %d,\n",parX,parY,parZ,mIndex(parX,parY,parZ));

  char title[250];
  CollieMasspoint* point = channel_->createMasspoint(parX, parY, parZ);
  point->Book((int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);
  point->BookData();
  //  printf("Internal rate2: %f, %f, %f\n",bkg[0]->Integral(), data->Integral(),sig[0]->Integral());

  ///fill signal dists
  map<string, CollieHistogram*> inputSig;
  map<string, double> inputAlphaS;
  for (uint i = 0; i < sig.size(); i++) {
    if (sig[i] == NULL)
      continue;
    const char* tmpName = channel_->getSignalName(i);
    applyCuts(sig[i]);

    CollieDistribution* dist = point->getSignalDistMutable(i);
    CollieHistogram* sigHist = new CollieHistogram();
    if (parZ > 0 && parY > 0)
      sprintf(title, "Signal %s Final Variable - %d,%d,%d", tmpName, parX, parY, parZ);
    else if (parY > 0)
      sprintf(title, "Signal %s Final Variable - %d,%d", tmpName, parX, parY);
    else
      sprintf(title, "Signal %s Final Variable - %d", tmpName, parX);
    sigHist->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

    if (!getCollie(sig[i], sigHist)) {
      printf("CollieIOFile::createMassPoint, signal histo error!\n");
      return;
    }
    inputSig[tmpName] = sigHist;

    if (alphaS.size() == sig.size())
      inputAlphaS[tmpName] = alphaS[i];
    else
      inputAlphaS[tmpName] = -1;

    if (!dist->fillFromHistogram(sigHist, rebinX_)) {
      printf("CollieIOFile::createMassPoint, Error extracting signal dist!\n");
      return;
    }
    dist->setModelXsec(1.0);
  }
  //  printf("information: %d, %d, %d, %d, %d, %d, %d, %d\n",inputSig.size(),sigAlpha_.size(),parX,parY,parZ,mIndex(parX,parY,parZ),inputAlphaS.size(),inputSig.size());
  //  nSig_ = inputSig.size();
  sigAlpha_[mIndex(parX, parY, parZ)] = inputAlphaS;
  sig_[mIndex(parX, parY, parZ)] = inputSig;

  //  printf("Internal rate3: %f, %f, %f\n",bkg[0]->Integral(), data->Integral(),sig[0]->Integral());

  ///fill bkgd dists
  CollieHistogram* abkgd = new CollieHistogram();
  if (parZ > 0 && parY > 0)
    sprintf(title, "All Bkgd Final Variable - %d,%d,%d", parX, parY, parZ);
  else if (parY > 0)
    sprintf(title, "All Bkgd Final Variable - %d,%d", parX, parY);
  else
    sprintf(title, "All Bkgd Final Variable - %d", parX);
  abkgd->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);
  //  printf("Internal rate4: %f, %f, %f\n",bkg[0]->Integral(), data->Integral(),sig[0]->Integral());

  map<string, CollieHistogram*> inputB;
  map<string, double> inputAlphaB;
  for (uint i = 0; i < bkg.size(); i++) {
    if (bkg[i] == NULL)
      continue;
    applyCuts(bkg[i]);
    const char* tmpName = channel_->getBackgroundName(i);
    CollieDistribution* dist = point->getBkgdDistMutable(i);

    CollieHistogram* hist = new CollieHistogram();
    if (parZ > 0 && parY > 0)
      sprintf(title, "%s Final Variable - %d,%d,%d", tmpName, parX, parY, parZ);
    else if (parY > 0)
      sprintf(title, "%s Final Variable - %d,%d", tmpName, parX, parY);
    else
      sprintf(title, "%s Final Variable - %d", tmpName, parX);
    hist->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

    if (!getCollie(bkg[i], hist)) {
      printf("CollieIOFile::createMassPoint, bkgd histo error!\n");
      return;
    }
    inputB[tmpName] = hist;

    if (alphaB.size() == bkg.size())
      inputAlphaB[tmpName] = alphaB[i];
    else
      inputAlphaB[tmpName] = -1;

    if (!dist->fillFromHistogram(hist, rebinX_)) {
      printf("CollieIOFile::createMassPoint, Error extracting bkgd dist %d!\n", i);
      return;
    }
    dist->setModelXsec(1.0);
    abkgd->Add(*hist);
  }

  bkgd_[mIndex(parX, parY, parZ)] = inputB;
  bkgdAlpha_[mIndex(parX, parY, parZ)] = inputAlphaB;
  allBkgd_[mIndex(parX, parY, parZ)] = abkgd;
  nBkgd_ = inputB.size();
  //  printf("Internal rate5: %f, %f, %f\n",bkg[0]->Integral(), data->Integral(),sig[0]->Integral());

  ///fill data dist
  applyCuts(data);
  CollieHistogram* mdata = new CollieHistogram();
  if (parZ > 0 && parY > 0)
    sprintf(title, "Data Final Variable - %d,%d,%d", parX, parY, parZ);
  else if (parY > 0)
    sprintf(title, "Data Final Variable - %d,%d", parX, parY);
  else
    sprintf(title, "Data Final Variable - %d", parX);
  mdata->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

  //  printf("Getting data: %f\n",data->Integral());

  if (!getCollie(data, mdata)) {
    printf("CollieIOFile::createMassPoint, data histo error!\n");
    return;
  }
  data_[mIndex(parX, parY, parZ)] = mdata;
  CollieDistribution* dist = point->getDataDistMutable();
  //  fillDist(mdata,dist);//WF 1/17/08
  if (!dist->fillFromHistogram(mdata, rebinX_)) {
    printf("CollieIOFile::createMassPoint, Error extracting data dist!\n");
    return;
  }

  //  printf("About to log: %f, %f\n",mdata->Integral(),abkgd->Integral());

  //log mass point
  point->logPoint();
  if (!logMassPoint(parX, parY, parZ)) {
    printf("CollieIOFile::createMassPoint, Error logging mass point!\n");
    return;
  }
  if (parZ > 0 && parY > 0)
    printf("==>Created mass point %d, %d, %d\n", parX, parY, parZ);
  else if (parY > 0)
    printf("==>Created mass point %d, %d\n", parX, parY);
  else
    printf("==>Created mass point %d\n", parX);

  return;
}

void CollieIOFile::createMassPoint2D(int mass, TH2D* data, TH2D* sig, double alphaS,
                                     vector<TH2D*>& bkg, vector<double>& alphaB) {
  vector<TH2D*> sigV;
  sigV.push_back(sig);
  vector<double> sigA;
  sigA.push_back(alphaS);
  return createMassPoint2D(mass, -1, -1, data, sigV, sigA, bkg, alphaB);
}

void CollieIOFile::createMassPoint2D(int massX, int massY, TH2D* data, TH2D* sig, double alphaS,
                                     vector<TH2D*>& bkg, vector<double>& alphaB) {
  vector<TH2D*> sigV;
  sigV.push_back(sig);
  vector<double> sigA;
  sigA.push_back(alphaS);
  return createMassPoint2D(massX, massY, -1, data, sigV, sigA, bkg, alphaB);
}

void CollieIOFile::createMassPoint2D(int massX, int massY, int massZ, TH2D* data, TH2D* sig,
                                     double alphaS, vector<TH2D*>& bkg, vector<double>& alphaB) {
  vector<TH2D*> sigV;
  sigV.push_back(sig);
  vector<double> sigA;
  sigA.push_back(alphaS);
  return createMassPoint2D(massX, massY, massZ, data, sigV, sigA, bkg, alphaB);
}

void CollieIOFile::createMassPoint2D(int mass, TH2D* data, vector<TH2D*>& sig,
                                     vector<double>& alphaS, vector<TH2D*>& bkg,
                                     vector<double>& alphaB) {
  return createMassPoint2D(mass, -1, -1, data, sig, alphaS, bkg, alphaB);
}

void CollieIOFile::createMassPoint2D(int massX, int massY, TH2D* data, vector<TH2D*>& sig,
                                     vector<double>& alphaS, vector<TH2D*>& bkg,
                                     vector<double>& alphaB) {
  return createMassPoint2D(massX, massY, -1, data, sig, alphaS, bkg, alphaB);
}

void CollieIOFile::createMassPoint2D(int massX, int massY, int massZ, TH2D* data,
                                     vector<TH2D*>& sig, vector<double>& alphaS, vector<TH2D*>& bkg,
                                     vector<double>& alphaB) {
  if (!check(true)) {
    printf("CollieIOFile::createMassPoint2D, Error: The file is not properly initialized!\n");
    return;
  }
  if (checkMassPoint(massX, massY, massZ)) {
    printf(
        "CollieIOFile::createMassPoint2D, Error: This mass point has already been created: %d, %d\n",
        massX, massY);
    return;
  }
  if (bkg.size() != nBkgd_) {
    printf("CollieIOFile::createMassPoint2D, Error: Mismatched bkgd sizes: %zu, %u\n", bkg.size(),
           nBkgd_);
    return;
  }

  checkBins2D(data, sig, bkg);

  char title[250];
  CollieMasspoint* point = channel_->createMasspoint(massX, massY, massZ);
  point->Book((int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_,
              (int)(histBinsY_ * 1.0 / rebinY_), histMinY_, histMaxY_);
  point->BookData();

  ///fill signal dists
  map<string, CollieHistogram*> inputSig;
  map<string, double> inputAlphaS;
  for (uint i = 0; i < sig.size(); i++) {
    if (sig[i] == NULL)
      continue;
    const char* tmpName = channel_->getSignalName(i);
    applyCuts(sig[i]);

    CollieDistribution* dist = point->getSignalDistMutable(i);
    CollieHistogram2d* sigHist = new CollieHistogram2d();
    if (massZ > 0 && massY > 0)
      sprintf(title, "Signal %s Final Variable - %d,%d,%d", tmpName, massX, massY, massZ);
    else if (massY > 0)
      sprintf(title, "Signal %s Final Variable - %d,%d", tmpName, massX, massY);
    else
      sprintf(title, "Signal %s Final Variable - %d", tmpName, massX);
    sigHist->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), (int)(histBinsY_ * 1.0 / rebinY_),
                  histMinX_, histMaxX_, histMinY_, histMaxY_);

    if (!getCollie2D(sig[i], sigHist)) {
      printf("CollieIOFile::createMassPoint2D, signal histo error!\n");
      return;
    }
    inputSig[tmpName] = sigHist;

    if (alphaS.size() == sig.size())
      inputAlphaS[tmpName] = alphaS[i];
    else
      inputAlphaS[tmpName] = -1;

    if (!dist->fillFromHistogram(sig[i], rebinX_)) {
      printf("CollieIOFile::createMassPoint, Error extracting signal dist!\n");
      return;
    }
    dist->setModelXsec(1.0);
  }
  //  nSig_ = inputSig.size();
  sigAlpha_[mIndex(massX, massY, massZ)] = inputAlphaS;
  sig_[mIndex(massX, massY, massZ)] = inputSig;

  ///fill bkgd dists
  CollieHistogram2d* abkgd = new CollieHistogram2d();
  if (massZ > 0)
    sprintf(title, "All Bkgd Final Variable - %d,%d,%d", massX, massY, massZ);
  else if (massY > 0)
    sprintf(title, "All Bkgd Final Variable - %d,%d", massX, massY);
  else
    sprintf(title, "All Bkgd Final Variable - %d", massX);
  abkgd->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), (int)(histBinsY_ * 1.0 / rebinY_),
              histMinX_, histMaxX_, histMinY_, histMaxY_);

  map<string, CollieHistogram*> input;
  map<string, double> inputAlpha;
  for (unsigned int i = 0; i < bkg.size(); i++) {
    if (bkg[i] == NULL)
      continue;
    const char* tmpName = channel_->getBackgroundName(i);
    applyCuts(bkg[i]);

    CollieHistogram2d* hist = new CollieHistogram2d();
    if (massZ > 0)
      sprintf(title, "%s Final Variable - %d,%d,%d", tmpName, massX, massY, massZ);
    else if (massY > 0)
      sprintf(title, "%s Final Variable - %d,%d", tmpName, massX, massY);
    else
      sprintf(title, "%s Final Variable - %d", tmpName, massX);
    hist->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), (int)(histBinsY_ * 1.0 / rebinY_),
               histMinX_, histMaxX_, histMinY_, histMaxY_);

    if (!getCollie2D(bkg[i], hist)) {
      printf("CollieIOFile::createMassPoint2D, signal histo error!\n");
      return;
    }
    input[tmpName] = hist;

    if (alphaB.size() == bkg.size())
      inputAlpha[tmpName] = alphaB[i];
    else
      inputAlpha[tmpName] = -1;

    CollieDistribution* dist = point->getBkgdDistMutable(i);

    if (!dist->fillFromHistogram(bkg[i], rebinX_, rebinY_)) {
      printf("CollieIOFile::createMassPoint, Error extracting bkgd dist %d!\n", i);
      return;
    }
    dist->setModelXsec(1.0);
    abkgd->Add(*hist);
  }
  bkgd_[mIndex(massX, massY, massZ)] = input;
  bkgdAlpha_[mIndex(massX, massY, massZ)] = inputAlpha;
  allBkgd_[mIndex(massX, massY, massZ)] = abkgd;
  nBkgd_ = input.size();

  ///fill data data
  applyCuts(data);
  CollieHistogram2d* mdata = new CollieHistogram2d();
  if (massZ > 0)
    sprintf(title, "Data Final Variable - %d,%d,%d", massX, massY, massZ);
  else if (massY > 0)
    sprintf(title, "Data Final Variable - %d,%d", massX, massY);
  else
    sprintf(title, "Data Final Variable - %d", massX);
  mdata->Book(title, (int)(histBinsX_ * 1.0 / rebinX_), (int)(histBinsY_ * 1.0 / rebinY_),
              histMinX_, histMaxX_, histMinY_, histMaxY_);

  getCollie2D(data, mdata);
  data_[mIndex(massX, massY, massZ)] = mdata;
  CollieDistribution* dist = point->getDataDistMutable();
  //  fillDist(mdata,dist);//WF 1/17/08
  if (!dist->fillFromHistogram(data, rebinX_, rebinY_)) {
    printf("CollieIOFile::createMassPoint, Error extracting data dist!\n");
    return;
  }

  ///log mass point
  point->logPoint();
  if (!logMassPoint(massX, massY, massZ)) {
    printf("CollieIOFile::createMassPoint2D, Error logging mass point!\n");
    return;
  }
  if (massZ > 0)
    printf("==>Created mass point %d, %d, %d\n", massX, massY, massZ);
  else if (massY > 0)
    printf("==>Created mass point %d, %d\n", massX, massY);
  else
    printf("==>Created mass point %d\n", massX);

  return;
}

////////////////////////////////////////////
//Systematic Uncertainty Creation methods
////////////////////////////////////////////
void CollieIOFile::createSigSystematic(string syst, TH1D* posFluct, TH1D* negFluct, int parX,
                                       int parY, int parZ) {
  createSigSystematic(0, syst, posFluct, negFluct, parX, parY, parZ);
}

void CollieIOFile::createSigSystematic(int sigIdx, string syst, TH1D* posFluct, TH1D* negFluct,
                                       int parX, int parY, int parZ) {
  if (!check() || !posFluct || !negFluct || sigIdx < 0 || sigIdx >= channel_->getNSignals()) {
    fprintf(stderr,
            "CollieIOFile::createSigSystematic: invalid channel, signal index, or histogram\n");
    return;
  }

  int mindex = mIndex(parX, parY, parZ);

  char title[250];
  snprintf(title, sizeof(title), "Signal %s systematic %d: %s, Pos",
           channel_->getSignalName(sigIdx), mindex, syst.c_str());
  TH1D* shPos = (TH1D*)posFluct->Clone(title);

  snprintf(title, sizeof(title), "Signal %s systematic %d: %s, Neg",
           channel_->getSignalName(sigIdx), mindex, syst.c_str());
  TH1D* shNeg = (TH1D*)negFluct->Clone(title);

  createSigSystematic_Intl(sigIdx, syst, shPos, shNeg, parX, parY, parZ);

  return;
}

void CollieIOFile::createBkgdSystematic(int bkgdIndex, string syst, TH1D* posFluct, TH1D* negFluct,
                                        int parX, int parY, int parZ) {
  if (!check() || !posFluct || !negFluct || bkgdIndex < 0 ||
      bkgdIndex >= channel_->getNBackgrounds()) {
    fprintf(
        stderr,
        "CollieIOFile::createBkgdSystematic: invalid channel, background index, or histogram\n");
    return;
  }
  int mindex = mIndex(parX, parY, parZ);

  char title[250];
  snprintf(title, sizeof(title), "%s systematic %d: %s, Pos",
           channel_->getBackgroundName(bkgdIndex), mindex, syst.c_str());
  TH1D* shPos = (TH1D*)posFluct->Clone(title);

  snprintf(title, sizeof(title), "%s systematic %d: %s, Neg",
           channel_->getBackgroundName(bkgdIndex), mindex, syst.c_str());
  TH1D* shNeg = (TH1D*)negFluct->Clone(title);

  createBkgdSystematic_Intl(bkgdIndex, syst, shPos, shNeg, parX, parY, parZ);

  return;
}

void CollieIOFile::createSigSystematic2D(string syst, TH2D* posFluct, TH2D* negFluct, int parX,
                                         int parY, int parZ) {
  return createSigSystematic2D(0, syst, posFluct, negFluct, parX, parY, parZ);
}

void CollieIOFile::createSigSystematic2D(int sigIndex, string syst, TH2D* posFluct, TH2D* negFluct,
                                         int parX, int parY, int parZ) {
  if (!check() || !posFluct || !negFluct || sigIndex < 0 || sigIndex >= channel_->getNSignals()) {
    fprintf(stderr,
            "CollieIOFile::createSigSystematic2D: invalid channel, signal index, or histogram\n");
    return;
  }
  int mindex = mIndex(parX, parY, parZ);

  char title[250];
  snprintf(title, sizeof(title), "Signal %s systematic %d: %s, Pos",
           channel_->getSignalName(sigIndex), mindex, syst.c_str());
  TH2D* shPos = (TH2D*)posFluct->Clone(title);

  snprintf(title, sizeof(title), "Signal %s systematic %d: %s, Neg",
           channel_->getSignalName(sigIndex), mindex, syst.c_str());
  TH2D* shNeg = (TH2D*)negFluct->Clone(title);

  createSigSystematic2D_Intl(sigIndex, syst, shPos, shNeg, parX, parY, parZ);

  return;
}

void CollieIOFile::createBkgdSystematic2D(int bkgdIndex, string syst, TH2D* posFluct,
                                          TH2D* negFluct, int parX, int parY, int parZ) {
  if (!check() || !posFluct || !negFluct || bkgdIndex < 0 ||
      bkgdIndex >= channel_->getNBackgrounds()) {
    fprintf(
        stderr,
        "CollieIOFile::createBkgdSystematic2D: invalid channel, background index, or histogram\n");
    return;
  }
  int mindex = mIndex(parX, parY, parZ);

  char title[250];
  snprintf(title, sizeof(title), "%s systematic %d: %s, Pos",
           channel_->getBackgroundName(bkgdIndex), mindex, syst.c_str());
  TH2D* shPos = (TH2D*)posFluct->Clone(title);

  snprintf(title, sizeof(title), "%s systematic %d: %s, Neg",
           channel_->getBackgroundName(bkgdIndex), mindex, syst.c_str());
  TH2D* shNeg = (TH2D*)negFluct->Clone(title);

  createBkgdSystematic2D_Intl(bkgdIndex, syst, shPos, shNeg, parX, parY, parZ);

  return;
}

void CollieIOFile::createFlatSigSystematic(string syst, double posV, double negV, int parX,
                                           int parY, int parZ) {
  return createFlatSigSystematic(0, syst, posV, negV, parX, parY, parZ);
}

void CollieIOFile::createFlatSigSystematic(int sigIdx, string syst, double posV, double negV,
                                           int parX, int parY, int parZ) {

  if (!check()) {
    printf("CollieIOFile::createFlatSigSystematic, Error: The file is not properly initialized!\n");
    return;
  }

  int mindex = mIndex(parX, parY, parZ);

  char title[250];
  sprintf(title, "Signal %s systematic %d: %s, Pos", channel_->getSignalName(sigIdx), mindex,
          syst.c_str());
  TH1D* pos = new TH1D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);
  sprintf(title, "Signal %s systematic %d: %s, Neg", channel_->getSignalName(sigIdx), mindex,
          syst.c_str());
  TH1D* neg = new TH1D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

  for (int b = 1; b <= pos->GetNbinsX(); b++) {
    pos->SetBinContent(b, posV);
    neg->SetBinContent(b, negV);
  }

  createSigSystematic_Intl(sigIdx, syst, pos, neg, parX, parY, parZ);

  return;
}

void CollieIOFile::createFlatBkgdSystematic(int bkgdIndex, string syst, double posV, double negV,
                                            int parX, int parY, int parZ) {

  if (!check()) {
    printf(
        "CollieIOFile::createFlatBkgdSystematic, Error: The file is not properly initialized!\n");
    return;
  }

  int mindex = mIndex(parX, parY, parZ);

  char title[250];
  sprintf(title, "%s systematic %d: %s, Pos", channel_->getBackgroundName(bkgdIndex), mindex,
          syst.c_str());
  TH1D* pos = new TH1D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);
  sprintf(title, "%s systematic %d: %s, Neg", channel_->getBackgroundName(bkgdIndex), mindex,
          syst.c_str());
  TH1D* neg = new TH1D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

  for (int b = 1; b <= pos->GetNbinsX(); b++) {
    pos->SetBinContent(b, posV);
    neg->SetBinContent(b, negV);
  }
  createBkgdSystematic_Intl(bkgdIndex, syst, pos, neg, parX, parY, parZ);

  return;
}

void CollieIOFile::createFlatSigSystematic2D(string syst, double posV, double negV, int parX,
                                             int parY, int parZ) {
  return createFlatSigSystematic2D(0, syst, posV, negV, parX, parY, parZ);
}

void CollieIOFile::createFlatSigSystematic2D(int sigIdx, string syst, double posV, double negV,
                                             int parX, int parY, int parZ) {

  if (!check(true)) {
    printf(
        "CollieIOFile::createFlatSigSystematic2D, Error: The file is not properly initialized!\n");
    return;
  }
  int mindex = mIndex(parX, parY, parZ);

  char title[250];
  sprintf(title, "Signal systematic %d: %s, Pos", mindex, syst.c_str());
  TH2D* pos = new TH2D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_,
                       (int)(histBinsY_ * 1.0 / rebinY_), histMinY_, histMaxY_);
  sprintf(title, "Signal systematic %d: %s, Neg", mindex, syst.c_str());
  TH2D* neg = new TH2D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_,
                       (int)(histBinsY_ * 1.0 / rebinY_), histMinY_, histMaxY_);

  for (int bx = 1; bx <= pos->GetNbinsX(); bx++) {
    for (int by = 1; by <= pos->GetNbinsY(); by++) {
      pos->SetBinContent(bx, by, posV);
      neg->SetBinContent(bx, by, negV);
    }
  }
  createSigSystematic2D_Intl(sigIdx, syst, pos, neg, parX, parY, parZ);

  return;
}

void CollieIOFile::createFlatBkgdSystematic2D(int bkgdIndex, string syst, double posV, double negV,
                                              int parX, int parY, int parZ) {

  if (!check(true)) {
    printf(
        "CollieIOFile::createFlatBkgdSystematic2D, Error: The file is not properly initialized!\n");
    return;
  }

  int mindex = mIndex(parX, parY, parZ);

  char title[250];
  sprintf(title, "%s systematic %d: %s, Pos", channel_->getBackgroundName(bkgdIndex), mindex,
          syst.c_str());
  TH2D* pos = new TH2D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_,
                       (int)(histBinsY_ * 1.0 / rebinY_), histMinY_, histMaxY_);
  sprintf(title, "%s systematic %d: %s, Neg", channel_->getBackgroundName(bkgdIndex), mindex,
          syst.c_str());
  TH2D* neg = new TH2D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_,
                       (int)(histBinsY_ * 1.0 / rebinY_), histMinY_, histMaxY_);

  for (int bx = 1; bx <= pos->GetNbinsX(); bx++) {
    for (int by = 1; by <= pos->GetNbinsY(); by++) {
      pos->SetBinContent(bx, by, posV);
      neg->SetBinContent(bx, by, negV);
    }
  }
  createBkgdSystematic2D_Intl(bkgdIndex, syst, pos, neg, parX, parY, parZ);

  return;
}

void CollieIOFile::createSigSystematic_Intl(int sigIdx, string syst, TH1D* pos, TH1D* neg, int parX,
                                            int parY, int parZ) {

  if (!check()) {
    printf("CollieIOFile::createSigSystematic, Error: The file is not properly initialized!\n");
    return;
  }

  if (!checkROOTHisto(pos, true)) {
    printf("CollieIOFile::createSigSystematic, Histogram Error\n");
    return;
  }
  if (!checkROOTHisto(neg, true)) {
    printf("CollieIOFile::createSigSystematic, Histogram Error\n");
    return;
  }

  if (syst.size() < 2) {
    printf("CollieIOFile::createSigSystematic, Error: Invalid systematic name size: syst=%s\n",
           syst.c_str());
    return;
  }

  int nnz = 0;
  int nnan = 0;
  int ninf = 0;
  for (int i = 1; i <= pos->GetNbinsX(); i++) {
    if (fabs(pos->GetBinContent(i)) > 1e-5)
      nnz++;
    if (fabs(neg->GetBinContent(i)) > 1e-5)
      nnz++;
    if (std::isnan(pos->GetBinContent(i)) != 0)
      nnan++;
    if (std::isnan(neg->GetBinContent(i)) != 0)
      nnan++;
    if (std::isinf(pos->GetBinContent(i)) != 0)
      ninf++;
    if (std::isinf(neg->GetBinContent(i)) != 0)
      ninf++;
  }
  if (nnan != 0 || ninf != 0) {
    printf("CollieIOFile::createSigSystematic, Rejecting systematic %s with %d/%d nan/inf bins\n",
           syst.c_str(), nnan, ninf);
    return;
  }
  if (nnz == 0) {
    printf("CollieIOFile::createSigSystematic, Rejecting empty systematic: %s, %s\n", syst.c_str(),
           pos->GetName());
    return;
  }

  CollieMasspoint* mp = NULL;
  CollieDistribution* dist = NULL;

  if (parY == -1 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX);
  else if (parY >= 0 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX, parY);
  else
    mp = channel_->getMasspoint(parX, parY, parZ);
  if (mp)
    dist = mp->getSignalDistMutable(sigIdx);
  else
    printf(
        "CollieIOFile::createSigSystematic, Error: This parameter point does not exist: %d, %d, %d\n",
        parX, parY, parZ);

  if (dist) {
    if (dist->hasSystematic(syst)) {
      printf(
          "CollieIOFile::createSigSystematic, signal %d already has systematic %s.  Rejecting...\n",
          sigIdx, syst.c_str());
      return;
    }

    if (usingSystPriors_) {
      double scl = 1;
      for (int bb = 1; bb <= systPriors_->GetNbinsX(); bb++) {
        string test = systPriors_->GetXaxis()->GetBinLabel(bb);
        if (test == syst)
          scl = systPriors_->GetBinContent(bb);
      }
      //      printf("Scaler %s: %f\n",syst.c_str(),scl);
      pos->Scale(scl);
      neg->Scale(scl);
    }

    dist->addSystematic(syst.c_str(), pos, neg);
    string test = pos->GetName();
    if (systNames_.find(test) == systNames_.end()) {
      systNames_[test] = 1;
      systHists_.push_back(pos);
      systHists_.push_back(neg);
    }
  } else
    printf(
        "CollieIOFile::createSigSystematic, Error: The signal distribution does not exist: %d, %d, %d\n",
        parX, parY, parZ);

  return;
}

void CollieIOFile::createShapeSigSystematic(string syst, TH1D* pos, TH1D* neg, int parX, int parY,
                                            int parZ, double sF, bool norm, int flatten,
                                            int testShape) {
  return createShapeSigSystematic(0, syst, pos, neg, parX, parY, parZ, sF, norm, flatten,
                                  testShape);
}

void CollieIOFile::createShapeSigSystematic(int sigIndex, string syst, TH1D* ipos, TH1D* ineg,
                                            int parX, int parY, int parZ, double sF, bool norm,
                                            int flatten, int testShape) {

  if (ipos == NULL || ineg == NULL) {
    printf("CollieIOFile::createShapeSigSystematic, Error: NULL histograms!\n");
    return;
  }

  if (!check()) {
    printf(
        "CollieIOFile::createShapeSigSystematic, Error: The file is not properly initialized!\n");
    return;
  }

  TH1D* pos = 0;
  if (usingBinMap_)
    pos = getBinMapROOT(ipos, "tmpPos");
  else
    pos = (TH1D*)ipos->Clone("tmpPos");

  TH1D* neg = 0;
  if (usingBinMap_)
    neg = getBinMapROOT(ineg, "tmpNeg");
  else
    neg = (TH1D*)ineg->Clone("tmpNeg");

  if (!checkROOTHisto(pos, false)) {
    printf("CollieIOFile::createShapeSigSystematic, Histogram Error\n");
    return;
  }

  if (!checkROOTHisto(neg, false)) {
    printf("CollieIOFile::createShapeSigSystematic, Histogram Error\n");
    return;
  }

  if (syst.size() < 2) {
    printf("CollieIOFile::createShapeSigSystematic, Error: Invalid name format: syst=%s\n",
           syst.c_str());
    return;
  }

  //Apply input cuts to systematics
  applyCuts(pos);
  applyCuts(neg);

  int nnan = 0;
  int ninf = 0;
  for (int i = 1; i <= pos->GetNbinsX(); i++) {
    if (std::isnan(pos->GetBinContent(i)) != 0)
      nnan++;
    if (std::isnan(neg->GetBinContent(i)) != 0)
      nnan++;
    if (std::isinf(pos->GetBinContent(i)) != 0)
      ninf++;
    if (std::isinf(neg->GetBinContent(i)) != 0)
      ninf++;
  }
  if (nnan != 0 || ninf != 0) {
    printf(
        "CollieIOFile::createShapeSigSystematic, Rejecting systematic %s with %d/%d nan/inf bins\n",
        syst.c_str(), nnan, ninf);
    return;
  }

  int mindex = mIndex(parX, parY, parZ);

  map<string, CollieHistogram*> sigMap = sig_[mindex];
  if (sigMap.find(channel_->getSignalName(sigIndex)) == sigMap.end()) {
    printf(
        "CollieIOFile::createShapeSigSystematic, Error: Cannot find signal %s (parameters: %d,%d,%d))",
        syst.c_str(), parX, parY, parZ);
    return;
  }
  CollieHistogram* refSig = sigMap[channel_->getSignalName(sigIndex)];

  if (refSig == NULL) {
    printf("CollieIOFile::createShapeSigSystematic, NULL reference distribution\n");
    return;
  }

  if (norm) {
    pos->Scale(refSig->Integral() / pos->Integral());
    neg->Scale(refSig->Integral() / neg->Integral());
  }

  //Check to see if the use has given us symmetric (identical) histograms
  bool symm = true;
  if (ipos == ineg) {
    printf(
        "CollieIOFile::createShapeSigSystematic, Warning: Using same input histogram for %s.\n  ===>Will make symmetric.",
        syst.c_str());
    symm = true;
  }
  for (int i = 1; i <= pos->GetNbinsX(); i++) {
    if (pos->GetBinContent(i) != neg->GetBinContent(i)) {
      symm = false;
      break;
    }
    if (symm)
      printf(
          "CollieIOFile::createShapeSigSystematic, Warning: Using same input histogram for %s.\n  ===>Will make symmetric.",
          syst.c_str());
  }

  char title[512];
  sprintf(title, "Shape Signal %s systematic %d: %s, Pos", channel_->getSignalName(sigIndex),
          mindex, syst.c_str());
  TH1D* shPos = new TH1D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);
  sprintf(title, "Shape Signal %s systematic %d: %s, Neg", channel_->getSignalName(sigIndex),
          mindex, syst.c_str());
  TH1D* shNeg = new TH1D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

  //Filter for coherent shape information
  if (flatten == 0 && !systOvrd_)
    flatten = checkForShape(refSig, pos, neg, syst);

  //  assert(flatten==false);

  if (flatten == 1 && !systOvrd_) {
    double vPos = 0.0;
    double vNeg = 0.0;
    double ePos = 0.0;
    double eNeg = 0.0;
    double pI = pos->Integral();
    double nI = neg->Integral();
    double rI = refSig->Integral();
    double pE = totErr(pos, 1, pos->GetNbinsX());
    double nE = totErr(neg, 1, neg->GetNbinsX());
    double rE = totErr(refSig, 1, pos->GetNbinsX());

    if (rI > 0) {
      vPos = (pI - rI) / rI;
      vNeg = (nI - rI) / rI;

      if (pI > 0) {
        ePos = getSystErr(rI, pI, rE, pE);
      }
      if (nI > 0) {
        eNeg = getSystErr(rI, nI, rE, nE);
      }
    }

    for (int b = 1; b <= shPos->GetNbinsX(); b++) {
      shPos->SetBinContent(b, vPos);
      shNeg->SetBinContent(b, vNeg);
      shPos->SetBinError(b, ePos);
      shNeg->SetBinError(b, eNeg);
    }
  } else if (testShape == 1 && !systOvrd_) {
    double co = 0.035;
    if (flatten == -1)
      co = 0.015;
    testShapeSystematics(refSig, pos, shPos, norm, syst, co);
    testShapeSystematics(refSig, neg, shNeg, norm, syst, co);
  } else if (testShape == 2 && !systOvrd_) {
    EqualProbDifference(refSig, pos, shPos, syst, norm);
    EqualProbDifference(refSig, neg, shNeg, syst, norm);
  } else {
    for (int b = 1; b <= shPos->GetNbinsX(); b++) {
      double vPos = 0.0;
      double vNeg = 0.0;
      double ePos = 0.0;
      double eNeg = 0.0;
      double pI = pos->GetBinContent(b);
      double nI = neg->GetBinContent(b);
      double rI = refSig->InBin(b - 1);
      if (rI > 0) {
        vPos = (pI - rI) / rI;
        vNeg = (nI - rI) / rI;

        if (pI > 0) {
          ePos = 1.0 / rI + 1.0 / pI;
        }
        if (nI > 0) {
          eNeg = 1.0 / rI + 1.0 / nI;
        }
        ePos = fabs(vPos * sqrt(ePos));
        eNeg = fabs(vNeg * sqrt(eNeg));
      }

      shPos->SetBinContent(b, vPos);
      shNeg->SetBinContent(b, vNeg);
      shPos->SetBinError(b, ePos);
      shNeg->SetBinError(b, eNeg);
    }
  }

  //if the input histos are NOT the same, flip the negative histo to get the sign right
  if (!symm)
    shNeg->Scale(-1.0);
  else {
    for (int ii = 1; ii <= shPos->GetNbinsX(); ii++)
      shNeg->SetBinContent(ii, shPos->GetBinContent(ii));
  }

  shPos->Scale(sF);
  shNeg->Scale(sF);

  //Now that we're done, check if we have any wacky rates.
  checkRates(refSig, shPos, shNeg, syst);

  createSigSystematic_Intl(sigIndex, syst, shPos, shNeg, parX, parY, parZ);

  delete pos;
  pos = NULL;
  delete neg;
  neg = NULL;

  return;
}

void CollieIOFile::createBkgdSystematic_Intl(int bkgdIndex, string syst, TH1D* pos, TH1D* neg,
                                             int parX, int parY, int parZ) {

  if (bkgdIndex >= (int)nBkgd_) {
    printf("CollieIOFile::createBkgdSystematic, Error: The bkgd index is out of range: index=%d\n",
           bkgdIndex);
    return;
  }

  if (!check()) {
    printf("CollieIOFile::createBkgdSystematic, Error: The file is not properly initialized!\n");
    return;
  }

  if (!checkROOTHisto(pos, true)) {
    printf("CollieIOFile::createBkgdSystematic, Histogram Error\n");
    return;
  }
  if (!checkROOTHisto(neg, true)) {
    printf("CollieIOFile::createBkgdSystematic, Histogram Error\n");
    return;
  }

  if (syst.size() < 2) {
    printf("CollieIOFile::createBkgdSystematic, Error: Invalid name format: syst=%s\n",
           syst.c_str());
    return;
  }

  int nnz = 0;
  int nnan = 0;
  int ninf = 0;
  for (int i = 1; i <= pos->GetNbinsX(); i++) {
    if (fabs(pos->GetBinContent(i)) > 1e-5)
      nnz++;
    if (fabs(neg->GetBinContent(i)) > 1e-5)
      nnz++;
    if (std::isnan(pos->GetBinContent(i)) != 0)
      nnan++;
    if (std::isnan(neg->GetBinContent(i)) != 0)
      nnan++;
    if (std::isinf(pos->GetBinContent(i)) != 0)
      ninf++;
    if (std::isinf(neg->GetBinContent(i)) != 0)
      ninf++;
  }
  if (nnan != 0 || ninf != 0) {
    printf("CollieIOFile::createBkgdSystematic, Rejecting systematic %s with %d/%d nan/inf bins\n",
           syst.c_str(), nnan, ninf);
    return;
  }
  if (nnz == 0) {
    printf("CollieIOFile::createBkgdSystematic, Rejecting empty systematic: %s, %s\n", syst.c_str(),
           pos->GetName());
    return;
  }

  CollieMasspoint* mp = NULL;
  CollieDistribution* dist = NULL;

  if (parY == -1 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX);
  else if (parY >= 0 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX, parY);
  else
    mp = channel_->getMasspoint(parX, parY, parZ);

  if (mp) {
    if (bkgdIndex >= 0) {
      dist = mp->getBkgdDistMutable(bkgdIndex);
      if (dist) {
        if (dist->hasSystematic(syst)) {
          printf(
              "CollieIOFile::createBkgdSystematic, background %d already has systematic %s.  Rejecting...\n",
              bkgdIndex, syst.c_str());
          return;
        }

        if (usingSystPriors_) {
          double scl = 1;
          for (int bb = 1; bb <= systPriors_->GetNbinsX(); bb++) {
            string test = systPriors_->GetXaxis()->GetBinLabel(bb);
            if (test == syst)
              scl = systPriors_->GetBinContent(bb);
          }
          //    printf("Scaler %s: %f\n",syst.c_str(),scl);

          pos->Scale(scl);
          neg->Scale(scl);
        }

        dist->addSystematic(syst.c_str(), pos, neg);
        string test = pos->GetName();
        if (systNames_.find(test) == systNames_.end()) {
          systNames_[test] = 1;
          systHists_.push_back(pos);
          systHists_.push_back(neg);
        }
      } else
        printf(
            "CollieIOFile::createBkgdSystematic, Error: This bkgd distribution does not exist: %d, %d, %d\n",
            parX, parY, parZ);
    } else if (bkgdIndex == -1) {
      for (uint b = 0; b < nBkgd_; b++) {
        dist = mp->getBkgdDistMutable(b);
        if (dist) {
          if (dist->hasSystematic(syst)) {
            printf(
                "CollieIOFile::createBkgdSystematic, background %d already has systematic %s.  Rejecting...\n",
                b, syst.c_str());
            continue;
          }
          dist->addSystematic(syst.c_str(), pos, neg);
          string test = pos->GetName();
          if (systNames_.find(test) == systNames_.end()) {
            systNames_[test] = 1;
            systHists_.push_back(pos);
            systHists_.push_back(neg);
          }
        } else
          printf(
              "CollieIOFile::createBkgdSystematic, Error: This bkgd distribution does not exist: %d, %d, %d\n",
              parX, parY, parZ);
      }
    }
  } else
    printf(
        "CollieIOFile::createBkgdSystematic, Error: This mass point does not exist: %d, %d, %d\n",
        parX, parY, parZ);

  return;
}

void CollieIOFile::createShapeBkgdSystematic(int bkgdIndex, string syst, TH1D* ipos, TH1D* ineg,
                                             int parX, int parY, int parZ, double sF, bool norm,
                                             int flatten, int testShape) {

  if (ipos == NULL || ineg == NULL) {
    printf("CollieIOFile::createShapeBkgdSystematic, Error: NULL histograms!\n");
    return;
  }

  if (bkgdIndex >= (int)nBkgd_) {
    printf(
        "CollieIOFile::createShapeBkgdSystematic, Error: The bkgd index is out of range: index=%d\n",
        bkgdIndex);
    return;
  }

  if (!check()) {
    printf(
        "CollieIOFile::createShapeBkgdSystematic, Error: The file is not properly initialized!\n");
    return;
  }

  TH1D* pos = 0;
  if (usingBinMap_)
    pos = getBinMapROOT(ipos, "tmpPos");
  else
    pos = (TH1D*)ipos->Clone("tmpPos");

  TH1D* neg = 0;
  if (usingBinMap_)
    neg = getBinMapROOT(ineg, "tmpNeg");
  else
    neg = (TH1D*)ineg->Clone("tmpNeg");

  if (!checkROOTHisto(pos, false)) {
    printf("CollieIOFile::createShapeBkgdSystematic, Histogram Error\n");
    return;
  }
  if (!checkROOTHisto(neg, false)) {
    printf("CollieIOFile::createShapeBkgdSystematic, Histogram Error\n");
    return;
  }

  if (syst.size() < 2) {
    printf("CollieIOFile::createShapeBkgdSystematic, Error: Invalid name format: syst=%s\n",
           syst.c_str());
    return;
  }

  //Apply input cuts to systematics
  applyCuts(pos);
  applyCuts(neg);

  int nnan = 0;
  int ninf = 0;
  for (int i = 1; i <= pos->GetNbinsX(); i++) {
    if (std::isnan(pos->GetBinContent(i)) != 0)
      nnan++;
    if (std::isnan(neg->GetBinContent(i)) != 0)
      nnan++;
    if (std::isinf(pos->GetBinContent(i)) != 0)
      ninf++;
    if (std::isinf(neg->GetBinContent(i)) != 0)
      ninf++;
  }

  if (nnan != 0 || ninf != 0) {
    printf(
        "CollieIOFile::createShapeBkgdSystematic, Rejecting systematic %s for Bkgd idx %d with %d/%d nan/inf bins\n",
        syst.c_str(), bkgdIndex, nnan, ninf);
    return;
  }

  int mindex = mIndex(parX, parY, parZ);

  map<string, CollieHistogram*> bkgMap = bkgd_[mindex];
  if (bkgMap.find(channel_->getBackgroundName(bkgdIndex)) == bkgMap.end()) {
    printf(
        "CollieIOFile::createShapeBkgdSystematic, Error: Cannot find background %s (parameters: %d,%d,%d))",
        syst.c_str(), parX, parY, parZ);
    return;
  }
  CollieHistogram* refBkgd = bkgMap[channel_->getBackgroundName(bkgdIndex)];

  if (refBkgd == NULL) {
    printf("CollieIOFile::createShapeBkgdSystematic, NULL reference distribution\n");
    return;
  }

  if (norm) {
    pos->Scale(refBkgd->Integral() / pos->Integral());
    neg->Scale(refBkgd->Integral() / neg->Integral());
  }

  char title[250];
  sprintf(title, "Shape %s systematic %d: %s, Pos", channel_->getBackgroundName(bkgdIndex), mindex,
          syst.c_str());
  TH1D* shPos = new TH1D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);
  sprintf(title, "Shape %s systematic %d: %s, Neg", channel_->getBackgroundName(bkgdIndex), mindex,
          syst.c_str());
  TH1D* shNeg = new TH1D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_);

  //Check to see if the user has given us symmetric (identical) histos
  bool symm = true;
  if (ipos == ineg) {
    printf(
        "CollieIOFile::createShapeSigSystematic, Warning: Using same input histogram for %s.\n  ===>Will make symmetric.",
        syst.c_str());
    symm = true;
  }
  for (int i = 1; i <= pos->GetNbinsX(); i++) {
    if (pos->GetBinContent(i) != neg->GetBinContent(i)) {
      symm = false;
      break;
    }
    if (symm)
      printf(
          "CollieIOFile::createShapeSigSystematic, Warning: Using same input histogram for %s.\n  ===>Will make symmetric.",
          syst.c_str());
  }

  //Filter for coherent shape information
  if (flatten == 0 && !systOvrd_)
    flatten = checkForShape(refBkgd, pos, neg, syst);

  if (flatten == 1 && !systOvrd_) {
    double vPos = 0.0;
    double vNeg = 0.0;
    double ePos = 0.0;
    double eNeg = 0.0;
    double pI = pos->Integral();
    double nI = neg->Integral();
    double rI = refBkgd->Integral();
    if (refBkgd->Integral() > 0) {
      vPos = (pI - rI) / rI;
      vNeg = (nI - rI) / rI;

      if (pI > 0) {
        ePos = 1.0 / rI + 1.0 / pI;
      }
      if (nI > 0) {
        eNeg = 1.0 / rI + 1.0 / nI;
      }
      ePos = vPos * sqrt(ePos);
      eNeg = vNeg * sqrt(eNeg);
    }
    for (int b = 1; b <= shPos->GetNbinsX(); b++) {
      shPos->SetBinContent(b, vPos);
      shNeg->SetBinContent(b, vNeg);
      shPos->SetBinError(b, ePos);
      shNeg->SetBinError(b, eNeg);
    }
  } else if (testShape && !systOvrd_) {
    double co = 0.035;
    if (flatten == -1)
      co = 0.015;
    testShapeSystematics(refBkgd, pos, shPos, norm, syst, co);
    testShapeSystematics(refBkgd, neg, shNeg, norm, syst, co);
  } else {
    for (int b = 1; b <= shPos->GetNbinsX(); b++) {
      double vPos = 0.0;
      double vNeg = 0.0;
      double ePos = 0.0;
      double eNeg = 0.0;
      double pI = pos->GetBinContent(b);
      double nI = neg->GetBinContent(b);
      double rI = refBkgd->InBin(b - 1);
      double pE = pos->GetBinError(b);
      double nE = neg->GetBinError(b);
      double rE = refBkgd->BinErr(b - 1);
      if (rI > 0) {
        vPos = (pI - rI) / rI;
        vNeg = (nI - rI) / rI;

        if (pI > 0) {
          ePos = getSystErr(rI, pI, rE, pE);
        }
        if (nI > 0) {
          eNeg = getSystErr(rI, nI, rE, nE);
        }
      }

      shPos->SetBinContent(b, vPos);
      shNeg->SetBinContent(b, vNeg);
      shPos->SetBinError(b, ePos);
      shNeg->SetBinError(b, eNeg);
    }
  }

  //if the input histos are NOT the same, flip the negative histo to get the sign right
  if (!symm)
    shNeg->Scale(-1.0);
  else {
    for (int ii = 1; ii <= shPos->GetNbinsX(); ii++) {
      shNeg->SetBinContent(ii, shPos->GetBinContent(ii));
    }
  }

  shPos->Scale(sF);
  shNeg->Scale(sF);

  //Now that we're all done, check to see if any rates are crazy.
  checkRates(refBkgd, shPos, shNeg, syst);

  createBkgdSystematic_Intl(bkgdIndex, syst, shPos, shNeg, parX, parY, parZ);

  delete pos;
  pos = NULL;
  delete neg;
  neg = NULL;

  return;
}
void CollieIOFile::createShapeSigSystematic2D(string syst, TH2D* pos, TH2D* neg, int parX, int parY,
                                              int parZ, double sF, bool norm, int flatten,
                                              int testShape) {
  createShapeSigSystematic2D(0, syst, pos, neg, parX, parY, parZ, sF, norm, flatten, testShape);
}

void CollieIOFile::createShapeSigSystematic2D(int sigIndex, string syst, TH2D* pos, TH2D* neg,
                                              int parX, int parY, int parZ, double sF, bool norm,
                                              int flatten, int testShape) {
  if (pos == NULL || neg == NULL) {
    printf("CollieIOFile::createShapeSigSystematic, Error: NULL histograms!\n");
    return;
  }
  if (!check(true)) {
    printf(
        "CollieIOFile::createShapeSigSystematic2D, Error: The file is not properly initialized!\n");
    return;
  }

  if (!checkROOTHisto(pos, false)) {
    printf("CollieIOFile::createShapeSigSystematic2D, Histogram Error\n");
    return;
  }
  if (!checkROOTHisto(neg, false)) {
    printf("CollieIOFile::createShapeSigSystematic2D, Histogram Error\n");
    return;
  }

  if (syst.size() < 2) {
    printf("CollieIOFile::createShapeSigSystematic2D, Error: Invalid name format: syst=%s\n",
           syst.c_str());
    return;
  }

  //Apply input cuts to systematics
  applyCuts(pos);
  applyCuts(neg);

  int nnan = 0;
  int ninf = 0;
  for (int i = 1; i <= pos->GetNbinsX(); i++) {
    for (int j = 1; j <= pos->GetNbinsY(); j++) {
      if (std::isnan(pos->GetBinContent(i, j)) != 0)
        nnan++;
      if (std::isnan(neg->GetBinContent(i, j)) != 0)
        nnan++;
      if (std::isinf(pos->GetBinContent(i, j)) != 0)
        ninf++;
      if (std::isinf(neg->GetBinContent(i, j)) != 0)
        ninf++;
    }
  }
  if (nnan != 0 || ninf != 0) {
    printf(
        "CollieIOFile::createShapeSigSystematic2D, Rejecting systematic %s with %d/%d nan/inf bins\n",
        syst.c_str(), nnan, ninf);
    return;
  }

  int mindex = mIndex(parX, parY, parZ);
  map<string, CollieHistogram*> sigMap = sig_[mindex];
  if (sigMap.find(channel_->getSignalName(sigIndex)) == sigMap.end()) {
    printf(
        "CollieIOFile::createShapeSigSystematic2D, Error: Cannot find signal %s (parameters: %d,%d,%d)",
        syst.c_str(), parX, parY, parZ);
    return;
  }

  CollieHistogram2d* refSig = (CollieHistogram2d*)sigMap[channel_->getSignalName(sigIndex)];

  if (refSig == NULL) {
    printf("CollieIOFile::createShapeSigSystematic2D, NULL reference distribution\n");
    return;
  }

  if (norm) {
    pos->Scale(refSig->Integral() / pos->Integral());
    neg->Scale(refSig->Integral() / neg->Integral());
  }

  char title[250];
  sprintf(title, "Shape Signal systematic %d: %s, Pos", mindex, syst.c_str());
  TH2D* shPos = new TH2D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_,
                         (int)(histBinsY_ * 1.0 / rebinY_), histMinY_, histMaxY_);

  sprintf(title, "Shape Signal systematic %d: %s, Neg", mindex, syst.c_str());
  TH2D* shNeg = new TH2D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_,
                         (int)(histBinsY_ * 1.0 / rebinY_), histMinY_, histMaxY_);

  //Check to see if the use has given us symmetric (identical) histograms
  bool symm = true;
  if (pos == neg) {
    printf(
        "CollieIOFile::createShapeSigSystematic, Warning: Using same input histogram for %s.\n  ===>Will make symmetric.",
        syst.c_str());
    symm = true;
  }
  for (int i = 1; i <= pos->GetNbinsX(); i++) {
    for (int j = 1; j <= pos->GetNbinsY(); j++) {
      if (pos->GetBinContent(i, j) != neg->GetBinContent(i, j)) {
        symm = false;
        break;
      }
    }
    if (symm)
      printf(
          "CollieIOFile::createShapeSigSystematic, Warning: Using same input histogram for %s.\n  ===>Will make symmetric.",
          syst.c_str());
  }

  //Filter for coherent shape information
  if (flatten == 0 && !systOvrd_)
    flatten = checkForShape(refSig, pos, neg, syst);

  if (flatten == 1 && !systOvrd_) {
    double vPos = 0.0;
    double vNeg = 0.0;
    double ePos = 0.0;
    double eNeg = 0.0;
    double pI = pos->Integral();
    double nI = neg->Integral();
    double rI = refSig->Integral();
    double pE = totErr(pos, 1, pos->GetNbinsX());
    double nE = totErr(neg, 1, neg->GetNbinsX());
    double rE = totErr(refSig, 1, pos->GetNbinsX());

    if (rI > 0) {
      vPos = (pI - rI) / rI;
      vNeg = (nI - rI) / rI;

      if (pI > 0) {
        ePos = getSystErr(rI, pI, rE, pE);
      }
      if (nI > 0) {
        eNeg = getSystErr(rI, nI, rE, nE);
      }
    }

    for (int bx = 0; bx < shPos->GetNbinsX(); bx++) {
      for (int by = 0; by < shPos->GetNbinsY(); by++) {
        shPos->SetBinContent(bx + 1, by + 1, vPos);
        shNeg->SetBinContent(bx + 1, by + 1, vNeg);
        shPos->SetBinError(bx + 1, by + 1, ePos);
        shNeg->SetBinError(bx + 1, by + 1, eNeg);
      }
    }
  } else {
    for (int bx = 0; bx < shPos->GetNbinsX(); bx++) {
      for (int by = 0; by < shPos->GetNbinsY(); by++) {
        double refBin = refSig->InBin(bx, by);
        if (refBin > 0) {
          double p = (pos->GetBinContent(bx + 1, by + 1) - refBin) / refBin;
          double n = (neg->GetBinContent(bx + 1, by + 1) - refBin) / refBin;
          if (p > 0.5)
            p = 0.5;
          if (n > 0.5)
            n = 0.5;
          if (p < -0.5)
            p = -0.5;
          if (n < -0.5)
            n = -0.5;
          shPos->SetBinContent(bx + 1, by + 1, p);
          shNeg->SetBinContent(bx + 1, by + 1, n);
        }
      }
    }
  }

  //if the input histos are NOT the same, flip the negative histo to get the sign right
  if (!symm)
    shNeg->Scale(-1.0);
  else {
    for (int ii = 1; ii <= shPos->GetNbinsX(); ii++)
      shNeg->SetBinContent(ii, shPos->GetBinContent(ii));
  }
  shPos->Scale(sF);
  shNeg->Scale(sF);

  //Now that we're done, check for crazy rates
  checkRates(refSig, shPos, shNeg, syst);

  createSigSystematic2D_Intl(sigIndex, syst, shPos, shNeg, parX, parY, parZ);

  return;
}

void CollieIOFile::createShapeBkgdSystematic2D(int bkgdIndex, string syst, TH2D* pos, TH2D* neg,
                                               int parX, int parY, int parZ, double sF, bool norm,
                                               int flatten, int testShape) {
  if (pos == NULL || neg == NULL) {
    printf("CollieIOFile::createShapeBkgdSystematic2D, Error: NULL histograms!\n");
    return;
  }
  if (bkgdIndex >= (int)nBkgd_) {
    printf(
        "CollieIOFile::createShapeBkgdSystematic2D, Error: The bkgd index is out of range: index=%d\n",
        bkgdIndex);
    return;
  }

  if (!check(true)) {
    printf(
        "CollieIOFile::createShapeBkgdSystematic2D, Error: The file is not properly initialized!\n");
    return;
  }

  if (!checkROOTHisto(pos, false)) {
    printf("CollieIOFile::createShapeBkgdSystematic2D, Histogram Error\n");
    return;
  }
  if (!checkROOTHisto(neg, false)) {
    printf("CollieIOFile::createShapeBkgdSystematic2D, Histogram Error\n");
    return;
  }

  if (syst.size() < 2) {
    printf("CollieIOFile::createShapeBkgdSystematic2D, Error: Invalid name format: syst=%s\n",
           syst.c_str());
    return;
  }

  //Apply input cuts to systematics
  applyCuts(pos);
  applyCuts(neg);

  int nnan = 0;
  int ninf = 0;
  for (int i = 1; i <= pos->GetNbinsX(); i++) {
    for (int j = 1; j <= pos->GetNbinsY(); j++) {
      if (std::isnan(pos->GetBinContent(i, j)) != 0)
        nnan++;
      if (std::isnan(neg->GetBinContent(i, j)) != 0)
        nnan++;
      if (std::isinf(pos->GetBinContent(i, j)) != 0)
        ninf++;
      if (std::isinf(neg->GetBinContent(i, j)) != 0)
        ninf++;
    }
  }
  if (nnan != 0 || ninf != 0) {
    printf(
        "CollieIOFile::createShapeBkgdSystematic2D, Rejecting systematic %s for Bkgd idx %d with %d/%d nan/inf bins\n",
        syst.c_str(), bkgdIndex, nnan, ninf);
    return;
  }

  int mindex = mIndex(parX, parY, parZ);
  map<string, CollieHistogram*> bkgMap = bkgd_[mindex];
  if (bkgMap.find(channel_->getBackgroundName(bkgdIndex)) == bkgMap.end()) {
    printf(
        "CollieIOFile::createShapeBkgdSystematic2D, Error: Cannot find background %s (parameters: %d,%d,%d)",
        syst.c_str(), parX, parY, parZ);
    return;
  }
  CollieHistogram2d* refBkgd = (CollieHistogram2d*)bkgMap[channel_->getBackgroundName(bkgdIndex)];

  if (refBkgd == NULL) {
    printf("CollieIOFile::createShapeBkgdSystematic2D, NULL reference distribution\n");
    return;
  }

  if (norm) {
    pos->Scale(refBkgd->Integral() / pos->Integral());
    neg->Scale(refBkgd->Integral() / neg->Integral());
  }

  char title[250];
  sprintf(title, "Shape Bkgd systematic %d: %s, Pos", mindex, syst.c_str());
  TH2D* shPos = new TH2D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_,
                         (int)(histBinsY_ * 1.0 / rebinY_), histMinY_, histMaxY_);

  sprintf(title, "Shape Bkgd systematic %d: %s, Neg", mindex, syst.c_str());
  TH2D* shNeg = new TH2D(title, title, (int)(histBinsX_ * 1.0 / rebinX_), histMinX_, histMaxX_,
                         (int)(histBinsY_ * 1.0 / rebinY_), histMinY_, histMaxY_);

  //Check to see if the use has given us symmetric (identical) histograms
  bool symm = true;
  if (pos == neg) {
    printf(
        "CollieIOFile::createShapeBkgdSystematic2D, Warning: Using same input histogram for %s.\n  ===>Will make symmetric.",
        syst.c_str());
    symm = true;
  }
  for (int i = 1; i <= pos->GetNbinsX(); i++) {
    for (int j = 1; j <= pos->GetNbinsY(); j++) {
      if (pos->GetBinContent(i, j) != neg->GetBinContent(i, j)) {
        symm = false;
        break;
      }
    }
    if (symm)
      printf(
          "CollieIOFile::createShapeSigSystematic, Warning: Using same input histogram for %s.\n  ===>Will make symmetric.",
          syst.c_str());
  }

  //Filter for coherent shape dependence...
  if (flatten == 0)
    flatten = checkForShape(refBkgd, pos, neg, syst);

  if (flatten == 1) {
    double vPos = 0.0;
    double vNeg = 0.0;
    double ePos = 0.0;
    double eNeg = 0.0;
    double pI = pos->Integral();
    double nI = neg->Integral();
    double rI = refBkgd->Integral();
    if (refBkgd->Integral() > 0) {
      vPos = (pI - rI) / rI;
      vNeg = (nI - rI) / rI;

      if (pI > 0) {
        ePos = 1.0 / rI + 1.0 / pI;
      }
      if (nI > 0) {
        eNeg = 1.0 / rI + 1.0 / nI;
      }
      ePos = fabs(vPos * sqrt(ePos));
      eNeg = fabs(vNeg * sqrt(eNeg));
    }
    for (int bx = 0; bx < shPos->GetNbinsX(); bx++) {
      for (int by = 0; by < shPos->GetNbinsY(); by++) {
        shPos->SetBinContent(bx + 1, by + 1, vPos);
        shNeg->SetBinContent(bx + 1, by + 1, vNeg);
        shPos->SetBinError(bx + 1, by + 1, ePos);
        shNeg->SetBinError(bx + 1, by + 1, eNeg);
      }
    }
  } else {
    for (int bx = 0; bx < shPos->GetNbinsX(); bx++) {
      for (int by = 0; by < shPos->GetNbinsY(); by++) {
        double refBin = refBkgd->InBin(bx, by);
        if (refBin > 0) {
          double p = (pos->GetBinContent(bx + 1, by + 1) - refBin) / refBin;
          double n = (neg->GetBinContent(bx + 1, by + 1) - refBin) / refBin;
          if (p > 0.5)
            p = 0.5;
          if (n > 0.5)
            n = 0.5;
          if (p < -0.5)
            p = -0.5;
          if (n < -0.5)
            n = -0.5;
          shPos->SetBinContent(bx + 1, by + 1, p);
          shNeg->SetBinContent(bx + 1, by + 1, n);
        }
      }
    }
  }

  //if the input histos are NOT the same, flip the negative histo to get the sign right
  if (!symm)
    shNeg->Scale(-1.0);
  else {
    for (int ii = 1; ii <= shPos->GetNbinsX(); ii++)
      shNeg->SetBinContent(ii, shPos->GetBinContent(ii));
  }

  shPos->Scale(sF);
  shNeg->Scale(sF);

  //Now that we're done, check for crazy rates
  checkRates(refBkgd, shPos, shNeg, syst);

  createBkgdSystematic2D_Intl(bkgdIndex, syst, shPos, shNeg, parX, parY, parZ);

  return;
}

void CollieIOFile::createSigSystematic2D_Intl(int sigIdx, string syst, TH2D* pos, TH2D* neg,
                                              int parX, int parY, int parZ) {

  if (!check(true)) {
    printf("CollieIOFile::createSigSystematic2D, Error: The file is not properly initialized!\n");
    return;
  }

  if (!checkROOTHisto(pos, true)) {
    printf("CollieIOFile::createSigSystematic2D, Histogram Error\n");
    return;
  }
  if (!checkROOTHisto(neg, true)) {
    printf("CollieIOFile::createSigSystematic2D, Histogram Error\n");
    return;
  }

  if (syst.size() < 2) {
    printf("CollieIOFile::createSigSystematic2D, Error: Invalid name format: syst=%s\n",
           syst.c_str());
    return;
  }

  //Apply input cuts to systematics
  applyCuts(pos);
  applyCuts(neg);

  int nnz = 0;
  int nnan = 0;
  int ninf = 0;
  for (int i = 1; i <= pos->GetNbinsX(); i++) {
    for (int j = 1; j <= pos->GetNbinsY(); j++) {
      if (fabs(pos->GetBinContent(i, j)) > 1e-5)
        nnz++;
      if (fabs(neg->GetBinContent(i, j)) > 1e-5)
        nnz++;
      if (std::isnan(pos->GetBinContent(i, j)) != 0)
        nnan++;
      if (std::isnan(neg->GetBinContent(i, j)) != 0)
        nnan++;
      if (std::isinf(pos->GetBinContent(i, j)) != 0)
        ninf++;
      if (std::isinf(neg->GetBinContent(i, j)) != 0)
        ninf++;
    }
  }
  if (nnan != 0 || ninf != 0) {
    printf("CollieIOFile::createSigSystematic2D, Rejecting systematic %s with %d/%d nan/inf bins\n",
           syst.c_str(), nnan, ninf);
    return;
  }
  if (nnz == 0) {
    printf("CollieIOFile::createSigSystematic2D, Rejecting empty systematic: %s, %s\n",
           syst.c_str(), pos->GetName());
    return;
  }

  CollieMasspoint* mp = NULL;
  CollieDistribution* dist = NULL;

  if (parY == -1 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX);
  else if (parY >= 0 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX, parY);
  else
    mp = channel_->getMasspoint(parX, parY, parZ);

  if (mp)
    dist = mp->getSignalDistMutable(sigIdx);
  else
    printf(
        "CollieIOFile::createSigSystematic2D, Error: This mass point does not exist: %d, %d, %d\n",
        parX, parY, parZ);

  if (dist) {
    if (dist->hasSystematic(syst)) {
      printf(
          "CollieIOFile::createSigSystematic2D, signal %d already has systematic %s.  Rejecting...\n",
          sigIdx, syst.c_str());
      return;
    }
    dist->addSystematic2D(syst.c_str(), pos, neg);
    string test = pos->GetName();
    if (systNames_.find(test) == systNames_.end()) {
      systNames_[test] = 1;
      systHists_.push_back(pos);
      systHists_.push_back(neg);
    }
  } else
    printf(
        "CollieIOFile::createSigSystematic2D, Error: This signal distribution does not exist: %d, %d, %d\n",
        parX, parY, parZ);
  return;
}

void CollieIOFile::createBkgdSystematic2D_Intl(int bkgdIndex, string syst, TH2D* pos, TH2D* neg,
                                               int parX, int parY, int parZ) {

  if (bkgdIndex >= (int)nBkgd_) {
    printf("CollieIOFile::createBkgdSystematic, Error: The bkgd index is out of range: index=%d\n",
           bkgdIndex);
    return;
  }

  if (!check()) {
    printf("CollieIOFile::createBkgdSystematic2D, Error: The file is not properly initialized!\n");
    return;
  }

  if (!checkROOTHisto(pos, true)) {
    printf("CollieIOFile::createSigSystematic2D, Histogram Error\n");
    return;
  }
  if (!checkROOTHisto(neg, true)) {
    printf("CollieIOFile::createSigSystematic2D, Histogram Error\n");
    return;
  }

  if (syst.size() < 2) {
    printf("CollieIOFile::createBkgdSystematic2D, Error: Invalid name format: syst=%s\n",
           syst.c_str());
    return;
  }

  //Apply input cuts to systematics
  applyCuts(pos);
  applyCuts(neg);

  int nnz = 0;
  int nnan = 0;
  int ninf = 0;
  for (int i = 1; i <= pos->GetNbinsX(); i++) {
    for (int j = 1; j <= pos->GetNbinsX(); j++) {
      if (fabs(pos->GetBinContent(i, j)) > 1e-5)
        nnz++;
      if (fabs(neg->GetBinContent(i, j)) > 1e-5)
        nnz++;
      if (std::isnan(pos->GetBinContent(i, j)) != 0)
        nnan++;
      if (std::isnan(neg->GetBinContent(i, j)) != 0)
        nnan++;
      if (std::isinf(pos->GetBinContent(i, j)) != 0)
        ninf++;
      if (std::isinf(neg->GetBinContent(i, j)) != 0)
        ninf++;
    }
  }
  if (nnan != 0 || ninf != 0) {
    printf(
        "CollieIOFile::createBkgdSystematic2D, Rejecting systematic %s with %d/%d nan/inf bins\n",
        syst.c_str(), nnan, ninf);
    return;
  }
  if (nnz == 0) {
    printf("CollieIOFile::createBkgdSystematic2D, Rejecting empty systematic: %s, %s\n",
           syst.c_str(), pos->GetName());
    return;
  }

  CollieMasspoint* mp = NULL;
  CollieDistribution* dist = NULL;

  if (parY == -1 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX);
  else if (parY >= 0 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX, parY);
  else
    mp = channel_->getMasspoint(parX, parY, parZ);
  if (mp) {
    if (bkgdIndex >= 0) {
      dist = mp->getBkgdDistMutable(bkgdIndex);
      if (dist) {
        if (dist->hasSystematic(syst)) {
          printf(
              "CollieIOFile::createBkgdSystematic2D, background %d already has systematic %s.  Rejecting...\n",
              bkgdIndex, syst.c_str());
          return;
        }
        dist->addSystematic2D(syst.c_str(), pos, neg);
        string test = pos->GetName();
        if (systNames_.find(test) == systNames_.end()) {
          systNames_[test] = 1;
          systHists_.push_back(pos);
          systHists_.push_back(neg);
        }
      } else
        printf(
            "CollieIOFile::createBkgdSystematic2D, Error: This bkgd distribution does not exist: %d, %d, %d\n",
            parX, parY, parZ);
    } else if (bkgdIndex == -1) {
      for (uint b = 0; b < nBkgd_; b++) {
        dist = mp->getBkgdDistMutable(b);
        if (dist) {
          if (dist->hasSystematic(syst)) {
            printf(
                "CollieIOFile::createBkgdSystematic2D, background %d already has systematic %s.  Rejecting...\n",
                b, syst.c_str());
            continue;
          }
          dist->addSystematic2D(syst.c_str(), pos, neg);
          string test = pos->GetName();
          if (systNames_.find(test) == systNames_.end()) {
            systNames_[test] = 1;
            systHists_.push_back(pos);
            systHists_.push_back(neg);
          }
        } else
          printf(
              "CollieIOFile::createBkgdSystematic2D, Error: This bkgd distribution does not exist: %d, %d, %d\n",
              parX, parY, parZ);
      }
    }
  } else
    printf(
        "CollieIOFile::createBkgdSystematic2D, Error: This mass point does not exist: %d, %d, %d\n",
        parX, parY, parZ);

  return;
}

////////////////////////
/// Utility Methods
//////////////////////////
void CollieIOFile::setSigFloatFlag(string syst, bool floatIt, int parX, int parY, int parZ) {

  if (!check()) {
    printf("CollieIOFile::setFloatFlag, Error: The file is not properly initialized!\n");
    return;
  }
  CollieMasspoint* mp = NULL;
  CollieDistribution* dist = NULL;

  if (parY == -1 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX);
  else if (parY >= 0 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX, parY);
  else
    mp = channel_->getMasspoint(parX, parY, parZ);
  for (uint s = 0; s < sig_.size(); s++) {
    if (mp) {
      dist = mp->getSignalDistMutable(s);
      if (dist)
        dist->setFloatFlag(syst.c_str(), floatIt);
      else
        printf(
            "CollieIOFile::setFloatFlag, Error: This sig distribution does not exist: %d, %d, %d\n",
            parX, parY, parZ);
    } else
      printf("CollieIOFile::setFloatFlag, Error: This mass point does not exist: %d, %d, %d\n",
             parX, parY, parZ);
  }
  return;
}

void CollieIOFile::setBkgdFloatFlag(int bkgdIndex, string syst, bool floatIt, int parX, int parY,
                                    int parZ) {

  if (!check()) {
    printf("CollieIOFile::setFloatFlag, Error: The file is not properly initialized!\n");
    return;
  }
  if (bkgdIndex >= (int)nBkgd_) {
    printf("CollieIOFile::setFloatFlag, Error: The bkgd index is out of range: index=%d\n",
           bkgdIndex);
    return;
  }

  CollieMasspoint* mp = NULL;
  CollieDistribution* dist = NULL;

  if (parY == -1 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX);
  else if (parY >= 0 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX, parY);
  else
    mp = channel_->getMasspoint(parX, parY, parZ);

  if (mp) {
    if (bkgdIndex >= 0) {
      dist = mp->getBkgdDistMutable(bkgdIndex);
      if (dist)
        dist->setFloatFlag(syst.c_str(), floatIt);
      else
        printf(
            "CollieIOFile::setFloatFlag, Error: This bkgd distribution does not exist: %d, %d, %d\n",
            parX, parY, parZ);
    }
  } else
    printf("CollieIOFile::setFloatFlag, Error: This mass point does not exist: %d, %d, %d\n", parX,
           parY, parZ);

  return;
}

void CollieIOFile::setLogNormalFlag(string syst, bool isNorm, int parX, int parY, int parZ) {
  if (!check()) {
    printf("CollieIOFile::setLogNormalFlag, Error: The file is not properly initialized!\n");
    return;
  }

  CollieMasspoint* mp = NULL;
  CollieDistribution* dist = NULL;

  if (parY == -1 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX);
  else if (parY >= 0 && parZ == -1 && parX >= 0)
    mp = channel_->getMasspoint(parX, parY);
  else
    mp = channel_->getMasspoint(parX, parY, parZ);
  if (mp) {
    for (uint s = 0; s < nSig_; s++) {
      dist = mp->getSignalDistMutable(s);
      if (dist)
        dist->setLogNormalFlag(syst.c_str(), isNorm);
      else
        printf(
            "CollieIOFile::setLogNormalFlag, Error: This sig distribution does not exist: %d, %d, %d\n",
            parX, parY, parZ);
    }
    for (uint bg = 0; bg < nBkgd_; bg++) {
      dist = mp->getBkgdDistMutable(bg);
      if (dist)
        dist->setLogNormalFlag(syst.c_str(), isNorm);
      else
        printf(
            "CollieIOFile::setLogNormalFlag, Error: This bkgd distribution does not exist: %d, %d, %d\n",
            parX, parY, parZ);
    }
  } else
    printf("CollieIOFile::setLogNormalFlag, Error: This mass point does not exist: %d, %d, %d\n",
           parX, parY, parZ);

  return;
}

bool CollieIOFile::getSmoothedROOT(double alpha, const TH1D* in, TH1D* out) {
  if (!check()) {
    printf("CollieIOFile::getSmoothedROOT, The file is not properly initialized!\n");
    return false;
  }

  CollieHistogram* c1 = new CollieHistogram();
  c1->Book("c1", in->GetNbinsX(), in->GetXaxis()->GetXmin(), in->GetXaxis()->GetXmax());

  if (!getSmoothedCollie(alpha, in, c1)) {
    printf("CollieIOFile::getSmoothedROOT, Error extracting smoothed histo\n");
    return false;
  }
  if (!getROOT(c1, out)) {
    printf("CollieIOFile::getSmoothedROOT, Error extracting output histo\n");
    return false;
  }
  delete c1;
  c1 = NULL;
  return false;
}

bool CollieIOFile::getSmoothedROOT(double alpha, const CollieHistogram* in, TH1D* out) {
  if (!check()) {
    printf("CollieIOFile::getSmoothedROOT, The file is not properly initialized!\n");
    return false;
  }

  CollieHistogram tmp(*in);
  if (!getSmoothedCollie(alpha, in, &tmp)) {
    printf("CollieIOFile::getSmoothedROOT, Error extracting smoothed histo\n");
    return false;
  }
  if (!getROOT(&tmp, out)) {
    printf("CollieIOFile::getSmoothedROOT, Error extracting output histo\n");
    return false;
  }

  return false;
}

bool CollieIOFile::getSmoothedCollie(double alpha, const TH1D* in, CollieHistogram* out) {
  if (!check()) {
    printf("CollieIOFile::getSmoothedCollie, Error: The file is not properly initialized!\n");
    return false;
  }
  if (in == NULL) {
    printf("CollieIOFile::getSmoothedCollie, Error: NULL histogram(s)!!\n");
    return false;
  }
  if (in->Integral() == 0) {
    printf("CollieIOFile::getSmoothedCollie, Error: Empty histogram(s)!!\n");
    return false;
  }

  if (noviceFlag_) {
    printf(
        "CollieIOFile::getSmoothedCollie, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing smoothing\n");
    return false;
  }

  //vars for smoothed hist
  CollieQuickKEYS smoother;
  char name[128];
  sprintf(name, "smoothed hist %d", smIncr_);

  smoother.Book(name, in->GetNbinsX(), in->GetXaxis()->GetXmin(), in->GetXaxis()->GetXmax());

  int iter = (int)(histNorm_); //// in->Integral());
  double sm = 100.0;
  if (in->GetMean() != 0)
    sm = in->GetRMS() / in->GetMean();
  if (sm > 0.25)
    iter /= 5;

  //Fill the intermediate smoothing histo
  // and collect info on statistical uncertainties
  double val = 0;
  double into = 0;
  double statErrTot = 0;
  double sqrInt = 0;
  for (int i = 1; i <= in->GetNbinsX(); i++) {
    val = in->GetBinCenter(i);

    into = in->GetBinContent(i);
    statErrTot += in->GetBinError(i) * in->GetBinError(i);
    sqrInt += into * into;
    for (int n = 0; n < (into * iter); n++) {
      smoother.Fill(val, 1.0 / iter);
    }
  }
  statErrTot = sqrt(statErrTot / (sqrInt + 1e-6));

  //prepare histos for merging
  if (alpha > 0)
    smoother.Smooth(alpha);
  out->Clear();

  //merge back
  //  step = (smoother.MaxEdge() - smoother.MinEdge())/smoother.nBins();
  //  int bin = 0;
  for (int i = 0; i < smoother.nBins(); i++) {
    //    val = (i*step+step/10000.0+smoother.MinEdge());
    //    bin = out->FindBin(val);
    out->FillBin(i, smoother.InBin(i));
    out->SetBinError(i, statErrTot * smoother.InBin(i));
  }

  return true;
}

bool CollieIOFile::getSmoothedCollie(double alpha, const CollieHistogram* in,
                                     CollieHistogram* out) {
  if (!check()) {
    printf("CollieIOFile::getSmoothedCollie, Error: The file is not properly initialized!\n");
    return false;
  }

  if (noviceFlag_) {
    printf(
        "CollieIOFile::getSmoothedCollie, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing smoothing\n");
    return false;
  }

  TH1D* tmpHist = new TH1D("tmphist", "tmphist", in->nBins(), in->MinEdge(), in->MaxEdge());
  tmpHist->Sumw2();
  getROOT(in, tmpHist);
  getSmoothedCollie(alpha, tmpHist, out);
  tmpHist->Delete();
  return true;
}

bool CollieIOFile::getSmoothedCollie2D(double alpha, const TH2D* in, CollieHistogram2d* out) {
  if (!check(true)) {
    printf("CollieIOFile::getSmoothedCollie2D, Error: The file is not properly initialized!\n");
    return false;
  }
  if (in == NULL) {
    printf("CollieIOFile::getSmoothedCollie2D, Error: NULL histogram(s)!!\n");
    return false;
  }
  if (in->Integral() == 0)
    return false;

  printf("CollieIOFile::getSmoothedCollie2D, cannot smooth in 2 dimensions yet...\n");

  if (noviceFlag_) {
    printf(
        "CollieIOFile::getSmoothedCollie2D, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing smoothing\n");
    return false;
  }

  getCollie2D(in, out);

  return true;
}

bool CollieIOFile::getSmoothedCollie2D(double alpha, const CollieHistogram2d* in,
                                       CollieHistogram2d* out) {
  if (!check(true)) {
    printf("CollieIOFile::getSmoothedCollie2D, Error: The file is not properly initialized!\n");
    return false;
  }
  if (in == NULL) {
    printf("CollieIOFile::getSmoothedCollie2D, Error: NULL histogram(s)!!\n");
    return false;
  }
  if (in->Integral() == 0)
    return false;

  printf("CollieIOFile::getSmoothedCollie2D, cannot smooth in 2 dimensions yet...\n");

  if (noviceFlag_) {
    printf(
        "CollieIOFile::getSmoothedCollie2D, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing smoothing\n");
    return false;
  }
  TH2D* tmpHist = new TH2D("tmphist", "tmphist", in->nxBins(), in->LeftEdge(), in->RightEdge(),
                           in->nyBins(), in->TopEdge(), in->BottomEdge());
  getROOT2D(in, tmpHist);
  getSmoothedCollie2D(alpha, tmpHist, out);
  tmpHist->Delete();
  return true;
}

bool CollieIOFile::getGroupedBinCollie(double group, const TH1D* in, CollieHistogram* out) {
  if (!check()) {
    printf("CollieIOFile::getGroupedBinCollie, Error: The file is not properly initialized!\n");
    return false;
  }
  if (in == NULL) {
    printf("CollieIOFile::getGroupedBinCollie, Error: NULL histogram(s)!!\n");
    return false;
  }
  if (in->Integral() == 0) {
    printf("CollieIOFile::getGroupedBinCollie, Error: Empty histogram(s)!!\n");
    return false;
  }

  if (noviceFlag_) {
    printf(
        "CollieIOFile::getGroupedBinCollie, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing smoothing\n");
    return false;
  }

  int rbFactor = int(fabs(group));
  if (rbFactor < 2) {
    printf("CollieIOFile::getGroupedBinCollie, Error: Illegal rebinning factor!! (%d)\n", rbFactor);
    return false;
  }

  TH1D* rbClone = (TH1D*)in->Clone("tmpGroupHisto");
  rbClone->Rebin(rbFactor);

  out->Clear();
  for (int i = 1; i <= in->GetNbinsX(); i++) {
    int inbin = int((i - 1) / rbFactor) + 1;
    out->FillBin(i - 1, rbClone->GetBinContent(inbin) / rbFactor);

    float inerr = rbClone->GetBinError(inbin);
    inerr = sqrt(inerr * inerr / rbFactor);
    out->SetBinError(i - 1, inerr);
  }
  delete rbClone;
  rbClone = NULL;

  return true;
}

bool CollieIOFile::getGroupedBinCollie(double group, const CollieHistogram* in,
                                       CollieHistogram* out) {
  if (!check()) {
    printf("CollieIOFile::getGroupedBinCollie, Error: The file is not properly initialized!\n");
    return false;
  }

  if (noviceFlag_) {
    printf(
        "CollieIOFile::getGroupedBinCollie, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing smoothing\n");
    return false;
  }

  TH1D* tmpHist =
      new TH1D("tmphistgroup", "tmphistgroup", in->nBins(), in->MinEdge(), in->MaxEdge());
  tmpHist->Sumw2();
  getROOT(in, tmpHist);
  getGroupedBinCollie(group, tmpHist, out);
  tmpHist->Delete();
  return true;
}

bool CollieIOFile::getGroupedBinCollie2D(double group, const TH2D* in, CollieHistogram2d* out) {
  if (!check(true)) {
    printf("CollieIOFile::getGroupedBinCollie2D, Error: The file is not properly initialized!\n");
    return false;
  }
  if (in == NULL) {
    printf("CollieIOFile::getGroupedBinCollie2D, Error: NULL histogram(s)!!\n");
    return false;
  }
  if (in->Integral() == 0)
    return false;

  if (noviceFlag_) {
    printf(
        "CollieIOFile::getGroupedBinCollie2D, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing smoothing\n");
    return false;
  }

  int rbFactor = int(fabs(group));
  if (rbFactor < 2) {
    printf("CollieIOFile::getGroupedBinCollie2D, Error: Illegal rebinning factor!! (%d)\n",
           rbFactor);
    return false;
  }

  double rbX = sqrt(rbFactor);
  double rbY = 1;
  if (rbX > 2) {
    rbX = (int)rbX;
    rbY = (int)(rbFactor / rbX);
  } else {
    rbX = rbFactor;
  }

  TH2D* rbClone = (TH2D*)in->Clone("tmpGroupHisto");
  rbClone->RebinX((int)rbX);
  rbClone->RebinY((int)rbY);

  out->Clear();
  for (int bx = 1; bx <= in->GetNbinsX(); bx++) {
    for (int by = 1; by <= in->GetNbinsX(); by++) {
      int inbinX = int((bx - 1) / rbX) + 1;
      int inbinY = int((by - 1) / rbY) + 1;
      out->Fill(bx - 1, by - 1, rbClone->GetBinContent(inbinX, inbinY) / (rbX * rbY));

      float inerr = rbClone->GetBinError(inbinX, inbinY);
      inerr = sqrt(inerr * inerr / rbFactor);
      out->SetBinError(bx - 1, by - 1, inerr);
    }
  }

  delete rbClone;
  rbClone = NULL;

  return true;
}

bool CollieIOFile::getGroupedBinCollie2D(double group, const CollieHistogram2d* in,
                                         CollieHistogram2d* out) {
  if (!check(true)) {
    printf("CollieIOFile::getGroupedBinCollie2D, Error: The file is not properly initialized!\n");
    return false;
  }

  if (in == NULL) {
    printf("CollieIOFile::getGroupedBinCollie2D, Error: NULL histogram(s)!!\n");
    return false;
  }

  if (in->Integral() == 0)
    return false;

  printf("CollieIOFile::getGroupedBinCollie2D, cannot smooth in 2 dimensions yet...\n");

  if (noviceFlag_) {
    printf(
        "CollieIOFile::getGroupedBinCollie2D, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing smoothing\n");
    return false;
  }
  TH2D* tmpHist = new TH2D("tmphistgroup", "tmphistgroup", in->nxBins(), in->LeftEdge(),
                           in->RightEdge(), in->nyBins(), in->TopEdge(), in->BottomEdge());
  getROOT2D(in, tmpHist);
  getGroupedBinCollie2D(group, tmpHist, out);
  tmpHist->Delete();
  return true;
}

bool CollieIOFile::getROOT(const CollieHistogram* in, TH1D* out) {
  if (!check()) {
    printf("CollieIOFile::getROOT, Error: The file is not properly initialized!\n");
    return false;
  }

  if (in == NULL || out == NULL) {
    printf("CollieIOFile::getROOT, Error: NULL histos!\n");
    return false;
  }

  if (!checkROOTHisto(out, false)) {
    printf("CollieIOFile::getROOT, Error: NULL histos!\n");
    return false;
  }

  out->Scale(0);

  for (int i = 0; i < in->nBins(); i++) {
    out->SetBinContent(i + 1, in->InBin(i));
    out->SetBinError(i + 1, in->BinErr(i));
  }
  return true;
}

bool CollieIOFile::getROOT2D(const CollieHistogram2d* in, TH2D* out) {
  if (!check()) {
    printf("CollieIOFile::getROOT2D, Error: The file is not properly initialized!\n");
    return false;
  }

  out->Scale(0);
  for (int i = 0; i < in->nxBins(); i++) {
    for (int j = 0; j < in->nxBins(); j++) {
      out->SetBinContent(i + 1, j + 1, in->InBin(i, j));
      out->SetBinError(i + 1, j + 1, in->BinErr(i, j));
    }
  }

  return true;
}

TH1D* CollieIOFile::getBinMapROOT(const TH1D* in, const char* name) {

  if (in == NULL) {
    printf("CollieIOFile::getRebinnedROOT, Error: NULL histograms!\n");
    return NULL;
  }

  if (!usingBinMap_) {
    printf("CollieIOFile::getRebinnedROOT, Error: No bin map generated!\n");
    return NULL;
  }

  TH1D* htmp = new TH1D(name, name, histBinsX_, &binEdges_[0]);
  for (int i = 1; i <= in->GetNbinsX(); i++) {
    int b = htmp->FindBin(in->GetBinCenter(i));

    double err = htmp->GetBinError(b);
    err = err * err + in->GetBinError(i) * in->GetBinError(i);

    htmp->AddBinContent(b, in->GetBinContent(i));
    htmp->SetBinError(b, sqrt(err));
  }

  return htmp;
}

bool CollieIOFile::getCollie(const TH1D* histo, CollieHistogram* out) {
  if (!check()) {
    printf("CollieIOFile::getCollie, Error: The file is not properly initialized!\n");
    return false;
  }

  if (histo == NULL || out == NULL) {
    printf("CollieIOFile::getCollie, Error: NULL histograms!\n");
    return false;
  }

  TH1D* htmp = 0;
  if (usingBinMap_)
    htmp = getBinMapROOT(histo, "temp1d");
  else
    htmp = (TH1D*)histo->Clone("temp1d");

  if (!checkROOTHisto(htmp, false)) {
    printf("CollieIOFile::getCollie, Histogram Error\n");
    return false;
  }

  out->Clear();
  for (int i = 0; i < htmp->GetNbinsX(); i++) {
    out->FillBin(i, htmp->GetBinContent(i + 1));
    out->SetBinError(i, htmp->GetBinError(i + 1));
  }

  delete htmp;
  htmp = NULL;

  return true;
}

bool CollieIOFile::getCollie2D(const TH2D* in, CollieHistogram2d* out) {
  if (!check(true)) {
    printf("CollieIOFile::getCollie, Error: The file is not properly initialized!\n");
    return false;
  }
  if (in == NULL || out == NULL) {
    printf("CollieIOFile::getCollie, Error: NULL histograms!\n");
    return false;
  }

  out->Clear();

  TH2D* htmp2 = (TH2D*)in->Clone("temp2d");
  if (!checkROOTHisto(htmp2, false)) {
    printf("CollieIOFile::getCollie2D, Histogram Error\n");
    return false;
  }

  for (int bx = 1; bx <= htmp2->GetNbinsX(); bx++)
    for (int by = 1; by <= htmp2->GetNbinsY(); by++) {
      out->Fill(bx - 1, by - 1, htmp2->GetBinContent(bx, by));
      out->SetBinError(bx - 1, by - 1, htmp2->GetBinError(bx, by));
    }

  delete htmp2;
  htmp2 = NULL;

  return true;
}

void CollieIOFile::fillDist(CollieHistogram* hist, CollieDistribution* dist) {
  double totaleff = hist->Integral();
  dist->setEfficiency(totaleff);
  if (totaleff > 0)
    totaleff = 1.0 / totaleff;
  string test("CollieHistogram2d");
  if (hist->getType() == test) {
    CollieHistogram2d* h2d = ((CollieHistogram2d*)hist);
    for (int i = 0; i < h2d->nxBins(); i++) {
      for (int j = 0; j < h2d->nyBins(); j++) {
        dist->setNormalizedBinValue(h2d->InBin(i, j) * totaleff, i, j);
      }
    }
  } else {
    for (int i = 0; i < hist->nBins(); i++) {
      dist->setNormalizedBinValue(hist->InBin(i) * totaleff, i);
    }
  }
  return;
}

//////////////////////////
//Interpolation Routines
//////////////////////////
void CollieIOFile::interpolate(TH1D* h1, double p1, double xsec1, TH1D* h2, double p2, double xsec2,
                               TH1D* hf, double pf, double xsecF) {

  if (h1 == NULL || h2 == NULL || hf == NULL) {
    printf("CollieIOFile::interpolate, Error: NULL histos!\n");
    return;
  }

  if (noviceFlag_) {
    printf(
        "CollieIOFile::interpolate, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing interpolation\n");
    return;
  }

  double offset = 0;
  for (int i = 1; i <= h1->GetNbinsX(); i++) {
    if (h1->GetBinContent(i) < offset)
      offset = h1->GetBinContent(i);
    if (h2->GetBinContent(i) < offset)
      offset = h2->GetBinContent(i);
  }

  if (offset < 0) {
    for (int i = 0; i <= h1->GetNbinsX() + 1; i++) {
      h1->AddBinContent(i, -1.0 * offset);
      h2->AddBinContent(i, -1.0 * offset);
    }
  }

  h1->Scale(1.0 / xsec1);
  h2->Scale(1.0 / xsec2);

  CollieHistogram* c1 = new CollieHistogram();
  CollieHistogram* c2 = new CollieHistogram();
  CollieHistogram* cf = new CollieHistogram();

  c1->Book("c1", h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
  c2->Book("c2", h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
  cf->Book("cf", h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());

  if (!getCollie(h1, c1)) {
    printf("CollieIOFile::interpolate, Error: Interpolation error histo 1\n");
    return;
  }
  if (!getCollie(h2, c2)) {
    printf("CollieIOFile::interpolate, Error: Interpolation error histo 2\n");
    return;
  }
  if (!getCollie(hf, cf)) {
    printf("CollieIOFile::interpolate, Error: Interpolation error histo 3\n");
    return;
  }

  cf->Interpolate(*c1, p1, *c2, p2, pf);

  if (!getROOT(cf, hf)) {
    printf("CollieIOFile::interpolate, Error: Interpolation error histo 4\n");
    return;
  }

  h1->Scale(xsec1);
  h2->Scale(xsec2);
  hf->Scale(xsecF);

  if (hf->Integral() < 0) {
    printf("CollieIOFile::interpolateR, Warning: Negative interpolation!\n");
    //hf->Scale(0);
  }

  if (offset < 0) {
    for (int i = 0; i <= h1->GetNbinsX() + 1; i++) {
      h1->AddBinContent(i, offset);
      h2->AddBinContent(i, offset);
      hf->AddBinContent(i, offset);
    }
  }

  delete c1;
  c1 = NULL;
  delete c2;
  c2 = NULL;
  delete cf;
  cf = NULL;

  return;
}

void CollieIOFile::interpolate(CollieHistogram* h1, double p1, double xsec1, CollieHistogram* h2,
                               double p2, double xsec2, CollieHistogram* hf, double pf,
                               double xsecF) {
  if (noviceFlag_) {
    printf(
        "CollieIOFile::interpolate, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing interpolation\n");
    return;
  }

  h1->Scale(1.0 / xsec1);
  h2->Scale(1.0 / xsec2);
  hf->Interpolate(*h1, p1, *h2, p2, pf);
  h1->Scale(xsec1);
  h2->Scale(xsec2);
  hf->Scale(xsecF);
  //  if(hf->Integral()<0) hf->Scale(0);
  return;
}

void CollieIOFile::interpolate2D(int axis, CollieHistogram2d* h1, double p1, double xsec1,
                                 CollieHistogram2d* h2, double p2, double xsec2,
                                 CollieHistogram2d* hf, double pf, double xsecF) {
  if (noviceFlag_) {
    printf(
        "CollieIOFile::interpolate2D, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing interpolation\n");
    return;
  }

  h1->Scale(1.0 / xsec1);
  h2->Scale(1.0 / xsec2);
  hf->Scale(0);
  //Scan across X-Axis to get profiles to be interpolated...
  CollieHistogram s1, s2, sf;
  if (axis == 0) { //x-axis is changing
    s1.Book("profile 1", h1->nxBins(), h1->LeftEdge(), h1->RightEdge());
    s2.Book("profile 2", h2->nxBins(), h2->LeftEdge(), h2->RightEdge());
    sf.Book("profile final", hf->nxBins(), hf->LeftEdge(), hf->RightEdge());
    for (int y = 0; y < h1->nyBins(); y++) {
      s1.Scale(0);
      s2.Scale(0);
      sf.Scale(0);
      for (int x = 0; x < h1->nxBins(); x++) {
        s1.Fill(x, h1->InBin(x, y));
        s2.Fill(x, h2->InBin(x, y));
      }
      sf.Interpolate(s1, p1, s2, p2, pf);
      for (int x = 0; x < h1->nxBins(); x++) {
        hf->Fill(x, y, sf.InBin(x));
      }
    }
  } else if (axis == 1) { //y-axis is changing
    s1.Book("profile 1", h1->nyBins(), h1->BottomEdge(), h1->TopEdge());
    s2.Book("profile 2", h2->nyBins(), h2->BottomEdge(), h2->TopEdge());
    sf.Book("profile final", hf->nyBins(), hf->BottomEdge(), hf->TopEdge());
    for (int x = 0; x < h1->nxBins(); x++) {
      s1.Scale(0);
      s2.Scale(0);
      sf.Scale(0);
      for (int y = 0; y < h1->nyBins(); y++) {
        s1.Fill(y, h1->InBin(x, y));
        s2.Fill(y, h2->InBin(x, y));
      }
      sf.Interpolate(s1, p1, s2, p2, pf);
      for (int y = 0; y < h1->nyBins(); y++) {
        hf->Fill(x, y, sf.InBin(y));
      }
    }
  }
  h1->Scale(xsec1);
  h2->Scale(xsec2);
  hf->Scale(xsecF);
  if (hf->Integral() < 0)
    hf->Scale(0);
  return;
}

void CollieIOFile::interpolateMassGrid(int dataFlag, int grid, int m1, int m1p, int m2, int m2p) {

  printf("Iterpolating doesn't work for now...:( \n");
  return;

  if (m2 == -1 && m2p == -1) {
    if (!check()) {
      printf(
          "CollieIOFile::interpolateMassGrid, Error: The file is not properly initialized(1)!\n");
      return;
    }
  } else {
    if (!check(true)) {
      printf(
          "CollieIOFile::interpolateMassGrid, Error: The file is not properly initialized(2)!\n");
      return;
    }
  }

  if (noviceFlag_) {
    printf(
        "CollieIOFile::interpolateMassGrid, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing interpolation\n");
    return;
  }

  if (verb_)
    printf(
        "CollieIOFile::interpolateMassGrid, Interpolating mass points with grid granularity %d\n",
        grid);
  if (sig_.size() <= 1)
    return;
  if (histBinsY_ > 0) {
    printf("CollieIOFile::interpolateMassGrid, cannot interpolate 2D histos...\n");
    return; ///can't interpolate on 2D hists for now...
  }

  list<int> massList;
  map<int, int>::iterator mIt, mIt2;
  if (m2 == -1 && m2p == -1) { ///a one-D interpolation...
    for (mIt = parList_.begin(); mIt != parList_.end(); mIt++)
      massList.push_back(mIt->first);
    massList.sort();
    interpolateMassList(dataFlag, grid, m1, m1p, massList);
  } else {
    printf("Two-D interpolation doesn't work!\n");

    //     //do the 2D grid...piecewise
    //     printf("Y grid...\n");
    //     //first fill in 1-step holes in y-grid...
    //     for(mIt = parList_.begin(); mIt != parList_.end(); mIt++){
    //       mIt->second.sort();
    //       oneStepPass(grid,m2,m2p,mIt->second,mIt->first);
    //     }
    //     printf("X grid...\n");
    //     //now fill in 1-step holes in x-grid...
    //     for(mIt = parList_.begin(); mIt != parList_.end(); mIt++){
    //       for(int a=m1; a<m1p; a+=grid){
    //   massList.clear();
    //   if(checkMassPoint(mIt->first,a)){
    //     massList.push_back(mIt->first);
    //     for(mIt2 = mIt; mIt2 != parList_.end(); mIt2++){
    //       if(mIt==mIt2) continue;
    //       if(checkMassPoint(mIt2->first,a)){
    //         massList.push_back(mIt2->first);
    //       }
    //     }
    //   }
    //   massList.sort();
    //   if(massList.size()>1) oneStepPass(grid,m1,m1p,massList, 0, a);
    //       }
    //     }
  }

  return;
}

void CollieIOFile::oneStepPass(int grid, int massLo, int massHi, list<int> masses, int massX,
                               int massY) {

  list<int> shortList;
  list<int>::iterator p1, p2;

  //   for(p1 = masses.begin(); p1!=masses.end(); p1++){
  //     for(p2 = p1; p2!=masses.end(); p2++){
  //       if(p1==p2) continue;
  //       shortList.clear();
  //       shortList.push_back(mIndex(*p1,massY)); shortList.push_back(mIndex(*p2,massY));
  //       if(massX>0){
  //   shortList.clear();
  //   shortList.push_back(mIndex(massX,*p1)); shortList.push_back(mIndex(massX,*p2));
  //       }
  //       if(abs(*p1-*p2)<=grid){ ///two adjacent points...thus get the two on either side...
  //   if((*p1-grid)<massLo || (*p2+grid)>massHi) continue;
  //   interpolateMassList(grid,*p1-grid,*p2+grid,shortList,massX,massY);
  //       }
  //       else if(abs(*p1-*p2)<=2.0*grid){ //two points separated by a grid spacing...fill it in...
  //   interpolateMassList(grid,*p1,*p2,shortList,massX,massY);
  //       }
  //     }
  //   }
}

void CollieIOFile::interpolateMassList(int dataFlag, int grid, int startMass, int stopMass,
                                       list<int> masses, int massX, int massY) {
  printf("CollieIOFile::interpolateMassList  Oops, I'm commented out!!!\n");
  return;
  /*
  if(masses.size()<2) {
    printf("CollieIOFile::interpolateMassList, Error: Not enough points to interpolate...\n");
    return;}
  ///masses should be mIndex'd
  if(massX!=0 && massY!=0) {
    printf("CollieIOFile::interpolateMassList, Error: cannot interpolate in 2D with both dimensions fixed...\n");
    return;  ///need to know which is fixed...can't be both
  }
  if(histBinsY_ > 0){
    printf("CollieIOFile::interpolateMassList, cannot interpolate 2D histos...\n");
    return; ///can't interpolate on 2D hists for now...
  }

  if(noviceFlag_){
    printf("CollieIOFile::interpolateMassList, Novice flag is set.  Are you sure you know what you're doing?\n==>Not performing interpolation\n");
    return;
  }

  char name[256];
  int first = masses.front(); int last = masses.back();
  list<int>::iterator p1, p2;
  int mX,mY,mZ;

  if(massY!=0){  ///fix y-mass, search in x-dir
    mLookUp(first,mX,mY,mZ);
    first = mX;
    mLookUp(last,mX,mY,mZ);
    last = mX;
  }
  else if(massX!=0){ //fix x-mass, search in y-dir
    mLookUp(first,mX,mY,mZ);
    first = mY;
    mLookUp(last,mX,mY,mZ);
    last = mY;
  }
  else {//if they're both zero, we have a 1D grid and search along it...
        //thus the indices are fine the way they are.
  }

  int mSearch=0; int mIgnore = 0;
  for(int imass = startMass; imass<=stopMass; imass+=grid){
    if(imass<first) {p1 = masses.begin(); p2=p1; p2++;}
    else if(imass>last) {p1 = masses.end(); p1--; p2=p1; p1--;}
    else{
      p1 = masses.begin();
      for(list<int>::iterator iter = masses.begin(); iter!=masses.end(); iter++){
  mSearch = *iter;
  if(massY!=0) mLookUp(*iter,mSearch,mIgnore,mZ);
  else if(massX!=0) mLookUp(*iter,mIgnore,mSearch,mZ);
  if(mSearch<imass) p1 = iter;
  p2 = p1; p2++;
      }
    }
    if(p2==masses.end()) continue;

    if(massX>0){
      if(checkMassPoint(massX,imass)) continue;
      sprintf(name,"Signal Final Variable %d, %d",massX,imass);
    }
    else if(massY>0){
      if(checkMassPoint(imass,massY)) continue;
      sprintf(name,"Signal Final Variable %d, %d",imass,massY);
    }
    else{
      if(checkMassPoint(imass)) continue;
      sprintf(name,"Signal Final Variable %d",imass);
    }

    CollieHistogram* ihist = new CollieHistogram();
    ihist->Book(name,(int)(histBinsX_/rebinX_),histMinX_,histMaxX_);

    double mH = imass*10.0;
    CollieHistogram* lo=NULL, *hi=NULL;
    int imassLo=0,imassHi=0,idxL=0, idxH=0;
    if(massX>0){
      mLookUp(*p1,mIgnore,mSearch,mZ);
      idxL = mIndex(mIgnore,mSearch);
      lo = sig_[mIndex(massX,mSearch)];
      imassLo = mSearch;
      mLookUp(*p2,mIgnore,mSearch,mZ);
      idxH = mIndex(mIgnore,mSearch);
      hi = sig_[mIndex(massX,mSearch)];
      imassHi = mSearch;
    }
    else if(massY>0){
      mLookUp(*p1,mSearch,mIgnore,mZ);
      idxL = mIndex(mIgnore,mSearch);
      lo = sig_[mIndex(mSearch,massY)];
      imassLo=mSearch;
      mLookUp(*p2,mSearch,mIgnore,mZ);
      idxH = mIndex(mIgnore,mSearch);
      hi = sig_[mIndex(mSearch,massY)];
      imassHi = mSearch;
    }
    else{ lo = sig_[*p1]; imassLo=*p1; idxL=imassLo; hi = sig_[*p2]; imassHi=*p2; idxH=imassHi;}

    interpolate(lo,imassLo*10.0,1,
    hi,imassHi*10.0,1,
    ihist, mH, 1);

    CollieMasspoint* point;
    if(massX>0) point = channel_->createMasspoint(massX,imass);
    else if(massY>0) point = channel_->createMasspoint(imass,massY);
    else point = channel_->createMasspoint(imass);

    point->Book((int)(histBinsX_/rebinX_), histMinX_, histMaxX_);
    point->BookData();
    //fill signal dist
    CollieDistribution* dist = point->getSignalDistMutable(0);
    fillDist(ihist,dist);  //Interpolated histogram, fill from histo
    dist->setModelXsec(1.0);

    if(massX>0) sig_[mIndex(massX,imass)] = ihist;
    else if(massY>0) sig_[mIndex(imass,massY)] = ihist;
    else sig_[imass] = ihist;

    map<string,CollieHistogram*> mapLo = bkgd_[idxL];
    map<string,CollieHistogram*> mapHi = bkgd_[idxH];
    map<string,CollieHistogram*> input;

    if(massX>0) sprintf(name,"All Bkgd Final Variable - %d,%d",massX,imass);
    if(massY>0) sprintf(name,"All Bkgd Final Variable - %d,%d",imass,massY);
    else sprintf(name,"All Bkgd Final Variable - %d",imass);
    CollieHistogram* allBkgd = new CollieHistogram();
    allBkgd->Book(name,(int)(histBinsX_*1.0/rebinX_),histMinX_,histMaxX_);

    for(map<string,CollieHistogram*>::iterator loI=mapLo.begin(); loI!=mapLo.end(); loI++){
      map<string,CollieHistogram*>::iterator hiI=mapHi.find(loI->first);
      assert(hiI!=mapHi.end());
      const char* t1 = loI->first.c_str();
      int idx=0;
      for(uint n=0; n<nBkgd_; n++){
  const char* t2 = channel_->getBackgroundName(n);
  if(strcmp(t1,t2)==0){idx=n; n=nBkgd_+10;}
      }
      const char* tmpName = channel_->getBackgroundName(idx);

      CollieHistogram* hist = new CollieHistogram();
      if(massX>0) sprintf(name,"%s Final Variable - %d,%d",tmpName,massX,imass);
      else if(massY>0) sprintf(name,"%s Final Variable - %d,%d",tmpName,imass,massY);
      else sprintf(name,"%s Final Variable - %d",tmpName,imass);
      hist->Book(name,(int)(histBinsX_*1.0/rebinX_),histMinX_,histMaxX_);

      interpolate(loI->second,imassLo*10.0,1,
      hiI->second,imassHi*10.0,1,
      hist, mH, 1);

      input[tmpName] = hist;
      CollieDistribution* dist = point->getBkgdDistMutable(idx);
      fillDist(hist,dist);//Interpolated, fill from dist...

      dist->setModelXsec(1.0);
      allBkgd->Add(*hist);
    }
    //9/18/08 WF comment next two lines.  Correct thing to do??
    //    bkgd_[mIndex(massX,massY)] = input;
    //    allBkgd_[mIndex(massX,massY)] = allBkgd;
    if(massX>0){ bkgd_[mIndex(massX,imass)] = input; allBkgd_[mIndex(massX,imass)] = allBkgd;}
    else if(massY>0){ bkgd_[mIndex(imass,massY)] = input; allBkgd_[mIndex(imass,massY)] = allBkgd;}
    else { bkgd_[imass] = input; allBkgd_[imass] = allBkgd;}

    const CollieEventList* evtList0 = 0;
    const CollieEventList* evtList2 = 0;

    ///Check to see if we have data event lists.  If so, linearly interpolate each event.  dataFlag==2
    // Else, interpolate the data histograms like the MC.                                  dataFlag==1
    // If not, arbitrarily choose the data point with the smallest chi2 to the bkgd.       dataFlag==0
    if(dataFlag==2){
      if(massX>0){
  evtList0 = getDataEventList(massX, imassLo);
  evtList2 = getDataEventList(massX, imassHi);
      }
      else if(massY>0){
  evtList0 = getDataEventList(imassLo,massY);
  evtList2 = getDataEventList(imassHi,massY);
      }
      else{
  evtList0 = getDataEventList(imassLo);
  evtList2 = getDataEventList(imassHi);
      }

      if(evtList0 && evtList2){

  int nDvars = evtList0->nDVars(); int nIvars = evtList0->nIVars();
  int nDims = evtList0->nDimensions();
  if((evtList2->nDVars()!=nDvars) ||  (evtList2->nIVars()!=nIvars) || (evtList2->nDimensions()!=nDims)){
    printf("CollieIOFile::interpolateMassList, Error: Data Event list parameter matching error\n");
    return;
  }

  point->ConfigEventList(nDvars,nIvars);
  CollieEventList* iList = point->getDataEventListMutable();

  for(int d=nDims; d<nDvars; d++){
    iList->setDVarName(d,evtList0->getDVarName(d));
  }
  for(int d=0; d<nIvars; d++){
    iList->setIVarName(d,evtList0->getIVarName(d));
  }

  double* dvars = new double[nDvars-nDims];
  int* ivars = new int[nIvars];

  bool validMatch = true;
  double x0=0, y0=0, x2=0, y2=0, ix=0, iy=0;

  for(int e0 = 0; e0<evtList0->getNEvents(); e0++){
    x0 = evtList0->getX(e0);
    if(nDims==2) y0 = evtList0->getY(e0);
    for(int d=nDims; d<nDvars; d++) dvars[d-nDims] = evtList0->getDVar(e0,d);
    for(int d=0; d<nIvars; d++) ivars[d] = evtList0->getIVar(e0,d);

    for(int e2 = 0; e2<evtList2->getNEvents(); e2++){
      validMatch = true;
      for(int d=nDims; d<nDvars; d++) if(evtList2->getDVar(e2,d) != dvars[d-nDims]) validMatch = false;
      for(int d=0; d<nIvars; d++) if(evtList2->getIVar(e2,d) != ivars[d]) validMatch = false;
      if(validMatch){
        x2 = evtList2->getX(e2);
        ix = x0 + (mH-imassLo*10.0)*(x2-x0)/(imassHi*10.0-imassLo*10.0);
        if(nDims==1) iList->addEvent(ix,dvars,ivars);
        else{
    y2 = evtList2->getY(e2);
    iy = y0 + (mH-imassLo*10.0)*(y2-y0)/(imassHi*10.0-imassLo*10.0);
    iList->addEvent(ix,iy,dvars,ivars);
        }
      }
    }
  }
      }
    }
    else{
      CollieHistogram* dataLo = data_[idxL];
      CollieHistogram* dataHi = data_[idxH];
      CollieHistogram* mdata= new CollieHistogram();
      if(massX>0) sprintf(name,"Data Final Variable - %d,%d",massX,imass);
      else if(massY>0) sprintf(name,"Data Final Variable - %d,%d",imass,massY);
      else sprintf(name,"Data Final Variable - %d",imass);
      mdata->Book(name,(int)(histBinsX_*1.0/rebinX_),histMinX_,histMaxX_);

      if(dataFlag==1){
  //linearly interpolate data histogram
  interpolate(dataLo,imassLo*10.0,1,
        dataHi,imassHi*10.0,1,
        mdata, mH, 1);
      }
      else{
  //Don't have evt lists.  Grab data with smallest chi2 to new bkgd
  double chi2_1 = dataLo->Chi2Test(allBkgd);
  double chi2_2 = dataHi->Chi2Test(allBkgd);

  if(chi2_1<chi2_2) mdata = dataLo;
  else mdata = dataHi;
      }

      //fill data dist
      dist = point->getDataDistMutable();
      fillDist(mdata,dist); //No parent TH*, fill from empty CollieHistogram

      if(massX>0) data_[mIndex(massX,imass)] = mdata;
      else if(massY>0) data_[mIndex(imass,massY)] = mdata;
      else data_[imass] = mdata;
    }

    //log mass point
    point->logPoint();
    if(massX>0){
      if(!logMassPoint(massX,imass)){ printf("CollieIOFile::interpolateMassList, Error logging mass point!\n"); return; }
      else printf("Created interpolated mass point %d, %d\n",massX,imass);
    }
    else if(massY>0){
      if(!logMassPoint(imass,massY)){ printf("CollieIOFile::interpolateMassList, Error logging mass point!\n"); return; }
      else printf("==>Created interpolated mass point %d, %d\n",imass,massY);
    }
    else {
      if(!logMassPoint(imass)){ printf("CollieIOFile::interpolateMassList, Error logging mass point!\n"); return; }
      else printf("==>Created interpolated mass point %d\n",imass);
    }

    if(evtList0 && evtList2){
      if(massX>0) finalizeDataEventList(massX,imass);
      else if(massY>0) finalizeDataEventList(imass,massY);
      else finalizeDataEventList(imass);
    }

    //If there are systematics available, interpolate those too!
    //But only if they are synchronized amongst the points.  Otherwise
    //the assumption of same structure isn't valid.

    dist = point->getSignalDistMutable(0);

    mLookUp(idxL,mX,mY,mZ);
    CollieMasspoint* pLo = channel_->getMasspoint(mX,mY,mZ);
    mLookUp(idxH,mX,mY,mZ);
    CollieMasspoint* pHi = channel_->getMasspoint(mX,mY,mZ);

    //Signal systematics...
    CollieDistribution* dLo = pLo->getSignalDistMutable(0);
    CollieDistribution* dHi = pHi->getSignalDistMutable(0);

    TH1D* hLo = 0;    TH1D* hHi = 0;
    for(int syst=0; syst<dLo->getNsystematics(); syst++){
      string sLo = dLo->getSystName(syst);
      string sHi = dHi->getSystName(syst);
      if(sLo == sHi){
  printf("Interpolating signal systematic %s\n",sLo.c_str());
  hLo = (TH1D*) dLo->getPositiveSyst(syst);
  hHi = (TH1D*) dHi->getPositiveSyst(syst);
  if(hLo==NULL || hHi==NULL) continue;

  if(massX>0) sprintf(name,"%s interpolated signal systematic, Pos - %d,%d",sLo.c_str(),massX,imass);
  else if(massY>0) sprintf(name,"%s interpolated signal systematic, Pos - %d,%d",sLo.c_str(),imass,massY);
  else sprintf(name,"%s interpolated signal systematic, Pos - %d",sLo.c_str(),imass);

  TH1D* hP = new TH1D(name,name,hLo->GetNbinsX(),hLo->GetXaxis()->GetXmin(),hLo->GetXaxis()->GetXmax());

  interpolate(hLo,imassLo*10.0,1,
        hHi,imassHi*10.0,1,
        hP, mH,1);

  hLo = (TH1D*) dLo->getNegativeSyst(syst);
  hHi = (TH1D*) dHi->getNegativeSyst(syst);
  if(hLo==NULL || hHi==NULL) continue;

  if(massX>0) sprintf(name,"%s interpolated signal systematic, Neg - %d,%d",sLo.c_str(),massX,imass);
  else if(massY>0) sprintf(name,"%s interpolated signal systematic, Neg - %d,%d",sLo.c_str(),imass,massY);
  else sprintf(name,"%s interpolated signal systematic, Neg - %d",sLo.c_str(),imass);

  TH1D* hN = new TH1D(name,name,hLo->GetNbinsX(),hLo->GetXaxis()->GetXmin(),hLo->GetXaxis()->GetXmax());

  interpolate(hLo,imassLo*10.0,1,
        hHi,imassHi*10.0,1,
        hN, mH,1);

  if(massX>0) createSigSystematic(sLo.c_str(),hP, hN,massX,imass);
  else if(massY>0) createSigSystematic(sLo.c_str(),hP, hN,imass,massY);
  else createSigSystematic(sLo.c_str(),hP, hN,imass);
      }
    }

    for(int nb = 0; nb<pLo->getNBkgdDists(); nb++){
      dLo = pLo->getBkgdDistMutable(nb);
      dHi = pHi->getBkgdDistMutable(nb);
      if(dLo==NULL || dHi==NULL) continue;

      for(int syst=0; syst<dLo->getNsystematics(); syst++){
  string sLo = dLo->getSystName(syst);
  string sHi = dHi->getSystName(syst);
  if(sLo == sHi){
    printf("Interpolating bkgd systematic %s\n",sLo.c_str());
    hLo = (TH1D*) dLo->getPositiveSyst(syst);
    hHi = (TH1D*) dHi->getPositiveSyst(syst);

    if(hLo==NULL || hHi==NULL) continue;

    if(massX>0) sprintf(name,"%s interpolated bkgd %d systematic, Pos - %d,%d",sLo.c_str(),nb,massX,imass);
    else if(massY>0) sprintf(name,"%s interpolated bkgd %d systematic, Pos - %d,%d",sLo.c_str(),nb, imass,massY);
    else sprintf(name,"%s interpolated bkgd %d systematic, Pos - %d",sLo.c_str(),nb,imass);

    TH1D* hP = new TH1D(name,name,hLo->GetNbinsX(),hLo->GetXaxis()->GetXmin(),hLo->GetXaxis()->GetXmax());

    interpolate(hLo,imassLo*10.0,1,
          hHi,imassHi*10.0,1,
          hP, mH,1);

    hLo = (TH1D*) dLo->getNegativeSyst(syst);
    hHi = (TH1D*) dHi->getNegativeSyst(syst);
    if(hLo==NULL || hHi==NULL) continue;

    if(massX>0) sprintf(name,"%s interpolated bkgd %d systematic, Neg - %d,%d",sLo.c_str(),nb,massX,imass);
    else if(massY>0) sprintf(name,"%s interpolated bkgd %d systematic, Neg - %d,%d",sLo.c_str(),nb,imass,massY);
    else sprintf(name,"%s interpolated bkgd %d systematic, Neg - %d",sLo.c_str(),nb,imass);

    TH1D* hN = new TH1D(name,name,hLo->GetNbinsX(),hLo->GetXaxis()->GetXmin(),hLo->GetXaxis()->GetXmax());

    interpolate(hLo,imassLo*10.0,1,
          hHi,imassHi*10.0,1,
          hN, mH,1);

    if(massX>0) createBkgdSystematic(nb,sLo.c_str(),hP, hN,massX,imass);
    else if(massY>0) createBkgdSystematic(nb,sLo.c_str(),hP, hN,imass,massY);
    else createBkgdSystematic(nb,sLo.c_str(),hP, hN,imass);
  }
      }
    }
  }
  */
  return;
}

//////////////////////////
//Sanity Check Methods...
//////////////////////////
bool CollieIOFile::check(bool twoD) {
  if (histMinX_ == UNDEF)
    return false;
  if (histMaxX_ == UNDEF)
    return false;
  if (histBinsX_ == 0)
    return false;
  if (twoD) {
    if (histMinY_ == UNDEF)
      return false;
    if (histMaxY_ == UNDEF)
      return false;
    if (histBinsY_ == 0)
      return false;
  }
  if (rebinX_ < 1)
    rebinX_ = 1;
  if (rebinY_ < 1)
    rebinY_ = 1;
  if (!init_)
    return false;
  return true;
}

bool CollieIOFile::checkMassPoint(int pX, int pY, int pZ) {
  map<int, int>::iterator iter = parList_.find(mIndex(pX, pY, pZ));
  if (iter == parList_.end())
    return false;
  return true;
}

///Systematics-creator helper-function
double CollieIOFile::getSystErr(double nv, double sv, double nerr, double serr) {
  if (serr == 0)
    serr = sv + 1;
  if (nerr == 0)
    nerr = nv + 1;
  if (nv == 0) {
    nv = 1e-6;
    nerr = 1;
  }

  if (sv == 0) {
    sv = 2e-6;
    serr = 1;
  }

  double ferr = serr * serr / (nv * nv);
  ferr += nerr * nerr * sv * sv / (nv * nv * nv * nv);
  double out = sqrt(ferr);
  if (sv == nv)
    out = 0.3;
  return out;
}

double CollieIOFile::totErr(CollieHistogram* hist, int i, int j) {
  double totErr = 0;
  for (int b = i; b <= j; b++) {
    totErr += hist->BinErr(b - 1) * hist->BinErr(b - 1);
  }
  return sqrt(totErr);
}

double CollieIOFile::totErr(TH1* hist, int i, int j) {
  double totErr = 0;
  for (int b = i; b <= j; b++) {
    totErr += hist->GetBinError(b) * hist->GetBinError(b);
  }
  return sqrt(totErr);
}

//Diagnostics for coherent effects in shape-dependent systematics
bool CollieIOFile::checkForShape(CollieHistogram* ref, TH1* pos, TH1* neg, string syst) {

  if (pos == NULL || ref == NULL || neg == NULL) {
    printf("CollieIOFile::checkForShape, NULL histograms\n");
    return true;
  }
  double pI = pos->Integral();
  double nI = neg->Integral();
  double rI = ref->Integral();

  if (rI <= 0)
    return true;

  bool twoD = false;
  CollieHistogram2d* ref2d = 0;
  if (pos->GetNbinsY() > 1) {
    twoD = true;
    ref2d = (CollieHistogram2d*)ref;
    if (pos->GetNbinsY() != ref2d->nyBins()) {
      printf("CollieIOFile::checkForShape, Binning mismatch  PX...\n");
      return true;
    }
    if (neg->GetNbinsY() != ref2d->nyBins()) {
      printf("CollieIOFile::checkForShape, Binning mismatch  NX...\n");
      return true;
    }
    if (pos->GetNbinsX() != ref2d->nxBins()) {
      printf("CollieIOFile::checkForShape, Binning mismatch  PY...\n");
      return true;
    }
    if (neg->GetNbinsX() != ref2d->nxBins()) {
      printf("CollieIOFile::checkForShape, Binning mismatch  NY...\n");
      return true;
    }
  } else {
    if (pos->GetNbinsX() != ref->nBins()) {
      printf("CollieIOFile::checkForShape, Binning mismatch  PX (%d, %d)...\n", pos->GetNbinsX(),
             ref->nBins());
      return true;
    }
    if (neg->GetNbinsX() != ref->nBins()) {
      printf("CollieIOFile::checkForShape, Binning mismatch  NX (%d, %d)...\n", neg->GetNbinsX(),
             ref->nBins());
      return true;
    }
  }

  double pDiff = (pI - rI) / rI;
  double nDiff = (nI - rI) / rI;

  ///Perform nodes/bins comparison
  //  Flatten any shape systematic that has more nodes than 40% of the possibilities.
  //  Protect regions of coherent shape deviation (10% of total bins, semi-arbitrary).
  //  Remove any overall normalization difference to protect real global shape deviations.
  if (ref->nBins() > 5) {
    double vLastP = 0;
    double vLastN = 0;
    int nNodesP = 0;
    int nNodesN = 0;
    int nConsecP = 0;
    int nConsecN = 0;
    int consecP = 0;
    int consecN = 0;
    if (!twoD) {
      for (int b = 0; b < ref->nBins(); b++) {
        double vThisP = 0;
        double vThisN = 0;
        if (ref->InBin(b) > 0) {
          vThisP = (pos->GetBinContent(b + 1) - ref->InBin(b)) / ref->InBin(b);
          vThisP -= pDiff;
          if ((vThisP < 0 && vLastP > 0) || (vThisP > 0 && vLastP < 0)) {
            nNodesP++;
            consecP = 0;
          } else
            consecP++;

          if (consecP > nConsecP)
            nConsecP = consecP;
          vLastP = vThisP;

          vThisN = (neg->GetBinContent(b + 1) - ref->InBin(b)) / ref->InBin(b);
          vThisN -= nDiff;
          if ((vThisN < 0 && vLastN > 0) || (vThisN > 0 && vLastN < 0)) {
            nNodesN++;
            consecN = 0;
          } else
            consecN++;

          if (consecN > nConsecN)
            nConsecN = consecN;
          vLastN = vThisN;
        }
      }
    } else {
      for (int bx = 0; bx < ref2d->nxBins(); bx++) {
        for (int by = 0; by < ref2d->nyBins(); by++) {
          double vThisP = 0;
          double vThisN = 0;
          if (ref2d->InBin(bx, by) > 0) {
            vThisP =
                (pos->GetBinContent(bx + 1, by + 1) - ref2d->InBin(bx, by)) / ref2d->InBin(bx, by);
            vThisP -= pDiff;
            if ((vThisP < 0 && vLastP > 0) || (vThisP > 0 && vLastP < 0)) {
              nNodesP++;
              consecP = 0;
            } else
              consecP++;

            if (consecP > nConsecP)
              nConsecP = consecP;
            vLastP = vThisP;

            vThisN =
                (neg->GetBinContent(bx + 1, by + 1) - ref2d->InBin(bx, by)) / ref2d->InBin(bx, by);
            vThisN -= nDiff;
            if ((vThisN < 0 && vLastN > 0) || (vThisN > 0 && vLastN < 0)) {
              nNodesN++;
              consecN = 0;
            } else
              consecN++;

            if (consecN > nConsecN)
              nConsecN = consecN;
            vLastN = vThisN;
          }
        }
      }
    }

    double nb = 0.40 * (pos->GetNbinsX() - 1);
    if (twoD)
      nb = nb * (pos->GetNbinsY() - 1);

    double cs = 0.10 * (pos->GetNbinsX());
    if (0.10 * (pos->GetNbinsY()) > cs)
      cs = 0.10 * (pos->GetNbinsY());

    if (nNodesP > nb && nConsecP < cs) {
      printf(
          "CollieIOFile::checkForShape, The positive fractional shape for %s has %d nodes in %d bins.\n          Only %d coherent bins.  Flattening...\n",
          syst.c_str(), nNodesP, pos->GetNbinsX(), nConsecP);
      return true;
    }
    if (nNodesN > nb && nConsecN < cs) {
      printf(
          "CollieIOFile::checkForShape, The negative fractional shape for %s has %d nodes in %d bins.\n          Only %d coherent bins.  Flattening...\n",
          syst.c_str(), nNodesN, neg->GetNbinsX(), nConsecN);
      return true;
    }
  }

  return false;
}

bool CollieIOFile::testShapeSystematics(CollieHistogram* nom, TH1* fluct, TH1* out, bool norm,
                                        string syst, double co) {

  if (nom == NULL || fluct == NULL || out == NULL) {
    printf("CollieIOFile::testShapeSystematics, NULL histograms! exiting...\n");
    return false;
  }

  int nbinsX = nom->nBins();
  if (fluct->GetNbinsX() != nbinsX || out->GetNbinsX() != nbinsX) {
    printf("CollieIOFile::testShapeSystematics, Difference in histogram binning! exiting...\n");
    return false;
  }

  //if desired, normalize integrals
  if (norm)
    fluct->Scale(nom->Integral() / fluct->Integral());

  //Clear output histo
  out->Scale(0.0);

  //Any syst unc less than 2% is preserved...should this be calculated in situ?  Can it?
  double systv = 1;
  double stat = 1;
  double nv = 1;
  double sv = 1;
  double nerr = 1;
  double serr = 1;
  double val = 1;
  bool found = 1;

  for (int i = 1; i <= nbinsX; i++) {

    nv = nom->InBin(i - 1);
    sv = fluct->GetBinContent(i);
    nerr = nom->BinErr(i - 1);
    serr = nerr;
    if (fluct->GetSumw2N() > 0)
      serr = fluct->GetBinError(i);

    //Fault protection circuitry
    if (nv == 0) {
      nv = 1e-6;
      nerr = nv;
    }
    if (sv == 0) {
      sv = 2e-6;
      serr = sv;
    }
    if (serr == 0)
      serr = sv + 1;
    if (nerr == 0)
      nerr = nv + 1;
    if (sv == nv)
      sv = nv + 1e-6;

    val = (sv - nv) / (nv + 1e-9);
    stat = getSystErr(nv, sv, nerr, serr);
    systv = fabs(val);

    if (getSystErr(nv, sv, nerr, serr) < 1e-6)
      stat *= 1.1;

    ///Is the stat uncertainty larger than the error we're trying assign??
    if (0.75 * stat > systv && fabs(val) > co) {
      found = false;
      for (int j = i + 1; j <= nbinsX && !found; j++) {
        nv = nom->Integral(i - 1, j - 1);
        sv = fluct->Integral(i, j);
        nerr = totErr(nom, i, j);
        serr = nerr;
        if (fluct->GetSumw2N() > 0)
          serr = totErr(fluct, i, j);
        //Fault protection circuitry
        if (nv == 0) {
          nv = 1e-6;
          nerr = nv;
        }
        if (sv == 0) {
          sv = 2e-6;
          serr = sv;
        }
        if (serr == 0)
          serr = sv + 1;
        if (nerr == 0)
          nerr = nv + 1;
        if (sv == nv)
          sv = nv + 1e-6;

        val = (sv - nv) / (nv + 1e-9);

        stat = getSystErr(nv, sv, nerr, serr);

        systv = fabs(val);

        if (0.75 * stat < systv || j == nbinsX || fabs(val) < co) {
          found = true;

          for (int b = i; b <= j; b++) {
            out->SetBinContent(b, val);
            out->SetBinError(b, getSystErr(nv, sv, nerr, serr) * sqrt(j - i));
          }
          i = j;
        }
      }
    } else {
      out->SetBinContent(i, val);
      out->SetBinError(i, getSystErr(nv, sv, nerr, serr));
    }
  }

  return true;
}

void CollieIOFile::EqualProbDifference(TH1* nom, TH1* fluct, TH1* out, string syst, bool norm,
                                       double prob) {

  if (nom == NULL || fluct == NULL || out == NULL) {
    printf("CollieIOFile::EqualProbDifference, NULL histograms! exiting...\n");
    return;
  }

  int nbinsX = nom->GetNbinsX();
  if (fluct->GetNbinsX() != nbinsX || out->GetNbinsX() != nbinsX) {
    printf("CollieIOFile::EqualProbDifference, Difference in histogram binning! exiting...\n");
    return;
  }

  //if desired, normalize integrals
  if (norm)
    fluct->Scale(nom->Integral() / fluct->Integral());

  double max = nom->GetMaximum();
  vector<int> binStart;
  vector<int> binEnd;

  //loop over bins to find regions with semi-equalized probability
  //depends on prob parameter (2.0 means everything within 2x of max)
  for (int i = 1; i <= nbinsX; i++) {
    if (nom->GetBinContent(i) < max / prob) {
      for (int j = i; j <= nbinsX; j++) {

        ///Tail catcher...
        if (nom->Integral(i, j) >= max / prob) {
          if (nom->Integral(j, nbinsX) < max / prob && nom->Integral(i, nbinsX) < max) {
            j = nom->GetNbinsX();
          }

          binStart.push_back(i);
          binEnd.push_back(j);
          i = j;

          break;
        }

        if (j == nbinsX) {
          binStart.push_back(i);
          binEnd.push_back(j);
          i = j;

          break;
        }
      }
    } else {
      binStart.push_back(i);
      binEnd.push_back(i);
    }
  }

  TH1D* rebinHist = new TH1D("rbh", "rbh", binStart.size(), 0, 1);

  //Formulate fractional per-bin fluctuations
  for (uint i = 0; i < binStart.size(); i++) {
    double val = fluct->Integral(binStart[i], binEnd[i]);
    double valNom = nom->Integral(binStart[i], binEnd[i]);
    double err = 1.0;

    if (valNom > 0) {
      double ap = 0, bp = 0;
      for (int j = binStart[i]; j <= binEnd[i]; j++) {
        ap += fluct->GetBinError(j) * fluct->GetBinError(j);
        bp += nom->GetBinError(j) * nom->GetBinError(j);
      }
      err = val / (valNom * valNom);
      err = err * err * bp;
      err += ap / (valNom * valNom);

      err = sqrt(err);
      val -= valNom;
      val /= valNom;
    } else {
      val = 0;
      if (val > 0)
        val = 1.0;
      err = 1.0;
    }

    if (std::isnan(val))
      printf("\nNAN: %f, %f\n", fluct->Integral(binStart[i] + 1, binEnd[i] + 1),
             nom->Integral(binStart[i], binEnd[i]));
    if (std::isinf(val))
      printf("\nINF: %f, %f\n", fluct->Integral(binStart[i] + 1, binEnd[i] + 1),
             nom->Integral(binStart[i], binEnd[i]));

    rebinHist->SetBinContent(i + 1, val);
    rebinHist->SetBinError(i + 1, err);
  }

  //Translate back to original binning
  out->Scale(0);
  for (uint i = 0; i < binStart.size(); i++) {
    for (int j = binStart[i]; j <= binEnd[i]; j++) {
      out->SetBinContent(j, rebinHist->GetBinContent(i + 1));
      out->SetBinError(j, rebinHist->GetBinError(i + 1));
    }
  }

  rebinHist->Delete();

  return;
}

void CollieIOFile::EqualProbDifference(CollieHistogram* nom, TH1* fluct, TH1* out, string syst,
                                       bool norm, double prob) {

  if (nom == NULL || fluct == NULL || out == NULL) {
    printf("CollieIOFile::EqualProbDifference, NULL histograms! exiting...\n");
    return;
  }

  int nbinsX = nom->nBins();
  if (fluct->GetNbinsX() != nbinsX || out->GetNbinsX() != nbinsX) {
    printf("CollieIOFile::EqualProbDifference, Difference in histogram binning! exiting...\n");
    return;
  }

  //if desired, normalize integrals
  if (norm)
    fluct->Scale(nom->Integral() / fluct->Integral());

  double max = nom->GetMaximum();
  vector<int> binStart;
  vector<int> binEnd;

  //loop over bins to find regions with semi-equalized probability
  //depends on prob parameter (2.0 means everything within 2x of max)
  for (int i = 0; i < nbinsX; i++) {
    if (nom->InBin(i) < max / prob) {
      for (int j = i; j < nbinsX; j++) {

        ///Tail catcher...
        if (nom->Integral(i, j) >= max / prob) {
          if (nom->Integral(j, nbinsX) < max / prob && nom->Integral(i, nbinsX) < max) {
            j = (nbinsX - 1);
          }

          binStart.push_back(i);
          binEnd.push_back(j);
          i = j;

          break;
        }

        if (j == (nbinsX - 1)) {
          binStart.push_back(i);
          binEnd.push_back(j);
          i = j;
          break;
        }
      }
    } else {
      binStart.push_back(i);
      binEnd.push_back(i);
    }
  }

  TH1D* rebinHist = new TH1D("rbh", "rbh", binStart.size(), 0, 1);

  //Formulate fractional per-bin fluctuations
  for (uint i = 0; i < binStart.size(); i++) {
    double val = fluct->Integral(binStart[i] + 1, binEnd[i] + 1);
    double valNom = nom->Integral(binStart[i], binEnd[i]);
    double err = 1.0;

    if (valNom > 0) {
      double ap = 0, bp = 0;
      for (int j = binStart[i]; j <= binEnd[i]; j++) {
        ap += fluct->GetBinError(j + 1) * fluct->GetBinError(j + 1);
        bp += nom->BinErr(j) * nom->BinErr(j);
      }
      if (ap == 0)
        ap = valNom * valNom;

      err = val / (valNom * valNom);
      err = err * err * bp;
      err += ap / (valNom * valNom);

      err = sqrt(err);
      val -= valNom;
      val /= valNom;
    } else {
      val = 0;
      if (val > 0)
        val = 1.0;
      err = 1.0;
    }

    if (std::isnan(val))
      printf("\nNAN: %f, %f\n", fluct->Integral(binStart[i] + 1, binEnd[i] + 1),
             nom->Integral(binStart[i], binEnd[i]));
    if (std::isinf(val))
      printf("\nINF: %f, %f\n", fluct->Integral(binStart[i] + 1, binEnd[i] + 1),
             nom->Integral(binStart[i], binEnd[i]));
    rebinHist->SetBinContent(i + 1, val);
    rebinHist->SetBinError(i + 1, err);
  }

  //Translate back to original binning
  out->Scale(0.0);
  for (uint i = 0; i < binStart.size(); i++) {
    for (int j = binStart[i]; j <= binEnd[i]; j++) {
      out->SetBinContent(j + 1, rebinHist->GetBinContent(i + 1));
      out->SetBinError(j + 1, rebinHist->GetBinError(i + 1));

      if (rebinHist->GetBinContent(i + 1) > 1.0) {
        printf(
            "**************************CollieIOFile::EqualProbDifference********************************\n");
        printf("Your systematics distribution indicate a change of %f percent in\n",
               100 * rebinHist->GetBinContent(i + 1));
        printf("bin %d of systematic: %s, %s\n", j, fluct->GetTitle(), fluct->GetName());
        out->SetBinContent(j + 1, 0.025);
        printf(
            "*******************************************************************************************\n");
      }
    }
  }

  rebinHist->Delete();

  return;
}

void CollieIOFile::checkRates(CollieHistogram* ref, TH1* pos, TH1* neg, string syst) {

  if (ref == NULL || pos == NULL || neg == NULL) {
    printf("CollieIOFile::checkRates, Error: NULL histos!\n");
    return;
  }

  double max = 0;
  double min = 0;
  double totP = 0;
  double totN = 0;
  for (int i = 0; i < ref->nBins(); i++) {
    if (pos->GetBinContent(i + 1) > max) {
      max = pos->GetBinContent(i + 1);
    }
    if (pos->GetBinContent(i + 1) < min) {
      min = pos->GetBinContent(i + 1);
    }
    if (neg->GetBinContent(i + 1) > max) {
      max = neg->GetBinContent(i + 1);
    }
    if (neg->GetBinContent(i + 1) < min) {
      min = neg->GetBinContent(i + 1);
    }

    totP += (1.0 + pos->GetBinContent(i + 1)) * ref->InBin(i);
    totN += (1.0 + neg->GetBinContent(i + 1)) * ref->InBin(i);
  }

  double in = ref->Integral();
  double ip = totP;
  double im = totN;

  if ((fabs(in - ip) > 0.50 * in) || (fabs(in - im) > 0.50 * in)) {
    printf("**************************CollieIOFile::checkRates********************************\n");
    printf(
        "These 1sigma histograms differ by %.2f%% (+1s), %.2f%% (-1s) in total rate (nom: %.2f, +1s: %.2f, -1s: %.2f)\n",
        100 * fabs(in - ip) / in, 100 * fabs(in - im) / in, in, ip, im);
    printf("from the nominal.  ID: %s, %s, %s.\n", syst.c_str(), pos->GetTitle(), pos->GetName());
    printf("Are you sure you used the right shape systematic method?  Please read\n");
    printf("collie/io/include/CollieIOFile.hh for appropriate instructions.\n");
    printf("\nIf you believe you have a systematic changing rates more than 30%% please make\n");
    printf(
        "sure you are using an appropriate systematics PDF (ie, choose Gaussian or log-Normal)\n");
    printf("**********************************************************************************\n");
  }

  if ((fabs(max) > 0.50) || (fabs(min) > 0.50)) {
    printf("**************************CollieIOFile::checkRates********************************\n");
    printf(
        "These 1sigma histograms contain systematic differences of up to %f%% (%f%%) for the maximum (minimum) fluctuations\n",
        100 * max, 100 * min);
    printf("from the nominal.  ID: %s, %s, %s.\n", syst.c_str(), pos->GetTitle(), pos->GetName());
    printf("Are you sure you used the right shape systematic method?  Please read\n");
    printf("collie/io/include/CollieIOFile.hh for appropriate instructions.\n");
    printf("\nIf you believe you have a systematic changing rates more than 30%% please make\n");
    printf(
        "sure you are using an appropriate systematics PDF (ie, choose Gaussian or log-Normal)\n");
    printf("**********************************************************************************\n");
  }

  return;
}

bool CollieIOFile::checkROOTHisto(TH1* histo, bool norm) {
  if (histo == NULL) {
    printf("CollieIOFile::checkROOTHisto, Error: NULL histos!\n");
    return false;
  }

  if (histo->GetNbinsY() == 1) {
    if (histo->GetNbinsX() != (int)(histBinsX_ * 1.0 / rebinX_)) {
      if (histo->GetNbinsX() == histBinsX_) {
        histo->Rebin(rebinX_);
        if (norm)
          histo->Scale(1.0 / rebinX_);
      } else {
        printf("CollieIOFile::checkROOTHisto, Error: Binning mismatch: %d vs %d\n",
               histo->GetNbinsX(), (int)histBinsX_ / rebinX_);
        return false;
      }
    }
  } else {
    if (histo->GetNbinsX() != (int)(histBinsX_ * 1.0 / rebinX_)) {
      if (histo->GetNbinsX() == histBinsX_) {
        ((TH2*)histo)->RebinX(rebinX_);
        if (norm)
          histo->Scale(1.0 / rebinX_);
      } else {
        printf("CollieIOFile::checkROOTHisto, Error: Binning mismatch: %d vs %d\n",
               histo->GetNbinsX(), (int)histBinsX_ / rebinX_);
        return false;
      }
    }
    if (histo->GetNbinsY() != (int)(histBinsY_ * 1.0 / rebinY_)) {
      if (histo->GetNbinsY() == histBinsY_) {
        ((TH2*)histo)->RebinY(rebinY_);
        if (norm)
          histo->Scale(1.0 / rebinY_);
      } else {
        printf("CollieIOFile::checkROOTHisto, Error: Binning mismatch: %d vs %d\n",
               histo->GetNbinsY(), (int)histBinsY_ / rebinY_);
        return false;
      }
    }
  }

  return true;
}

void CollieIOFile::checkBins(TH1D* data, vector<TH1D*>& sig, vector<TH1D*>& bkgd) {

  if (!data) {
    printf("CollieIOFile::checkBins, error: NULL data histogram\n");
    return;
  }

  bool of = 0;
  bool uf = 0;
  if (data->GetBinContent(0) > 0)
    uf = 1;
  if (data->GetBinContent(data->GetNbinsX() + 1) > 0)
    of = 1;
  for (uint b = 0; b < bkgd.size(); b++) {
    if (!bkgd[b]) {
      printf("CollieIOFile::checkBins, error: NULL bkgd histogram (number %d)\n", b);
      return;
    }
    if (bkgd[b]->GetBinContent(0) > 0)
      uf = 1;
    if (bkgd[b]->GetBinContent(data->GetNbinsX() + 1) > 0)
      of = 1;

    if (of) {
      printf("\nCollieIOFile::checkBins, warning: Detected non-zero overflow contents (%s).\n",
             bkgd[b]->GetTitle());
      printf("CollieIOFile::checkBins, warning: Collie ignores over and underflow bins.\n");
    }
    if (uf) {
      printf("\nCollieIOFile::checkBins, warning: Detected non-zero underflow contents (%s).\n",
             bkgd[b]->GetTitle());
      printf("CollieIOFile::checkBins, warning: Collie ignores over and underflow bins.\n");
    }
    of = 0;
    uf = 0;
  }
  for (uint s = 0; s < sig.size(); s++) {
    if (!sig[s]) {
      printf("CollieIOFile::checkBins, error: NULL sig histogram (number %d)\n", s);
      return;
    }
    if (sig[s]->GetBinContent(0) > 0)
      uf = 1;
    if (sig[s]->GetBinContent(data->GetNbinsX() + 1) > 0)
      of = 1;
    if (of) {
      printf("\nCollieIOFile::checkBins, warning: Detected non-zero overflow contents (%s).\n",
             sig[s]->GetTitle());
      printf("CollieIOFile::checkBins, warning: Collie ignores over and underflow bins.\n");
    }
    if (uf) {
      printf("\nCollieIOFile::checkBins, warning: Detected non-zero underflow contents (%s).\n",
             sig[s]->GetTitle());
      printf("CollieIOFile::checkBins, warning: Collie ignores over and underflow bins.\n");
    }
    of = 0;
    uf = 0;
  }

  int nprobBins = 0;
  for (int ibin = 1; ibin <= data->GetNbinsX(); ibin++) {
    double sumBkgd = 0;
    double sumSig = 0;
    for (uint b = 0; b < bkgd.size(); b++) {
      sumBkgd += bkgd[b]->GetBinContent(ibin);
    }
    for (uint s = 0; s < sig.size(); s++) {
      sumSig += sig[s]->GetBinContent(ibin);
    }

    if (sumSig < 0 || sumBkgd < 0 || data->GetBinContent(ibin) < 0) {
      printf(
          "CollieIOFile::checkBins, error: Negative bin contents (bin %d): data=%.0f, bkgd=%.1f, sig=%.1f\n",
          ibin, data->GetBinContent(ibin), sumSig, sumBkgd);
      return;
    }

    if (sumBkgd == 0 && sumSig > 0) {
      nprobBins++;

      //Migrate to neighboring bin with highest s/b
      double sbPlus = 0;
      double sbMinus = 0;
      double ssPlus = 0;
      double ssMinus = 0;
      //start by looking for non-zero neighbors
      int distance = 1;
      int pBin = ibin + distance;
      int mBin = ibin - distance;
      while (sbPlus == 0 && sbMinus == 0 && pBin < data->GetNbinsX() && mBin > 0) {
        pBin = ibin + distance;
        mBin = ibin - distance;
        for (uint b = 0; b < bkgd.size(); b++) {
          sbPlus += bkgd[b]->GetBinContent(pBin);
          sbMinus += bkgd[b]->GetBinContent(mBin);
        }
        for (uint s = 0; s < sig.size(); s++) {
          ssPlus += sig[s]->GetBinContent(pBin);
          ssMinus += sig[s]->GetBinContent(mBin);
        }

        if (sbPlus > 0) {
          sbPlus = 1.0 / sbPlus;
          sbPlus *= ssPlus;
        }
        if (sbMinus > 0) {
          sbMinus = 1.0 / sbMinus;
          sbMinus *= ssMinus;
        }

        //if we didn't find it, take one more step away.
        distance++;
      }

      int theBin = mBin;
      if (sbPlus > sbMinus) {
        theBin = pBin;
      }

      data->AddBinContent(theBin, data->GetBinContent(ibin));
      data->SetBinContent(ibin, 0);
      for (uint s = 0; s < sig.size(); s++) {
        sig[s]->AddBinContent(theBin, sig[s]->GetBinContent(ibin));
        sig[s]->SetBinContent(ibin, 0);
      }
    }
  }

  if (nprobBins > 0)
    printf("CollieIOFile, Fixed %d bin(s) with 0 bkgd and non-zero signal\n", nprobBins);

  return;
}

void CollieIOFile::checkBins2D(TH2D* data, vector<TH2D*>& sig, vector<TH2D*>& bkgd) {

  int nprobBins = 0;
  for (int by = 1; by <= data->GetNbinsY(); by++) {
    for (int bx = 1; bx <= data->GetNbinsX(); bx++) {
      double sumBkgd = 0;
      for (uint b = 0; b < bkgd.size(); b++) {
        sumBkgd += bkgd[b]->GetBinContent(bx, by);
      }
      double sumSig = 0;
      for (uint s = 0; s < sig.size(); s++) {
        sumSig += sig[s]->GetBinContent(bx, by);
      }
      if (sumSig < 0 || sumBkgd < 0 || data->GetBinContent(bx, by) < 0) {
        printf(
            "CollieIOFile::checkBins, error: Negative bin contents (bin %d/%d): data=%.0f, bkgd=%.1f, sig=%.1f\n",
            bx, by, data->GetBinContent(bx, by), sumSig, sumBkgd);
        return;
      }
      if (sumBkgd == 0 && sumSig > 0) {
        nprobBins++;

        //Migrate to neighboring bin with highest s/b
        double sbPlus = 0;
        double sbMinus = 0;
        double ssPlus = 0;
        double ssMinus = 0;
        //start by looking for non-zero neighbor
        int distance = 1;
        int pBinX = bx;
        int mBinX = bx;
        int pBinY = by;
        int mBinY = by;
        bool doX = true;
        bool doY = false;
        while (sbPlus == 0 && sbMinus == 0) {
          if (doX) {
            pBinX = bx + distance;
            mBinX = bx - distance;
            pBinY = by;
            mBinY = by;
          } else {
            pBinX = bx;
            mBinX = bx;
            pBinY = by + distance;
            mBinY = by - distance;
          }
          for (uint b = 0; b < bkgd.size(); b++) {
            sbPlus += bkgd[b]->GetBinContent(pBinX, pBinY);
            sbMinus += bkgd[b]->GetBinContent(mBinX, mBinY);
          }
          for (uint s = 0; s < sig.size(); s++) {
            ssPlus += sig[s]->GetBinContent(pBinX, pBinY);
            ssMinus += sig[s]->GetBinContent(mBinX, mBinY);
          }

          if (sbPlus > 0) {
            sbPlus = 1.0 / sbPlus;
            sbPlus *= ssPlus;
          }
          if (sbMinus > 0) {
            sbMinus = 1.0 / sbMinus;
            sbMinus *= ssMinus;
          }

          //if we didn't find it, take one more step away.
          if (doX && !doY) {
            doX = false;
            doY = true;
          } else {
            doY = false;
            doX = true;
            distance++;
          }
        }

        int theXbin = mBinX;
        int theYbin = mBinY;
        if (sbPlus > sbMinus) {
          theXbin = pBinX;
          theYbin = pBinY;
        }

        data->SetBinContent(theXbin, theYbin,
                            data->GetBinContent(bx, by) + data->GetBinContent(theXbin, theYbin));
        data->SetBinContent(bx, by, 0);

        for (uint s = 0; s < sig.size(); s++) {
          sig[s]->SetBinContent(theXbin, theYbin,
                                sig[s]->GetBinContent(bx, by) +
                                    sig[s]->GetBinContent(theXbin, theYbin));
          sig[s]->SetBinContent(bx, by, 0);
        }
      }
    }
  }

  if (nprobBins > 0)
    printf("CollieIOFile, Fixed %d bin(s) with 0 bkgd and non-zero signal\n", nprobBins);
  return;
}

void CollieIOFile::applyCuts(TH1* hist) {
  if (cutLowX_ == UNDEF || cutHighX_ == UNDEF)
    return;
  //  printf("apply cuts: %f, %f\n",cutLowX_,cutHighX_);
  if (hist->GetNbinsY() > 1) { //twoD!
    if (cutHighY_ == UNDEF || cutLowY_ == UNDEF)
      return;
    TAxis* aX = hist->GetXaxis();
    TAxis* aY = hist->GetYaxis();
    for (int x = 1; x <= hist->GetNbinsX(); x++) {
      for (int y = 1; y <= hist->GetNbinsX(); y++) {
        double cX = aX->GetBinCenter(x);
        double cY = aY->GetBinCenter(y);
        if (cX < cutLowX_ || cX > cutHighX_) {
          hist->SetBinContent(x, y, 0);
          hist->SetBinError(x, y, 0);
        }
        if (cY < cutLowY_ || cY > cutHighY_) {
          hist->SetBinContent(x, y, 0);
          hist->SetBinError(x, y, 0);
        }
      }
    }
  } else {
    for (int x = 1; x <= hist->GetNbinsX(); x++) {
      if (hist->GetBinCenter(x) < cutLowX_ || hist->GetBinCenter(x) > cutHighX_) {
        hist->SetBinContent(x, 0);
        hist->SetBinError(x, 0);
      }
    }
  }
}

void CollieIOFile::applyCuts(CollieHistogram* hist) {
  if (cutLowX_ == UNDEF || cutHighX_ == UNDEF)
    return;

  //  printf("apply cuts: %f, %f\n",cutLowX_,cutHighX_);
  //  if(hist->GetBinCenter(x)<cutLowX_ || hist->GetBinCenter(x)>cutHighX_){

  string test("CollieHistogram2d");
  if (hist->getType() == test) {
    CollieHistogram2d* h2d = ((CollieHistogram2d*)hist);
    if (cutHighY_ == UNDEF || cutLowY_ == UNDEF)
      return;
    for (int x = 0; x < h2d->nxBins(); x++) {
      for (int y = 0; y < h2d->nyBins(); y++) {
        double cX = h2d->GetBinCenterX(x);
        double cY = h2d->GetBinCenterY(y);
        if (cX < cutLowX_ || cX > cutHighX_) {
          h2d->SetBinContent(x, y, 0);
          h2d->SetBinError(x, y, 0);
        }
        if (cY < cutLowY_ || cY > cutHighY_) {
          h2d->SetBinContent(x, y, 0);
          h2d->SetBinError(x, y, 0);
        }
      }
    }
  } else {
    for (int x = 0; x < hist->nBins(); x++) {
      if (hist->GetBinCenter(x) < cutLowX_ || hist->GetBinCenter(x) > cutHighX_) {
        hist->SetBinContent(x, 0);
        hist->SetBinError(x, 0);
      }
    }
  }
  return;
}

int CollieIOFile::generateBinMap(TH1D* bkgd, TH1D* sig, string style, double highS, double midS,
                                 double minSB, double minB, double minStatSB, double minStatB) {

  if (bkgd == NULL || sig == NULL) {
    printf("CollieIOFile::generateBinMap, NULL Input Histograms\n");
    usingBinMap_ = false;
    return -1;
  }

  if (histBinsY_ > 0) {
    printf("CollieIOFile::generateBinMap, Cannot yet rebin in 2D\n");
    usingBinMap_ = false;
    return -1;
  }

  if (bkgd->GetNbinsX() == histBinsX_) {
    printf(
        "CollieIOFile::generateBinMap, Bins in (%d) equals bins requested (%d), no rebinning necessary.\n",
        bkgd->GetNbinsX(), histBinsX_);
    usingBinMap_ = false;
    return -1;
  }

  if (style == "MVA")
    return generateMVABinMap(bkgd, sig, highS, midS, minSB, minB, minStatSB, minStatB);
  if (style == "InvMass")
    return generateInvMassBinMap(bkgd, sig, highS, midS, minSB, minB, minStatSB, minStatB);
  else {
    printf("CollieIOFile::generateBinMap, Bin style %s is not recognized\n", style.c_str());
    return -1;
  }
}

int CollieIOFile::generateInvMassBinMap(TH1D* bkgd, TH1D* sig, double highS, double midS,
                                        double minSB, double minB, double minStatSB,
                                        double minStatB) {
  printf("Still Working on the InvMass bin mapping...sorry.\n");
  return -1;
}

int CollieIOFile::generateMVABinMap(TH1D* bkgd, TH1D* sig, double highS, double midS, double minSB,
                                    double minB, double minStatSB, double minStatB) {

  binMap_.clear();
  vector<double> tmpEdges;
  rebinX_ = 1;
  usingBinMap_ = true;

  int nbins = 1;
  int sbin = bkgd->GetNbinsX();
  int ebin = sbin;

  int hBins = (int)(histBinsX_ * highS);
  int mBins = (int)(histBinsX_ * midS);
  int lBins = histBinsX_ - mBins - hBins;

  int nH = 0;
  int nM = 0;
  int nL = 0;
  double minMed = 0;
  double minLo = 0;

  double oVal = 0;
  double b2 = 0;
  double slope1 = -1;
  double slope2 = -1;

  bool lowQuad = true;

  tmpEdges.push_back(histMaxX_);
  for (int i = bkgd->GetNbinsX(); i > 0; i--) {
    double sumSB = 0;
    double sumB = 0;
    double sumerrSB = 0;
    double sumerrB = 0;

    for (int b = ebin; b <= sbin; b++) {
      sumB += bkgd->GetBinContent(b);
      sumSB = sumB + sig->GetBinContent(b);

      sumerrB += bkgd->GetBinError(b) * bkgd->GetBinError(b);
      sumerrSB +=
          bkgd->GetBinError(b) * bkgd->GetBinError(b) + sig->GetBinError(b) * sig->GetBinError(b);
    }
    sumerrB = sqrt(sumerrB) / (sumB + 1e-10);
    sumerrSB = sqrt(sumerrSB) / (sumSB + 1e-10);

    if (nH < hBins) {
      if (sumSB > minSB && sumB > minB && sumerrSB < minStatSB && sumerrB < minStatB) {
        binMap_[i] = nbins;
        tmpEdges.push_back(bkgd->GetBinLowEdge(ebin));
        printf("Adding high bins: %f\n", bkgd->GetBinLowEdge(ebin));
        nH++;
        nbins++;
        sbin = i - 1;
        ebin = i - 1;
      } else {
        binMap_[i] = nbins;
        ebin--;
      }
    } else if (nM < mBins) {
      if (slope1 == -1) {
        slope1 = getSlope(bkgd->Integral(1, ebin) + sig->Integral(1, ebin), minSB, mBins + lBins);
        if (slope1 > minSB)
          slope1 = minSB;
      }
      minMed = minSB * 1.5 + (nbins - nH - 1) * slope1;
      if (sumSB > minMed && sumB > minB && sumerrSB < minStatSB && sumerrB < minStatB) {
        binMap_[i] = nbins;
        tmpEdges.push_back(bkgd->GetBinLowEdge(ebin));
        printf("Adding mid bins: %f\n", bkgd->GetBinLowEdge(ebin));
        nM++;
        nbins++;
        sbin = i - 1;
        ebin = i - 1;
      } else {
        binMap_[i] = nbins;
        ebin--;
      }
    } else {
      if (slope2 == -1) {
        oVal = bkgd->Integral(1, ebin) + sig->Integral(1, ebin);

        slope2 = getSlope(oVal, minMed, lBins);
        if (slope2 < 3.0 * minSB)
          lowQuad = false;

        if (lowQuad) {
          for (int s = 0; s <= (histBinsX_ - nM - nH); s++)
            b2 += s * s;
          double sub = lBins * minMed * 1.5;
          if (sub < oVal)
            oVal -= sub;
          else
            lowQuad = false;
        }
      }

      int ibin = nbins - nH - nM - 1;
      if (lowQuad)
        minLo = minMed * 1.5 + ibin * ibin * oVal / b2;
      else
        minLo = minMed * 1.5 + ibin * slope2;

      if (sumSB > minLo && sumB > minB && nbins < histBinsX_) {
        binMap_[i] = nbins;
        tmpEdges.push_back(bkgd->GetBinLowEdge(ebin));
        printf("Adding low bins: %f\n", bkgd->GetBinLowEdge(ebin));
        nL++;
        nbins++;
        sbin = i - 1;
        ebin = i - 1;
      } else {
        if (i > 1 || (i == 1 && sbin != ebin)) {
          binMap_[i] = nbins;
          ebin--;
        } else {
          binMap_[i] = nbins - 1;
        }
      }
    }
  }
  if (nbins == histBinsX_) {
    tmpEdges.push_back(histMinX_);
    printf("Done adding bins: %f\n", histMinX_);
    nL++;
  } else if (nbins < histBinsX_) {
    double last = tmpEdges[nbins - 1];
    for (int b = (histBinsX_ - nbins); b >= 0; b--) {
      tmpEdges.push_back(histMinX_ + b * (last - histMinX_) / (histBinsX_ - nbins + 1));
      printf("Adding empty bins: %f\n",
             histMinX_ + b * (last - histMinX_) / (histBinsX_ - nbins + 1));
      nL++;
    }
  } else
    printf("Too many bins??\n");

  binEdges_.clear();
  printf("\nGenerated %zu edges:\n", tmpEdges.size());
  for (uint i = 1; i <= tmpEdges.size(); i++) {
    binEdges_.push_back(tmpEdges[tmpEdges.size() - i]);
    printf("Bin Low Edge %d: %f\n", i, tmpEdges[tmpEdges.size() - i]);
  }

  printf("//*************************************************//\n");
  printf("CollieIOFile::generateBinMap Report:\n");
  printf("N original bins: %d\n", bkgd->GetNbinsX());
  printf("N final bins   : %d\n", nbins);
  printf("NbinsHi: %d, NbinsMid: %d, NbinsLo: %d\n", nH, nM, nL);
  printf("//*************************************************//\n");

  return nbins;
}

TH1D* CollieIOFile::h1FMT(TH1D* inhist) {
  if (rebinX_ <= 1)
    return inhist;
  if (inhist == NULL)
    return NULL;

  inhist = (TH1D*)(inhist->Rebin(rebinX_));
  inhist->Scale(1.0 / rebinX_);
  return inhist;
}

void CollieIOFile::setSystPriors(TH1D* inHist) {
  if (inHist == NULL)
    return;

  systPriors_ = inHist;
  usingSystPriors_ = true;
}
