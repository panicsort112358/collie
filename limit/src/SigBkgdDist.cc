#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <iomanip>
#include "SigBkgdDist.hh"

using namespace std;

///Standard constructor
SigBkgdDist::SigBkgdDist(int var1, int var2, int var3)
    : m_delBinsLookup(),
      m_signalExclusionSumsReady(false),
      m_bkgdExclusionSumsReady(false) {
  setup(string(""), 0, 0, 1, var1, var2, var3);
  m_signal = NULL;
  m_bkgd = NULL;

  m_signalParent = NULL;
  m_bkgdParent = NULL;
  m_signalErr = NULL;
  m_bkgdErr = NULL;
  m_data = NULL;
  m_rebinBins = NULL;
  gl_systRand.clear(); // USED TO BE  gl_systrand = NULL;
  gl_systFloat = NULL;
  gl_systFit = NULL;
  gl_bkgdFit = NULL;
  gl_sigFit = NULL;
  gl_systCenter = NULL;
  gl_exclChan = "";
  gl_fitAllSources = true;

  m_fitValueHist = NULL;
  m_errMatHist = NULL;

  m_sigScale.clear();
  m_sigScalesVary = false;
  gl_cSigScale = 1.0;
  m_bkgdScale.clear();
  m_bkgdScale.push_back(1.0);
  m_dataScale = 1.00;
  n_bins = 0;
  n_trueBins = 0;
  m_delBins.clear();
  m_deleteArrays = true;
  m_rebinned = false;
  m_sigFluct = true;
  m_varied = false;
  m_init = false;
  m_adjSyst = false;
  m_useStat = false;
  m_condensed = false;
  m_verbose = false;

  m_sbhypo = false;
  m_bohypo = false;

  m_systScale = 1;
  m_systOffset = 0;

  n_cacheData = 0;

  fillArrays(true, false);
}

///copy constructor
SigBkgdDist::SigBkgdDist(const SigBkgdDist& sbd)
    : m_delBinsLookup(sbd.m_delBinsLookup),
      m_signalExclusionSumsReady(sbd.m_signalExclusionSumsReady),
      m_bkgdExclusionSumsReady(sbd.m_bkgdExclusionSumsReady) {
  m_var1 = sbd.m_var1;
  m_var2 = sbd.m_var2;
  m_var3 = sbd.m_var3;
  m_min = sbd.m_min;
  m_max = sbd.m_max;
  n_bins = sbd.n_bins;
  n_trueBins = 0;
  m_delBins.clear();
  m_chan = sbd.m_chan;

  m_sigScale = sbd.m_sigScale;
  m_sigScalesVary = sbd.m_sigScalesVary;
  gl_cSigScale = sbd.gl_cSigScale;

  m_bkgdScale = sbd.m_bkgdScale;
  m_dataScale = sbd.m_dataScale;
  m_useStat = sbd.m_useStat;
  n_cacheData = 0;
  m_fitValueHist = NULL;
  m_errMatHist = NULL;

  m_sbhypo = sbd.m_sbhypo;
  m_bohypo = sbd.m_bohypo;

  const map<string, map<string, CollieDistribution*>>& tmpS = sbd.m_SigDist;
  for (map<string, map<string, CollieDistribution*>>::const_iterator s = tmpS.begin();
       s != tmpS.end(); ++s) { //loop over channels
    const map<string, CollieDistribution*>& tmpSS = s->second;
    for (map<string, CollieDistribution*>::const_iterator d = tmpSS.begin(); d != tmpSS.end();
         ++d) { //loop over signals
      addSigDist(s->first, d->first, d->second);
    }
  }

  const map<string, vector<CollieDistribution*>>& tmpB = sbd.m_BkgdDist;
  for (map<string, vector<CollieDistribution*>>::const_iterator b = tmpB.begin(); b != tmpB.end();
       ++b)
    for (vector<CollieDistribution*>::const_iterator d = b->second.begin(); d != b->second.end();
         ++d)
      addBkgdDist(b->first, *d);

  const map<string, CollieDistribution*>& tmpD = sbd.m_DataDist;
  for (map<string, CollieDistribution*>::const_iterator d = tmpD.begin(); d != tmpD.end(); ++d)
    addDataDist(d->first, d->second);

  if (n_bins > 0) {
    m_signal = new double[n_bins];
    m_bkgd = new double[n_bins];

    m_signalParent = new double[n_bins];
    m_bkgdParent = new double[n_bins];
    m_signalErr = new double[n_bins];
    m_bkgdErr = new double[n_bins];
    m_data = new double[n_bins];

    memcpy(m_signal, sbd.m_signal, sizeof(double) * n_bins);
    memcpy(m_bkgd, sbd.m_bkgd, sizeof(double) * n_bins);

    memcpy(m_signalParent, sbd.m_signalParent, sizeof(double) * n_bins);
    memcpy(m_bkgdParent, sbd.m_bkgdParent, sizeof(double) * n_bins);
    memcpy(m_signalErr, sbd.m_signalErr, sizeof(double) * n_bins);
    memcpy(m_bkgdErr, sbd.m_bkgdErr, sizeof(double) * n_bins);
    memcpy(m_data, sbd.m_data, sizeof(double) * n_bins);

  } else {
    m_signal = NULL;
    m_bkgd = NULL;
    m_signalParent = NULL;
    m_bkgdParent = NULL;
    m_signalErr = NULL;
    m_bkgdErr = NULL;
    m_data = NULL;
  }

  m_deleteArrays = true;
  m_rebinned = sbd.m_rebinned;
  m_sigFluct = sbd.m_sigFluct;

  if (m_rebinned) {
    n_trueBins = sbd.n_trueBins;
    m_rebinBins = new int[n_trueBins];
    memcpy(m_rebinBins, sbd.m_rebinBins, sizeof(int) * n_trueBins);
  } else
    m_rebinBins = NULL;

  gl_fluctMap = map<string, double>(sbd.gl_fluctMap);
  gl_systCuts = map<string, double>(sbd.gl_systCuts);
  gl_systNames = vector<string>(sbd.gl_systNames);
  gl_bkgdNames = vector<string>(sbd.gl_bkgdNames);
  gl_sigNames = vector<string>(sbd.gl_sigNames);
  gl_fitAllSources = sbd.gl_fitAllSources;

  gl_systRand = std::vector<double>(gl_systNames.size(), 0.0);

  gl_systFloat = new bool[gl_systNames.size()];
  memset(gl_systFloat, 0, sizeof(bool) * gl_systNames.size());

  gl_systFit = new bool[gl_systNames.size()];
  memset(gl_systFit, 1, sizeof(bool) * gl_systNames.size());

  gl_bkgdFit = new bool[gl_bkgdNames.size()];
  memset(gl_bkgdFit, 1, sizeof(bool) * gl_bkgdNames.size());

  gl_sigFit = new bool[gl_sigNames.size()];
  memset(gl_sigFit, 1, sizeof(bool) * gl_sigNames.size());

  gl_systCenter = new double[gl_systNames.size()];
  memset(gl_systCenter, 0, sizeof(double) * gl_systNames.size());

  for (uint i = 0; i < gl_systNames.size(); ++i) {
    if (sbd.getFloatFlag(i))
      gl_systFloat[i] = 1;
    if (!sbd.getSystFitFlag(i))
      gl_systFit[i] = 0;
  }

  for (uint i = 0; i < gl_bkgdNames.size(); ++i) {
    if (!sbd.getBkgdFitFlag(i))
      gl_bkgdFit[i] = 0;
  }

  for (uint i = 0; i < gl_sigNames.size(); ++i) {
    if (!sbd.getSigFitFlag(i))
      gl_sigFit[i] = 0;
  }

  gl_exclChan = sbd.gl_exclChan;
  m_varied = false;
  m_init = false;
  m_condensed = false;
  m_verbose = false;
  m_adjSyst = false;
  m_systScale = 1;
  m_systOffset = 0;

  this->setBaselineModel(true);
}

//detailed constructor
SigBkgdDist::SigBkgdDist(string chan, int bins, double min, double max, int var1, int var2,
                         int var3)
    : m_signalExclusionSumsReady(false),
      m_bkgdExclusionSumsReady(false) {
  setup(chan, bins, min, max, var1, var2, var3);
  m_signal = new double[bins];
  m_bkgd = new double[bins];

  m_signalParent = new double[bins];
  m_bkgdParent = new double[bins];
  m_signalErr = new double[bins];
  m_bkgdErr = new double[bins];
  m_data = new double[bins];

  memset(m_signal, 0, sizeof(double) * bins);
  memset(m_bkgd, 0, sizeof(double) * bins);

  memset(m_signalParent, 0, sizeof(double) * bins);
  memset(m_bkgdParent, 0, sizeof(double) * bins);
  memset(m_signalErr, 0, sizeof(double) * bins);
  memset(m_bkgdErr, 0, sizeof(double) * bins);
  memset(m_data, 0, sizeof(double) * bins);

  m_sigScale.clear();
  m_sigScalesVary = false;
  gl_cSigScale = 1.0;
  m_bkgdScale.clear();
  m_bkgdScale.push_back(1.0);
  m_dataScale = 1.0;
  m_deleteArrays = true;

  m_rebinned = false;
  n_trueBins = 0;
  m_condensed = false;
  m_verbose = false;
  m_useStat = false;
  m_sigFluct = true;
  m_sbhypo = false;
  m_bohypo = false;
  m_rebinBins = NULL;

  m_fitValueHist = NULL;
  m_errMatHist = NULL;

  gl_systRand.clear(); // USED TO BE: gl_systRand=NULL;
  gl_systFloat = NULL;
  gl_systFit = NULL;
  gl_bkgdFit = NULL;
  gl_sigFit = NULL;
  gl_systCenter = NULL;
  gl_fluctMap.clear();
  gl_systCuts.clear();
  gl_fitAllSources = true;

  m_varied = false;
  m_init = false;
  m_adjSyst = false;
  m_systScale = 1;
  m_systOffset = 0;
  n_cacheData = 0;
  m_delBins.clear();

  // WARNING - THIS STEP VIOLATES THE USUAL C++ COPY SEMANTICS

  fillArrays(true, false);
  // WARNING - THIS STEP MAY VIOLATE THE USUAL C++ COPY SEMANTICS
}

//really detailed constructor
SigBkgdDist::SigBkgdDist(string chan, int bins, double min, double max, double* sigarray,
                         double* bkgdarray, double* dataarray, int var1, int var2, int var3)
    : m_signalExclusionSumsReady(false),
      m_bkgdExclusionSumsReady(false) {
  setup(chan, bins, min, max, var1, var2, var3);
  m_signal = sigarray;
  m_bkgd = bkgdarray;

  m_signalParent = sigarray;
  m_bkgdParent = bkgdarray;
  m_signalErr = new double[bins];
  m_bkgdErr = new double[bins];
  memset(m_signalErr, 0, sizeof(double) * bins);
  memset(m_bkgdErr, 0, sizeof(double) * bins);
  m_data = dataarray;
  m_sigScale.clear();
  m_sigScalesVary = false;
  gl_cSigScale = 1.0;
  m_bkgdScale.clear();
  m_bkgdScale.push_back(1.0);
  m_dataScale = 1.0;
  m_deleteArrays = false;
  m_sigFluct = true;
  m_useStat = false;
  m_condensed = false;
  m_verbose = false;
  m_rebinned = false;
  m_sbhypo = false;
  m_bohypo = false;
  n_trueBins = 0;
  m_rebinBins = NULL;
  m_fitValueHist = NULL;
  m_errMatHist = NULL;

  gl_systRand.clear(); // USED TO BE: gl_systRand=NULL;
  gl_systFloat = NULL;
  gl_systFit = NULL;
  gl_bkgdFit = NULL;
  gl_sigFit = NULL;
  gl_systCenter = NULL;
  gl_fluctMap.clear();
  gl_systCuts.clear();
  gl_fitAllSources = true;

  m_varied = false;
  m_init = false;
  m_adjSyst = false;
  m_systScale = 1;
  m_systOffset = 0;
  n_cacheData = 0;
  m_delBins.clear();

  fillArrays(true, false);
}

//destructor
SigBkgdDist::~SigBkgdDist() {

  m_SigDist.clear();
  m_BkgdDist.clear();
  m_DataDist.clear();

  if (m_deleteArrays) {
    delete[] m_signal;
    delete[] m_bkgd;
    delete[] m_signalParent;
    delete[] m_bkgdParent;
    delete[] m_data;
    delete[] m_rebinBins;
  }
  delete[] m_signalErr;
  delete[] m_bkgdErr;
  delete[] gl_systFloat;
  delete[] gl_systFit;
  delete[] gl_bkgdFit;
  delete[] gl_sigFit;
  delete[] gl_systCenter;
  setup(string(""), 0, 0, 0, 0, 0, 0);
}

///std setup method
void SigBkgdDist::setup(string chan, int bins, double min, double max, int var1, int var2,
                        int var3) {

  m_var1 = var1;
  m_var2 = var2;
  m_var3 = var3;
  m_min = min;
  m_max = max;
  n_bins = bins;
  m_delBinsLookup.resize(n_bins);
  std::fill(m_delBinsLookup.begin(), m_delBinsLookup.end(), 0);
  m_chan = chan;
  gl_fluctMap.clear();
  gl_systCuts.clear();

  return;
}

void SigBkgdDist::addSigDist(string chan, string sig, CollieDistribution* dist) {
  map<string, map<string, CollieDistribution*>>::iterator mItS;
  mItS = m_SigDist.find(chan);

  if (mItS == m_SigDist.end()) {
    //We don't even have the channel!
    map<string, CollieDistribution*> inmap;
    inmap[sig] = dist;
    m_SigDist[chan] = inmap;

    map<string, double> smap;
    smap[sig] = 1.0;
    m_sigScale[chan] = smap;
  } else {
    //We have the channel, but do we have this signal?
    map<string, CollieDistribution*>::iterator mItD;
    mItD = mItS->second.find(sig);
    if (mItD == mItS->second.end()) {
      //OK, don't have the signal...so we add it.
      mItS->second[sig] = dist;

      map<string, map<string, double>>::iterator smap = m_sigScale.find(chan);
      smap->second[sig] = 1.0;
    } else {
      printf("SigBkgdDist::addSigDist, Already have signal %s for channel %s\n", sig.c_str(),
             chan.c_str());
    }
  }
}

void SigBkgdDist::addBkgdDist(string s, CollieDistribution* dist) {

  map<string, vector<CollieDistribution*>>::iterator iter = m_BkgdDist.find(s);
  if (iter == m_BkgdDist.end()) {
    vector<CollieDistribution*> invect;
    invect.push_back(dist);
    m_BkgdDist[s] = invect;
  } else {
    iter->second.push_back(dist);
    m_bkgdScale.push_back(1.0);
  }
}

void SigBkgdDist::addDataDist(string s, CollieDistribution* dist) {
  m_DataDist[s] = dist;
}

const CollieDistribution* SigBkgdDist::getSigDist(string chan, string sig) const {
  map<string, map<string, CollieDistribution*>>::const_iterator channel = m_SigDist.find(chan);
  if (channel != m_SigDist.end()) {
    map<string, CollieDistribution*>::const_iterator distribution = channel->second.find(sig);
    if (distribution != channel->second.end())
      return distribution->second;
  }
  return NULL;
}

const CollieDistribution* SigBkgdDist::getSigDist(string chan, unsigned int i) const {
  map<string, map<string, CollieDistribution*>>::const_iterator channel = m_SigDist.find(chan);
  if (channel != m_SigDist.end()) {
    if (i < channel->second.size()) {
      map<string, CollieDistribution*>::const_iterator distribution = channel->second.begin();
      for (uint ss = 0; ss < i; ss++)
        ++distribution;
      return distribution->second;
    }
  }
  return NULL;
}

const CollieDistribution* SigBkgdDist::getBkgdDist(string s, unsigned int i) const {
  map<string, vector<CollieDistribution*>>::const_iterator channel = m_BkgdDist.find(s);
  if (channel != m_BkgdDist.end()) {
    if (i < channel->second.size()) {
      return channel->second[i];
    }
  }
  return NULL;
}
const CollieDistribution* SigBkgdDist::getDataDist(string s) const {
  map<string, CollieDistribution*>::const_iterator distribution = m_DataDist.find(s);
  if (distribution != m_DataDist.end()) {
    return distribution->second;
  }
  return NULL;
}

vector<string> SigBkgdDist::getChannelNames() const {
  vector<string> out;
  for (map<string, CollieDistribution*>::const_iterator distribution = m_DataDist.begin();
       distribution != m_DataDist.end(); ++distribution)
    out.push_back(distribution->first);
  return out;
}

unsigned int SigBkgdDist::getNsigDist(string s) const {
  map<string, map<string, CollieDistribution*>>::const_iterator channel = m_SigDist.find(s);
  if (channel != m_SigDist.end()) {
    return channel->second.size();
  }
  return 0;
}

unsigned int SigBkgdDist::getNbkgdDist(string s) const {
  map<string, vector<CollieDistribution*>>::const_iterator channel = m_BkgdDist.find(s);
  if (channel != m_BkgdDist.end()) {
    return channel->second.size();
  }
  return 0;
}

//Read in fluctuations from a histogram
void SigBkgdDist::fluctuate(TH1D* fluctMap) {

  if (fluctMap == NULL) {
    printf("SigBkgdDist::fluctuate, NULL Histogram!!\n");
    return;
  }

  if (fluctMap->GetNbinsX() != getNsyst()) {
    printf("SigBkgdDist::fluctuate, Mismatch in fit param number!\n");
    return;
  }

  std::vector<double> f_params(fluctMap->GetNbinsX(), 0.0);

  for (int i = 0; i < fluctMap->GetNbinsX(); ++i) {
    f_params[i] = fluctMap->GetBinContent(i);
  }

  fluctuate(f_params);

  return;
}

//Read in fluctuations from a previous fit output file
void SigBkgdDist::fluctuate(TFile* fitFile, string fitType) {
  if (fitFile == NULL) {
    printf("SigBkgdDist::fluctuate, NULL File!!\n");
    return;
  }

  TH1D* fluctMap = NULL;
  if (fitType == string("Null Fit Params")) {
    fluctMap = (TH1D*)fitFile->Get(fitType.c_str());
  } else if (fitType == string("Test Fit Params")) {
    fluctMap = (TH1D*)fitFile->Get(fitType.c_str());
  }

  if (fluctMap == NULL) {
    printf("SigBkgdDist::fluctuate, NULL Histogram!!\n");
    return;
  } else
    fluctuate(fluctMap);

  return;
}

void SigBkgdDist::fillArrays(bool fillData, bool varySyst, std::vector<double> const& fluctMap,
                             TrialPoint& lastMinuitTrialPoint) {

  if (m_rebinned)
    return fillRebinnedArrays(fillData, varySyst, fluctMap.data());

  bool isStdFit = gl_fitAllSources;
  if (gl_exclChan != "")
    isStdFit = false;

  m_init = true;
  m_varied = varySyst;
  CollieDistribution* dist = 0;
  int bins = n_bins;
  if (n_cacheData != m_DataDist.size()) {
    bins = 0;
    n_cacheData = m_DataDist.size();
    for (iterD = m_DataDist.begin(); iterD != m_DataDist.end(); ++iterD) {
      dist = iterD->second;
      bins += (dist->getNYbins() < 1) ? dist->getNXbins() : (dist->getNYbins() * dist->getNXbins());
    }
  }
  if (bins == 0)
    return;

  if (bins != n_bins) {
    double* tempS = new double[bins];
    double* tempB = new double[bins];

    double* tempSP = new double[bins];
    double* tempBP = new double[bins];
    double* tempSErr = new double[bins];
    double* tempBErr = new double[bins];
    double* tempD = new double[bins];
    memset(tempS, 0, sizeof(double) * bins);
    memset(tempB, 0, sizeof(double) * bins);
    memset(tempSP, 0, sizeof(double) * bins);
    memset(tempBP, 0, sizeof(double) * bins);
    memset(tempSErr, 0, sizeof(double) * bins);
    memset(tempBErr, 0, sizeof(double) * bins);
    memset(tempD, 0, sizeof(double) * bins);

    if (m_signal != NULL)
      delete[] m_signal;
    if (m_bkgd != NULL)
      delete[] m_bkgd;
    if (m_signalParent != NULL)
      delete[] m_signalParent;
    if (m_bkgdParent != NULL)
      delete[] m_bkgdParent;
    if (m_signalErr != NULL)
      delete[] m_signalErr;
    if (m_bkgdErr != NULL)
      delete[] m_bkgdErr;
    if (m_data != NULL)
      delete[] m_data;
    m_signal = tempS;
    m_bkgd = tempB;

    m_signalParent = tempSP;
    m_bkgdParent = tempBP;
    m_signalErr = tempSErr;
    m_bkgdErr = tempBErr;
    m_data = tempD;
    n_bins = bins;
    fillData = true;
  } else {
    memset(m_signal, 0, sizeof(double) * n_bins);
    memset(m_bkgd, 0, sizeof(double) * n_bins);
    std::fill(m_bkgdParent, m_bkgdParent + n_bins, 0.0);     // (mf) fix issue
                                                             // of never clearing
                                                             // m_bkgdParent
    std::fill(m_signalParent, m_signalParent + n_bins, 0.0); // (mf) fix issue
                                                             // of never clearing
                                                             // m_signalParent
    std::fill(m_bkgdErr, m_bkgdErr + n_bins, 0.0);           // (mf) fix issue
                                                             // of never clearing
                                                             // m_bkgdErr
    std::fill(m_signalErr, m_signalErr + n_bins, 0.0);       // (mf) fix issue
                                                             // of never clearing
                                                             // m_signalErr

    if (fillData)
      memset(m_data, 0, sizeof(double) * n_bins);
  }

  int doneBin = 0;
  int lin = 0;
  int tBin = 0;
  int maxY = 0;
  int maxX = 0;
  int binY = 0;
  int binX = 0;
  double inval = 0;
  bool doVary = false;
  bool doStat = false;
  double scale = 1.0;

  std::vector<double>::const_iterator onlyNonBaseS = lastMinuitTrialPoint.diff(fluctMap);

  if (varySyst && m_sigFluct && isStdFit && !m_useStat)
  // Optimization is applied only if fitting
  // for a PE.  Warning:  External code
  // guarantees that in this case,
  // m_useStat is false.  If it is not, this
  // optimization will not faithfully
  // reproduce original behavior.
  {
    if (onlyNonBaseS == static_cast<std::vector<double>::const_iterator>(fluctMap.end())) {
      setNewPointSignalExpectations(&fluctMap[0], m_signal);
    } else {
      setPerturbedSignalExpectations(fluctMap, onlyNonBaseS, m_signal);
    }
  } else { // the not  ( varySyst && m_sigFluct ) case: original method used

    for (iterS = m_SigDist.begin(); iterS != m_SigDist.end(); ++iterS) { //loop over channels

      //query the data for the number of bins for this channel...
      dist = m_DataDist[iterS->first];
      maxY = (dist->getNYbins() < 1) ? 1 : dist->getNYbins();
      maxX = dist->getNXbins();

      for (iterD = iterS->second.begin(); iterD != iterS->second.end();
           ++iterD) { //loop over channels
        dist = iterD->second;
        lin = 0;
        inval = 0;

        //Collect information
        doVary = varySyst && m_sigFluct && (iterS->first != gl_exclChan) &&
                 getSigFitFlag(dist->GetName());
        doStat = m_useStat && varySyst;
        scale = getSignalScale(iterS->first, iterD->first);

        for (binY = 0; binY < maxY; ++binY)
          for (binX = 0; binX < maxX; ++binX) {

            inval = dist->getEfficiency(binX, binY) * scale;

            if (doStat) {
              double fstat = -1 * inval * 4;
              while ((fstat + inval) < 0) {
                fstat = m_randgaus.fire() * dist->getBinStatErr(binX, binY) * scale;
              }
              inval += fstat;
            }

            if (doVary) {
              double nom = dist->getEfficiency(binX, binY) * scale -
                           dist->getEfficiencyVaried(binX, binY, &fluctMap[0]) * scale;
              inval += nom;
            }

            /*
      if(doVary){
        inval = dist->getEfficiencyVaried(binX,binY,&fluctMap[0])*scale;
      }
      else{
        inval = dist->getEfficiency(binX,binY)*scale;
      }

      if(doStat){
       inval += m_randgaus.fire()*dist->getBinStatErr(binX,binY)*scale;
      }
      */

            if (inval < 0.0)
              inval = 1e-9;

            tBin = lin + doneBin;

            signal(tBin) += inval;

            if (fillData)
              m_signalErr[tBin] +=
                  dist->getBinErrValue(binX, binY) * dist->getBinErrValue(binX, binY);

            m_signalParent[tBin] += dist->getEfficiency(binX, binY) * scale;

            ++lin;
          }
      }
      doneBin += lin;
    }
  } // Terminate loop with original signal computation (if/else on varyStst && )

  doneBin = 0;
  if (varySyst && m_sigFluct && isStdFit && !m_useStat) {
    // Optimization is applied only if fitting
    if (onlyNonBaseS == static_cast<std::vector<double>::const_iterator>(fluctMap.end())) {
      setNewPointBkgdExpectations(&fluctMap[0], m_bkgd);
    } else {
      setPerturbedBkgdExpectations(fluctMap, onlyNonBaseS, m_bkgd);
    }
  } else { // the not  ( varySyst && m_sigFluct ) case: original method
    for (iterB = m_BkgdDist.begin(); iterB != m_BkgdDist.end(); ++iterB) {

      //query the data for the number of bins for this channel...
      dist = m_DataDist[iterB->first];
      maxY = (dist->getNYbins() < 1) ? 1 : dist->getNYbins();
      maxX = dist->getNXbins();

      gl_thisIdx = 0;
      for (iterV = iterB->second.begin(); iterV != iterB->second.end(); ++iterV) {
        dist = (*iterV);
        lin = 0;
        inval = 0;

        //Collect booleans
        doVary = varySyst && (iterB->first != gl_exclChan) && getBkgdFitFlag(dist->GetName());
        doStat = m_useStat && varySyst;
        scale = m_bkgdScale[gl_thisIdx];

        for (binY = 0; binY < maxY; ++binY)
          for (binX = 0; binX < maxX; ++binX) {

            if (doVary) {
              inval = (*iterV)->getEfficiencyVaried(binX, binY, &fluctMap[0]) * scale;
            } else {
              inval = (*iterV)->getEfficiency(binX, binY) * scale;
            }

            if (doStat) {
              inval += m_randgaus.fire() * dist->getBinStatErr(binX, binY) * scale;
            }
            if (inval < 0.0)
              inval = 1e-9;

            tBin = lin + doneBin;

            bkgd(tBin) += inval;

            if (fillData)
              m_bkgdErr[tBin] +=
                  dist->getBinErrValue(binX, binY) * dist->getBinErrValue(binX, binY);

            m_bkgdParent[tBin] += dist->getEfficiency(binX, binY) * scale;

            ++lin;
          }
        ++gl_thisIdx;
      }
      doneBin += lin;
    }
  } // Terminates if else for ( varySyst && m_sigFluct )

  doneBin = 0;

  if (fillData) {
    doneBin = 0;
    for (iterD = m_DataDist.begin(); iterD != m_DataDist.end(); ++iterD) {
      dist = iterD->second;
      lin = 0;
      maxY = (dist->getNYbins() < 1) ? 1 : dist->getNYbins();
      maxX = dist->getNXbins();
      for (binY = 0; binY < maxY; ++binY)
        for (binX = 0; binX < maxX; ++binX) {
          data(lin + doneBin) += dist->getEfficiency(binX, binY) * m_dataScale;
          ++lin;
        }
      doneBin += lin;
    }

    for (int b = 0; b < n_bins; ++b) {
      m_bkgdErr[b] = sqrt(m_bkgdErr[b]);
      m_signalErr[b] = sqrt(m_signalErr[b]);
      if (m_sbhypo)
        m_data[b] = m_bkgdParent[b] + m_signalParent[b];
      if (m_bohypo)
        m_data[b] = m_bkgdParent[b];
    }
  }

  return;
}

void SigBkgdDist::fillRebinnedArrays(bool fillData, bool varySyst, const double* fluctMap) {

  m_varied = varySyst;
  m_init = true;
  memset(m_signal, 0, sizeof(double) * n_bins);
  memset(m_bkgd, 0, sizeof(double) * n_bins);
  if (fillData) {
    memset(m_signalParent, 0, sizeof(double) * n_bins);
    memset(m_bkgdParent, 0, sizeof(double) * n_bins);
    memset(m_data, 0, sizeof(double) * n_bins);
  }
  CollieDistribution* dist;
  int doneBin = 0;
  int lin = 0;
  int maxJ = 0;
  double inval = 0;
  int j, i;
  for (iterS = m_SigDist.begin(); iterS != m_SigDist.end(); ++iterS) { //loop over channels
    lin = 0;
    for (iterD = iterS->second.begin(); iterD != iterS->second.end(); ++iterD) { //loop over signals
      lin = 0;
      dist = iterD->second;
      maxJ = (dist->getNYbins() < 1) ? 1 : dist->getNYbins();
      inval = 0;
      for (j = 0; j < maxJ; ++j)
        for (i = 0; i < dist->getNXbins(); ++i) {
          if (varySyst && m_sigFluct) {
            inval = dist->getEfficiencyVaried(i, j, fluctMap) *
                    getSignalScale(iterS->first, iterD->first);
          } else {
            inval = dist->getEfficiency(i, j) * getSignalScale(iterS->first, iterD->first);
          }

          if (m_useStat && varySyst) {
            inval += m_randgaus.fire() * dist->getBinStatErr(i, j) *
                     getSignalScale(iterS->first, iterD->first);
          }

          if (inval < 0)
            inval = 0;
          signal(m_rebinBins[lin + doneBin]) += inval;
          ++lin;
        }
    }
    doneBin += lin;
  }

  doneBin = 0;
  for (iterB = m_BkgdDist.begin(); iterB != m_BkgdDist.end(); ++iterB) {
    lin = 0;
    gl_thisIdx = 0;
    for (iterV = iterB->second.begin(); iterV != iterB->second.end(); ++iterV) {
      lin = 0;
      dist = (*iterV);
      maxJ = (dist->getNYbins() < 1) ? 1 : dist->getNYbins();
      inval = 0;
      for (j = 0; j < maxJ; ++j)
        for (i = 0; i < dist->getNXbins(); ++i) {
          // USED TO BE (and will evolve): if(varySyst) inval=dist->getEfficiencyVaried(i,j,fluctMap)*m_bkgdScale;
          if (varySyst) {
            inval = dist->getEfficiencyVaried(i, j, fluctMap) * m_bkgdScale[gl_thisIdx];
          } else {
            inval = dist->getEfficiency(i, j) * m_bkgdScale[gl_thisIdx];
          }

          if (m_useStat && varySyst) {
            inval += m_randgaus.fire() * dist->getBinStatErr(i, j) * m_bkgdScale[gl_thisIdx];
          }

          bkgd(m_rebinBins[lin + doneBin]) += inval;
          ++lin;
        }
      ++gl_thisIdx;
    }
    doneBin += lin;
  }

  if (fillData) {
    doneBin = 0;
    for (iterD = m_DataDist.begin(); iterD != m_DataDist.end(); ++iterD) {
      dist = iterD->second;
      int lin = 0;
      int maxJ = (dist->getNYbins() < 1) ? 1 : dist->getNYbins();
      for (int j = 0; j < maxJ; ++j)
        for (int i = 0; i < dist->getNXbins(); ++i) {
          data(m_rebinBins[lin + doneBin]) += dist->getEfficiency(i, j) * m_dataScale;
          ++lin;
        }
      doneBin += lin;
    }
  }

  return;
}

///integration methods
double SigBkgdDist::totSignal() const {
  double tot = 0;
  for (int i = 0; i < n_bins; ++i)
    tot += m_signal[i];

  return tot;
}

double SigBkgdDist::totBkgd() const {
  double tot = 0;
  for (int i = 0; i < n_bins; ++i)
    tot += m_bkgd[i];
  return tot;
}

double SigBkgdDist::totData() const {
  double tot = 0;
  for (int i = 0; i < n_bins; ++i)
    tot += m_data[i];
  return tot;
}

// fill data with S+B prediction
void SigBkgdDist::fillDataSB(bool sb) {
  m_sbhypo = sb;
  if (sb) {
    m_bohypo = false;
    for (int i = 0; i < n_bins; ++i)
      m_data[i] = m_bkgdParent[i] + m_signalParent[i];
  }
}
// fill data with B-only prediction
void SigBkgdDist::fillDataBOnly(bool bo) {
  m_bohypo = bo;
  if (bo) {
    m_bohypo = false;
    for (int i = 0; i < n_bins; ++i)
      m_data[i] = m_bkgdParent[i];
  }
}

//scaling methods
void SigBkgdDist::scaleSignal(double f) {
  gl_cSigScale *= f;

  for (iterSclO = m_sigScale.begin(); iterSclO != m_sigScale.end(); ++iterSclO) {
    for (iterSclI = iterSclO->second.begin(); iterSclI != iterSclO->second.end(); ++iterSclI) {
      if (!getSigFitFlag(iterSclI->first))
        continue;
      iterSclI->second *= f;
    }
  }

  for (int i = 0; i < n_bins; ++i)
    m_signal[i] *= f;

  return;
}

void SigBkgdDist::fillSignalArrays(std::vector<double> const& fluctMap) {

  bool varySyst = (fluctMap.size() > 0);
  int lin = 0;
  int tBin = 0;
  int doneBin = 0;
  double inval = 0;
  int binY = 0;
  int binX = 0;
  int maxY = 0;
  int maxX = 0;
  double scale = 1.0;

  for (int b = 0; b < n_bins; ++b)
    signal(b) = 0;

  CollieDistribution* dist = 0;
  for (iterS = m_SigDist.begin(); iterS != m_SigDist.end(); ++iterS) { //loop over channels

    //query the data for the number of bins for this channel...
    dist = m_DataDist[iterS->first];
    maxY = (dist->getNYbins() < 1) ? 1 : dist->getNYbins();
    maxX = dist->getNXbins();

    for (iterD = iterS->second.begin(); iterD != iterS->second.end(); ++iterD) { //loop over signals
      lin = 0;
      inval = 0;
      dist = iterD->second;
      scale = getSignalScale(iterS->first, iterD->first);

      for (binY = 0; binY < maxY; ++binY)
        for (binX = 0; binX < maxX; ++binX) {
          if (varySyst && m_sigFluct) {
            inval = dist->getEfficiencyVaried(binX, binY, &fluctMap[0]) * scale;
          } else
            inval = dist->getEfficiency(binX, binY) * scale;

          tBin = lin + doneBin;

          m_signalParent[tBin] += dist->getEfficiency(binX, binY) * scale;

          if (inval < 0.0)
            inval = 1e-9;

          signal(tBin) += inval;

          m_signalParent[tBin] += dist->getEfficiency(binX, binY) * scale;

          ++lin;
        }
    }
    doneBin += lin;
  }
  return;
}

void SigBkgdDist::scaleSignalChan(string chan, double f, std::vector<double> const& fluctMap) {
  m_sigScalesVary = true;

  iterSclO = m_sigScale.find(chan);

  if (iterSclO != m_sigScale.end()) {
    for (iterSclI = iterSclO->second.begin(); iterSclI != iterSclO->second.end(); ++iterSclI) {
      if (!getSigFitFlag(iterSclI->first))
        continue;
      iterSclI->second *= f;
    }
  } else {
    printf("SigBkgdDist::scaleSignal, This channel doesn't exist! %s\n", chan.c_str());
    return;
  }

  fillSignalArrays(fluctMap);

  return;
}

void SigBkgdDist::scaleSignal(string sig, double f, std::vector<double> const& fluctMap) {

  if (!getSigFitFlag(sig))
    return;

  m_sigScalesVary = true;

  for (iterSclO = m_sigScale.begin(); iterSclO != m_sigScale.end(); ++iterSclO) {
    iterSclI = iterSclO->second.find(sig);
    if (iterSclI != iterSclO->second.end()) {
      iterSclI->second = iterSclI->second * f;
    }
  }

  fillSignalArrays(fluctMap);

  return;
}

void SigBkgdDist::scaleSignal(string chan, string sig, double f,
                              std::vector<double> const& fluctMap) {
  m_sigScalesVary = true;

  iterSclO = m_sigScale.find(chan);

  if (iterSclO != m_sigScale.end()) {
    iterSclI = iterSclO->second.find(sig);
    if (iterSclI != iterSclO->second.end()) {
      iterSclI->second *= f;
    } else {
      printf("SigBkgdDist::scaleSignal, This signal doesn't exist! %s\n", sig.c_str());
      return;
    }
  } else {
    printf("SigBkgdDist::scaleSignal, This channel doesn't exist! %s\n", chan.c_str());
    return;
  }

  fillSignalArrays(fluctMap);

  return;
}

void SigBkgdDist::scaleBackground(double f) {
  for (uint i = 0; i < m_bkgdScale.size(); ++i)
    m_bkgdScale[i] *= f;

  for (int i = 0; i < n_bins; ++i)
    m_bkgd[i] *= f;

  return;
}

void SigBkgdDist::scaleBackground(uint i, double f, std::vector<double> const& fluctMap) {
  if (i >= 0 && i < m_bkgdScale.size())
    m_bkgdScale[i] *= f;
  else
    return;

  bool varySyst = (fluctMap.size() > 0);
  int lin = 0;
  int doneBin = 0;
  double inval = 0;
  int binY = 0;
  int binX = 0;
  int maxY = 0;
  int maxX = 0;

  for (int b = 0; b < n_bins; ++b)
    bkgd(b) = 0;

  for (iterB = m_BkgdDist.begin(); iterB != m_BkgdDist.end(); ++iterB) {
    gl_thisIdx = 0;
    for (iterV = iterB->second.begin(); iterV != iterB->second.end(); ++iterV) {
      lin = 0;
      inval = 0;
      maxY = ((*iterV)->getNYbins() < 1) ? 1 : (*iterV)->getNYbins();
      maxX = (*iterV)->getNXbins();
      for (binY = 0; binY < maxY; ++binY)
        for (binX = 0; binX < maxX; ++binX) {
          if (varySyst) {
            inval =
                (*iterV)->getEfficiencyVaried(binX, binY, &fluctMap[0]) * m_bkgdScale[gl_thisIdx];
          } else
            inval = (*iterV)->getEfficiency(binX, binY) * m_bkgdScale[gl_thisIdx];

          m_bkgdParent[lin + doneBin] +=
              (*iterV)->getEfficiency(binX, binY) * m_bkgdScale[gl_thisIdx];

          if (inval < 0.0)
            inval = 1e-9;

          bkgd(lin + doneBin) += inval;
          m_bkgdParent[lin + doneBin] +=
              (*iterV)->getEfficiency(binX, binY) * m_bkgdScale[gl_thisIdx];
          ++lin;
        }
      ++gl_thisIdx;
    }
    doneBin += lin;
  }

  return;
}

void SigBkgdDist::scaleData(double f) {
  m_dataScale *= f;

  for (int i = 0; i < n_bins; ++i)
    m_data[i] *= f;

  return;
}

double SigBkgdDist::getSignalScale(string chan, string sig) const {

  if (!m_sigScalesVary) {
    return gl_cSigScale;
  }

  map<string, map<string, double>>::const_iterator iterO = m_sigScale.find(chan);

  if (iterO == m_sigScale.end()) {
    printf("SigBkgdDist::getSignalScale, index out of range\n");
    return 1;
  }

  map<string, double>::const_iterator iterI = iterO->second.find(sig);
  if (iterI == iterO->second.end()) {
    printf("SigBkgdDist::getSignalScale, index out of range\n");
    return 1;
  }

  return iterI->second;
}

double SigBkgdDist::getBackgroundScale(int i) const {
  if (i < 0 || static_cast<size_t>(i) >= m_bkgdScale.size()) {
    printf("SigBkgdDist::getBackgroundScale, index out of range\n");
    return 1;
  }
  return m_bkgdScale[i];
}

void SigBkgdDist::setSignalScale(double f) {
  gl_cSigScale = f;

  for (iterSclO = m_sigScale.begin(); iterSclO != m_sigScale.end(); ++iterSclO) {
    for (iterSclI = iterSclO->second.begin(); iterSclI != iterSclO->second.end(); ++iterSclI) {
      iterSclI->second = f;
    }
  }

  m_sigScalesVary = false;

  fillArrays(true);
  return;
}

void SigBkgdDist::setSignalScaleChan(string chan, double f) {

  iterSclO = m_sigScale.find(chan);

  if (iterSclO != m_sigScale.end()) {
    for (iterSclI = iterSclO->second.begin(); iterSclI != iterSclO->second.end(); ++iterSclI) {
      iterSclI->second = f;
    }
  } else {
    printf("SigBkgdDist::setSignalScale, This channel doesn't exist! %s\n", chan.c_str());
    return;
  }

  m_sigScalesVary = true;

  fillArrays(true);
  return;
}

void SigBkgdDist::setSignalScale(int sigIndex, double f) {

  m_sigScalesVary = true;

  if (sigIndex < 0 || static_cast<size_t>(sigIndex) >= gl_sigNames.size()) {
    printf("SigBkgdDist::setSignalScale, Invalid vector index %d", sigIndex);
    return;
  }

  string sig = gl_sigNames[sigIndex];
  //  printf("Index %d, Sig: %s\n",sigIndex,sig.c_str());

  for (iterSclO = m_sigScale.begin(); iterSclO != m_sigScale.end(); ++iterSclO) {
    iterSclI = iterSclO->second.find(sig);
    if (iterSclI != iterSclO->second.end()) {
      iterSclI->second = f;
    }
  }

  m_sigScalesVary = true;

  fillArrays(true);
  return;
}

void SigBkgdDist::setSignalScale(string sig, double f) {

  m_sigScalesVary = true;

  for (iterSclO = m_sigScale.begin(); iterSclO != m_sigScale.end(); ++iterSclO) {
    iterSclI = iterSclO->second.find(sig);
    if (iterSclI != iterSclO->second.end()) {
      iterSclI->second = f;
    }
  }

  m_sigScalesVary = true;

  fillArrays(true);
  return;
}

void SigBkgdDist::setSignalScale(string chan, string sig, double f) {
  iterSclO = m_sigScale.find(chan);

  if (iterSclO != m_sigScale.end()) {
    iterSclI = iterSclO->second.find(sig);
    if (iterSclI != iterSclO->second.end()) {
      iterSclI->second = f;
    } else {
      printf("SigBkgdDist::scaleSignal, This signal doesn't exist! %s\n", sig.c_str());
      return;
    }
  } else {
    printf("SigBkgdDist::scaleSignal, This channel doesn't exist! %s\n", chan.c_str());
    return;
  }

  m_sigScalesVary = true;

  fillArrays(true);
  return;
}

void SigBkgdDist::setBackgroundScale(double f) {
  for (uint i = 0; i < m_bkgdScale.size(); ++i)
    m_bkgdScale[i] = f;

  fillArrays(true);
  return;
}

void SigBkgdDist::setDataScale(double f) {
  m_dataScale = f;

  fillArrays(true);
  return;
}

///appending methods
int SigBkgdDist::append(const SigBkgdDist& asbd) {

  if (m_rebinned) {
    printf("SigBkgdDist::append, Undoing the rebinning process for appending\n");
    m_rebinned = false;
    fillArrays(true); //reset to standard array sizes...
  }

  const map<string, map<string, CollieDistribution*>>& tmpS = asbd.m_SigDist;
  for (map<string, map<string, CollieDistribution*>>::const_iterator s = tmpS.begin();
       s != tmpS.end(); ++s)
    for (map<string, CollieDistribution*>::const_iterator d = s->second.begin();
         d != s->second.end(); ++d)
      addSigDist(s->first, d->first, d->second);

  const map<string, vector<CollieDistribution*>>& tmpB = asbd.m_BkgdDist;
  for (map<string, vector<CollieDistribution*>>::const_iterator b = tmpB.begin(); b != tmpB.end();
       ++b)
    for (vector<CollieDistribution*>::const_iterator d = b->second.begin(); d != b->second.end();
         ++d)
      addBkgdDist(b->first, *d);

  const map<string, CollieDistribution*>& tmpD = asbd.m_DataDist;
  for (map<string, CollieDistribution*>::const_iterator d = tmpD.begin(); d != tmpD.end(); ++d)
    addDataDist(d->first, d->second);

  fillCondensedMap();
  fillArrays(true);

  m_delBinsLookup.resize(n_bins);
  std::fill(m_delBinsLookup.begin(), m_delBinsLookup.end(), 0);

  return n_bins;
}

void SigBkgdDist::varySystematics() {
  //Condense systematics names.  Do this because each correlated source
  //   is fluctuated only once...  Only do it once to save time...
  if (gl_fluctMap.size() == 0) {
    fillCondensedMap();
  }

  //Fill systematics maps, every time anew...
  const uint ns = gl_systNames.size();
  for (uint s = 0; s < ns; ++s) {
    if (gl_systFit[s])
      gl_systRand[s] = m_randgaus.fire();

    /*
    gl_systRand[s] += gl_systCenter[s];
    if(m_adjSyst){
      gl_systRand[s] *= m_systScale;
      gl_systRand[s] += m_systOffset;
    }
    */
  }

  return;
}

void SigBkgdDist::clearSystFluctValues(bool zeroSyst) {
  assert((gl_systRand.size() == gl_systNames.size()) || (gl_systRand.size() == 0));
  std::fill(gl_systRand.begin(), gl_systRand.end(), 0.0);

  if (zeroSyst) {
    if (gl_systCenter != NULL)
      memset(gl_systCenter, 0, sizeof(double) * gl_systNames.size());
  } else {
    for (unsigned int i = 0; i < gl_systNames.size(); ++i)
      gl_systRand[i] += gl_systCenter[i];
  }
  return;
}

void SigBkgdDist::setBaselineModel(bool copy) {
  if (gl_fluctMap.size() == 0) {
    fillCondensedMap(copy);
  }
  clearSystFluctValues(true);
  fillArrays(true);
  return;
}

void SigBkgdDist::zeroModel() {
  if (gl_fluctMap.size() == 0) {
    fillCondensedMap();
  }
  clearSystFluctValues(false);
  fillArrays(true, true, gl_systRand);
  return;
}

int SigBkgdDist::getSOBbin(double sob) {
  assert(m_rebinned);
  assert(m_max > m_min);

  int b = int((sob - m_min) / (m_max - m_min) * double(n_bins));
  if (b < 0)
    b = 0;
  if (b >= n_bins)
    b = n_bins - 1;

  return b;
}

int SigBkgdDist::rebin_sob(double min_l10_sob, double max_l10_sob, int ibins) {
  if (min_l10_sob >= max_l10_sob)
    return false;
  if (ibins <= 0)
    return false;

  double f_min = 1000;
  double f_max = -1000;

  delete[] m_rebinBins;
  m_rebinBins = NULL;
  m_rebinBins = new int[n_bins];
  memset(m_rebinBins, 0, sizeof(int) * n_bins);

  n_trueBins = n_bins;
  m_rebinned = true;
  n_bins = ibins;
  m_min = min_l10_sob;
  m_max = max_l10_sob;

  double* os = m_signal;
  m_signal = new double[ibins];
  double* ob = m_bkgd;
  m_bkgd = new double[ibins];
  double* od = m_data;
  m_data = new double[ibins];

  delete[] m_signalParent;
  m_signalParent = new double[ibins];
  delete[] m_bkgdParent;
  m_bkgdParent = new double[ibins];
  memset(m_signal, 0, sizeof(double) * ibins);
  memset(m_bkgd, 0, sizeof(double) * ibins);
  memset(m_signalParent, 0, sizeof(double) * ibins);
  memset(m_bkgdParent, 0, sizeof(double) * ibins);
  memset(m_data, 0, sizeof(double) * ibins);

  char title[256];
  sprintf(title, "S/B Histogram, Signal Var: %d", m_var1);
  TH1D* sob_sig = new TH1D(title, title, ibins, min_l10_sob, max_l10_sob);
  sprintf(title, "S/B Histogram, Bkgd Var: %d", m_var1);
  TH1D* sob_bkg = new TH1D(title, title, ibins, min_l10_sob, max_l10_sob);
  sprintf(title, "S/B Histogram, Data Var: %d", m_var1);
  TH1D* sob_dat = new TH1D(title, title, ibins, min_l10_sob, max_l10_sob);

  double sob = 0;
  for (int i = 0; i < n_trueBins; ++i) {
    if (ob[i] == 0) {
      if (os[i] != 0) {
        printf("SigBkgdDist::rebin_sob, zero background non-zero signal!!\n");
        m_signal[ibins - 1] += os[i];
        m_signalParent[ibins - 1] += os[i];
        m_data[ibins - 1] += od[i];
        sob_sig->AddBinContent(ibins - 1, os[i]);
        sob_dat->AddBinContent(ibins - 1, od[i]);
        m_rebinBins[i] = ibins - 1;
      }
      continue; // skip zero-sig/bkgd bins
    }
    sob = os[i] / ob[i];
    sob = log10(sob);
    if (os[i] > 0) {
      if (sob < f_min)
        f_min = sob;
      else if (sob > f_max)
        f_max = sob;
    }
    int sb = getSOBbin(sob);
    m_rebinBins[i] = sb;
    m_signal[sb] += os[i];
    m_signalParent[sb] += os[i];
    m_bkgd[sb] += ob[i];
    m_bkgdParent[sb] += ob[i];
    m_data[sb] += od[i];
    sob_sig->AddBinContent(sb, os[i]);
    sob_bkg->AddBinContent(sb, ob[i]);
    sob_dat->AddBinContent(sb, od[i]);
  }

  //  printf("SOB Max: %f, SOB Min: %f\n",f_max,f_min);
  delete[] os;
  os = NULL;
  delete[] ob;
  ob = NULL;
  delete[] od;
  od = NULL;

  return true;
}

int SigBkgdDist::rebin_sob_minimized(double min_l10_sob, double max_l10_sob, int bins) {

  if (min_l10_sob > max_l10_sob)
    return false;
  if (bins == 0)
    return false;

  double f_min = 1000;
  double f_max = -1000;
  sobMinMax_minimized(f_min, f_max);

  if (min_l10_sob < f_min) {
    printf("Adapting SOB Min from %.2f to %.2f\n", min_l10_sob, f_min - (0.01 * fabs(f_min)));
    min_l10_sob = f_min - (0.01 * fabs(f_min));
  }
  if (max_l10_sob > f_max) {
    printf("Adapting SOB Max from %.2f to %.2f\n", max_l10_sob, f_max + (0.01 * fabs(f_max)));
    max_l10_sob = f_max + (0.01 * fabs(f_max));
  }

  return rebin_sob(min_l10_sob, max_l10_sob, bins);
}

int SigBkgdDist::rebin_sob_adaptive(double min_l10_sob, double max_l10_sob, int bins) {

  if (min_l10_sob > max_l10_sob)
    return false;
  if (bins == 0)
    return false;

  double f_min = 1000;
  double f_max = -1000;
  sobMinMax(f_min, f_max);

  if (min_l10_sob < f_min) {
    printf("Adapting SOB Min from %.2f to %.2f\n", min_l10_sob, f_min - (0.01 * fabs(f_min)));
    min_l10_sob = f_min - (0.01 * fabs(f_min));
  }
  if (max_l10_sob > f_max) {
    printf("Adapting SOB Max from %.2f to %.2f\n", max_l10_sob, f_max + (0.01 * fabs(f_max)));
    max_l10_sob = f_max + (0.01 * fabs(f_max));
  }

  return rebin_sob(min_l10_sob, max_l10_sob, bins);
}

void SigBkgdDist::sobMinMax(double& min, double& max) {
  ////find boundary max/min s/b values
  for (int i = 0; i < n_bins; ++i) {
    if (m_bkgd[i] == 0 || m_signal[i] == 0)
      continue; // skip zero-background bins
    double sob = m_signal[i] / m_bkgd[i];
    sob = log10(sob);
    if (sob < min)
      min = sob;
    else if (sob > max)
      max = sob;
  }

  return;
}

void SigBkgdDist::sobMinMax_minimized(double& min, double& max) {
  ////find boundary max/min s/b values
  for (int i = 0; i < n_bins; ++i) {
    if (m_bkgd[i] == 0 || m_signal[i] == 0)
      continue; // skip zero-background bins
    double sob = m_signal[i] / m_bkgd[i];
    sob = log10(sob);
    if (sob < min)
      min = sob;
    else if (sob > max)
      max = sob;
  }

  //find min s/b enclosing 99% of signal...
  double sumSig = 0;
  double sigFrac = 1.0;
  while (sigFrac > 0.99) {
    sumSig = 0;
    for (int i = 0; i < n_bins; ++i) {
      if (m_bkgd[i] == 0 || m_signal[i] == 0)
        continue; // skip zero-background bins
      double sob = m_signal[i] / m_bkgd[i];
      sob = log10(sob);
      if (sob >= min)
        sumSig += m_signal[i];
    }
    sigFrac = sumSig / totSignal();
    if (sigFrac > 0.99)
      min += fabs(min / 1000.0);
  }

  return;
}

double SigBkgdDist::calculateDeltaLLR() const {
  double llrtotal = 0;
  double llr;
  for (int i = 0; i < n_bins; ++i) {
    if (m_bkgd[i] == 0 && m_signal[i] == 0) {
      continue;
    } //skip zero bins...
    if (m_bkgd[i] < 1e-6)
      llr = m_signal[i] * log(1.0 + m_signal[i] * 1e6);
    else
      llr = m_signal[i] * log(1.0 + m_signal[i] / m_bkgd[i]);
    llrtotal += llr;
  }
  return 2.0 * llrtotal;
}

double SigBkgdDist::calculateDeltaLLRobs() const {
  double llrtotal = 0;
  double llr;

  for (int i = 0; i < n_bins; ++i) {
    if (m_bkgd[i] == 0 && m_signal[i] == 0) {
      continue;
    } //skip zero bins...
    if (m_bkgd[i] < 1e-6)
      llr = -m_signal[i] + m_data[i] * log(1.0 + m_signal[i] * 1e6);
    else
      llr = -m_signal[i] + (m_signal[i] + m_bkgd[i]) * log(1.0 + m_signal[i] / m_bkgd[i]);
    llrtotal += llr;

    if (m_bkgd[i] < 1e-6)
      llr = -m_signal[i] + m_data[i] * log(1.0 + m_signal[i] * 1e6);
    else
      llr = -m_signal[i] + m_data[i] * log(1.0 + m_signal[i] / m_bkgd[i]);
    llrtotal -= llr;
  }
  return 2.0 * fabs(llrtotal);
}

double SigBkgdDist::calculateLLR() const {
  double llrtotal = 0;
  double llr;
  for (int i = 0; i < n_bins; ++i) {
    if (m_bkgd[i] == 0 && m_signal[i] == 0) {
      continue;
    } //skip zero bins...
    if (m_bkgd[i] < 1e-6)
      llr = -m_signal[i] + m_data[i] * log(1.0 + m_signal[i] * 1e6);
    else
      llr = -m_signal[i] + m_data[i] * log(1.0 + m_signal[i] / m_bkgd[i]);
    llrtotal += llr;
  }
  return -2.0 * llrtotal;
}

double SigBkgdDist::calculateLLR(int i) const {
  double llr;
  if (m_bkgd[i] == 0 && m_signal[i] == 0) {
    return 0;
  } //skip zero bins...
  if (m_bkgd[i] < 1e-6)
    llr = -m_signal[i] + m_data[i] * log(1.0 + m_signal[i] * 1e6);
  else
    llr = -m_signal[i] + m_data[i] * log(1.0 + m_signal[i] / m_bkgd[i]);
  return llr;
}

double SigBkgdDist::calculateLLRmax() const {
  double llrmax = 0;
  double llr;
  for (int i = 0; i < n_bins; ++i) {
    if (m_bkgd[i] == 0 && m_signal[i] == 0) {
      continue;
    } //skip zero bins...
    if (m_bkgd[i] < 1e-6)
      llr = -m_signal[i] + m_data[i] * log(1.0 + m_signal[i] * 1e6);
    else
      llr = -m_signal[i] + m_data[i] * log(1.0 + m_signal[i] / m_bkgd[i]);
    if (llr > llrmax)
      llrmax = llr;
  }
  return llrmax;
}

double SigBkgdDist::calculateLLRw() const {
  double llrtotal = 0;
  double llr;
  printf("\nCalc LLR: %d bins\n", n_bins);
  int bad = 0;
  for (int i = 0; i < n_bins; ++i) {
    if (m_bkgd[i] == 0 && m_signal[i] == 0)
      continue; //skip zero bins...
    if (m_bkgd[i] < 1e-6) {
      llr = -m_signal[i] + m_data[i] * log(1 + m_signal[i] * 1e6);
      ++bad;
    } else
      llr = -m_signal[i] + m_data[i] * log(1 + m_signal[i] / m_bkgd[i]);
    llrtotal += -llr;
  }
  printf("llr: %f, bad bins: %d\n", 2.0 * llrtotal, bad);
  return 2.0 * llrtotal;
}

double SigBkgdDist::calculateLLRb() const {
  double llrtotal = 0;
  double llr;
  for (int i = 0; i < n_bins; ++i) {
    if (m_bkgd[i] == 0 && m_signal[i] == 0)
      continue; //skip zero bins...
    if (m_bkgd[i] < 1e-6)
      llr = -m_signal[i] + m_bkgd[i] * log(1 + m_signal[i] * 1e6);
    else
      llr = -m_signal[i] + m_bkgd[i] * log(1 + m_signal[i] / m_bkgd[i]);
    llrtotal += -llr;
  }
  return 2.0 * llrtotal;
}

double SigBkgdDist::calculateLLRsb() const {
  double llrtotal = 0;
  double llr;
  for (int i = 0; i < n_bins; ++i) {
    if (m_bkgd[i] == 0 && m_signal[i] == 0)
      continue; //skip zero bins...
    if (m_bkgd[i] < 1e-6)
      llr = -m_signal[i] + (m_signal[i] + m_bkgd[i]) * log(1 + m_signal[i] * 1e6);
    else
      llr = -m_signal[i] + (m_signal[i] + m_bkgd[i]) * log(1 + m_signal[i] / m_bkgd[i]);
    llrtotal += -llr;
  }
  return 2.0 * llrtotal;
}

double SigBkgdDist::calculateFmax() const {
  double llrtotal = 0;
  double llr;
  for (int i = 0; i < n_bins; ++i) {
    if (m_bkgd[i] == 0 && m_signal[i] == 0)
      continue; //skip zero bins...
    if (m_bkgd[i] < 1e-6)
      llr = (m_signal[i] + m_bkgd[i]) * log(1 + m_signal[i] * 1e6);
    else
      llr = (m_signal[i] + m_bkgd[i]) * log(1 + m_signal[i] / m_bkgd[i]);
    llrtotal += llr;
  }
  return 5.0 * llrtotal;
}

double SigBkgdDist::calculateF(int i) const {
  double llrtotal = 0;
  double llr;
  if (i >= 0 && i < n_bins) {
    if (m_bkgd[i] < 1e-6)
      llr = log(1 + m_signal[i] * 1e6);
    else
      llr = log(1 + m_signal[i] / m_bkgd[i]);
    llrtotal += llr;
  }
  return llrtotal;
}

// The following is the orignal (pre-optimization) calculateChi2LLR
/*
double SigBkgdDist::calculateChi2LLR(bool fitSig, double sigLLR) const{

  double llrtotal=0;
  double mod; double dat;
  double sigLog = 0;
  double A=0; double B=0;

  for (int i=0; i<n_bins; ++i) {
    sigLog=10;
    if(m_bkgdParent[i]>0) sigLog = log10(1.0+m_signalParent[i]/m_bkgdParent[i]);
    if(sigLog>sigLLR && !fitSig) continue;

    for(uint j=0; j<m_delBins.size(); ++j) if(i==m_delBins[j]){ i++;  break; }

    mod = m_bkgd[i];
    if(fitSig) mod += m_signal[i];

    dat = m_data[i];
    if(mod==0 && dat==0) continue;

    A = mod-dat; B = 0;
    if(dat>0 && mod>0) B = dat * log(mod/dat);

    llrtotal += (A-B);
  }

  return 2.0*llrtotal;
}
*/
/*
double SigBkgdDist::calculateChi2LLRderivative(int ipar, bool fitSig, double sigLLR, double* fluctMap) {

  double dllrtotal=0;
  double mod; double dat; double dmod;
  double sigLog = 0;

  for (int i=0; i<n_bins; ++i) {
    sigLog=1;
    if(m_bkgdParent[i]>0) sigLog = log10(1.0+m_signalParent[i]/m_bkgdParent[i]);
    if(sigLog>sigLLR && !fitSig) continue;

    mod  = m_bkgd[i];
    dmod = d_bkgd[ipar][i]*m_bkgdParent[i];
    if(fitSig){
       mod += m_signal[i];
      dmod += d_signal[ipar][i]*m_signalParent[i];
    }
    dat = m_data[i];

    if(mod!=0) dllrtotal += dmod*(1.0-dat/mod);
  }
  return 2.0*dllrtotal;
}
*/

int SigBkgdDist::getSystIndex(string name) {
  for (uint s = 0; s < gl_systNames.size(); ++s) {
    if (gl_systNames[s] == name)
      return s;
  }
  return -1;
}

bool SigBkgdDist::getFloatFlag(string name) {
  for (uint s = 0; s < gl_systNames.size(); ++s) {
    if (gl_systNames[s] == name)
      return gl_systFloat[s];
  }
  return false;
}

bool SigBkgdDist::getSystFitFlag(string name) {
  for (uint s = 0; s < gl_systNames.size(); ++s) {
    if (gl_systNames[s] == name)
      return gl_systFit[s];
  }
  return false;
}

bool SigBkgdDist::getBkgdFitFlag(string name) {
  for (uint s = 0; s < gl_bkgdNames.size(); ++s) {
    if (gl_bkgdNames[s] == name)
      return gl_bkgdFit[s];
  }
  return false;
}

bool SigBkgdDist::getSigFitFlag(string name) {
  for (uint s = 0; s < gl_sigNames.size(); ++s) {
    if (gl_sigNames[s] == name)
      return gl_sigFit[s];
  }
  return false;
}

double SigBkgdDist::calculateSystDiff(std::vector<double> const& syst) const {
  //  return 0;
  double tot = 0;
  double diff = 0;
  const uint ns = gl_systNames.size();
  for (uint s = 0; s < ns; ++s) {
    if (gl_systFloat[s]) {
      //      printf("Skipping %d, %s\n",s,getSystName(s).c_str());
      continue;
    }
    if (std::isinf(syst[s]) || std::isnan(syst[s]))
      diff = 1.0;
    else
      diff = gl_systRand[s] - syst[s];
    tot += diff * diff;
  }
  return tot;
}

/*
double SigBkgdDist::calculateSystDiff(std::vector<double> const & syst) const{
  double tot = 0; double diff = 0;
  const uint ns = gl_systNames.size();
  for(uint s=0; s<ns; ++s){
    if(gl_systFloat[s]) continue;
    if(std::isinf(syst[s]) || std::isnan(syst[s])) diff = 1.0;
    else diff = gl_systRand[s] - syst[s];
    tot += diff*diff;
  }
  return tot;
}
*/

void SigBkgdDist::fillCondensedMap(bool copy) {

  gl_fluctMap.clear();
  gl_bkgdMap.clear();
  gl_sigMap.clear();
  gl_systNames.clear();
  gl_bkgdNames.clear();
  gl_sigNames.clear();
  gl_systCuts.clear();

  CollieDistribution* dist;
  int bins = 0;
  for (iterD = m_DataDist.begin(); iterD != m_DataDist.end(); ++iterD) {
    dist = iterD->second;
    bins += (dist->getNYbins() < 1) ? dist->getNXbins() : (dist->getNYbins() * dist->getNXbins());
  }

  uint bidx = 0;
  if (!copy)
    m_bkgdScale.clear();
  for (iterB = m_BkgdDist.begin(); iterB != m_BkgdDist.end(); ++iterB) {
    bidx = 0;
    for (iterV = iterB->second.begin(); iterV != iterB->second.end(); ++iterV) {
      dist = (*iterV);

      if (bidx >= m_bkgdScale.size() && !copy)
        m_bkgdScale.push_back(1.0);

      if (gl_bkgdMap.find(dist->GetName()) == gl_bkgdMap.end()) {
        gl_bkgdMap[dist->GetName()] = 1;
        gl_bkgdNames.push_back(dist->GetName());
      }
      ++bidx;

      for (int s = 0; s < dist->getNsystematics(); ++s) {
        //  double cut = dist->getSystLimit(s);
        double cut = 5;
        if (gl_fluctMap.find(dist->getSystName(s)) != gl_fluctMap.end()) {
          if (gl_systCuts[dist->getSystName(s)] > cut)
            gl_systCuts[dist->getSystName(s)] = cut;
          continue;
        } else {
          gl_systCuts[dist->getSystName(s)] = cut;
          gl_fluctMap[dist->getSystName(s)] = 0.0;
          gl_systNames.push_back(dist->getSystName(s));
        }
      }
    }
  }

  if (!copy)
    m_sigScale.clear();
  for (iterS = m_SigDist.begin(); iterS != m_SigDist.end(); ++iterS) {
    iterSclO = m_sigScale.find(iterS->first);
    for (iterD = iterS->second.begin(); iterD != iterS->second.end(); ++iterD) {
      dist = iterD->second;

      if (!copy) {
        if (iterSclO != m_sigScale.end()) { ///we already have this channel!
          iterSclI = iterSclO->second.find(iterD->first);
          if (iterSclI == iterSclO->second.end()) { //we don't have this signal
            iterSclO->second[iterD->first] = 1.0;
          } else
            continue; //we already have this signal for this channel??
        } else {      //We don't have this channel yet.
          map<string, double> inM;
          inM[iterD->first] = 1.0;
          m_sigScale[iterS->first] = inM;
        }
      }

      if (gl_sigMap.find(dist->GetName()) == gl_sigMap.end()) {
        gl_sigNames.push_back(dist->GetName());
        gl_sigMap[dist->GetName()] = 1;
      }

      for (int s = 0; s < dist->getNsystematics(); ++s) {
        if (gl_fluctMap.find(dist->getSystName(s)) != gl_fluctMap.end())
          continue;
        else {
          gl_fluctMap[dist->getSystName(s)] = 0.0;
          gl_systNames.push_back(dist->getSystName(s));
        }
      }
    }
  }

  //  for(map<string,int>::iterator iter = gl_bkgdMap.begin(); iter!=gl_bkgdMap.end(); ++iter)
  //    gl_bkgdNames.push_back(iter->first);

  //  for(map<string,int>::iterator iter = gl_sigMap.begin(); iter!=gl_sigMap.end(); ++iter)
  //    gl_sigNames.push_back(iter->first);

  gl_systRand = std::vector<double>(gl_systNames.size(), 0.0);

  delete[] gl_systFloat;
  gl_systFloat = NULL;
  gl_systFloat = new bool[gl_systNames.size()];
  memset(gl_systFloat, 0, sizeof(bool) * gl_systNames.size());

  delete[] gl_systFit;
  gl_systFit = NULL;
  gl_systFit = new bool[gl_systNames.size()];
  memset(gl_systFit, 1, sizeof(bool) * gl_systNames.size());

  delete[] gl_bkgdFit;
  gl_bkgdFit = NULL;
  gl_bkgdFit = new bool[gl_bkgdNames.size()];
  memset(gl_bkgdFit, 1, sizeof(bool) * gl_bkgdNames.size());

  delete[] gl_sigFit;
  gl_sigFit = NULL;
  gl_sigFit = new bool[gl_sigNames.size()];
  memset(gl_sigFit, 1, sizeof(bool) * gl_sigNames.size());

  delete[] gl_systCenter;
  gl_systCenter = NULL;
  gl_systCenter = new double[gl_systNames.size()];
  memset(gl_systCenter, 0, sizeof(double) * gl_systNames.size());

  /*
  d_signal = new double*[gl_systNames.size()];
  d_bkgd = new double*[gl_systNames.size()];
  for(uint i=0; i<gl_systNames.size(); ++i){
    d_signal[i] = new double[bins];
    d_bkgd[i] = new double[bins];
    memset(d_signal[i],0,sizeof(double)*bins);
    memset(d_bkgd[i],0,sizeof(double)*bins);
  }
  */
  for (iterS = m_SigDist.begin(); iterS != m_SigDist.end(); ++iterS) {
    for (iterD = iterS->second.begin(); iterD != iterS->second.end(); ++iterD) {
      iterD->second->linearize(gl_systNames);
      iterD->second->getFloatFlagList(gl_systNames, gl_systFloat);
    }
  }

  for (iterB = m_BkgdDist.begin(); iterB != m_BkgdDist.end(); ++iterB) {
    for (iterV = iterB->second.begin(); iterV != iterB->second.end(); ++iterV) {
      (*iterV)->linearize(gl_systNames);
      (*iterV)->getFloatFlagList(gl_systNames, gl_systFloat);
    }
  }

  for (iterD = m_DataDist.begin(); iterD != m_DataDist.end(); ++iterD)
    iterD->second->linearize(gl_systNames);

  m_condensed = true;

  return;
}

bool SigBkgdDist::setSystFluctValue(unsigned int i, double value) {
  if (i < 0 || i >= gl_systNames.size())
    return false;
  gl_systRand[i] = value;
  return true;
}

double SigBkgdDist::getSystFluctValue(string name) {
  for (uint i = 0; i < gl_systNames.size(); ++i) {
    if (gl_systNames[i] == name) {
      return gl_systRand[i];
    }
  }
  return -1e6;
}

bool SigBkgdDist::setSystFluctValue(string name, double value) {
  for (uint i = 0; i < gl_systNames.size(); ++i) {
    if (gl_systNames[i] == name) {
      gl_systRand[i] = value;
      return true;
    }
  }
  return false;
}

bool SigBkgdDist::setFloatFlag(unsigned int i, bool floatit) {
  if (!isInit()) {
    fluctuate();
    setBaselineModel();
  }
  if (i < 0 || i >= gl_systNames.size())
    return false;
  gl_systFloat[i] = floatit;
  return true;
}

bool SigBkgdDist::setFloatFlag(string name, bool floatit) {
  if (!isInit()) {
    fluctuate();
    setBaselineModel();
  }
  for (uint i = 0; i < gl_systNames.size(); ++i) {
    if (gl_systNames[i] == name) {
      gl_systFloat[i] = floatit;
      return true;
    }
  }
  return false;
}

bool SigBkgdDist::setSystFitFlag(unsigned int i, bool fitit) {
  if (!isInit()) {
    fluctuate();
    setBaselineModel();
  }
  if (i < 0 || i >= gl_systNames.size())
    return false;
  gl_systFit[i] = fitit;
  return true;
}

bool SigBkgdDist::setSystFitFlag(string name, bool fitit) {
  if (!isInit()) {
    fluctuate();
    setBaselineModel();
  }

  for (uint i = 0; i < gl_systNames.size(); ++i) {
    if (gl_systNames[i] == name) {
      gl_systFit[i] = fitit;
      return true;
    }
  }
  return false;
}

bool SigBkgdDist::setBkgdFitFlag(unsigned int i, bool fitit) {
  if (!isInit()) {
    fluctuate();
    setBaselineModel();
  }

  if (i < 0 || i >= gl_bkgdNames.size())
    return false;
  gl_bkgdFit[i] = fitit;

  gl_fitAllSources = true;
  for (uint i = 0; i < gl_bkgdNames.size(); ++i)
    if (gl_bkgdFit[i] == false)
      gl_fitAllSources = false;

  return true;
}

bool SigBkgdDist::setBkgdFitFlag(string name, bool fitit) {
  if (!isInit()) {
    fluctuate();
    setBaselineModel();
  }

  for (uint i = 0; i < gl_bkgdNames.size(); ++i) {
    if (gl_bkgdNames[i] == name) {
      gl_bkgdFit[i] = fitit;
      return true;
    }
  }

  gl_fitAllSources = true;
  for (uint i = 0; i < gl_bkgdNames.size(); ++i)
    if (gl_bkgdFit[i] == false)
      gl_fitAllSources = false;

  return false;
}

bool SigBkgdDist::setSigFitFlag(unsigned int i, bool fitit) {
  if (!isInit()) {
    fluctuate();
    setBaselineModel();
  }
  if (i < 0 || i >= gl_sigNames.size())
    return false;
  gl_sigFit[i] = fitit;

  gl_fitAllSources = true;
  for (uint i = 0; i < gl_sigNames.size(); ++i)
    if (gl_sigFit[i] == false)
      gl_fitAllSources = false;

  return true;
}

bool SigBkgdDist::setSigFitFlag(string name, bool fitit) {
  if (!isInit()) {
    fluctuate();
    setBaselineModel();
  }

  for (uint i = 0; i < gl_sigNames.size(); ++i) {
    if (gl_sigNames[i] == name) {
      gl_sigFit[i] = fitit;
      return true;
    }
  }

  gl_fitAllSources = true;
  for (uint i = 0; i < gl_sigNames.size(); ++i)
    if (gl_sigFit[i] == false)
      gl_fitAllSources = false;

  return false;
}

void SigBkgdDist::setSystCentralValue(string name, double par) {
  if (!isInit()) {
    fluctuate();
    setBaselineModel();
  }
  for (uint i = 0; i < gl_systNames.size(); ++i) {
    if (gl_systNames[i] == name)
      gl_systCenter[i] = par;
  }
  return;
}

void SigBkgdDist::setSystCentralValue(uint i, double par) {
  if (!isInit()) {
    fluctuate();
    setBaselineModel();
  }
  if (i < 0 || i > gl_systNames.size())
    return;
  else
    gl_systCenter[i] = par;
  return;
}

void SigBkgdDist::excludeChannel(string chan) {
  gl_exclChan = "";
  iterD = m_DataDist.find(chan);
  if (iterD != m_DataDist.end())
    gl_exclChan = chan;

  //Don't fit any systematics that are exclusive to this channel!

  /// Make a copy of what we're fitting or not
  /// Turn off all systematics
  vector<int> systCopy;
  for (int s = 0; s < getNsyst(); ++s) {
    systCopy.push_back(getSystFitFlag(s));
    setSystFitFlag(s, false);
  }

  //Turn on any systematics that are in the remaining channels
  for (iterB = m_BkgdDist.begin(); iterB != m_BkgdDist.end(); ++iterB) {
    if (iterB->first == gl_exclChan)
      continue;
    for (iterV = iterB->second.begin(); iterV != iterB->second.end(); ++iterV) {
      CollieDistribution* dist = (*iterV);
      for (int s = 0; s < dist->getNsystematics(); ++s) {
        setSystFitFlag(getSystIndex(dist->getSystName(s)), true);
      }
    }
  }

  for (iterS = m_SigDist.begin(); iterS != m_SigDist.end(); ++iterS) {
    for (iterD = iterS->second.begin(); iterD != iterS->second.end(); ++iterD) {
      if (iterS->first == gl_exclChan)
        continue;
      CollieDistribution* dist = iterD->second;

      for (int s = 0; s < dist->getNsystematics(); ++s) {
        setSystFitFlag(getSystIndex(dist->getSystName(s)), true);
      }
    }
  }

  //Turn off any systematics that were originally off
  for (int s = 0; s < getNsyst(); ++s) {
    if (!systCopy[s])
      setSystFitFlag(s, false);
  }

  if (gl_exclChan == "") {
    for (int s = 0; s < getNsyst(); ++s)
      setSystFitFlag(s, true);
  }

  return;
}

void SigBkgdDist::excludeChannel(int chan) {
  gl_exclChan = "";

  if (chan >= 0) {
    if ((uint)chan > m_DataDist.size())
      return;
    iterD = m_DataDist.begin();
    for (uint i = 0; i < (uint)chan; ++i)
      ++iterD;
    if (iterD == m_DataDist.end())
      return;
    gl_exclChan = iterD->first;

    //Don't fit any systematics that are exclusive to this channel!

    /// Make a copy of what we're fitting or not
    /// Turn off all systematics
    vector<int> systCopy;
    for (int s = 0; s < getNsyst(); ++s) {
      systCopy.push_back(getSystFitFlag(s));
      setSystFitFlag(s, false);
    }

    //Turn on any systematics that are in the remaining channels
    for (iterB = m_BkgdDist.begin(); iterB != m_BkgdDist.end(); ++iterB) {
      if (iterB->first == gl_exclChan)
        continue;
      for (iterV = iterB->second.begin(); iterV != iterB->second.end(); ++iterV) {
        CollieDistribution* dist = (*iterV);
        for (int s = 0; s < dist->getNsystematics(); ++s) {
          setSystFitFlag(getSystIndex(dist->getSystName(s)), true);
        }
      }
    }

    for (iterS = m_SigDist.begin(); iterS != m_SigDist.end(); ++iterS) {
      for (iterD = iterS->second.begin(); iterD != iterS->second.end(); ++iterD) {
        if (iterS->first == gl_exclChan)
          continue;
        CollieDistribution* dist = iterD->second;

        for (int s = 0; s < dist->getNsystematics(); ++s) {
          setSystFitFlag(getSystIndex(dist->getSystName(s)), true);
        }
      }
    }

    //Turn off any systematics that were originally off
    for (int s = 0; s < getNsyst(); ++s) {
      if (!systCopy[s])
        setSystFitFlag(s, false);
    }
  } else if (chan == -1) {
    for (int s = 0; s < getNsyst(); ++s)
      setSystFitFlag(s, true);
  }

  return;
}

string SigBkgdDist::getChannelName(uint chan) const {
  if (chan >= m_DataDist.size())
    return "";
  map<string, CollieDistribution*>::const_iterator distribution = m_DataDist.begin();
  for (uint i = 0; i < chan; ++i)
    ++distribution;
  return distribution->first;
}

void SigBkgdDist::calculateBkgdUncertainty(double& total, double& average, double& sigWeighted) {
  setBaselineModel();
  SigBkgdDist nom(*this);
  double varied[n_bins];

  for (int b = 0; b < n_bins; ++b)
    varied[b] = 0;
  total = 0;
  average = 0;
  sigWeighted = 0;

  for (int i = 0; i < 1e4; ++i) {
    fluctuate();
    for (int b = 0; b < n_bins; ++b)
      varied[b] += (nom.bkgd(b) - bkgd(b)) * (nom.bkgd(b) - bkgd(b)) / 1e4;
    total += (nom.totBkgd() - totBkgd()) * (nom.totBkgd() - totBkgd()) / 1e4;
  }
  total = sqrt(total) / nom.totBkgd();

  double sigBkgd = 0;
  double bkgdBins = 0;
  for (int b = 0; b < n_bins; ++b) {
    if (nom.bkgd(b) > 0)
      ++bkgdBins;
    average += varied[b];

    sigWeighted += varied[b] * nom.signal(b);
    if (nom.signal(b) > 0) {
      sigBkgd += nom.bkgd(b) * nom.signal(b);
    }
  }

  sigWeighted /= nom.totSignal();
  sigBkgd /= nom.totSignal();
  sigWeighted = sqrt(sigWeighted) / sigBkgd;

  average /= bkgdBins;
  average = sqrt(average) / (nom.totBkgd() / bkgdBins);

  setBaselineModel();
  return;
}

///ROOT Accessor functions
TH1D* SigBkgdDist::getSigHistogram() {
  char title[100];
  sprintf(title, "%s: Signal Dist", m_chan.c_str());
  TH1D* outHist = new TH1D(title, title, n_bins, m_min, m_max);

  for (int b = 0; b < n_bins; ++b)
    outHist->SetBinContent(b + 1, m_signal[b]);

  return outHist;
}

TH1D* SigBkgdDist::getBkgdHistogram() {
  char title[100];
  sprintf(title, "%s: Bkgd Dist", m_chan.c_str());
  TH1D* outHist = new TH1D(title, title, n_bins, m_min, m_max);

  for (int b = 0; b < n_bins; ++b)
    outHist->SetBinContent(b + 1, m_bkgd[b]);

  return outHist;
}

TH1D* SigBkgdDist::getDataHistogram() {
  char title[100];
  sprintf(title, "%s: Data Dist", m_chan.c_str());
  TH1D* outHist = new TH1D(title, title, n_bins, m_min, m_max);

  for (int b = 0; b < n_bins; ++b)
    outHist->SetBinContent(b + 1, m_data[b]);

  return outHist;
}

void SigBkgdDist::drawPlots() {

  TH1D* bkgd = getBkgdHistogram();
  TH1D* sig = getSigHistogram();
  TH1D* data = getDataHistogram();

  sig->SetFillColor(4);
  bkgd->SetLineWidth(3);
  bkgd->SetLineColor(2);
  data->SetLineWidth(2);

  return;
}

void SigBkgdDist::setNewPointSignalExpectations(const double* fluctMap, double* sig) {
  // We are relying on the maxX and maxY values obtained from each distribution
  // to match the space in the sig array.

  // Outer loops can still be on distributions, since if there is a commonality
  // in an s-related computation, it will already be amortized over all the bins
  // computed for.

  gl_addedBins = 0;
  gl_thisBin = 0;
  CollieDistribution* dist = 0;
  for (iterS = m_SigDist.begin(); iterS != m_SigDist.end(); ++iterS) {
    /*
      if(iterB->first==gl_exclChan){
      iterV=iterB->second.begin();
      CollieDistribution* dist = (*iterV);
      gl_addedBins += dist->getNXbins()*dist->getNYbins();
      gl_thisBin=0;
      continue;
      }
    */
    for (iterD = iterS->second.begin(); iterD != iterS->second.end(); ++iterD) {
      dist = iterD->second;

      // Now we loop over s.  This loop is dependent on info
      // owned by dist.  So we will set up the bins and then hand-off to dist.
      // But dist already knows about maxX and maxY so all we need to give it
      // is m_sigScale and where we want the results put.
      dist->addBinEfficiencies(fluctMap, getSignalScale(iterS->first, iterD->first), sig,
                               gl_addedBins);
      if (gl_thisBin == 0)
        gl_thisBin = dist->getNXbins() * dist->getNYbins();
    }
    gl_addedBins += gl_thisBin;
    gl_thisBin = 0;
  }

  m_signalExclusionSumsReady = false;
} // setNewPointSignalExpectations

void SigBkgdDist::setNewPointBkgdExpectations(const double* fluctMap, double* bkg) {

  gl_addedBins = 0;
  gl_thisBin = 0;
  gl_thisIdx = 0;
  for (iterB = m_BkgdDist.begin(); iterB != m_BkgdDist.end(); ++iterB) {
    /*
    if(iterB->first==gl_exclChan){
      iterV=iterB->second.begin();
      CollieDistribution* dist = (*iterV);
      gl_addedBins += dist->getNXbins()*dist->getNYbins();
      gl_thisBin=0;
      continue;
    }
    */
    for (iterV = iterB->second.begin(); iterV != iterB->second.end(); ++iterV) {
      (*iterV)->addBinEfficiencies(fluctMap, m_bkgdScale[gl_thisIdx++], bkg, gl_addedBins);
      if (gl_thisBin == 0)
        gl_thisBin = (*iterV)->getNXbins() * (*iterV)->getNYbins();
    }
    gl_addedBins += gl_thisBin;
    gl_thisBin = 0;
    gl_thisIdx = 0;
  }

  m_bkgdExclusionSumsReady = false;
} // setNewPointBkgdExpectations

void SigBkgdDist::setPerturbedSignalExpectations(std::vector<double> const& fluctMap,
                                                 std::vector<double>::const_iterator s,
                                                 double* sig) {
  if (!m_signalExclusionSumsReady) {
    prepareSignalExclusionSums();
  }
  int perturbed_s = 0;
  gl_addedBins = 0;
  gl_thisBin = 0;
  CollieDistribution* dist = 0;
  for (iterS = m_SigDist.begin(); iterS != m_SigDist.end(); ++iterS) {
    /*
      if(iterB->first==gl_exclChan){
      iterV=iterB->second.begin();
      CollieDistribution* dist = (*iterV);
      gl_addedBins += dist->getNXbins()*dist->getNYbins();
      gl_thisBin=0;
      continue;
      }
    */
    for (iterD = iterS->second.begin(); iterD != iterS->second.end(); ++iterD) {
      dist = iterD->second;
      perturbed_s = s - fluctMap.begin();
      dist->addBinEfficiencies(perturbed_s, fluctMap[perturbed_s],
                               getSignalScale(iterS->first, iterD->first), sig, gl_addedBins);
      if (gl_thisBin == 0)
        gl_thisBin = dist->getNXbins() * dist->getNYbins();
    } // end of loop over N dists
    gl_addedBins += gl_thisBin;
    gl_thisBin = 0;
  } // end of loop over collection of channels
}

void SigBkgdDist::setPerturbedBkgdExpectations(std::vector<double> const& fluctMap,
                                               std::vector<double>::const_iterator s, double* sig) {

  if (!m_bkgdExclusionSumsReady) {
    prepareBkgdExclusionSums();
  }
  int perturbed_s = 0;
  gl_addedBins = 0;
  gl_thisBin = 0;
  gl_thisIdx = 0;
  for (iterB = m_BkgdDist.begin(); iterB != m_BkgdDist.end(); ++iterB) {
    /*
    if(iterB->first==gl_exclChan){
      iterV=iterB->second.begin();
      CollieDistribution* dist = (*iterV);
      gl_addedBins += dist->getNXbins()*dist->getNYbins();
      gl_thisBin=0;
      continue;
    }
    */
    for (iterV = iterB->second.begin(); iterV != iterB->second.end(); ++iterV) {
      perturbed_s = s - fluctMap.begin();
      (*iterV)->addBinEfficiencies(perturbed_s, fluctMap[perturbed_s], m_bkgdScale[gl_thisIdx++],
                                   sig, gl_addedBins);
      if (gl_thisBin == 0)
        gl_thisBin = (*iterV)->getNXbins() * (*iterV)->getNYbins();
    } // end of loop over dist
    gl_addedBins += gl_thisBin;
    gl_thisBin = 0;
    gl_thisIdx = 0;
  } // end of loop over collection of dists
}

void SigBkgdDist::prepareSignalExclusionSums() {
  for (iterS = m_SigDist.begin(); iterS != m_SigDist.end(); ++iterS) {
    //    if(iterB->first==gl_exclChan) continue;
    for (iterD = iterS->second.begin(); iterD != iterS->second.end(); ++iterD) {
      iterD->second->prepareExclusionSums();
    }
  }
  m_signalExclusionSumsReady = true;
}

void SigBkgdDist::prepareBkgdExclusionSums() {
  for (iterB = m_BkgdDist.begin(); iterB != m_BkgdDist.end(); ++iterB) {
    //    if(iterB->first==gl_exclChan) continue;
    for (iterV = iterB->second.begin(); iterV != iterB->second.end(); ++iterV) {
      CollieDistribution* dist = (*iterV);
      dist->prepareExclusionSums();
    }
  }
  m_bkgdExclusionSumsReady = true;
}

double SigBkgdDist::calculateChi2LLR(bool fitSig, double sigLLR) const {
  (void)sigLLR;
  double llrtotal = 0;
  double mod;
  double dat;
  double llr;
  const uint nb = n_bins;
  const uint db = m_delBins.size();
  for (uint i = 0; i < nb; ++i) {
    //    if(m_data[i]==0 && m_bkgd[i]==0) continue;
    //    if(m_bkgd[i]<1e-5 && m_signal[i]>0) cout << "Bin barf!! " << i << " " << m_data[i] << " " << m_bkgd[i] << " " << m_signal[i] << " " << m_signal[i]/m_bkgd[i] << endl;

    if ((db != 0) && (m_delBinsLookup[i] != 0))
      continue;

    /*
    if(!fitSig){
      double bkp = m_bkgdParent[i];
      if ( bkp > 0 ) {
        if (m_signalParent[i] > bkp * sigLLRcut) continue;
      } else {
        if (sigLLR<10.0) continue;
      }
    }
    */

    mod = m_bkgd[i] + 1e-9;
    if (fitSig)
      mod += m_signal[i];

    dat = m_data[i];

    if (std::isinf(mod) || std::isnan(mod) || mod < 0) {
      cout << "nan barf!!" << endl;
      mod = 1e-9;
    }

    llr = mod - dat;
    if (dat > 0 && mod > 0)
      llr -= dat * log(mod / dat);

    llrtotal += llr;
  }

  return 2.0 * llrtotal;
} // calculateChi2LLR

void SigBkgdDist::excludeBins(vector<uint> delBins) {
  m_delBins = delBins;

  std::fill(m_delBinsLookup.begin(), m_delBinsLookup.end(), 0);
  for (uint i = 0; i != delBins.size(); ++i) {
    m_delBinsLookup[delBins[i]] = 1;
  }
}

void SigBkgdDist::clearExludedBins() {
  m_delBins.clear();
  std::fill(m_delBinsLookup.begin(), m_delBinsLookup.end(), 0);
}

void SigBkgdDist::generateFitMap(TH1* fitValuesIn, TH2* errMatIn, string type) {

  if (m_fitValueHist != NULL) {
    delete m_fitValueHist;
    m_fitValueHist = NULL;
    delete m_errMatHist;
    m_errMatHist = NULL;
  }

  char name[256];
  sprintf(name, "Intl Fit Param Histo %s", type.c_str());
  m_fitValueHist = new TH1D(name, name, getNsyst(), 0, 1);
  sprintf(name, "Intl Error Matrix Histo %s", type.c_str());
  m_errMatHist = new TH2D(name, name, getNsyst(), 0, 1, getNsyst(), 0, 1);

  TAxis* inAx = fitValuesIn->GetXaxis();
  TAxis* outAx = m_fitValueHist->GetXaxis();

  TAxis* twoDAx = m_errMatHist->GetXaxis();
  TAxis* twoDAy = m_errMatHist->GetYaxis();

  int testBins = fitValuesIn->GetNbinsX();

  for (int ss = 0; ss < getNsyst(); ss++) {
    outAx->SetBinLabel(ss + 1, getSystName(ss).c_str());
    twoDAx->SetBinLabel(ss + 1, getSystName(ss).c_str());
    twoDAy->SetBinLabel(ss + 1, getSystName(ss).c_str());

    TString iName(getSystName(ss).c_str());

    bool found1 = false;

    for (int bb = 1; bb <= testBins; bb++) {
      TString oName(inAx->GetBinLabel(bb));

      if (iName == oName) {
        m_fitValueHist->SetBinContent(ss + 1, fitValuesIn->GetBinContent(bb));
        found1 = true;

        for (int sp = 0; sp < getNsyst(); sp++) {
          TString ipName(getSystName(sp).c_str());
          bool found2 = false;

          for (int bp = 1; bp <= testBins; bp++) {

            TString opName(inAx->GetBinLabel(bp));

            if ((ipName == opName) && !found2) {
              //    printf("%s, %s, %s, %s\n",iName.Data(),oName.Data(),ipName.Data(),opName.Data());
              //    printf("%d, %d, %d, %d\n",ss+1,bb,sp+1,bp);
              m_errMatHist->SetBinContent(ss + 1, sp + 1, errMatIn->GetBinContent(bb, bp));
              found2 = true;
            }
          }
          if (!found2) {
            double mv = 0;
            if (ss == sp)
              mv = 1;
            m_errMatHist->SetBinContent(ss + 1, sp + 1, mv);
          }
        }
      }
    }
    if (!found1) {
      m_fitValueHist->SetBinContent(ss + 1, 0);
    }
  }

  m_fitValueHist->Write();
  m_errMatHist->Write();

  return;
}

void SigBkgdDist::generateFitHistos(TH1* signalSF, TH1* fitParams, TH2* errMat, string type) {
  setBaselineModel();

  if ((type != string("TEST")) && (type != string("NULL"))) {
    printf("SigBkgdDist::generateFitHistos, must specify TEST or NULL fit type!\n");
    return;
  }

  if (errMat == NULL || fitParams == NULL || signalSF == NULL) {
    printf("SigBkgdDist::generateFitHistos, NULL input histogram!\n");
    return;
  }

  //  printf("syst sizes: %d, %d\n",fitParams->GetNbinsX(),getNsyst());
  generateFitMap(fitParams, errMat, type);

  signalSF->Write();

  char title[256];

  std::vector<double> fitValues(m_fitValueHist->GetNbinsX(), 0.0);
  for (int i = 1; i <= m_fitValueHist->GetNbinsX(); ++i) {
    //    printf("%s/%s, %f\n",m_fitValueHist->GetXaxis()->GetBinLabel(i),fitParams->GetXaxis()->GetBinLabel(i),m_fitValueHist->GetBinContent(i)/(1e-9+fitParams->GetBinContent(i)));
    fitValues[i - 1] = m_fitValueHist->GetBinContent(i);
  }
  //fit templates
  CollieDistribution* dist = 0;
  int maxY = 0;
  int maxX = 0;
  int binY = 0;
  int binX = 0;

  int nBinsCommon = 0;
  double rminCommon = 0;
  double rmaxCommon = 0;
  bool commonBins = true;

  for (iterD = m_DataDist.begin(); iterD != m_DataDist.end(); ++iterD) { //loop over channels
    dist = iterD->second;
    maxY = (dist->getNYbins() < 1) ? 1 : dist->getNYbins();
    maxX = dist->getNXbins();

    if (maxY > 1)
      commonBins = false;
    if (nBinsCommon == 0) {
      nBinsCommon = maxX;
      rminCommon = dist->getMinX();
      rmaxCommon = dist->getMaxX();
    } else if (maxX != nBinsCommon)
      commonBins = false;

    sprintf(title, "%s Fit, Data %s", type.c_str(), iterD->first.c_str());
    TH1D* data = new TH1D(title, title, maxX * maxY, dist->getMinX(), dist->getMaxX());
    data->Sumw2();
    for (binY = 0; binY < maxY; ++binY)
      for (binX = 0; binX < maxX; ++binX) {
        data->SetBinContent((binX + maxX * binY) + 1, dist->getEfficiency(binX, binY));
        data->SetBinError((binX + maxX * binY) + 1,
                          sqrt(data->GetBinContent(binX + maxX * binY + 1)));
      }
  }
  if (nBinsCommon == 0)
    nBinsCommon = 1;

  for (iterS = m_SigDist.begin(); iterS != m_SigDist.end(); ++iterS) { //loop over channels

    //query the data for the number of bins for this channel...
    dist = m_DataDist[iterS->first];
    maxY = (dist->getNYbins() < 1) ? 1 : dist->getNYbins();
    maxX = dist->getNXbins();

    for (iterD = iterS->second.begin(); iterD != iterS->second.end(); ++iterD) { //loop over signals
      dist = iterD->second;
      sprintf(title, "%s Fit, Sig %s %s", type.c_str(), iterD->first.c_str(), iterS->first.c_str());
      TH1D* sig = new TH1D(title, title, maxX * maxY, dist->getMinX(), dist->getMaxX());
      sig->Sumw2();
      for (binY = 0; binY < maxY; ++binY)
        for (binX = 0; binX < maxX; ++binX) {
          double ss = dist->getEfficiencyVaried(binX, binY, &fitValues[0]);
          if (ss < 0)
            ss = 0;
          sig->SetBinContent((binX + maxX * binY) + 1, ss);
          sig->SetBinError((binX + maxX * binY) + 1, dist->getBinStatErr(binX, binY));
        }
      if (signalSF->GetNbinsX() > 1) {
        int ibin = signalSF->GetXaxis()->FindBin(iterD->first.c_str());
        sig->Scale(signalSF->GetBinContent(ibin));
      } else
        sig->Scale(signalSF->GetBinContent(1));
    }
  }

  for (iterB = m_BkgdDist.begin(); iterB != m_BkgdDist.end(); ++iterB) { //loop over channels

    //query the data for the number of bins for this channel...
    dist = m_DataDist[iterB->first];
    maxY = (dist->getNYbins() < 1) ? 1 : dist->getNYbins();
    maxX = dist->getNXbins();

    sprintf(title, "%s Fit, Sum Bkgd %s", type.c_str(), iterB->first.c_str());
    TH1D* bkgT = new TH1D(title, title, maxX * maxY, dist->getMinX(), dist->getMaxX());
    bkgT->Sumw2();
    int idx = 0;
    for (iterV = iterB->second.begin(); iterV != iterB->second.end(); ++iterV) {
      dist = (*iterV);
      sprintf(title, "%s Fit, Bkgd %d %s", type.c_str(), idx, iterB->first.c_str());
      TH1D* bkg = new TH1D(title, title, maxX * maxY, dist->getMinX(), dist->getMaxX());
      bkg->Sumw2();

      for (binY = 0; binY < maxY; ++binY)
        for (binX = 0; binX < maxX; ++binX) {
          double bb = dist->getEfficiencyVaried(binX, binY, &fitValues[0]);
          if (bb < 0)
            bb = 0;
          bkg->SetBinContent((binX + maxX * binY) + 1, bb);
          bkg->SetBinError((binX + maxX * binY) + 1, dist->getBinStatErr(binX, binY));
        }
      bkgT->Add(bkg);
      idx++;
    }
  }

  if (type == string("TEST")) {

    TAxis* ax = m_fitValueHist->GetXaxis();
    TAxis* emx = m_errMatHist->GetXaxis();
    TAxis* emy = m_errMatHist->GetYaxis();

    //Calculate total uncertainty based on the covariance matrix from MINUIT
    const int nSystI = fitValues.size();

    //Total uncertainties with and without error matrix applied
    std::vector<double> params(nSystI, 0.0); // USED TO BE: double[..]
    int dbins = 0;

    double nomSystPTS_emat[nSystI][nBinsCommon];
    double nomSystMTS_emat[nSystI][nBinsCommon];
    double nomSystPTS[nSystI][nBinsCommon];
    double nomSystMTS[nSystI][nBinsCommon];

    for (int i = 0; i < nSystI; i++) {
      for (int j = 0; j < nBinsCommon; j++) {
        nomSystPTS_emat[i][j] = 0;
        nomSystMTS_emat[i][j] = 0;
        nomSystPTS[i][j] = 0;
        nomSystMTS[i][j] = 0;
      }
    }

    const vector<string> channelNames = getChannelNames();
    for (uint chanIdx = 0; chanIdx < channelNames.size(); chanIdx++) {
      const string& channelName = channelNames[chanIdx];
      const CollieDistribution* cs = getSigDist(channelName, 0);
      if (cs == NULL) {
        printf("SigBkgdDist::generateFitHistos, this CollieDistribution is NULL!\n");
        return;
      }

      int tbins = cs->getNXbins() * cs->getNYbins();

      //Do the total for the channel first...
      double nomSystPT[nSystI][tbins];
      double nomSystMT[nSystI][tbins];

      double nomSystPT_emat[nSystI][tbins];
      double nomSystMT_emat[nSystI][tbins];

      int sIdx = 0;

      for (int s = 0; s < nSystI; s++) {
        //SBD index of the sth systematic in the error matrix
        sIdx = getSystIndex(ax->GetBinLabel(s + 1));
        if (sIdx < 0) {
          printf("SigBkgdDist::generateFitHistos  Could not identify systematic %d: %s\n", s + 1,
                 ax->GetBinLabel(s + 1));
          continue;
        }
        for (int st = 0; st < nSystI; st++)
          params[st] = 0;

        //Get upward fluctuated bin values for each systematic
        setBaselineModel();
        if (getSystFitFlag(sIdx))
          params[sIdx] = 1.0;
        fluctuate(params);
        for (int bx = 0; bx < cs->getNXbins(); bx++) {
          for (int by = 0; by < cs->getNYbins(); by++) {
            int bb = bx + cs->getNXbins() * by;
            nomSystPT[sIdx][bb] = signal(bb + dbins) + bkgd(bb + dbins);
            if (commonBins)
              nomSystPTS[sIdx][bb] += signal(bb + dbins) + bkgd(bb + dbins);
          }
        }

        //Get downward fluctuated bin values for each systematic
        setBaselineModel();
        params[sIdx] = -1.0;
        fluctuate(params);

        for (int bx = 0; bx < cs->getNXbins(); bx++) {
          for (int by = 0; by < cs->getNYbins(); by++) {
            int bb = bx + cs->getNXbins() * by;
            nomSystMT[sIdx][bb] = signal(bb + dbins) + bkgd(bb + dbins);
            if (commonBins)
              nomSystMTS[sIdx][bb] += signal(bb + dbins) + bkgd(bb + dbins);
          }
        }

        //Subtract baseline background
        params[sIdx] = 0;
        setBaselineModel();
        fluctuate(params);
        setBaselineModel();
        for (int bx = 0; bx < cs->getNXbins(); bx++) {
          for (int by = 0; by < cs->getNYbins(); by++) {
            int bb = bx + cs->getNXbins() * by;
            nomSystPT[sIdx][bb] -= (signal(bb + dbins) + bkgd(bb + dbins));
            nomSystMT[sIdx][bb] -= (signal(bb + dbins) + bkgd(bb + dbins));

            if (commonBins)
              nomSystPTS[sIdx][bb] -= (signal(bb + dbins) + bkgd(bb + dbins));
            if (commonBins)
              nomSystMTS[sIdx][bb] -= (signal(bb + dbins) + bkgd(bb + dbins));

            nomSystPT_emat[sIdx][bb] = 0;
            nomSystMT_emat[sIdx][bb] = 0;
          }
        }
      }
      dbins += cs->getNXbins() * cs->getNYbins();

      //factor in error matrix & generate sigma^2 for each systematic
      int sIdx1 = 0;
      int sIdx2 = 0;
      for (int s1 = 0; s1 < nSystI; s1++) { //rows
        sIdx1 = getSystIndex(ax->GetBinLabel(s1 + 1));
        for (int s2 = 0; s2 < nSystI; s2++) { //columns
          sIdx2 = getSystIndex(ax->GetBinLabel(s2 + 1));
          if (sIdx1 < 0 || sIdx2 < 0) {
            //      printf("Failed to find syst pair: %s/%s\n",ax->GetBinLabel(s1+1),ax->GetBinLabel(s2+1));
            continue;
          }

          for (int bx = 0; bx < cs->getNXbins(); bx++) {
            for (int by = 0; by < cs->getNYbins(); by++) {
              int bb = bx + cs->getNXbins() * by;

              int b1 = emx->FindBin(ax->GetBinLabel(s1 + 1));
              int b2 = emy->FindBin(ax->GetBinLabel(s2 + 1));
              double emv = m_errMatHist->GetBinContent(b1, b2);
              if (b1 > emx->GetNbins() || b1 < 1 || b2 > emy->GetNbins() || b2 < 1) {
                if (s1 == s2)
                  emv = 1;
                else
                  emv = 0;
              }
              nomSystPT_emat[sIdx1][bb] += emv * nomSystPT[sIdx1][bb] * nomSystPT[sIdx2][bb];
              nomSystMT_emat[sIdx1][bb] += emv * nomSystMT[sIdx1][bb] * nomSystMT[sIdx2][bb];
            }
          }
        }
      }

      sprintf(title, "Positive Error Matrix Systematics: %s, Total", channelName.c_str());
      TH1D* hpt = new TH1D(title, title, tbins, cs->getMinX(), cs->getMaxX());
      sprintf(title, "Negative Error Matrix Systematics: %s, Total", channelName.c_str());
      TH1D* hmt = new TH1D(title, title, tbins, cs->getMinX(), cs->getMaxX());

      //sum sigma^2 for all systematics
      for (int s = 0; s < nSystI; s++) {
        for (int b = 0; b < tbins; b++) {
          hpt->AddBinContent(b + 1, nomSystPT_emat[s][b]);
          hmt->AddBinContent(b + 1, nomSystMT_emat[s][b]);
        }
      }

      //take sqrt for each bin
      for (int b = 0; b < tbins; b++) {
        //  assert(hpt->GetBinContent(b+1)>=0);
        if (hpt->GetBinContent(b + 1) < 0) {
          printf("Negative Err Matrix Propagation! (Bin %d content = %f)\n", b + 1,
                 hpt->GetBinContent(b + 1));
        }
        if (hmt->GetBinContent(b + 1) < 0) {
          printf("Negative Err Matrix Propagation! (Bin %d content = %f)\n", b + 1,
                 hmt->GetBinContent(b + 1));
        }
        hpt->SetBinContent(b + 1, sqrt(fabs(hpt->GetBinContent(b + 1))));
        hmt->SetBinContent(b + 1, sqrt(fabs(hmt->GetBinContent(b + 1))));
      }

      //Now do the individual sources...
      uint nsigd = getNsigDist(channelName);
      uint nbkgd = getNbkgdDist(channelName);
      for (uint d = 0; d < (nsigd + nbkgd); d++) {

        const CollieDistribution* cd =
            (d < nsigd) ? getSigDist(channelName, d) : getBkgdDist(channelName, d - nsigd);
        if (cd == NULL) {
          printf("SigBkgdDist::generateFitHistos NULL Distribution!\n");
          continue;
        }

        int tbins = cd->getNXbins() * cd->getNYbins();

        //Because the error matrix is indexed differently than the SBD, we need to
        // make sure we index correctly.
        double nomSystP[nSystI][tbins];
        double nomSystM[nSystI][tbins];
        double nomSystP_emat[nSystI][tbins];
        double nomSystM_emat[nSystI][tbins];

        int sIdx = 0;
        for (int s = 0; s < nSystI; s++) {

          //SBD index of the sth systematic in the error matrix
          sIdx = getSystIndex(ax->GetBinLabel(s + 1));
          if (sIdx < 0)
            continue;

          for (int st = 0; st < nSystI; st++)
            params[st] = 0;

          //Get upward fluctuated bin values for each systematic
          setBaselineModel();
          if (getSystFitFlag(sIdx))
            params[sIdx] = 1.0;
          fluctuate(params);
          for (int bx = 0; bx < cs->getNXbins(); bx++) {
            for (int by = 0; by < cs->getNYbins(); by++) {
              int bb = bx + cs->getNXbins() * by;
              nomSystP[sIdx][bb] = cd->getEfficiencyVaried(bx, by, &params[0]);
            }
          }

          //Get downward fluctuated bin values for each systematic
          setBaselineModel();
          params[sIdx] = -1.0;
          fluctuate(params);
          for (int bx = 0; bx < cs->getNXbins(); bx++) {
            for (int by = 0; by < cs->getNYbins(); by++) {
              int bb = bx + cs->getNXbins() * by;
              nomSystM[sIdx][bb] = cd->getEfficiencyVaried(bx, by, &params[0]);
            }
          }

          //Subtract baseline background
          params[sIdx] = 0;
          setBaselineModel();
          setBaselineModel();
          for (int bx = 0; bx < cs->getNXbins(); bx++) {
            for (int by = 0; by < cs->getNYbins(); by++) {
              int bb = bx + cs->getNXbins() * by;
              nomSystP[sIdx][bb] -= cd->getEfficiencyVaried(bx, by, &params[0]);
              nomSystM[sIdx][bb] -= cd->getEfficiencyVaried(bx, by, &params[0]);
              nomSystP_emat[sIdx][bb] = 0;
              nomSystM_emat[sIdx][bb] = 0;
            }
          }
        }

        //factor in error matrix & generate sigma^2 for each systematic
        for (int s1 = 0; s1 < nSystI; s1++) { //rows
          int sIdx1 = getSystIndex(ax->GetBinLabel(s1 + 1));

          for (int s2 = 0; s2 < nSystI; s2++) { //columns
            int sIdx2 = getSystIndex(ax->GetBinLabel(s2 + 1));

            if (sIdx1 < 0 || sIdx2 < 0) {
              //        printf("SigBkgdDist::generateFitHistos  Negative syst index: %d, %d\n",sIdx1,sIdx2);
              continue;
            }

            for (int bx = 0; bx < cs->getNXbins(); bx++) {
              for (int by = 0; by < cs->getNYbins(); by++) {
                int bb = bx + cs->getNXbins() * by;
                int b1 = emx->FindBin(ax->GetBinLabel(s1 + 1));
                int b2 = emy->FindBin(ax->GetBinLabel(s2 + 1));
                double emv = m_errMatHist->GetBinContent(b1, b2);

                if (b1 > emx->GetNbins() || b1 < 1 || b2 > emy->GetNbins() || b2 < 1) {
                  if (s1 == s2)
                    emv = 1;
                  else
                    emv = 0;
                }
                nomSystP_emat[sIdx1][bb] += emv * nomSystP[sIdx1][bb] * nomSystP[sIdx2][bb];
                nomSystM_emat[sIdx1][bb] += emv * nomSystM[sIdx1][bb] * nomSystM[sIdx2][bb];
              }
            }
          }
        }

        if (d < nsigd)
          sprintf(title, "Positive Error Matrix Systematics: %s, Sig %u", channelName.c_str(), d);
        else
          sprintf(title, "Positive Error Matrix Systematics: %s, Bkgd %u", channelName.c_str(),
                  d - nsigd);
        TH1D* hp = new TH1D(title, title, tbins, cd->getMinX(), cd->getMaxX());
        if (d < nsigd)
          sprintf(title, "Negative Error Matrix Systematics: %s, Sig %u", channelName.c_str(), d);
        else
          sprintf(title, "Negative Error Matrix Systematics: %s, Bkgd %u", channelName.c_str(),
                  d - nsigd);
        TH1D* hm = new TH1D(title, title, tbins, cd->getMinX(), cd->getMaxX());

        //sum sigma^2 for all systematics
        for (int s = 0; s < nSystI; s++) {
          for (int b = 0; b < tbins; b++) {
            hp->AddBinContent(b + 1, nomSystP_emat[s][b]);
            hm->AddBinContent(b + 1, nomSystM_emat[s][b]);
          }
        }

        //take sqrt for each bin
        for (int b = 0; b < tbins; b++) {
          hp->SetBinContent(b + 1, sqrt(fabs(hp->GetBinContent(b + 1))));
          hm->SetBinContent(b + 1, sqrt(fabs(hm->GetBinContent(b + 1))));
        }
      }
    }

    if (commonBins) {
      int sIdx1 = 0;
      int sIdx2 = 0;
      for (int s1 = 0; s1 < nSystI; s1++) { //rows
        sIdx1 = getSystIndex(ax->GetBinLabel(s1 + 1));
        for (int s2 = 0; s2 < nSystI; s2++) { //columns
          sIdx2 = getSystIndex(ax->GetBinLabel(s2 + 1));
          if (sIdx1 < 0 || sIdx2 < 0)
            continue;

          for (int bb = 0; bb < nBinsCommon; bb++) {
            int b1 = emx->FindBin(ax->GetBinLabel(s1 + 1));
            int b2 = emy->FindBin(ax->GetBinLabel(s2 + 1));
            double emv = m_errMatHist->GetBinContent(b1, b2);
            if (b1 > emx->GetNbins() || b1 < 1 || b2 > emy->GetNbins() || b2 < 1) {
              if (s1 == s2)
                emv = 1;
              else
                emv = 0;
            }
            //      printf("Total: %f\n",emv);
            nomSystPTS_emat[sIdx1][bb] += emv * nomSystPTS[sIdx1][bb] * nomSystPTS[sIdx2][bb];
            nomSystMTS_emat[sIdx1][bb] += emv * nomSystMTS[sIdx1][bb] * nomSystMTS[sIdx2][bb];
          }
        }
      }

      sprintf(title, "Positive Error Matrix Systematics, Channel Sum");
      TH1D* hpts = new TH1D(title, title, nBinsCommon, rminCommon, rmaxCommon);
      sprintf(title, "Negative Error Matrix Systematics, Channel Sum");
      TH1D* hmts = new TH1D(title, title, nBinsCommon, rminCommon, rmaxCommon);

      sprintf(title, "Positive Pre-Fit Systematics, Channel Sum");
      TH1D* hptst = new TH1D(title, title, nBinsCommon, rminCommon, rmaxCommon);
      sprintf(title, "Negative Pre-Fit Systematics, Channel Sum");
      TH1D* hmtst = new TH1D(title, title, nBinsCommon, rminCommon, rmaxCommon);

      //sum sigma^2 for all systematics
      for (int s = 0; s < nSystI; s++) {
        for (int b = 0; b < nBinsCommon; b++) {
          hpts->AddBinContent(b + 1, nomSystPTS_emat[s][b]);
          hmts->AddBinContent(b + 1, nomSystMTS_emat[s][b]);
          hptst->AddBinContent(b + 1, nomSystPTS[s][b] * nomSystPTS[s][b]);
          hmtst->AddBinContent(b + 1, nomSystMTS[s][b] * nomSystMTS[s][b]);
        }
      }

      //take sqrt for each bin
      for (int b = 0; b < nBinsCommon; b++) {
        //  assert(hpt->GetBinContent(b+1)>=0);
        if (hpts->GetBinContent(b + 1) < 0)
          printf("Negative Err Matrix Propagation! (Bin %d content = %f)\n", b + 1,
                 hpts->GetBinContent(b + 1));
        if (hmts->GetBinContent(b + 1) < 0)
          printf("Negative Err Matrix Propagation! (Bin %d content = %f)\n", b + 1,
                 hmts->GetBinContent(b + 1));
        hpts->SetBinContent(b + 1, sqrt(fabs(hpts->GetBinContent(b + 1))));
        hmts->SetBinContent(b + 1, sqrt(fabs(hmts->GetBinContent(b + 1))));

        hptst->SetBinContent(b + 1, sqrt(fabs(hptst->GetBinContent(b + 1))));
        hmtst->SetBinContent(b + 1, sqrt(fabs(hmtst->GetBinContent(b + 1))));
      }
    }
  }
}
