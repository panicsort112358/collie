#include "ProfileLH.hh"
#include "ProfileLikelihoodEngine.hh"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory>
#include <stdexcept>
#include <vector>
#include "TrialPoint.hh"

///Standard constructor
ProfileLH::ProfileLH() {
  setup();
}

///Constructor with SBD initializer
ProfileLH::ProfileLH(SigBkgdDist* asbd) {
  setup();
  if (asbd == NULL) {
    printf("ProfileLH::ProfileLH, the input SBD is NULL!\n");
    return;
  }
  setModel(asbd);
}

ProfileLH::~ProfileLH() = default;

void ProfileLH::setup() {
  pf_sbd = NULL;
  pf_params.clear(); // USED TO BE:  pf_params = NULL;
  pf_index.clear();
  pf_paramErrs.clear();
  pf_nParm = 0;
  pf_iter = 0;
  pf_fitSig = true;
  pf_sigFloat = false;
  pf_sigLLR = 1e6;
  pf_scl = 1.0;

  m_init = false;
  m_fitTest = false;
  m_fitSig = true;
  m_verbose = false;
  m_sigFloat = false;
  m_sigScale = 1.0;
  m_sigScaleNom = 1.0;
  m_sigScaleErr = 0.0;
  m_sigScaleErrP = 0.0;
  m_sigScaleErrM = 0.0;
  m_sigScaleParab = 0.0;
  m_sigScaleGlobCC = 0.0;
  m_sigLLR = 1e6;
  m_chi2min = -1.0;
  m_edm = -1.0;
  m_lun = 99;
  m_status = -1;
  m_maxIter = 250000;

  m_systNamesFitMap.clear();
  m_systNamesFitBkgdMap.clear();
  m_systNamesFit.clear();
  m_systNamesFitBkgd.clear();
  m_systNames.clear();
  m_systValues.clear();
  m_systErrors.clear();
  m_systErrorsUp.clear();
  m_systErrorsDn.clear();
  m_systErrorsParab.clear();
  m_systErrorsGlobCC.clear();

  m_condI.sig = 0;
  m_condI.bkgd = 0;
  m_condI.data = 0;
  m_condI.chi2 = 0;
  m_condI.chi2_bins = 0;
  m_condI.chi2_syst = 0;

  m_condF.sig = 0;
  m_condF.bkgd = 0;
  m_condF.data = 0;
  m_condF.chi2 = 0;
  m_condF.chi2_bins = 0;
  m_condF.chi2_syst = 0;
}

void ProfileLH::setModel(SigBkgdDist* asbd) {
  pf_sbd = asbd;
  pf_fitSig = m_fitSig;
  pf_sigLLR = m_sigLLR;
  fillSyst();
}

int ProfileLH::getNiterations() {
  return pf_iter;
}

void ProfileLH::getMNerrs() {
  if (!m_init || !m_engine) {
    printf("ProfileLH::getMNerrs, not initialized\n");
    return;
  }

  m_systErrorsUp.clear();
  m_systErrorsDn.clear();
  m_systErrorsParab.clear();
  m_systErrorsGlobCC.clear();

  double errP = 0.0;
  double errM = 0.0;
  double eparab = 0.0;
  double globcc = 0.0;
  for (uint num = 0; num < pf_nParm; ++num) {
    if (!m_engine->minosError(num, errP, errM, eparab, globcc)) {
      eparab = num < m_systErrors.size() ? m_systErrors[num] : 0.0;
      errP = eparab;
      errM = -eparab;
      globcc = 0.0;
    }
    m_systErrorsUp.push_back(errP);
    m_systErrorsDn.push_back(errM);
    m_systErrorsParab.push_back(eparab);
    m_systErrorsGlobCC.push_back(globcc);
  }

  if (m_sigFloat) {
    if (!m_engine->minosError(pf_nParm, errP, errM, eparab, globcc)) {
      eparab = m_sigScaleErr;
      errP = eparab;
      errM = -eparab;
      globcc = 0.0;
    }
    m_sigScaleErrP = errP;
    m_sigScaleErrM = errM;
    m_sigScaleParab = eparab;
    m_sigScaleGlobCC = globcc;
  }

  return;
}

bool ProfileLH::getEMat(double* emat) {

  if (!m_init || !m_engine) {
    printf("ProfileLH::getEMat, not initialized\n");
    return false;
  }
  getMNerrs();
  return m_engine->covariance(emat, pf_nParm);
}

void ProfileLH::fillSyst() {
  m_systNames.clear();
  m_systNamesFit.clear();
  m_systNamesFitBkgd.clear();
  m_systNamesFitMap.clear();
  m_systNamesFitBkgdMap.clear();

  if (pf_sbd == NULL) {
    printf(
        "ProfileLH::fillSyst, You're trying to fit a profile with a NULL signal/background distribution!!\n");
    return;
  }

  //Make sure all the systematics are properly set
  if (!pf_sbd->isInit()) {
    pf_sbd->varySystematics();
    pf_sbd->setBaselineModel();
  }

  vector<string> chanNames = pf_sbd->getChannelNames();
  const CollieDistribution* dist;
  string systName = "";
  map<string, int> m_systMap;

  for (vector<string>::iterator its = chanNames.begin(); its != chanNames.end(); ++its) {
    for (uint b = 0; b < pf_sbd->getNbkgdDist(*its); ++b) {

      dist = pf_sbd->getBkgdDist(*its, b);
      if (dist != NULL) {

        for (int s = 0; s < dist->getNsystematics(); ++s) {
          systName = dist->getSystName(s);
          if (m_systMap.find(systName) != m_systMap.end())
            continue;
          else {
            m_systNames.push_back(systName);
            if (pf_sbd->getSystFitFlag(systName)) {
              m_systNamesFitMap[systName] = pf_sbd->getSystIndex(systName);
              m_systNamesFitBkgdMap[systName] = pf_sbd->getSystIndex(systName);
              m_systNamesFit.push_back(systName);
              m_systNamesFitBkgd.push_back(systName);
            } else {
              if (m_verbose)
                printf("ProfileLH:: Not Fitting %s\n", dist->getSystName(s).c_str());
            }
            m_systMap[systName] = pf_sbd->getSystIndex(systName);
          }
        }
      }
    }
  }

  for (vector<string>::iterator its = chanNames.begin(); its != chanNames.end(); ++its) {
    for (uint b = 0; b < pf_sbd->getNsigDist(*its); ++b) {
      //      dist = pf_sbd->getSigDist(*its,pf_sbd->getSigName(b));
      dist = pf_sbd->getSigDist(*its, b);
      if (dist != NULL) {
        for (int s = 0; s < dist->getNsystematics(); ++s) {
          systName = dist->getSystName(s);

          if (m_systMap.find(systName) != m_systMap.end())
            continue;
          else {
            m_systNames.push_back(systName);
            if (pf_sbd->getSystFitFlag(systName) && pf_fitSig) {
              m_systNamesFitMap[systName] = pf_sbd->getSystIndex(systName);
              m_systNamesFit.push_back(systName);
            } else {
              if (m_verbose)
                printf("ProfileLH:: Not Fitting %s\n", systName.c_str());
            }
            m_systMap[systName] = pf_sbd->getSystIndex(systName);
          }
        }
      } // else printf("This sig dist is NULL???\n");
    }
  }

  //our adjustable copy of the syst params
  pf_params = std::vector<double>(m_systNames.size(), 0.0);
  pf_paramErrs = std::vector<double>(m_systNames.size(), 0.0);
  //  m_systErrors.clear();
  //  m_systValues.clear();

  //collect existing parameters
  for (unsigned int i = 0; i < m_systNames.size(); ++i)
    pf_params[i] = pf_sbd->getSystFluctValue(i);

  initFitParams();

  return;
}

void ProfileLH::initFitParams() {
  //Make a storage container with only the values we're fitting...
  //   either background-only or sig+bkgd
  //   for bkgd-only, make index so we know where it fits in master copy.
  //   Allow individual fit parameters to be excluded by the user.

  pf_index.clear();
  map<string, int>::iterator im;

  if (m_fitSig) {
    pf_nParm = m_systNamesFit.size();
    pf_index = std::vector<uint>(pf_nParm, 0);

    //now fill in the syst indices
    for (uint s = 0; s < m_systNamesFit.size(); ++s) {
      im = m_systNamesFitMap.find(m_systNamesFit[s]);
      pf_index[s] = im->second;
    }
  } else {
    pf_nParm = m_systNamesFitBkgd.size();
    pf_index = std::vector<uint>(pf_nParm, 0);

    //now fill in the syst indices
    for (uint s = 0; s < m_systNamesFitBkgd.size(); ++s) {
      im = m_systNamesFitBkgdMap.find(m_systNamesFitBkgd[s]);
      pf_index[s] = im->second;
    }
  }

  return;
}

bool ProfileLH::getMinosError(uint idx, double& errP, double& errM, double& parab, double& globCC) {
  if (m_systErrorsUp.size() == 0)
    getMNerrs();
  if (m_systErrorsUp.size() == 0) {
    printf("ProfileLH::getMinosError, getMinosError failed...\n");
    return false;
  } //failed..
  if (idx >= m_systErrorsUp.size() || idx >= m_systErrorsDn.size() ||
      idx >= m_systErrorsParab.size() || idx >= m_systErrorsGlobCC.size())
    return false;

  errP = m_systErrorsUp[idx];
  errM = m_systErrorsDn[idx];
  parab = m_systErrorsParab[idx];
  globCC = m_systErrorsGlobCC[idx];
  return true;
}

bool ProfileLH::getMinosErrXsec(double& errP, double& errM, double& parab, double& globCC) {
  if (m_systErrorsUp.size() == 0)
    getMNerrs();
  if (m_systErrorsUp.size() == 0) {
    printf("ProfileLH::getMinosErrXsec, getMinosErrXsec failed...\n");
    return false;
  } //failed..

  errP = m_sigScaleErrP;
  errM = m_sigScaleErrM;
  parab = m_sigScaleParab;
  globCC = m_sigScaleGlobCC;
  return true;
}

double ProfileLH::getParamVal(uint p) {
  if (p >= pf_params.size())
    return 0;
  return pf_params[p];
}

double ProfileLH::getParamError(uint p) {
  if (p >= pf_paramErrs.size())
    return 0;
  return pf_paramErrs[p];
}

int ProfileLH::getFitSystIndex(string name) {
  map<string, int>::iterator iter;
  if (m_fitSig) {
    for (iter = m_systNamesFitMap.begin(); iter != m_systNamesFitMap.end(); ++iter)
      if (iter->first == name)
        return iter->second;
  } else {
    for (iter = m_systNamesFitBkgdMap.begin(); iter != m_systNamesFitBkgdMap.end(); ++iter)
      if (iter->first == name)
        return iter->second;
  }
  return -1;
}

double ProfileLH::getFitParamErr(string name) const {

  if (m_fitSig) {
    for (uint s = 0; s < m_systNamesFit.size(); ++s) {
      if (m_systNamesFit[s] == name)
        return m_systErrors[s];
    }
  } else {
    for (uint s = 0; s < m_systNamesFitBkgd.size(); ++s) {
      if (m_systNamesFitBkgd[s] == name)
        return m_systErrors[s];
    }
  }
  return 1;
}

double ProfileLH::getFitParamVal(string name) const {
  if (m_fitSig) {
    for (uint s = 0; s < m_systNamesFit.size(); ++s) {
      if (m_systNamesFit[s] == name)
        return m_systValues[s];
    }
  } else {
    for (uint s = 0; s < m_systNamesFitBkgd.size(); ++s) {
      if (m_systNamesFitBkgd[s] == name)
        return m_systValues[s];
    }
  }
  return 0;
}

string ProfileLH::getFitSystName(uint idx) {

  if (idx < 0)
    return "";
  if (m_fitSig) {
    if (idx >= m_systNamesFit.size())
      return "";
    return m_systNamesFit[idx];
  } else {
    if (idx >= m_systNamesFitBkgd.size())
      return "";
    return m_systNamesFitBkgd[idx];
  }
  return "";
}

double ProfileLH::getSystSum2() {
  double tot = 0;
  for (uint i = 0; i < m_systNames.size(); ++i) {
    if (pf_sbd->getFloatFlag(i))
      continue;
    tot += pf_sbd->getSystFluctValue(i) * pf_sbd->getSystFluctValue(i);
  }
  return tot;
}

void ProfileLH::fitProfile() {

  if (pf_sbd == NULL) {
    printf(
        "ProfileLH::fitProfile, You're trying to fit a profile with a NULL signal/background distribution!!\n");
    return;
  }

  //Initialize our parameters
  pf_iter = 0;
  m_status = -1;
  m_chi2min = -1.0;
  m_edm = -1.0;
  m_sigScale = m_sigScaleNom;
  m_sigScaleErr = 0.0;

  //check to see if the fitting model has changed...
  if ((pf_fitSig != m_fitSig) || (m_sigFloat != pf_sigFloat)) {
    pf_fitSig = m_fitSig;
    pf_sigFloat = m_sigFloat;
    fillSyst();
    m_init = false;
  }

  //M Owen 16-3-08
  //Only want to vary the systematics in each fit iterations
  //do not want to vary with the stat uncertainties.
  bool histostat = pf_sbd->usingHistoStats();
  pf_sbd->useHistoStats(false);

  if (m_verbose)
    printf("ProfileLH: %zu systematics, %zu fit params\n", m_systNames.size(),
           m_systNamesFit.size());

  //collect existing parameters
  for (uint i = 0; i < m_systNames.size(); ++i) {
    pf_params[i] = pf_sbd->getSystFluctValue(i);
    pf_paramErrs[i] = 0;
  }

  //Record initial conditions...
  m_condI.sig = pf_sbd->totSignal();
  m_condI.bkgd = pf_sbd->totBkgd();
  m_condI.data = pf_sbd->totData();
  m_condI.chi2 =
      pf_sbd->calculateChi2LLR(pf_fitSig, pf_sigLLR) + pf_sbd->calculateSystDiff(pf_params);

  m_condI.chi2_bins = pf_sbd->calculateChi2LLR(pf_fitSig, pf_sigLLR);
  m_condI.chi2_syst = pf_sbd->calculateSystDiff(pf_params);

  // Establish the memory of the last Minuit trial point, so that
  // sbd->fillArrays() can operate more efficiently
  pf_lastMinuitTrialPoint.clear();

  std::vector<collie::ProfileParameter> parameters;
  parameters.reserve(pf_nParm + (m_sigFloat ? 1U : 0U));
  for (uint num = 0; num < pf_nParm; ++num) {
    collie::ProfileParameter parameter;
    parameter.name = "Par" + std::to_string(num + 1);
    parameter.initial = m_fitTest || m_sigFloat                     ? 0.0
                        : m_init && m_systValues.size() == pf_nParm ? m_systValues[num]
                                                                    : pf_params[pf_index[num]];
    parameter.step = 0.1;
    parameter.lower = -4.5;
    parameter.upper = 4.5;
    parameter.limited = true;
    parameters.push_back(parameter);
  }
  if (m_sigFloat) {
    collie::ProfileParameter parameter;
    parameter.name = "SigRate";
    parameter.initial = m_sigScaleNom;
    parameter.step = 0.1;
    parameter.lower = 1e-5;
    parameter.upper = 1e5;
    parameter.limited = true;
    parameters.push_back(parameter);
  }

  collie::ProfileFitResult result;
  if (parameters.empty()) {
    m_engine.reset();
    pf_iter = 1;
    m_chi2min = evaluate(nullptr);
    m_edm = 0.0;
    m_status = 3;
    m_init = true;
    result.converged = true;
  } else {
    try {
      m_engine = std::make_unique<collie::ProfileLikelihoodEngine>();
      result = m_engine->fit(
          parameters, [this](const double* values) { return evaluate(values); },
          static_cast<std::size_t>(m_maxIter), 1.0, m_verbose ? 1 : 0);
      pf_iter = static_cast<uint>(result.evaluations);
      m_chi2min = result.minimum;
      m_edm = result.edm;
      m_status = result.covarianceStatus;
      m_init = result.converged;
    } catch (const std::exception& error) {
      fprintf(stderr, "ProfileLH: Minuit2 setup failed: %s\n", error.what());
      m_status = -1;
      m_init = false;
      pf_sbd->useHistoStats(histostat);
      return;
    }
  }

  if (m_status != 3 && m_verbose)
    printf("ProfileLH: Warning! Minuit2 convergence failing (%d)!\n", m_status);
  if (m_verbose)
    printf("ProfileLH: %u Minuit2 evaluations\n", pf_iter);
  if (pf_iter > 0.95 * m_maxIter)
    printf("ProfileLH: Minuit2 evaluation warning: %u evaluations\n", pf_iter);

  m_systErrors.clear();
  m_systValues.clear();
  m_systErrorsUp.clear();
  m_systErrorsDn.clear();
  m_systErrorsParab.clear();
  m_systErrorsGlobCC.clear();

  for (uint num = 0; num < pf_nParm; ++num) {
    const double val = result.values[num];
    const double error = result.errors[num];
    m_systErrors.push_back(error);
    m_systValues.push_back(val);
    pf_params[pf_index[num]] = val;
    pf_paramErrs[pf_index[num]] = error;
  }

  if (m_sigFloat) {
    m_sigScale = result.values[pf_nParm];
    m_sigScaleErr = result.errors[pf_nParm];
    pf_sbd->scaleSignal(m_sigScale);
  }

  pf_sbd->fluctuate(pf_params);

  //Record final conditions...
  m_condF.sig = pf_sbd->totSignal();
  m_condF.bkgd = pf_sbd->totBkgd();
  m_condF.data = pf_sbd->totData();
  m_condF.chi2 = m_chi2min;
  m_condF.chi2_bins = pf_sbd->calculateChi2LLR(pf_fitSig, pf_sigLLR);
  m_condF.chi2_syst = pf_sbd->calculateSystDiff(pf_params);
  if (m_sigFloat) {
    pf_sbd->scaleSignal(1.0 / m_sigScale);
    pf_sbd->fluctuate(pf_params);
  }

  /// Put Histo stats back (M Owen 16-3-08)
  pf_sbd->useHistoStats(histostat);
  //  if(doPrint) print();
  return;
}

// USED TO BE const double*
std::vector<double> const& ProfileLH::getFitParams() {
  return pf_params;
}

void ProfileLH::print() {

  printf("\n****************************************:\n");
  printf("\nProfileLH Fit Results:\n");
  printf("Fit signal: %d, Signal LLR cut: %f, Fit Status: %d\n", m_fitSig, m_sigLLR, m_status);
  printf("MINUIT chi2 min: %f, edm: %f\n", m_chi2min, m_edm);
  printf("Initial Conditions:\n");
  printf("==>Signal: %.3f\n", m_condI.sig);
  printf("==>Bkgd  : %.3f\n", m_condI.bkgd);
  printf("==>Data  : %.0f\n", m_condI.data);
  printf("==>Chi2  : %.3f (%.3f/%.3f)\n", m_condI.chi2, m_condI.chi2_bins, m_condI.chi2_syst);
  printf("Fitted Conditions:\n");
  printf("==>Signal: %.3f\n", m_condF.sig);
  printf("==>Bkgd  : %.3f\n", m_condF.bkgd);
  printf("==>Data  : %.0f\n", m_condF.data);
  printf("==>Chi2  : %.3f (%.3f/%.3f)\n", m_condF.chi2, m_condF.chi2_bins, m_condF.chi2_syst);
  for (uint i = 0; i < pf_nParm; ++i) {
    printf("==>Fit Param: %14s, Value: %.4f,\t Error: %.4f\n", m_systNames[pf_index[i]].c_str(),
           m_systValues[i], m_systErrors[i]);
  }
  if (m_sigFloat) {
    printf("==>Sig Scale: %f, Error: %f\n", m_sigScale, m_sigScaleErr);
  }
  printf("****************************************:\n");
  return;
}

double ProfileLH::evaluate(const double* parameters) {
  ++pf_iter;

  pf_scl = 1.0;
  if (pf_sigFloat) {
    pf_scl = parameters[pf_nParm];
    if (pf_scl < 1e-5 || !std::isfinite(pf_scl))
      pf_scl = 1e-5;
  }

  for (uint index = 0; index < pf_nParm; ++index)
    pf_params[pf_index[index]] = parameters[index];

  pf_sbd->fluctuate(pf_params, pf_lastMinuitTrialPoint);
  if (pf_sigFloat)
    pf_sbd->scaleSignal(pf_scl);
  const double value =
      pf_sbd->calculateChi2LLR(pf_fitSig, pf_sigLLR) + pf_sbd->calculateSystDiff(pf_params);
  if (pf_sigFloat)
    pf_sbd->scaleSignal(1.0 / pf_scl);
  return value;
}
