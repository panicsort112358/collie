#include "ProfileLH_2D.hh"
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
ProfileLH_2D::ProfileLH_2D() {
  setup();
}

///Constructor with SBD initializer
ProfileLH_2D::ProfileLH_2D(SigBkgdDist* asbd) {
  setup();
  if (asbd == NULL) {
    printf("ProfileLH_2D::ProfileLH_2D, the input SBD is NULL!\n");
    return;
  }
  setModel(asbd);
}

ProfileLH_2D::~ProfileLH_2D() = default;

void ProfileLH_2D::setup() {
  pf2d_sbd = NULL;
  pf2d_params.clear(); // USED TO BE:  pf2d_params = NULL;
  pf2d_index.clear();
  pf2d_paramErrs.clear();
  pf2d_nParm = 0;
  pf2d_iter = 0;
  pf2d_fitSig = true;
  pf2d_sigFloat = false;
  pf2d_sigLLR = 1e6;
  pf2d_scl1 = 1.0;
  pf2d_scl2 = 1.0;

  m_init = false;
  m_fitTest = false;
  m_fitSig = true;
  m_verbose = false;
  m_sigFloat = false;

  m_sigScale1 = 1.0;
  m_sigScale2 = 1.0;

  m_sigScaleNom1 = 1.0;
  m_sigScaleNom2 = 1.0;

  m_sigScaleErr2 = 0.0;
  m_sigScaleErr1 = 0.0;

  m_sigLLR = 1e6;
  m_chi2min = -1.0;
  m_edm = -1.0;
  m_lun = 99;
  m_status = -1;
  m_maxIter = 55000;

  m_systNamesFit.clear();
  m_systNamesFitBkgd.clear();
  m_systNamesFitMap.clear();
  m_systNamesFitBkgdMap.clear();
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

void ProfileLH_2D::setModel(SigBkgdDist* asbd) {
  pf2d_sbd = asbd;
  pf2d_fitSig = m_fitSig;
  pf2d_sigLLR = m_sigLLR;
  fillSyst();
}

int ProfileLH_2D::getNiterations() {
  return pf2d_iter;
}

void ProfileLH_2D::getMNerrs() {
  if (!m_init || !m_engine) {
    printf("ProfileLH_2D::getMNerrs, not initialized\n");
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
  for (uint num = 0; num < pf2d_nParm; ++num) {
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
    if (m_engine->minosError(pf2d_nParm, errP, errM, eparab, globcc))
      m_sigScaleErr1 = eparab;
    if (m_engine->minosError(pf2d_nParm + 1, errP, errM, eparab, globcc))
      m_sigScaleErr2 = eparab;
  }
}

bool ProfileLH_2D::getEMat(double* emat) {
  if (!m_init || !m_engine) {
    printf("ProfileLH_2D::getEMat, not initialized\n");
    return false;
  }
  getMNerrs();
  return m_engine->covariance(emat, pf2d_nParm);
}

void ProfileLH_2D::fillSyst() {
  m_systNames.clear();
  m_systNamesFit.clear();
  m_systNamesFitBkgd.clear();
  m_systNamesFitMap.clear();
  m_systNamesFitBkgdMap.clear();

  if (pf2d_sbd == NULL) {
    printf(
        "ProfileLH_2D::fillSyst, You're trying to fit a profile with a NULL signal/background distribution!!\n");
    return;
  }

  //Make sure all the systematics are properly set
  if (!pf2d_sbd->isInit()) {
    pf2d_sbd->varySystematics();
    pf2d_sbd->setBaselineModel();
  }

  vector<string> chanNames = pf2d_sbd->getChannelNames();
  const CollieDistribution* dist;
  string systName = "";
  map<string, int> m_systMap;

  for (vector<string>::iterator its = chanNames.begin(); its != chanNames.end(); ++its) {
    for (uint b = 0; b < pf2d_sbd->getNbkgdDist(*its); ++b) {

      dist = pf2d_sbd->getBkgdDist(*its, b);
      if (dist != NULL) {

        for (int s = 0; s < dist->getNsystematics(); ++s) {
          systName = dist->getSystName(s);
          if (m_systMap.find(systName) != m_systMap.end())
            continue;
          else {
            m_systNames.push_back(systName);
            if (pf2d_sbd->getSystFitFlag(systName)) {
              m_systNamesFitMap[systName] = pf2d_sbd->getSystIndex(systName);
              m_systNamesFitBkgdMap[systName] = pf2d_sbd->getSystIndex(systName);
              m_systNamesFit.push_back(systName);
              m_systNamesFitBkgd.push_back(systName);
            } else {
              if (m_verbose)
                printf("ProfileLH_2D:: Not Fitting %s\n", dist->getSystName(s).c_str());
            }
            m_systMap[systName] = pf2d_sbd->getSystIndex(systName);
          }
        }
      }
    }
  }

  for (vector<string>::iterator its = chanNames.begin(); its != chanNames.end(); ++its) {
    for (uint b = 0; b < pf2d_sbd->getNsigDist(*its); ++b) {

      dist = pf2d_sbd->getSigDist(*its, b);
      if (dist != NULL) {

        for (int s = 0; s < dist->getNsystematics(); ++s) {
          systName = dist->getSystName(s);

          if (m_systMap.find(systName) != m_systMap.end())
            continue;
          else {
            m_systNames.push_back(systName);
            if (pf2d_sbd->getSystFitFlag(systName) && pf2d_fitSig) {
              m_systNamesFitMap[systName] = pf2d_sbd->getSystIndex(systName);
              m_systNamesFit.push_back(systName);
            } else {
              if (m_verbose)
                printf("ProfileLH_2D:: Not Fitting %s\n", systName.c_str());
            }
            m_systMap[systName] = pf2d_sbd->getSystIndex(systName);
          }
        }
      }
    }
  }

  //our adjustable copy of the syst params
  pf2d_params = std::vector<double>(m_systNames.size(), 0.0);
  pf2d_paramErrs = std::vector<double>(m_systNames.size(), 0.0);
  //  m_systErrors.clear();
  //  m_systValues.clear();

  //collect existing parameters
  for (unsigned int i = 0; i < m_systNames.size(); ++i)
    pf2d_params[i] = pf2d_sbd->getSystFluctValue(i);

  initFitParams();

  return;
}

void ProfileLH_2D::initFitParams() {
  //Make a storage container with only the values we're fitting...
  //   either background-only or sig+bkgd
  //   for bkgd-only, make index so we know where it fits in master copy.
  //   Allow individual fit parameters to be excluded by the user.

  pf2d_index.clear();
  map<string, int>::iterator im;

  if (m_fitSig) {
    pf2d_nParm = m_systNamesFit.size();
    pf2d_index = std::vector<uint>(pf2d_nParm, 0);

    //now fill in the syst indices
    for (uint s = 0; s < m_systNamesFit.size(); ++s) {
      im = m_systNamesFitMap.find(m_systNamesFit[s]);
      pf2d_index[s] = im->second;
    }
  } else {
    pf2d_nParm = m_systNamesFitBkgd.size();
    pf2d_index = std::vector<uint>(pf2d_nParm, 0);

    //now fill in the syst indices
    for (uint s = 0; s < m_systNamesFitBkgd.size(); ++s) {
      im = m_systNamesFitBkgdMap.find(m_systNamesFitBkgd[s]);
      pf2d_index[s] = im->second;
    }
  }

  return;
}

bool ProfileLH_2D::getMinosError(uint idx, double& errP, double& errM, double& parab,
                                 double& globCC) {
  if (m_systErrorsUp.size() == 0)
    getMNerrs();
  if (m_systErrorsUp.size() == 0) {
    printf("ProfileLH_2D::getMinosError, getMinosError failed...\n");
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

double ProfileLH_2D::getParamVal(uint p) {
  if (p >= pf2d_params.size())
    return 0;
  return pf2d_params[p];
}

double ProfileLH_2D::getParamError(uint p) {
  if (p >= pf2d_paramErrs.size())
    return 0;
  return pf2d_paramErrs[p];
}

int ProfileLH_2D::getFitSystIndex(string name) {
  map<string, int>::iterator iter;

  if (name == "Sig1")
    return pf2d_nParm + 1;
  if (name == "Sig2")
    return pf2d_nParm + 2;

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

double ProfileLH_2D::getFitParamErr(string name) const {

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
  //  printf("ProfileLH2D::getFitParamErr, No such systematic: %s\n",name.c_str());
  return 0;
}

double ProfileLH_2D::getFitParamVal(string name) const {
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
  //  printf("ProfileLH2D::getFitParamValue, No such systematic: %s\n",name.c_str());
  return 0;
}

string ProfileLH_2D::getFitSystName(uint idx) {
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

double ProfileLH_2D::getSystSum2() {
  double tot = 0;
  for (uint i = 0; i < m_systNames.size(); ++i) {
    if (pf2d_sbd->getFloatFlag(i))
      continue;
    tot += pf2d_sbd->getSystFluctValue(i) * pf2d_sbd->getSystFluctValue(i);
  }
  return tot;
}

void ProfileLH_2D::fitProfile() {

  if (pf2d_sbd == NULL) {
    printf(
        "ProfileLH_2D::fitProfile, You're trying to fit a profile with a NULL signal/background distribution!!\n");
    return;
  }

  //Initialize our parameters
  pf2d_iter = 0;
  m_status = -1;
  m_chi2min = -1.0;
  m_edm = -1.0;
  m_sigScale1 = m_sigScaleNom1;
  m_sigScaleErr1 = 0.0;

  m_sigScale2 = m_sigScaleNom2;
  m_sigScaleErr2 = 0.0;

  //check to see if the fitting model has changed...
  if ((pf2d_fitSig != m_fitSig) || (m_sigFloat != pf2d_sigFloat)) {
    pf2d_fitSig = m_fitSig;
    pf2d_sigFloat = m_sigFloat;
    fillSyst();
    m_init = false;
  }

  //M Owen 16-3-08
  //Only want to vary the systematics in each fit iterations
  //do not want to vary with the stat uncertainties.
  bool histostat = pf2d_sbd->usingHistoStats();
  pf2d_sbd->useHistoStats(false);

  if (m_verbose)
    printf("ProfileLH_2D: %zu systematics, %zu fit params\n", m_systNames.size(),
           m_systNamesFit.size());

  //collect existing parameters
  for (uint i = 0; i < m_systNames.size(); ++i) {
    pf2d_params[i] = pf2d_sbd->getSystFluctValue(i);
    pf2d_paramErrs[i] = 0;
  }

  //Record initial conditions...
  m_condI.sig = pf2d_sbd->totSignal();
  m_condI.bkgd = pf2d_sbd->totBkgd();
  m_condI.data = pf2d_sbd->totData();
  m_condI.chi2 = pf2d_sbd->calculateChi2LLR(pf2d_fitSig, pf2d_sigLLR) +
                 pf2d_sbd->calculateSystDiff(pf2d_params);

  m_condI.chi2_bins = pf2d_sbd->calculateChi2LLR(pf2d_fitSig, pf2d_sigLLR);
  m_condI.chi2_syst = pf2d_sbd->calculateSystDiff(pf2d_params);

  // Establish the memory of the last Minuit trial point, so that
  // sbd->fillArrays() can operate more efficiently
  pf2d_lastMinuitTrialPoint.clear();

  std::vector<collie::ProfileParameter> parameters;
  parameters.reserve(pf2d_nParm + (m_sigFloat ? 2U : 0U));
  for (uint num = 0; num < pf2d_nParm; ++num) {
    collie::ProfileParameter parameter;
    parameter.name = "Par" + std::to_string(num + 1);
    parameter.initial = m_fitTest || m_sigFloat ? 0.0
                        : m_init && m_systValues.size() == pf2d_nParm
                            ? m_systValues[num]
                            : pf2d_params[pf2d_index[num]];
    parameter.step = 0.1;
    parameter.lower = -4.5;
    parameter.upper = 4.5;
    parameter.limited = true;
    parameters.push_back(parameter);
  }
  if (m_sigFloat) {
    collie::ProfileParameter signal1;
    signal1.name = "SigRate1";
    signal1.initial = m_sigScaleNom1;
    signal1.step = 0.1;
    signal1.lower = 1e-5;
    signal1.upper = 1e5;
    signal1.limited = true;
    parameters.push_back(signal1);

    collie::ProfileParameter signal2;
    signal2.name = "SigRate2";
    signal2.initial = m_sigScaleNom2;
    signal2.step = 0.1;
    signal2.lower = 1e-5;
    signal2.upper = 1e5;
    signal2.limited = true;
    parameters.push_back(signal2);
  }

  collie::ProfileFitResult result;
  if (parameters.empty()) {
    m_engine.reset();
    pf2d_iter = 1;
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
      pf2d_iter = static_cast<uint>(result.evaluations);
      m_chi2min = result.minimum;
      m_edm = result.edm;
      m_status = result.covarianceStatus;
      m_init = result.converged;
    } catch (const std::exception& error) {
      fprintf(stderr, "ProfileLH_2D: Minuit2 setup failed: %s\n", error.what());
      m_status = -1;
      m_init = false;
      pf2d_sbd->useHistoStats(histostat);
      return;
    }
  }

  if (m_status != 3 && m_verbose)
    printf("ProfileLH_2D: Warning! Minuit2 convergence failing (%d)!\n", m_status);
  if (m_verbose)
    printf("ProfileLH_2D: %u Minuit2 evaluations\n", pf2d_iter);
  if (pf2d_iter > 0.95 * m_maxIter)
    printf("ProfileLH_2D: Minuit2 evaluation warning: %u evaluations\n", pf2d_iter);

  m_systErrors.clear();
  m_systValues.clear();
  m_systErrorsUp.clear();
  m_systErrorsDn.clear();
  m_systErrorsParab.clear();
  m_systErrorsGlobCC.clear();

  for (uint num = 0; num < pf2d_nParm; ++num) {
    const double val = result.values[num];
    const double error = result.errors[num];
    m_systErrors.push_back(error);
    m_systValues.push_back(val);
    pf2d_params[pf2d_index[num]] = val;
    pf2d_paramErrs[pf2d_index[num]] = error;
  }

  if (m_sigFloat) {
    m_sigScale1 = result.values[pf2d_nParm];
    m_sigScaleErr1 = result.errors[pf2d_nParm];
    pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(0), m_sigScale1, pf2d_params);

    m_sigScale2 = result.values[pf2d_nParm + 1];
    m_sigScaleErr2 = result.errors[pf2d_nParm + 1];
    pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(1), m_sigScale2, pf2d_params);
  }

  pf2d_sbd->fluctuate(pf2d_params);

  //Record final conditions...
  m_condF.sig = pf2d_sbd->totSignal();
  m_condF.bkgd = pf2d_sbd->totBkgd();
  m_condF.data = pf2d_sbd->totData();
  m_condF.chi2 = m_chi2min;
  m_condF.chi2_bins = pf2d_sbd->calculateChi2LLR(pf2d_fitSig, pf2d_sigLLR);
  m_condF.chi2_syst = pf2d_sbd->calculateSystDiff(pf2d_params);
  if (m_sigFloat) {
    pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(0), 1.0 / m_sigScale1, pf2d_params);
    pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(1), 1.0 / m_sigScale2, pf2d_params);
    pf2d_sbd->fluctuate(pf2d_params);
  }

  /// Put Histo stats back (M Owen 16-3-08)
  pf2d_sbd->useHistoStats(histostat);
  //  if(doPrint) print();
  return;
}

// USED TO BE const double*
std::vector<double> const& ProfileLH_2D::getFitParams() {
  return pf2d_params;
}

void ProfileLH_2D::print() {

  printf("\n****************************************:\n");
  printf("\nProfileLH_2D Fit Results:\n");
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
  for (uint i = 0; i < pf2d_nParm; ++i) {
    printf("==>Fit Param: %14s, Value: %.4f,\t Error: %.4f\n", m_systNames[pf2d_index[i]].c_str(),
           m_systValues[i], m_systErrors[i]);
  }
  if (m_sigFloat) {
    printf("==>Sig1 Scale: %f, Error: %f\n", m_sigScale1, m_sigScaleErr1);
    printf("==>Sig2 Scale: %f, Error: %f\n", m_sigScale2, m_sigScaleErr2);
  }
  printf("****************************************:\n");
  return;
}

double ProfileLH_2D::evaluate(const double* parameters) {
  ++pf2d_iter;

  pf2d_scl1 = 1.0;
  pf2d_scl2 = 1.0;
  if (pf2d_sigFloat) {
    pf2d_scl1 = parameters[pf2d_nParm];
    if (pf2d_scl1 < 1e-5 || !std::isfinite(pf2d_scl1))
      pf2d_scl1 = 1e-5;
    pf2d_scl2 = parameters[pf2d_nParm + 1];
    if (pf2d_scl2 < 1e-5 || !std::isfinite(pf2d_scl2))
      pf2d_scl2 = 1e-5;
  }

  for (uint index = 0; index < pf2d_nParm; ++index)
    pf2d_params[pf2d_index[index]] = parameters[index];

  pf2d_sbd->fluctuate(pf2d_params, pf2d_lastMinuitTrialPoint);
  if (pf2d_sigFloat) {
    pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(0), pf2d_scl1, pf2d_params);
    pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(1), pf2d_scl2, pf2d_params);
  }
  const double value = pf2d_sbd->calculateChi2LLR(pf2d_fitSig, pf2d_sigLLR) +
                       pf2d_sbd->calculateSystDiff(pf2d_params);
  if (pf2d_sigFloat) {
    pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(0), 1.0 / pf2d_scl1, pf2d_params);
    pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(1), 1.0 / pf2d_scl2, pf2d_params);
  }
  return value;
}
