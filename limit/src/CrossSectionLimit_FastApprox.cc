#include "CrossSectionLimit_FastApprox.hh"
#include "CrossSectionLimit.hh"
#include "CLfit2.hh"
#include "CLfast.hh"

/* Calculate Cross-Section Branching Ratio Limits */
bool CrossSectionLimit_FastApprox::calculate(const SigBkgdDist& dist, CLpoint& CLs) {

  //Equivalant 2-sigma values for CLs overcoverage
  double nomSigmas[5];
  nomSigmas[4] = 1.671; //p2s
  nomSigmas[3] = 1.743;
  nomSigmas[2] = 1.978; //med
  nomSigmas[1] = 2.436;
  nomSigmas[0] = 3.090; //m2s

  //Gather some information on the channel
  calculateFast(dist, CLs);
  calculateFit(dist, CLs);

  //Calibrate our expectations and calculate limits
  double expGuess = sqrt(7.82 / m_dLLRnomExp);
  double glCf = m_expFact / expGuess;
  for (int i = 0; i < 5; i++)
    nomSigmas[i] *= glCf;
  printf("Fast: %f, Fit: %f\n", m_dLLRnomExp, m_dLLRfitExp);

  //Poisson-adjusted expected limit, postfit
  CLs.xsec_medfactor_m2s =
      sqrt(2.0 * (nomSigmas[0] - 2 * glCf) * (nomSigmas[0] - 2 * glCf) / m_dLLRfitExp);
  CLs.xsec_medfactor_m1s =
      sqrt(2.0 * (nomSigmas[1] - 1 * glCf) * (nomSigmas[1] - 1 * glCf) / m_dLLRfitExp);
  CLs.xsec_medfactor =
      sqrt(2.0 * (nomSigmas[2] + 0 * glCf) * (nomSigmas[2] + 0 * glCf) / m_dLLRfitExp);
  CLs.xsec_medfactor_p1s =
      sqrt(2.0 * (nomSigmas[3] + 1 * glCf) * (nomSigmas[3] + 1 * glCf) / m_dLLRfitExp);
  CLs.xsec_medfactor_p2s =
      sqrt(2.0 * (nomSigmas[4] + 2 * glCf) * (nomSigmas[4] + 2 * glCf) / m_dLLRfitExp);

  //Locate observed delta-llr value in units of chi2
  double obsChi2 = m_fitBOnlyDLLR[0];
  for (int i = 1; i < 5; i++) {
    if (m_fitLLRobs > m_fitBOnlyLLR[i]) {
      obsChi2 = m_fitBOnlyDLLR[i - 1];
      break;
    }
  }

  //Translate to N-sigmas and determine the reference value
  double obsSigma = m_dLLRfitObs / obsChi2;
  double obsSigmaRef = 0;
  if (obsSigma < -2)
    obsSigmaRef = nomSigmas[0] + (nomSigmas[0] - nomSigmas[1]) * (-2 - obsSigma);
  else if (obsSigma < -1)
    obsSigmaRef = nomSigmas[1] + (nomSigmas[0] - nomSigmas[1]) * (-1 - obsSigma);
  else if (obsSigma < 0)
    obsSigmaRef = nomSigmas[2] + (nomSigmas[1] - nomSigmas[2]) * (0 - obsSigma);
  else if (obsSigma < 1)
    obsSigmaRef = nomSigmas[3] + (nomSigmas[2] - nomSigmas[3]) * (1 - obsSigma);
  else if (obsSigma < 2)
    obsSigmaRef = nomSigmas[4] + (nomSigmas[3] - nomSigmas[4]) * (2 - obsSigma);
  else
    obsSigmaRef = nomSigmas[4] + (nomSigmas[3] - nomSigmas[4]) * (2 - obsSigma);

  obsSigma *= glCf;
  obsSigmaRef *= glCf;
  CLs.xsec_obsfactor =
      sqrt(2.0 * (obsSigmaRef + obsSigma) * (obsSigmaRef + obsSigma) / m_dLLRfitExp);

  myCLpoint.xsec_medfactor_m2s = CLs.xsec_medfactor_m2s;
  myCLpoint.xsec_medfactor_m1s = CLs.xsec_medfactor_m1s;
  myCLpoint.xsec_medfactor_p2s = CLs.xsec_medfactor_p2s;
  myCLpoint.xsec_medfactor_p1s = CLs.xsec_medfactor_p1s;
  myCLpoint.xsec_medfactor = CLs.xsec_medfactor;
  myCLpoint.xsec_obsfactor = CLs.xsec_obsfactor;

  return true;
}

void CrossSectionLimit_FastApprox::calculateFast(const SigBkgdDist& dist, CLpoint& CLs) {

  CLfast clcomputeFast;

  CrossSectionLimit csLim;
  csLim.setup(&clcomputeFast);
  csLim.setVerbose(m_verbose);
  csLim.setCLlevel(m_cllevel);
  csLim.setAccuracy(m_accuracy);
  csLim.setPrecision(m_precision);

  m_dLLRnomExp = dist.calculateDeltaLLR();

  m_seed = dist.calculateDeltaLLR();
  if (m_seed > 0)
    m_seed = sqrt(7.82 / m_seed);
  else
    m_seed = 1;
  csLim.setSearchSeed(m_seed);
  csLim.calculateExpected(1);
  csLim.calculateObserved(0);

  myCLpoint.reset(CLs.getVar1(), CLs.getVar2(), CLs.getVar2());
  csLim.calculate(dist, myCLpoint);
  m_expFact = myCLpoint.xsec_medfactor;

  return;
}

void CrossSectionLimit_FastApprox::calculateFit(const SigBkgdDist& dist, CLpoint& CLs) {

  CLfit2 clcomputeFit;

  clcomputeFit.calculateCLs(dist, myCLpoint, CLcompute::LEVEL_TEST);

  m_fitBOnlyLLR.clear();
  m_fitBOnlyDLLR.clear();
  m_fitBOnlyLLR.push_back(myCLpoint.llrb_m2s);
  m_fitBOnlyLLR.push_back(myCLpoint.llrb_m1s);
  m_fitBOnlyLLR.push_back(myCLpoint.llrb);
  m_fitBOnlyLLR.push_back(myCLpoint.llrb_p1s);
  m_fitBOnlyLLR.push_back(myCLpoint.llrb_p2s);
  m_fitLLRobs = myCLpoint.llrobs;

  for (int i = 0; i < 4; i++)
    m_fitBOnlyDLLR.push_back(m_fitBOnlyLLR[i] - m_fitBOnlyLLR[i + 1]);

  CLs.llrb_m2s = myCLpoint.llrb_m2s;
  CLs.llrb_m1s = myCLpoint.llrb_m1s;
  CLs.llrb = myCLpoint.llrb;
  CLs.llrb_p1s = myCLpoint.llrb_p1s;
  CLs.llrb_p2s = myCLpoint.llrb_p2s;
  CLs.llrsb = myCLpoint.llrsb;
  CLs.llrobs = myCLpoint.llrobs;

  m_dLLRfitExp = myCLpoint.llrb - myCLpoint.llrsb;
  m_dLLRfitObs = myCLpoint.llrb - myCLpoint.llrobs;

  return;
}

void CrossSectionLimit_FastApprox::print() {
  printf("\n****************************************:\n");
  printf("CrossSectionLimit Fast CLfit2 Approximation results::\n");
  printf("\nLLR Results:");
  myCLpoint.print();
  printf("\nLimit Calculation Results:\n");
  printf("==>CL level: %.4f\n", m_cllevel);
  printf("==>Accuracy: %.4f, Precision: %d\n", m_accuracy, m_precision);
  printf("==>Scaling factor Exp: %.4f, Obs: %.3f\n", myCLpoint.xsec_medfactor,
         myCLpoint.xsec_obsfactor);
  printf("==>N-sigma factors -2,-1,0,+1,+2: %.3f / %.3f / %.3f / %.3f / %.3f\n",
         myCLpoint.xsec_medfactor_m2s, myCLpoint.xsec_medfactor_m1s, myCLpoint.xsec_medfactor,
         myCLpoint.xsec_medfactor_p1s, myCLpoint.xsec_medfactor_p2s);
  printf("\n****************************************:\n");
  return;
}
