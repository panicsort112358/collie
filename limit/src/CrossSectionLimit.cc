#include <stdio.h>
#include "CLpoint.hh"
#include "CrossSectionLimit.hh"
#include <math.h>
#include <TH1D.h>
#include <TF1.h>

Double_t erfFunc(Double_t* x, Double_t* par) {
  Double_t fitval = TMath::Erf((x[0] - par[0]) / (par[1] * sqrt(2) + 1e-9));
  return fitval;
}

class CrossSectionLimitmedFunc {
public:
  CLcompute* p_cl;
  double func(double s, SigBkgdDist& dist, int level, double limit, int nSigma, double& err) {
    CLpoint test(100);
    dist.scaleSignal(s);
    p_cl->calculateCLs(dist, test, level);
    dist.scaleSignal(1.0 / s);

    double retval = test.cls_med - limit;
    double clbval = 0.50;
    if (nSigma == -2) {
      retval = test.cls_med_m2s - limit;
      clbval = (1.0 - 0.0455 / 2.0);
    }
    if (nSigma == -1) {
      retval = test.cls_med_m1s - limit;
      clbval = (1.0 - 0.3173 / 2.0);
    }
    if (nSigma == 1) {
      retval = test.cls_med_p1s - limit;
      clbval = (0.3173 / 2.0);
    }
    if (nSigma == 2) {
      retval = test.cls_med_p2s - limit;
      clbval = (0.0455 / 2.0);
    }

    double esb = test.clsb_med / (1.0 * test.nTrials_exp);
    double eb = clbval / (1.0 * test.nTrials_exp);
    err = sqrt(test.cls_med * test.cls_med *
               (esb / (test.clsb_med * test.clsb_med) + eb / (clbval * clbval)));
    if (err < 0.01)
      err = 0.01;
    if (std::isnan(err))
      err = 0.1;
    if (std::isinf(err))
      err = 0.1;

    return retval;
  }
};

class CrossSectionLimitobsFunc {
public:
  CLcompute* p_cl;
  double func(double s, SigBkgdDist& dist, int level, double limit, int nSigma, double& err) {
    CLpoint test(100);
    dist.scaleSignal(s);
    p_cl->calculateCLs(dist, test, level);
    dist.scaleSignal(1.0 / s);

    double esb = test.clsb_obs / (1.0 * test.nTrials_obs);
    double eb = test.clb_obs / (1.0 * test.nTrials_obs);
    err = sqrt(test.cls_obs * test.cls_obs *
               (esb / (test.clsb_obs * test.clsb_obs) + eb / (test.clb_obs * test.clb_obs)));
    if (err < 0.01)
      err = 0.01;
    if (std::isnan(err))
      err = 0.1;
    if (std::isinf(err))
      err = 0.1;

    return test.cls_obs - limit;
  }
};

//#define SIGN(a,b) ((b)>=0.0 ? ((a)<0.0?-(a):(a)):((a)<0.0?(a):-(a)))

template <class T>

class Seeker1d {

public:
  Seeker1d() {
    p_sbd = NULL;
  }

  SigBkgdDist* p_sbd;

  double seek(double seed = 1.0) {

    double divisor = seed;
    double xmax = -1;
    double xmin = 1e6;

    map<double, double> varMap;
    map<double, double> errMap;
    map<double, double>::iterator iter;
    map<double, double>::iterator iterE;

    int precis = CLcompute::LEVEL_VERYFAST;
    if (m_precision == 0)
      precis = CLcompute::LEVEL_VERYFAST;
    if (m_precision == 1)
      precis = CLcompute::LEVEL_FAST;
    if (m_precision == 2)
      precis = CLcompute::LEVEL_STANDARD;
    if (m_precision == 3)
      precis = CLcompute::LEVEL_FINE;
    if (m_precision == 4)
      precis = CLcompute::LEVEL_VERYFINE;
    if (m_verbose)
      printf("   CrossSectionLimit: Calculating at level %d\n", m_precision);

    //Test the first guess...
    double errval = 0;
    double test1 = m_functor.func(seed, *p_sbd, precis, m_cllevel, m_sigma, errval);
    while ((test1 + m_cllevel) < 0.50) {
      //If true, somehow we are way off!
      //Scan until we get above CLs=0.50
      seed *= 2.0;
      test1 = m_functor.func(seed, *p_sbd, precis, m_cllevel, m_sigma, errval);
    }

    if (m_verbose)
      printf("   CrossSectionLimit: Found delta CL=%.5f%% at X=%f\n", test1 * 100.0, seed);

    //Test the first point to see if we're lucky
    if (fabs(test1) <= m_accuracy) {
      return seed;
    } else {
      varMap[seed] = test1;
      errMap[seed] = errval;
    }
    //Build a range of values...
    //If we're low, check a little higher.  Or vice-versa.
    //If we're far away, boost the jump to converge faster
    double nval = seed * getBoost(test1, m_cllevel);
    test1 = m_functor.func(nval, *p_sbd, precis, m_cllevel, m_sigma, errval);
    varMap[nval] = test1;
    errMap[nval] = errval;
    if (m_verbose)
      printf("   CrossSectionLimit: Found delta CL=%.5f%% at X=%f\n", test1 * 100.0, nval);

    //If we were low, we should be high now.  Are we?
    //Recalculate our boost factor
    nval *= getBoost(test1, m_cllevel);
    test1 = m_functor.func(nval, *p_sbd, precis, m_cllevel, m_sigma, errval);
    varMap[nval] = test1;
    errMap[nval] = errval;
    if (m_verbose)
      printf("   CrossSectionLimit: Found delta CL=%.5f%% at X=%f\n", test1 * 100.0, nval);

    double minVal = 100;
    for (iter = varMap.begin(); iter != varMap.end(); iter++) {
      if ((iter->second + m_cllevel) < 0.55)
        continue;
      if (fabs(iter->second) < minVal) {
        divisor = iter->first;
        minVal = fabs(iter->second);
      }
      if (iter->first < xmin)
        xmin = iter->first;
      if (iter->first > xmax)
        xmax = iter->first;
    }
    xmin *= 0.7 / divisor;
    xmax *= 1.3 / divisor;

    TH1D* hf1 = new TH1D("hf1", "hf1", 10000, xmin, xmax);
    iterE = errMap.begin();
    for (iter = varMap.begin(); iter != varMap.end(); iter++, iterE++) {
      if ((iter->second + m_cllevel) < 0.55)
        continue;

      int ibin = hf1->FindBin(iter->first / divisor);

      while (hf1->GetBinContent(ibin) > 0) {
        ibin += 1;
      }

      hf1->SetBinContent(ibin, iter->second + m_cllevel);
      hf1->SetBinError(ibin, iterE->second);

      if (m_verbose)
        printf("       CrossSectionLimit: Adding point %f, %f\n", iter->first,
               iter->second + m_cllevel);
    }

    //Fit this function...
    if (m_verbose)
      printf("       CrossSectionLimit: Fit Range %f - %f\n", xmin * divisor, xmax * divisor);

    TF1* tf1 = new TF1("myFitErf1", erfFunc, xmin, xmax, 2);
    tf1->SetParameters(-0.1, 0.55);

    hf1->Fit(tf1, "QRN", "", xmin, xmax);
    nval = tf1->GetX(m_cllevel) * divisor;
    if (m_verbose)
      printf("     =>CrossSectionLimit: Fit = %f (%f/%f)\n", nval, tf1->GetParameter(0),
             tf1->GetParameter(1));

    test1 = m_functor.func(nval, *p_sbd, precis, m_cllevel, m_sigma, errval);
    if (m_verbose)
      printf("   CrossSectionLimit: Found delta CL=%.5f%% at X=%f\n\n", test1 * 100.0, nval);
    delete hf1;
    hf1 = NULL;
    delete tf1;
    tf1 = NULL;

    //Test the second guess...
    if (fabs(test1) <= m_accuracy)
      return nval;
    else {
      varMap[nval] = test1;
      errMap[nval] = errval;
    }
    //Didn't work, so try a fifth point!
    //Fit it again...

    minVal = 100;
    xmin = 1e6;
    xmax = -1;
    for (iter = varMap.begin(); iter != varMap.end(); iter++) {
      if ((iter->second + m_cllevel) < 0.55)
        continue;
      if (fabs(iter->second) < minVal) {
        divisor = iter->first;
        minVal = fabs(iter->second);
      }
      if (iter->first < xmin)
        xmin = iter->first;
      if (iter->first > xmax)
        xmax = iter->first;
    }
    xmin *= 0.7 / divisor;
    xmax *= 1.3 / divisor;

    TH1D* hf2 = new TH1D("hf2", "hf2", 10000, xmin, xmax);
    iterE = errMap.begin();
    for (iter = varMap.begin(); iter != varMap.end(); iter++, iterE++) {
      if ((iter->second + m_cllevel) < 0.55)
        continue;

      int ibin = hf2->FindBin(iter->first / divisor);

      while (hf2->GetBinContent(ibin) > 0) {
        ibin += 1;
      }

      hf2->SetBinContent(ibin, iter->second + m_cllevel);
      hf2->SetBinError(ibin, iterE->second);

      if (m_verbose)
        printf("       CrossSectionLimit: Adding point %f, %f\n", iter->first,
               iter->second + m_cllevel);
    }

    //Fit this function...
    if (m_verbose)
      printf("       CrossSectionLimit: Fit Range %f - %f\n", xmin * divisor, xmax * divisor);

    TF1* tf2 = new TF1("myFitErf2", erfFunc, xmin, xmax, 2);
    tf2->SetParameters(-0.1, 0.55);

    hf2->Fit(tf2, "QRN", "", xmin, xmax);

    nval = tf2->GetX(m_cllevel) * divisor;
    if (m_verbose)
      printf("     =>CrossSectionLimit: Fit = %f (%f/%f)\n", nval, tf2->GetParameter(0),
             tf2->GetParameter(1));

    delete tf2;
    tf2 = NULL;
    delete hf2;
    hf2 = NULL;

    test1 = m_functor.func(nval, *p_sbd, precis, m_cllevel, m_sigma, errval);
    if (m_verbose)
      printf("   CrossSectionLimit: Found delta CL=%.5f%% at X=%f\n\n", test1 * 100.0, nval);

    //Test the third guess...
    if (fabs(test1) <= m_accuracy)
      return nval;
    else {
      varMap[nval] = test1;
      errMap[nval] = errval;
    }

    //Didn't find it, try one more time!
    minVal = 100;
    xmin = 1e6;
    xmax = -1;
    for (iter = varMap.begin(); iter != varMap.end(); iter++) {
      if ((iter->second + m_cllevel) < 0.55)
        continue;
      if (fabs(iter->second) < minVal) {
        divisor = iter->first;
        minVal = fabs(iter->second);
      }
      if (iter->first < xmin)
        xmin = iter->first;
      if (iter->first > xmax)
        xmax = iter->first;
    }
    xmin *= 0.7 / divisor;
    xmax *= 1.3 / divisor;

    TH1D* hf3 = new TH1D("hf3", "hf3", 10000, xmin, xmax);
    iterE = errMap.begin();
    for (iter = varMap.begin(); iter != varMap.end(); iter++, iterE++) {
      if ((iter->second + m_cllevel) < 0.55)
        continue;

      int ibin = hf3->FindBin(iter->first / divisor);

      while (hf3->GetBinContent(ibin) > 0) {
        ibin += 1;
      }

      hf3->SetBinContent(ibin, iter->second + m_cllevel);
      hf3->SetBinError(ibin, iterE->second);

      if (m_verbose)
        printf("       CrossSectionLimit: Adding point %f, %f\n", iter->first,
               iter->second + m_cllevel);
    }

    //Fit this function...
    if (m_verbose)
      printf("       CrossSectionLimit: Fit Range %f - %f\n", xmin * divisor, xmax * divisor);

    TF1* tf3 = new TF1("myFitErf3", erfFunc, xmin, xmax, 2);
    tf3->SetParameters(-0.1, 0.55);

    hf3->Fit(tf3, "QRN", "", xmin, xmax);

    nval = tf3->GetX(m_cllevel) * divisor;
    if (m_verbose)
      printf("     =>CrossSectionLimit: Fit = %f (%f/%f)\n", nval, tf3->GetParameter(0),
             tf3->GetParameter(1));
    delete tf3;
    tf3 = NULL;
    delete hf3;
    hf3 = NULL;
    test1 = m_functor.func(nval, *p_sbd, precis, m_cllevel, m_sigma, errval);
    if (m_verbose)
      printf("   CrossSectionLimit: Found delta CL=%.5f%% at X=%f\n", test1 * 100.0, nval);

    return nval;
  }

  bool m_verbose;

  T m_functor;
  double m_cllevel;
  double m_accuracy;
  int m_precision;
  int m_sigma;

protected:
  double isign(double val) {
    if (val <= 0)
      return -1.0;
    else
      return 1.0;
  }

  double getBoost(double test, double goal) {

    double cls = test + goal;

    double nomVal = TMath::ErfInverse(goal) * sqrt(2) * 0.5 + 0.5;
    double iVal = TMath::ErfInverse(0.9999 * cls) * sqrt(2) * 0.5 + 0.5;

    double boost = nomVal / iVal;

    //There is some noise in the system, so don't over do it!
    if (fabs(boost - 1) > 0.50)
      boost = 1.0 + isign(boost) * 0.50;
    if (fabs(boost - 1) < 0.10)
      boost = 1.0 + isign(boost) * 0.10;

    return boost;
  }
};

/* Calculate Cross-Section Branching Ratio Limits */
bool CrossSectionLimit::calculate(const SigBkgdDist& dist, CLpoint& CLs) {
  if (dist.totSignal() < 1e-6) {
    printf("CrossSectionLimit::calculate, no signal present! %.3fx1e-6\n", dist.totSignal() * 1e6);
    return false;
  }

  if (!p_cl) {
    printf("CrossSectionLimit::calculate, Error: CLcompute class not specified!\n");
    return false;
  }

  SigBkgdDist scaled(dist);

  if (p_cl->getNoviceFlag())
    p_cl->noviceWarning();

  // signal may need a very-very small flat background to converge
  double tot = scaled.totSignal();
  tot = tot / scaled.nbins() * (1e-6);
  for (int i = 0; i < scaled.nbins(); i++)
    scaled.signal(i) += tot;

  // we have a distribution...
  Seeker1d<CrossSectionLimitmedFunc> medSeeker;
  Seeker1d<CrossSectionLimitobsFunc> obsSeeker;

  medSeeker.p_sbd = &scaled;
  medSeeker.m_functor.p_cl = p_cl;
  medSeeker.m_accuracy = m_accuracy;
  medSeeker.m_precision = m_precision;
  medSeeker.m_verbose = m_verbose;
  medSeeker.m_cllevel = m_cllevel;
  medSeeker.m_sigma = m_sigma;

  obsSeeker.p_sbd = &scaled;
  obsSeeker.m_functor.p_cl = p_cl;
  obsSeeker.m_accuracy = m_accuracy;
  obsSeeker.m_precision = m_precision;
  obsSeeker.m_verbose = m_verbose;
  obsSeeker.m_cllevel = m_cllevel;
  obsSeeker.m_sigma = m_sigma;

  if (m_verbose) {
    printf("==>CL target: %f \n", m_cllevel);
    printf("==>CL accuracy: %.3f%% \n", m_accuracy * 100);
  }

  double saveSeed = m_seed;

  if (saveSeed == 1.00) {
    m_seed = medSeeker.p_sbd->calculateDeltaLLR();
    if (m_seed > 0)
      m_seed = sqrt(7.82 / m_seed);
    else
      m_seed = 1;
    if (m_verbose)
      printf("Exp search seed set to %f\n", m_seed);
  }

  double factor_med = -1;
  if (m_calcExp)
    factor_med = medSeeker.seek(m_seed);

  if (saveSeed == 1.00) {
    m_seed = medSeeker.p_sbd->calculateDeltaLLRobs();
    if (m_seed > 0)
      m_seed = sqrt(7.82 / m_seed);
    else
      m_seed = 1;
    if (m_verbose)
      printf("Obs search seed set to %f\n", m_seed);
  }

  double factor_obs = -1;
  if (m_calcObs)
    factor_obs = obsSeeker.seek(m_seed);

  if (m_verbose)
    printf("Found exp (obs) %f (%f)\n", factor_med, factor_obs);

  CLs.xsec_cl = m_cllevel;
  CLs.xsec_obsfactor = factor_obs;
  CLs.xsec_medfactor = factor_med;

  m_expFact = factor_med;
  m_obsFact = factor_obs;

  return true;
}

void CrossSectionLimit::print() {
  printf("\n****************************************:\n");

  printf("CrossSectionLimit calculator results::\n");

  printf("==>CL level: %.4f, N sigma: %d\n", m_cllevel, m_sigma);
  printf("==>Accuracy: %.4f, Precision: %d\n", m_accuracy, m_precision);
  printf("==>Scaling factor Exp: %.6f, Obs: %.3f\n", m_expFact, m_obsFact);
  if (p_cl->getMedianExpected())
    printf("\n==>Expected limit calculated using median expected outcome\n");
  else
    printf("\n==>Expected limit calculated using Data==Bkgd condition\n");
  printf("****************************************:\n");
  return;
}
