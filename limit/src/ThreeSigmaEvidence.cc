#include <stdio.h>
#include "ThreeSigmaEvidence.hh"
#include <math.h>

class ThreeSigmaEvidencemedFunc {
public:
  CLcompute* p_cl;
  double func(double s, SigBkgdDist& dist, int level, double limit, int nSigma = 0) {
    CLpoint test(100);
    dist.scaleSignal(s);
    dist.scaleBackground(s);
    p_cl->calculateCLs(dist, test, level);
    dist.scaleSignal(1.0 / s);
    dist.scaleBackground(1.0 / s);

    float retval = test.clb_sb - limit;
    if (nSigma == -1)
      retval = test.clb_sb_m1s - limit;
    if (nSigma == -2)
      retval = test.clb_sb_m2s - limit;
    if (nSigma == 1)
      retval = test.clb_sb_p1s - limit;
    if (nSigma == 2)
      retval = test.clb_sb_p2s - limit;

    return retval;
  }
};
class ThreeSigmaEvidenceobsFunc {
public:
  CLcompute* p_cl;
  double func(double s, SigBkgdDist& dist, int level, double limit, int nSigma = 0) {
    CLpoint test(100);
    dist.scaleSignal(s);
    dist.scaleBackground(s);
    p_cl->calculateCLs(dist, test, level);
    dist.scaleSignal(1.0 / s);
    dist.scaleBackground(1.0 / s);
    return (1.0 - test.clb_obs) - limit;
  }
};

#define SIGN(a, b) ((b) >= 0.0 ? ((a) < 0.0 ? -(a) : (a)) : ((a) < 0.0 ? (a) : -(a)))

template <class T> class Seeker1d {
public:
  Seeker1d() {
    p_sbd = NULL;
  }
  SigBkgdDist* p_sbd;
  double seek(double seed = 1.0) {
    double x1 = seed;
    double x2 = seed * 4;
    double factor = 0; // final factor
    int its = 10;
    // first we seek a bracket
    do {
      bracket(x1, x2);
      if (m_verbose)
        printf("==>Limit bracketed at Xlow=%f and Xhigh=%f\n", x1, x2);
      // now we can use Ridder
      factor = zriddr(x1, x2);
      if (factor < 0)
        factor = zriddr(x1 / 2, x2 * 2);
      its--;
    } while (factor < 0 && its > 0);
    return factor;
  }
  bool m_verbose;

  T m_functor;
  double m_cllevel;
  double m_accuracy;
  int m_precision;
  int m_sigma;

protected:
  double zriddr(double x1, double x2) {
    int j;
    double ans, fh, fl, fm, fnew, s, xh, xl, xm;

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

    fl = m_functor.func(x1, *p_sbd, precis, m_cllevel, m_sigma);
    fh = m_functor.func(x2, *p_sbd, precis, m_cllevel, m_sigma);
    if ((fl > 0.0f && fh < 0.0f) || (fl < 0.0f && fh > 0.0f)) {
      xl = x1;
      xh = x2;
      ans = (-1.11e30);
      for (j = 0; j < 20; j++) {
        ///try midpoint of two scaling factors...
        xm = 0.5 * (xl + xh);
        fm = m_functor.func(xm, *p_sbd, precis, m_cllevel, m_sigma);
        if (m_verbose)
          printf("   ThreeSigmaEvidence: Found delta CL=%.5f%% at X=%f\n", fm * 100.0, xm);
        s = sqrt(fm * fm - fl * fh);
        if (s == 0.0)
          return xm;
        if (fabs(fm) <= m_accuracy)
          return xm;

        ans = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
        fnew = m_functor.func(ans, *p_sbd, precis, m_cllevel, m_sigma);
        if (m_verbose)
          printf("   ThreeSigmaEvidence: Found 2nd delta CL=%.5f%% at X=%f\n", fnew * 100.0, ans);
        if (fnew == 0.0)
          return ans;
        if (fabs(fnew) <= m_accuracy)
          return ans;

        if (SIGN(fm, fnew) != fm) {
          xl = xm;
          fl = fm;
          xh = ans;
          fh = fnew;
        } else if (SIGN(fl, fnew) != fl) {
          xh = ans;
          fh = fnew;
        } else if (SIGN(fh, fnew) != fh) {
          xl = ans;
          fl = fnew;
        }
      }
      return ans;
    }
    return -1;
  }

  void bracket(double& x1, double& x2) {

    int precis = CLcompute::LEVEL_VERYFAST;
    if (m_precision >= 2)
      precis = CLcompute::LEVEL_FAST;
    if (m_precision == 4)
      precis = CLcompute::LEVEL_STANDARD;

    for (int i = 0; i < 100; i++) {
      double value = m_functor.func(x2, *p_sbd, precis, m_cllevel, m_sigma);
      if (m_verbose)
        printf("ThreeSigmaEvidence, Found %f at %f\n", value, x2);
      if (x2 > 1.0) {
        if (value > 0.0)
          break;
        x1 = x2;
        x2 *= 2.0;
      } else {
        if (value < 0.0)
          break;
        x1 = x2;
        x2 /= 2.0;
      }
    }
    if (x1 > x2) {
      x1 *= 1.5;
      x2 /= 1.5;
    } else {
      x1 /= 1.5;
      x2 *= 1.5;
    }
    return;
  }
};

/* Calculate Luminosity Ratio Limits */
bool ThreeSigmaEvidence::calculate(const SigBkgdDist& dist, CLpoint& CLs) {
  if (dist.totSignal() < 1e-6)
    return false;

  SigBkgdDist scaled(dist);

  // signal may need a very-very small flat background
  double tot = scaled.totSignal();
  tot = tot / scaled.nbins() * (1e-8);
  for (int i = 0; i < scaled.nbins(); i++)
    scaled.signal(i) += tot;

  // we have a distribution...
  Seeker1d<ThreeSigmaEvidencemedFunc> medSeeker;
  Seeker1d<ThreeSigmaEvidenceobsFunc> obsSeeker;

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
    printf("==>Seed value: %.2f\n", m_seed);
  }
  double factor_med = -1;
  if (m_calcExp)
    factor_med = medSeeker.seek(m_seed);

  double factor_obs = -1;
  if (m_calcObs)
    factor_obs = obsSeeker.seek(m_seed);

  if (m_verbose)
    printf("%f %f\n", factor_med, factor_obs);

  CLs.xsec_cl = m_cllevel;
  CLs.xsec_obsfactor = factor_obs;
  CLs.xsec_medfactor = factor_med;

  m_expFact = factor_med;
  m_obsFact = factor_obs;

  return true;
}

void ThreeSigmaEvidence::print() {
  printf("\n****************************************:\n");
  printf("\nThreeSigmaEvidence calculator results:\n");
  printf("==>CL level: %.4f, N sigma: %d\n", m_cllevel, m_sigma);
  printf("==>Accuracy: %.4f, Precision: %d\n", m_accuracy, m_precision);
  printf("==>Lumi scaling factor Exp: %.3f, Obs: %.3f\n", m_expFact, m_obsFact);
  printf("****************************************:\n");
  return;
}
