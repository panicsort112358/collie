#ifndef FitTest_hh_included
#define FitTest_hh_included

#include "SigBkgdDist.hh"
#include "CLfast.hh"
#include "ProfileLH.hh"
#include <sys/time.h>

/** Profile Likelihood fit testing program. **/

/**
   This class performs a series of tests on the fit performed in
   either CLfit or CLfit2.  The tests are based on the systematics
   provided and their correlations.
**/

class FitTest {
public:
  /// constructor
  FitTest();
  //  ~FitTest();

  inline void setIterations(int niter) {
    m_iterations = niter;
  }
  inline int getIterations(void) {
    return m_iterations;
  }
  inline void testPE(bool choice) {
    m_fillPE = choice;
  }

  void runTest(SigBkgdDist* sbd, double llrcut, bool doPE = true, bool doNM1 = true);

  //  void runSystematicsTest(SigBkgdDist* sbd);

  //Benchmark timer for fit timing
  void benchMark(SigBkgdDist* sbd, int navg = 100, bool fitSig = true);

private:
  int m_iterations;
  bool m_fillPE;

  float m_chiSB;
  float m_chiB;
  int m_maxSyst;
  int m_maxSystSB;
  int m_maxSystB;

  void runPETest(ProfileLH& pflh, SigBkgdDist* sbd, double llrcut);

  void runNM1Test(ProfileLH& pflh, SigBkgdDist* sbd, double llrcut);
};

#endif
