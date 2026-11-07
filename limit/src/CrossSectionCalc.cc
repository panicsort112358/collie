#include <stdio.h>
#include "CrossSectionCalc.hh"
#include "ProfileLH.hh"
#include <math.h>
#include "CLfast.hh"

CLpoint CrossSectionCalc::calcValue(SigBkgdDist dist, double sf) {

  int precis = CLcompute::LEVEL_VERYVERYFAST;
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

  CLpoint test(100);
  dist.scaleSignal(sf);
  p_cl->calculateCLs(dist, test, precis);
  dist.scaleSignal(1.0 / sf);

  return test;
}

/* Calculate Cross-Section Branching Ratio Value */
bool CrossSectionCalc::calculate(const SigBkgdDist& dist, CLpoint& CLs, bool doPlots) {

  if (p_cl == NULL) {
    printf("CrossSectionCalc, you must first setup a CLcompute class!\n");
    return false;
  }

  double errP = 0;
  double errM = 0;
  double parab = 0;
  double globCC = 0;
  char title[256];

  if (m_runPE) {

    SigBkgdDist scaled(dist);
    int bins = (int)((m_maximum - m_seed) / m_granularity + m_granularity / 100.0) + 1;
    double low = m_seed - m_granularity / 2.0;
    double high = m_maximum + m_granularity / 2.0;
    TH1D* clsbobs = new TH1D("1-CLsb_obs", "1-CLsb_obs", bins, low, high);
    TH1D* clsbmed = new TH1D("1-CLsb_med", "1-CLsb_med", bins, low, high);
    TH1D* clbobs = new TH1D("1-CLb_obs", "1-CLb_obs", bins, low, high);
    TH1D* clbmed = new TH1D("1-CLb_med", "1-CLb_med", bins, low, high);
    TH1D* clsmed = new TH1D("CLs_med", "CLs_med", bins, low, high);
    TH1D* clsobs = new TH1D("CLs_obs", "CLs_obs", bins, low, high);
    TH1D* llrobs = new TH1D("LLRobs", "LLRobs", bins, low, high);

    for (int s = 1; s <= bins; s++) {
      double scaler = m_seed + (s - 1) * m_granularity;
      if (m_verbose)
        printf("CrossSectionCalc, scaling xsec at %.3f\n", scaler);

      CLpoint cl = calcValue(scaled, scaler);
      clsbobs->SetBinContent(s, 1.0 - cl.clsb_obs);
      clsbmed->SetBinContent(s, 1.0 - cl.clsb_med);
      clbobs->SetBinContent(s, 1.0 - cl.clb_obs);
      clbmed->SetBinContent(s, 1.0 - cl.clb_med);
      clsobs->SetBinContent(s, cl.cls_obs);
      clsmed->SetBinContent(s, cl.cls_med);
      llrobs->SetBinContent(s, cl.llrobs);
    }
  }

  SigBkgdDist nullfit(dist);
  if (m_sbhypo)
    nullfit.fillDataSB();
  if (m_bohypo)
    nullfit.fillDataBOnly();
  m_fitValues.signalI = nullfit.totSignal();
  m_fitValues.bkgdI = nullfit.totBkgd();
  m_fitValues.data = nullfit.totData();

  //nominal templates
  const vector<string> channelNames = nullfit.getChannelNames();
  for (uint chanIdx = 0; chanIdx < channelNames.size(); chanIdx++) {
    const string& channelName = channelNames[chanIdx];

    const CollieDistribution* cs = nullfit.getSigDist(channelName, 0);
    sprintf(title, "Nominal, Sig %s", channelName.c_str());

    if (cs == NULL) {
      printf("CrossSectionCalc::calculate(), This CollieDistribution is NULL! Idx: %d, %s\n",
             chanIdx, channelName.c_str());
      return false;
    }

    int tbins = cs->getNXbins() * cs->getNYbins();
    //    printf("Tbins: %d\n",tbins);
    TH1D* nsig = new TH1D(title, title, tbins, cs->getMinX(), cs->getMaxX());
    nsig->Sumw2();

    for (uint f = 0; f < nullfit.getNsigDist(channelName); f++) {
      const CollieDistribution* css = nullfit.getSigDist(channelName, f);
      if (css == NULL) {
        printf(
            "CrossSectionCalc::calculate(), This signal CollieDistribution is NULL! Idx: %d, %s\n",
            chanIdx, channelName.c_str());
        return false;
      }
      for (int bx = 0; bx < cs->getNXbins(); bx++) {
        for (int by = 0; by < cs->getNYbins(); by++) {
          nsig->SetBinContent((bx + cs->getNXbins() * by) + 1, css->getEfficiency(bx, by));
          nsig->SetBinError((bx + cs->getNXbins() * by) + 1, css->getBinStatErr(bx, by));
        }
      }
    }

    const CollieDistribution* cba = nullfit.getBkgdDist(channelName, 0);
    sprintf(title, "Nominal, Sum Bkgd %s", channelName.c_str());
    TH1D* nbkgT = new TH1D(title, title, tbins, cba->getMinX(), cba->getMaxX());
    nbkgT->Sumw2();

    for (uint f = 0; f < nullfit.getNbkgdDist(channelName); f++) {
      const CollieDistribution* cb = nullfit.getBkgdDist(channelName, f);
      if (cb == NULL) {
        printf("CrossSectionCalc::calculate(), This bkgd CollieDistribution is NULL! Idx: %d, %s\n",
               chanIdx, channelName.c_str());
        return false;
      }

      for (int bx = 0; bx < cb->getNXbins(); bx++) {
        for (int by = 0; by < cb->getNYbins(); by++) {
          nbkgT->AddBinContent((bx + cs->getNXbins() * by) + 1, cb->getEfficiency(bx, by));
          double e1 = nbkgT->GetBinError((bx + cs->getNXbins() * by) + 1);
          double e2 = cb->getBinStatErr(bx, by);
          nbkgT->SetBinError((bx + cs->getNXbins() * by) + 1, sqrt(e1 * e1 + e2 * e2));
        }
      }
    }
  }

  //Set up the fitter
  ProfileLH pfLH;

  pfLH.sigLLR(1e10);
  pfLH.setLUN(99);

  //Now do the S+B Fit
  SigBkgdDist testfit(dist);
  if (m_sbhypo)
    testfit.fillDataSB();
  if (m_bohypo)
    testfit.fillDataBOnly();
  pfLH.fitSignal(true);
  pfLH.floatSignal(true);
  pfLH.setModel(&testfit);
  testfit.setBaselineModel();
  pfLH.setModel(&testfit);
  testfit.setBaselineModel();
  pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy
  testfit.setBaselineModel();
  pfLH.fitProfile();
  testfit.setBaselineModel();
  pfLH.fitProfile();
  pfLH.fitProfile();

  //our adjustable copy of the best fit syst params
  std::vector<double> f_paramsTEST(testfit.getNsyst(), 0.0);
  TH1D* testParamHist = new TH1D("Test Fit Params", "Test Fit Params", testfit.getNsyst(), 0, 1);
  for (int i = 0; i < testfit.getNsyst(); i++) {
    double pv = 0;
    if (testfit.getSystFitFlag(i))
      pv = pfLH.getFitParamVal(testfit.getSystName(i));

    f_paramsTEST[i] = pv;
    testParamHist->SetBinContent(i + 1, pv);
    testParamHist->GetXaxis()->SetBinLabel(i + 1, testfit.getSystName(i).c_str());
  }

  m_fitValues.chi2_sb = pfLH.getFuncMin();
  m_fitValues.status_sb = pfLH.getStatus();
  m_fitValues.sigScale = pfLH.getSignalSF();
  m_fitValues.sigScaleErrTot = pfLH.getSignalSFerr();
  //  printf("SigScales: %.3f +/- %.3f\n", m_fitValues.sigScale, m_fitValues.sigScaleErrTot);
  pfLH.getMinosErrXsec(errP, errM, parab, globCC);
  //  printf("Signal scale factor: %f ( + %f / %f)\n",m_fitValues.sigScale,errP,errM);
  m_fitValues.sigScaleErrP = errP;
  m_fitValues.sigScaleErrM = errM;

  if (errP > m_fitValues.sigScale && fabs(errM) < 1e-4)
    errM = -1 * m_fitValues.sigScale;

  testfit.scaleSignal(m_fitValues.sigScale);
  m_fitValues.signalF_sb = testfit.totSignal();
  m_fitValues.bkgdF_sb = testfit.totBkgd();

  TH2D* errMatHistT = 0;
  TAxis* ax = 0;
  TAxis* ay = 0;

  if (doPlots)
    printf("Generating TEST fit error matrix...\n");
  double errMatrixT[pfLH.getNfitSyst()][pfLH.getNfitSyst()];
  if (doPlots)
    pfLH.getErrorMatrix(&errMatrixT[0][0]);
  //Fix any wandering in the extra error prop...
  testfit.fluctuate(f_paramsTEST);

  m_fitValues.params.clear();
  m_fitValues.sigmas.clear();
  m_fitValues.names.clear();

  errMatHistT = new TH2D("TEST Fit Error Matrix", "TEST Fit Error Matrix", pfLH.getNfitSyst(), 0, 1,
                         pfLH.getNfitSyst(), 0, 1);
  ax = errMatHistT->GetXaxis();
  ay = errMatHistT->GetYaxis();

  for (int i = 0; i < testfit.getNsyst(); i++) {
    TString csyst(testfit.getSystName(i));

    vector<double> row;
    bool found1 = false;

    for (int j = 0; j < testfit.getNsyst(); j++) {
      TString dsyst(testfit.getSystName(j));

      bool found2 = false;

      for (int k = 0; k < pfLH.getNfitSyst(); k++) {
        TString tsyst(pfLH.getFitSystName(k).c_str());

        if (tsyst == csyst) {
          if (!found1) {
            pfLH.getMinosError(i, errP, errM, parab, globCC);

            if (fabs(errP) < 1.1 && fabs(errM) < 1.1)
              m_fitValues.sigmas.push_back((fabs(errP) + fabs(errM)) / 2.0);
            else
              m_fitValues.sigmas.push_back(pfLH.getFitParamErr(i));

            m_fitValues.params.push_back(pfLH.getFitParamVal(i));
            m_fitValues.names.push_back(pfLH.getFitSystName(i));
            found1 = true;
          }
          for (int l = 0; l < pfLH.getNfitSyst(); l++) {
            TString psyst(pfLH.getFitSystName(l).c_str());

            if (dsyst == psyst) {
              row.push_back(errMatrixT[i][j]);
              errMatHistT->SetBinContent(i + 1, j + 1, errMatrixT[i][j]);
              ax->SetBinLabel(i + 1, csyst.Data());
              ay->SetBinLabel(j + 1, dsyst.Data());
              found2 = true;
            }
          }
        }
      }

      if (!found2) {
        //  printf("TEST Did not find %s/%s\n", csyst.Data(),dsyst.Data());
        double mv = 0;
        if (i == j)
          mv = 1;
        row.push_back(mv);
        errMatHistT->SetBinContent(i + 1, j + 1, mv);
        ax->SetBinLabel(i + 1, csyst.Data());
        ay->SetBinLabel(j + 1, dsyst.Data());
      }
    }
    if (found1) {
      m_fitValues.emat.push_back(row);
    } else {
      printf("CrossSectionCalc::calc, Error finding systematic %s\n", csyst.Data());
    }
  }

  TH1D* sigScaleValue = new TH1D("Signal Scale Factor", "Signal Scale Factor", 1, 0, 1);
  sigScaleValue->SetBinContent(1, m_fitValues.sigScale);
  sigScaleValue->GetXaxis()->SetBinLabel(1, testfit.getSigName(0).c_str());

  //Calculate statistics-only uncertainty
  for (int i = 0; i < testfit.getNsyst() - 1; i++)
    testfit.setSystFitFlag(i, false);
  pfLH.setModel(&testfit);
  testfit.fluctuate(f_paramsTEST);
  testfit.setBaselineModel();
  testfit.fluctuate(f_paramsTEST);
  pfLH.setModel(&testfit);
  testfit.fluctuate(f_paramsTEST);
  testfit.setBaselineModel();
  testfit.fluctuate(f_paramsTEST);

  for (int i = 0; i < testfit.getNsyst(); i++) {
    testfit.setSystFluctValue(i, f_paramsTEST[i]);
    testfit.setSystCentralValue(i, f_paramsTEST[i]);
  }
  testfit.fluctuate(f_paramsTEST);
  //  pfLH.fitProfile();
  //  pfLH.fitProfile();//fit twice on the first fit to make sure MINUIT is happy
  //  printf("totSig: %f\n",testfit.totSignal());

  double nom = testfit.calculateChi2LLR(true, 1e9);
  double start = (m_fitValues.sigScale + m_fitValues.sigScaleErrM) / m_fitValues.sigScale;
  double end = (m_fitValues.sigScale + m_fitValues.sigScaleErrP) / m_fitValues.sigScale;
  double step = (end - start) / 200.0;
  errP = -1e6;
  errM = 1e6;
  double mymin = 0;
  double myminc = 1e6;
  double last = 1e6;
  // printf("nom: %f, start: %f, Stop: %f, step: %f\n",nom, start, end, step);
  for (int v = 0; v < 200; v++) {
    double mv = start + 1.0 * v * step + 1e-6;
    //    printf("Scale: %f\n",mv);
    testfit.scaleSignal(mv);

    double test = testfit.calculateChi2LLR(true, 1e9);
    if ((test - nom) < 1 && (last - nom) > 1)
      errM = mv;
    if ((test - nom) > 1 && (last - nom) < 1)
      errP = mv;

    testfit.scaleSignal(1.0 / (mv));
    //    printf("chi2: %f, %f\n",test, test-nom);
    if (test < myminc) {
      myminc = test;
      mymin = mv;
    }

    last = test;
  }
  errP = errP * m_fitValues.sigScale - m_fitValues.sigScale;
  errM = errM * m_fitValues.sigScale - m_fitValues.sigScale;

  m_fitValues.sigScaleErrStatP = errP;
  m_fitValues.sigScaleErrStatM = errM;
  //  printf("Stat test signal scale factor: %f ( + %f / %f)\n",mymin*m_fitValues.sigScale,errP,errM);

  m_fitValues.sigScaleErrStat = (errP - errM) / 2.0;

  double psyst = m_fitValues.sigScaleErrTot * m_fitValues.sigScaleErrTot -
                 m_fitValues.sigScaleErrStat * m_fitValues.sigScaleErrStat;

  if (psyst > 0)
    m_fitValues.sigScaleErrSyst = sqrt(psyst);
  else
    m_fitValues.sigScaleErrSyst = -1.0;

  psyst = m_fitValues.sigScaleErrP * m_fitValues.sigScaleErrP -
          m_fitValues.sigScaleErrStatP * m_fitValues.sigScaleErrStatP;

  if (psyst > 0)
    m_fitValues.sigScaleErrSystP = sqrt(psyst);
  else
    m_fitValues.sigScaleErrSystP = -1.0;

  psyst = m_fitValues.sigScaleErrM * m_fitValues.sigScaleErrM -
          m_fitValues.sigScaleErrStatM * m_fitValues.sigScaleErrStatM;

  if (psyst > 0)
    m_fitValues.sigScaleErrSystM = sqrt(psyst);
  else
    m_fitValues.sigScaleErrSystM = -1.0;

  testfit.setBaselineModel();
  for (int i = 0; i < testfit.getNsyst(); i++)
    testfit.setSystFitFlag(i, true);
  pfLH.setModel(&testfit);
  testfit.setBaselineModel();
  pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy
  pfLH.setModel(&testfit);
  testfit.setBaselineModel();
  pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy

  //Check bkgd-only fit first...
  pfLH.fitSignal(false);
  pfLH.floatSignal(false);
  pfLH.setModel(&nullfit);
  nullfit.setBaselineModel(); //zeros nuisance param central values...
  pfLH.setModel(&nullfit);
  nullfit.setBaselineModel(); //zeros nuisance param central values...
  pfLH.fitProfile();          //fit twice on the first fit to make sure MINUIT is happy
  nullfit.setBaselineModel();
  pfLH.fitProfile();
  nullfit.setBaselineModel();
  pfLH.fitProfile();
  pfLH.fitProfile();

  //our adjustable copy of the syst params
  std::vector<double> f_paramsNULL(nullfit.getNsyst(), 0.0);
  TH1D* nullParamHist = new TH1D("Null Fit Params", "Null Fit Params", nullfit.getNsyst(), 0, 1);
  for (int i = 0; i < nullfit.getNsyst(); i++) {
    double pv = 0;
    if (nullfit.getSystFitFlag(i))
      pv = pfLH.getFitParamVal(nullfit.getSystName(i));

    f_paramsNULL[i] = pv;
    nullParamHist->SetBinContent(i + 1, pv);
    nullParamHist->GetXaxis()->SetBinLabel(i + 1, nullfit.getSystName(i).c_str());
  }

  TH2D* errMatHistN = 0;

  if (doPlots)
    printf("Generating NULL fit error matrix...\n");
  errMatHistN = new TH2D("NULL Fit Error Matrix", "NULL Fit Error Matrix", nullfit.getNsyst(), 0, 1,
                         nullfit.getNsyst(), 0, 1);
  ax = errMatHistN->GetXaxis();
  ay = errMatHistN->GetYaxis();

  double errMatrixN[pfLH.getNfitSyst()][pfLH.getNfitSyst()];
  if (doPlots)
    pfLH.getErrorMatrix(&errMatrixN[0][0]);
  //Fix any wandering in the extra error prop...
  nullfit.fluctuate(f_paramsNULL);

  for (int i = 0; i < nullfit.getNsyst(); i++) {
    TString csyst(nullfit.getSystName(i));

    for (int j = 0; j < nullfit.getNsyst(); j++) {
      TString dsyst(nullfit.getSystName(j));

      bool found = false;

      for (int k = 0; k < pfLH.getNfitSyst(); k++) {
        TString tsyst(pfLH.getFitSystName(k).c_str());

        if (tsyst == csyst) {

          for (int l = 0; l < pfLH.getNfitSyst(); l++) {
            TString psyst(pfLH.getFitSystName(l).c_str());

            if (dsyst == psyst) {
              errMatHistN->SetBinContent(i + 1, j + 1, errMatrixN[i][j]);
              ax->SetBinLabel(i + 1, csyst.Data());
              ay->SetBinLabel(j + 1, dsyst.Data());
              found = true;
            }
          }
        }
      }

      if (!found) {
        //  printf("NULL Did not find %s/%s\n", csyst.Data(),dsyst.Data());
        double mv = 0;
        if (i == j)
          mv = 1;
        errMatHistN->SetBinContent(i + 1, j + 1, mv);
        ax->SetBinLabel(i + 1, csyst.Data());
        ay->SetBinLabel(j + 1, dsyst.Data());
      }
    }
  }

  m_fitValues.signalF_b = nullfit.totSignal();
  m_fitValues.bkgdF_b = nullfit.totBkgd();
  m_fitValues.chi2_b = pfLH.getFuncMin();
  m_fitValues.status_b = pfLH.getStatus();

  print();
  if (doPlots) {
    printf("Generating TEST fit histograms...\n");
    testfit.generateFitHistos(sigScaleValue, testParamHist, errMatHistT, "TEST");

    printf("Generating NULL fit histograms...\n");
    nullfit.generateFitHistos(sigScaleValue, nullParamHist, errMatHistN, "NULL");
  }

  CLs.fit_sigScale = m_fitValues.sigScale;
  CLs.fit_sigScale_Err = m_fitValues.sigScaleErrTot;
  CLs.fit_sigScale_ErrStat = m_fitValues.sigScaleErrStat;
  CLs.fit_sigScale_ErrSyst = m_fitValues.sigScaleErrSyst;
  CLs.fit_sigScale_ErrP = m_fitValues.sigScaleErrP;
  CLs.fit_sigScale_ErrStatP = m_fitValues.sigScaleErrStatP;
  CLs.fit_sigScale_ErrSystP = m_fitValues.sigScaleErrSystP;
  CLs.fit_sigScale_ErrM = m_fitValues.sigScaleErrM;
  CLs.fit_sigScale_ErrStatM = m_fitValues.sigScaleErrStatM;
  CLs.fit_sigScale_ErrSystM = m_fitValues.sigScaleErrSystM;
  CLs.fit_chi2_b = m_fitValues.chi2_b;
  CLs.fit_chi2_sb = m_fitValues.chi2_sb;

  return true;
}

bool CrossSectionCalc::generateFitValues(const SigBkgdDist& dist, bool doEMAT) {

  printf("Generating Fit Values\n");

  SigBkgdDist nullfit(dist);
  if (m_sbhypo)
    nullfit.fillDataSB();
  if (m_bohypo)
    nullfit.fillDataBOnly();
  m_fitValues.signalI = nullfit.totSignal();
  m_fitValues.bkgdI = nullfit.totBkgd();
  m_fitValues.data = nullfit.totData();

  //Set up the fitter
  ProfileLH pfLH;
  pfLH.sigLLR(1e10);
  pfLH.setLUN(99);

  //First do the S+B Fit
  printf("\nTest Hypothesis Fit:\n");
  SigBkgdDist testfit(dist);
  if (m_sbhypo)
    testfit.fillDataSB();
  if (m_bohypo)
    testfit.fillDataBOnly();
  pfLH.fitSignal(true);
  pfLH.floatSignal(false);
  pfLH.setModel(&testfit);
  testfit.setBaselineModel();
  pfLH.setModel(&testfit);
  testfit.setBaselineModel();
  pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy
  testfit.setBaselineModel();
  pfLH.fitProfile();
  testfit.setBaselineModel();
  pfLH.fitProfile();
  pfLH.fitProfile();
  printf("Fit status: %d\n", pfLH.getStatus());

  //our adjustable copy of the best fit syst params
  std::vector<double> f_paramsTEST(testfit.getNsyst(), 0.0);
  TH1D* testParamHist = new TH1D("Test Fit Params", "Test Fit Params", testfit.getNsyst(), 0, 1);
  TH1D* testPostHist =
      new TH1D("Test Fit Posteriors", "Test Fit Posteriors", testfit.getNsyst(), 0, 1);
  for (int i = 0; i < testfit.getNsyst(); i++) {
    double pv = 0;
    double err = 1;
    if (testfit.getSystFitFlag(i)) {
      pv = pfLH.getFitParamVal(testfit.getSystName(i));
      err = pfLH.getFitParamErr(testfit.getSystName(i));
    }

    printf("Adding %s: %.3f\n", testfit.getSystName(i).c_str(), err);
    f_paramsTEST[i] = pv;
    testParamHist->SetBinContent(i + 1, pv);
    testParamHist->GetXaxis()->SetBinLabel(i + 1, testfit.getSystName(i).c_str());
    testPostHist->SetBinContent(i + 1, err);
    testPostHist->GetXaxis()->SetBinLabel(i + 1, testfit.getSystName(i).c_str());
  }

  printf("Generating TEST fit error matrix...\n");
  double errMatrixT[pfLH.getNfitSyst()][pfLH.getNfitSyst()];
  for (int i = 0; i < pfLH.getNfitSyst(); i++) {
    for (int j = 0; j < pfLH.getNfitSyst(); j++) {
      double val = 0;
      if (i == j)
        val = 1;
      errMatrixT[i][j] = val;
    }
  }
  if (doEMAT)
    pfLH.getErrorMatrix(&errMatrixT[0][0]);

  TH2D* errMatHistT = new TH2D("TEST Fit Error Matrix", "TEST Fit Error Matrix", pfLH.getNfitSyst(),
                               0, 1, pfLH.getNfitSyst(), 0, 1);
  TAxis* ax = errMatHistT->GetXaxis();
  TAxis* ay = errMatHistT->GetYaxis();

  for (int i = 0; i < testfit.getNsyst(); i++) {
    TString csyst(testfit.getSystName(i));

    vector<double> row;
    //    bool found1 = false;

    for (int j = 0; j < testfit.getNsyst(); j++) {
      TString dsyst(testfit.getSystName(j));

      bool found2 = false;

      for (int k = 0; k < pfLH.getNfitSyst(); k++) {
        TString tsyst(pfLH.getFitSystName(k).c_str());

        if (tsyst == csyst) {
          for (int l = 0; l < pfLH.getNfitSyst(); l++) {
            TString psyst(pfLH.getFitSystName(l).c_str());

            if (dsyst == psyst) {
              row.push_back(errMatrixT[i][j]);
              errMatHistT->SetBinContent(i + 1, j + 1, errMatrixT[i][j]);
              ax->SetBinLabel(i + 1, csyst.Data());
              ay->SetBinLabel(j + 1, dsyst.Data());
              found2 = true;
            }
          }
        }
      }

      if (!found2) {
        //  printf("TEST Did not find %s/%s\n", csyst.Data(),dsyst.Data());
        double mv = 0;
        if (i == j)
          mv = 1;
        row.push_back(mv);
        errMatHistT->SetBinContent(i + 1, j + 1, mv);
        ax->SetBinLabel(i + 1, csyst.Data());
        ay->SetBinLabel(j + 1, dsyst.Data());
      }
    }
  }

  testfit.setBaselineModel();
  for (int i = 0; i < testfit.getNsyst(); i++)
    testfit.setSystFitFlag(i, true);
  pfLH.setModel(&testfit);
  testfit.setBaselineModel();
  pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy
  pfLH.setModel(&testfit);
  testfit.setBaselineModel();
  pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy

  //Check bkgd-only fit
  printf("\nNull Hypothesis Fit:\n");
  pfLH.fitSignal(false);
  pfLH.floatSignal(false);
  pfLH.setModel(&nullfit);
  nullfit.setBaselineModel(); //zeros nuisance param central values...
  pfLH.setModel(&nullfit);
  nullfit.setBaselineModel(); //zeros nuisance param central values...
  pfLH.fitProfile();          //fit twice on the first fit to make sure MINUIT is happy
  nullfit.setBaselineModel();
  pfLH.fitProfile();
  nullfit.setBaselineModel();
  pfLH.fitProfile();
  pfLH.fitProfile();
  printf("Fit status: %d\n", pfLH.getStatus());

  //our adjustable copy of the syst params
  std::vector<double> f_paramsNULL(nullfit.getNsyst(), 0.0);
  TH1D* nullParamHist = new TH1D("Null Fit Params", "Null Fit Params", nullfit.getNsyst(), 0, 1);
  TH1D* nullPostHist =
      new TH1D("Null Fit Posteriors", "Null Fit Posteriors", nullfit.getNsyst(), 0, 1);

  for (int i = 0; i < nullfit.getNsyst(); i++) {
    double pv = 0;
    double err = 1;
    if (nullfit.getSystFitFlag(i)) {
      pv = pfLH.getFitParamVal(nullfit.getSystName(i));
      err = pfLH.getFitParamErr(nullfit.getSystName(i));
    }
    printf("Adding %s: %.3f\n", nullfit.getSystName(i).c_str(), err);

    f_paramsNULL[i] = pv;
    nullParamHist->SetBinContent(i + 1, pv);
    nullParamHist->GetXaxis()->SetBinLabel(i + 1, nullfit.getSystName(i).c_str());
    nullPostHist->SetBinContent(i + 1, err);
    nullPostHist->GetXaxis()->SetBinLabel(i + 1, nullfit.getSystName(i).c_str());
  }

  TH2D* errMatHistN = 0;

  printf("Generating NULL fit error matrix...\n");
  errMatHistN = new TH2D("NULL Fit Error Matrix", "NULL Fit Error Matrix", nullfit.getNsyst(), 0, 1,
                         nullfit.getNsyst(), 0, 1);
  ax = errMatHistN->GetXaxis();
  ay = errMatHistN->GetYaxis();

  double errMatrixN[pfLH.getNfitSyst()][pfLH.getNfitSyst()];
  for (int i = 0; i < pfLH.getNfitSyst(); i++) {
    for (int j = 0; j < pfLH.getNfitSyst(); j++) {
      double val = 0;
      if (i == j)
        val = 1;
      errMatrixN[i][j] = val;
    }
  }

  if (doEMAT)
    pfLH.getErrorMatrix(&errMatrixN[0][0]);

  for (int i = 0; i < nullfit.getNsyst(); i++) {
    TString csyst(nullfit.getSystName(i));

    for (int j = 0; j < nullfit.getNsyst(); j++) {
      TString dsyst(nullfit.getSystName(j));

      bool found = false;

      for (int k = 0; k < pfLH.getNfitSyst(); k++) {
        TString tsyst(pfLH.getFitSystName(k).c_str());

        if (tsyst == csyst) {

          for (int l = 0; l < pfLH.getNfitSyst(); l++) {
            TString psyst(pfLH.getFitSystName(l).c_str());

            if (dsyst == psyst) {
              errMatHistN->SetBinContent(i + 1, j + 1, errMatrixN[i][j]);
              ax->SetBinLabel(i + 1, csyst.Data());
              ay->SetBinLabel(j + 1, dsyst.Data());
              found = true;
            }
          }
        }
      }

      if (!found) {
        //  printf("NULL Did not find %s/%s\n", csyst.Data(),dsyst.Data());
        double mv = 0;
        if (i == j)
          mv = 1;
        errMatHistN->SetBinContent(i + 1, j + 1, mv);
        ax->SetBinLabel(i + 1, csyst.Data());
        ay->SetBinLabel(j + 1, dsyst.Data());
      }
    }
  }

  printf("\n......Done!\n");

  return true;
}

bool CrossSectionCalc::testFitPE(const SigBkgdDist& dist, double sigVal, int nPE) {
  collie::PoissonRandom poisson;

  SigBkgdDist pfit(dist);
  //Set up the fitter
  ProfileLH pfLH;
  pfLH.sigLLR(1e4);
  pfLH.setLUN(99);
  pfLH.fitSignal(true);
  pfLH.floatSignal(true);
  pfLH.setSignalSF(sigVal);
  pfLH.setModel(&pfit);
  pfit.setBaselineModel(); //zeros nuisance param central values...

  pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy
  pfit.setBaselineModel();
  pfLH.fitProfile();
  pfit.setBaselineModel();

  TH1D* vals = new TH1D("Fitted Xsec", "Fitted Xsec", 2000, 0, 8);
  TH1D* sigR = new TH1D("Signal Rate", "Signal Rate", 200, 0, 1000);
  TH1D* bkgR = new TH1D("Background Rate", "Background Rate", 1000, 0, 20000);
  TH1D* dataR = new TH1D("Data Rate", "Data Rate", 1000, 0, 20000);
  TH1D* errs = new TH1D("Fit Error", "Fit Error", 200, 0, 1.5);

  double pdata[pfit.nbins()];
  for (int i = 0; i < nPE; i++) {
    if ((i % 1000) == 0)
      printf("processed %d events\n", i);

    //generate pseudodata...
    pfit.setBaselineModel();
    pfit.fluctuate();
    sigR->Fill(pfit.totSignal() * sigVal);
    bkgR->Fill(pfit.totBkgd());

    //fire off random poisson...
    for (int j = 0; j < pfit.nbins(); j++) {
      pdata[j] = poisson.fire(pfit.bkgd(j) + pfit.signal(j) * sigVal);
    }

    pfit.setBaselineModel(); //zeros nuisance param central values...
    for (int j = 0; j < pfit.nbins(); j++)
      pfit.data(j) = pdata[j];
    dataR->Fill(pfit.totData());

    pfLH.fitProfile();
    if (pfLH.getStatus() != 3) {
      printf("CrossSectionCalc::testFitPE, Fit failure (%d,%d)! Skipping iteration...\n",
             pfLH.getStatus(), i);
      continue;
    }

    vals->Fill(pfLH.getSignalSF());
    errs->Fill(pfLH.getSignalSFerr());
  }

  return true;
}

TH1D* CrossSectionCalc::get1DContour(const SigBkgdDist& dist, double xMin, double xMax, int bins) {

  TH1D* h = new TH1D("1D Cross Section ChiSq Contour", "1D Cross Section ChiSq Contour", bins, xMin,
                     xMax);
  h->GetXaxis()->SetTitle("Signal Scale");

  //Set up the fitter
  SigBkgdDist pfit(dist);
  ProfileLH pfLH;
  pfLH.sigLLR(1e4);
  pfLH.setLUN(99);
  pfLH.fitSignal(true);
  pfLH.floatSignal(false);
  pfLH.setModel(&pfit);
  pfit.setBaselineModel(); //zeros nuisance param central values...

  pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy
  pfit.setBaselineModel();
  pfLH.fitProfile();
  pfit.setBaselineModel();

  for (int i = 0; i < bins; i++) {

    double sigVal = xMin + i * (xMax - xMin) / bins + 1e-5;

    pfit.scaleSignal(sigVal);

    pfit.fluctuate();
    pfit.setBaselineModel();
    pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy
    pfit.setBaselineModel();
    pfLH.fitProfile();

    h->SetBinContent(i + 1, pfLH.getFuncMin());
    printf("value: %f, %f\n", sigVal, pfLH.getFuncMin());
    pfit.scaleSignal(1 / sigVal);
  }

  return h;
}

TH2D* CrossSectionCalc::get2DContour(const SigBkgdDist& dist, string sname, double xMin,
                                     double xMax, int binsx, double yMin, double yMax, int binsy) {

  TH2D* h =
      new TH2D("2D Cross Section/Systematic ChiSq Contour",
               "2D Cross Section/Systematic ChiSq Contour", binsx, xMin, xMax, binsy, yMin, yMax);

  h->GetXaxis()->SetTitle("Signal Scale");
  h->GetYaxis()->SetTitle(sname.c_str());

  //Set up the fitter
  SigBkgdDist pfit(dist);
  pfit.setFloatFlag(sname, false);
  ProfileLH pfLH;
  pfLH.sigLLR(1e4);
  pfLH.setLUN(99);
  pfLH.fitSignal(true);
  pfLH.floatSignal(false);
  pfLH.setModel(&pfit);
  pfit.setBaselineModel(); //zeros nuisance param central values...

  pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy
  pfit.setBaselineModel();
  pfLH.fitProfile();
  pfit.setBaselineModel();

  std::vector<double> f_paramsTEST(pfit.getNsyst(), 0.0);

  int idx = pfit.getSystIndex(sname);

  //  printf("A: %f\n",pfit.totBkgd());
  pfit.setSystFitFlag(sname, false);
  pfLH.setModel(&pfit);
  //  printf("B: %f\n",pfit.totBkgd());

  for (int x = 0; x < binsx; x++) {
    for (int y = 0; y < binsy; y++) {

      double sigVal = xMin + x * (xMax - xMin) / binsx + 1e-5;
      double systval = yMin + y * (yMax - yMin) / binsy + 1e-5;

      pfit.scaleSignal(sigVal);

      f_paramsTEST[idx] = systval;
      //      printf("setting %d, %f\n",idx,systval);
      pfit.setBaselineModel();
      //      printf("1: %f\n",pfit.totBkgd());
      pfit.fluctuate(f_paramsTEST);
      // printf("2: %f\n",pfit.totBkgd());
      //      pfit.setBaselineModel();
      pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy
      // printf("3: %f\n",pfit.totBkgd());
      //  pfit.setBaselineModel();
      pfLH.fitProfile();
      // printf("4: %f\n",pfit.totBkgd());
      h->SetBinContent(x + 1, y + 1, pfLH.getFuncMin());
      // printf("value: %f, %f\n",systval,pfLH.getFuncMin());
      pfit.scaleSignal(1 / sigVal);

      pfit.setBaselineModel();
      //      printf("5: %f\n",pfit.totBkgd());
      pfit.fluctuate();
    }
  }
  pfit.setSystFitFlag(sname, true);

  return h;
}

void CrossSectionCalc::printEMAT() {
  printf("\n****************************************\n");
  printf("CrossSectionCalc error matrix:\n");
  printf("%15s", "");
  for (uint i = 0; i < m_fitValues.names.size(); i++)
    printf("%15s", m_fitValues.names[i].c_str());
  printf("\n");
  for (uint i = 0; i < m_fitValues.names.size(); i++) {
    printf("%15s", m_fitValues.names[i].c_str());
    for (uint j = 0; j < m_fitValues.names.size(); j++)
      printf("%15f", m_fitValues.emat[i][j]);
    printf("\n");
  }

  printf("****************************************\n");
}

void CrossSectionCalc::print() {

  printf("\n****************************************\n");
  printf("CrossSectionCalc results:\n");
  if (m_runPE) {
    printf("==>Precision: %d, Granularity: %.3f\n", m_precision, m_granularity);
    printf("==>Scaling factor Min: %.3f, Max: %.3f\n", m_seed, m_maximum);
  }
  printf("==>Data: %f\n", m_fitValues.data);
  printf("Prefit==> Signal: %.3f, Bkgd: %.3f\n", m_fitValues.signalI, m_fitValues.bkgdI);
  printf("B-Only Fit==> status: %d, Signal: ---, Bkgd: %.3f, Chi2: %.3f\n", m_fitValues.status_b,
         m_fitValues.bkgdF_b, m_fitValues.chi2_b);
  printf("S+B Fit=====> status: %d, Signal: %.3f, Bkgd: %.3f, Chi2: %.3f\n", m_fitValues.status_sb,
         m_fitValues.signalF_sb, m_fitValues.bkgdF_sb, m_fitValues.chi2_sb);

  for (uint i = 0; i < m_fitValues.names.size(); i++) {
    printf("==>S+B Fit Param: %14s, Value: %.4f,\t Error: %.4f\n", m_fitValues.names[i].c_str(),
           m_fitValues.params[i], m_fitValues.sigmas[i]);
  }
  printf("> Signal Scale Factor: %.3f, Err: %.3f [ %.3f (stat), %.3f (syst) ]\n",
         m_fitValues.sigScale, m_fitValues.sigScaleErrTot, m_fitValues.sigScaleErrStat,
         m_fitValues.sigScaleErrSyst);
  printf("> Assym Err: +%.3f [ %.3f (stat), %.3f (syst) ], %.3f [ %.3f (stat), -%.3f (syst) ]\n",
         m_fitValues.sigScaleErrP, m_fitValues.sigScaleErrStatP, m_fitValues.sigScaleErrSystP,
         m_fitValues.sigScaleErrM, m_fitValues.sigScaleErrStatM, m_fitValues.sigScaleErrSystM);
  printf("****************************************\n");

  return;
}
