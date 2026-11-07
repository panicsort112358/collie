#include <stdio.h>
#include "CrossSectionCalc2D.hh"
#include "ProfileLH_2D.hh"
#include <math.h>
#include "CLfast.hh"

CLpoint CrossSectionCalc2D::calcValue(SigBkgdDist dist, double sf) {

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
bool CrossSectionCalc2D::calculate(const SigBkgdDist& dist, CLpoint& CLs, bool doPlots) {

  if (p_cl == NULL) {
    printf("CrossSectionCalc2D, you must first setup a CLcompute class!\n");
    return false;
  }

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
        printf("CrossSectionCalc2D, scaling xsec at %.3f\n", scaler);

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

  double errP = 0;
  double errM = 0;
  double parab = 0;
  double globCC = 0;
  char title[256];

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

    const CollieDistribution* csa = nullfit.getSigDist(channelName, nullfit.getSigName(0));
    sprintf(title, "Nominal, Sig %s", channelName.c_str());
    int tbins = csa->getNXbins() * csa->getNYbins();
    TH1D* nsig = new TH1D(title, title, tbins, csa->getMinX(), csa->getMaxX());
    nsig->Sumw2();

    for (uint sig = 0; sig < nullfit.getNsigDist(channelName); sig++) {
      const CollieDistribution* cs = nullfit.getSigDist(channelName, nullfit.getSigName(sig));

      for (int bx = 0; bx < cs->getNXbins(); bx++) {
        for (int by = 0; by < cs->getNYbins(); by++) {
          double e1 = nsig->GetBinError((bx + cs->getNXbins() * by) + 1);
          e1 = e1 * e1 * nsig->GetBinContent((bx + cs->getNXbins() * by) + 1);
          double e2 = cs->getBinStatErr(bx, by);
          e2 = e2 * e2 * cs->getEfficiency(bx, by);
          double totv =
              cs->getEfficiency(bx, by) + nsig->GetBinContent((bx + cs->getNXbins() * by) + 1);
          nsig->AddBinContent((bx + cs->getNXbins() * by) + 1, cs->getEfficiency(bx, by));
          nsig->SetBinError((bx + cs->getNXbins() * by) + 1, sqrt((e1 + e2) / totv));
        }
      }
    }

    const CollieDistribution* cba = nullfit.getBkgdDist(channelName, 0);
    sprintf(title, "Nominal, Sum Bkgd %s", channelName.c_str());
    TH1D* nbkgT = new TH1D(title, title, tbins, cba->getMinX(), cba->getMaxX());
    nbkgT->Sumw2();
    for (uint f = 0; f < nullfit.getNbkgdDist(channelName); f++) {
      const CollieDistribution* cb = nullfit.getBkgdDist(channelName, f);

      for (int bx = 0; bx < cb->getNXbins(); bx++) {
        for (int by = 0; by < cb->getNYbins(); by++) {
          double e1 = nsig->GetBinError((bx + cb->getNXbins() * by) + 1);
          e1 = e1 * e1 * nbkgT->GetBinContent((bx + cb->getNXbins() * by) + 1);
          double e2 = cb->getBinStatErr(bx, by);
          e2 = e2 * e2 * cb->getEfficiency(bx, by);
          double totv =
              cb->getEfficiency(bx, by) + nbkgT->GetBinContent((bx + cb->getNXbins() * by) + 1);
          nbkgT->AddBinContent((bx + cb->getNXbins() * by) + 1, cb->getEfficiency(bx, by));
          nbkgT->SetBinError((bx + cb->getNXbins() * by) + 1, sqrt((e1 + e2) / totv));
        }
      }
    }
  }

  //Set up the fitter
  ProfileLH_2D pfLH;

  pfLH.sigLLR(1e10);
  pfLH.setLUN(99);

  //Check bkgd-only fit first...
  pfLH.fitSignal(false);
  pfLH.setModel(&nullfit);
  nullfit.setBaselineModel(); //zeros nuisance param central values...
  pfLH.setModel(&nullfit);
  nullfit.setBaselineModel(); //zeros nuisance param central values...
  pfLH.fitProfile();          //fit twice on the first fit to make sure MINUIT is happy
  nullfit.setBaselineModel();
  pfLH.fitProfile();
  nullfit.setBaselineModel();
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
    nullParamHist->GetXaxis()->SetBinLabel(i + 1, pfLH.getFitSystName(i).c_str());
  }

  m_fitValues.signalF_b = nullfit.totSignal();
  m_fitValues.bkgdF_b = nullfit.totBkgd();
  m_fitValues.chi2_b = pfLH.getFuncMin();
  m_fitValues.status_b = pfLH.getStatus();

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
  m_fitValues.sigScale1 = pfLH.getSignalSF1();
  m_fitValues.sigScaleErrTot1 = pfLH.getSignalSFerr1();
  m_fitValues.sigScale2 = pfLH.getSignalSF2();
  m_fitValues.sigScaleErrTot2 = pfLH.getSignalSFerr2();

  //  printf("SigScales: %.3f +/- %.3f; %.3f +/- %.3f\n", m_fitValues.sigScale1, m_fitValues.sigScaleErrTot1,m_fitValues.sigScale2,m_fitValues.sigScaleErrTot2);

  testfit.scaleSignal(testfit.getSigName(0), m_fitValues.sigScale1, f_paramsTEST);
  testfit.scaleSignal(testfit.getSigName(1), m_fitValues.sigScale2, f_paramsTEST);
  m_fitValues.signalF_sb = testfit.totSignal();
  m_fitValues.bkgdF_sb = testfit.totBkgd();

  double errMatrix[pfLH.getNfitSyst()][pfLH.getNfitSyst()];
  if (doPlots)
    pfLH.getErrorMatrix(&errMatrix[0][0]);
  //Fix any wandering in the extra error prop...
  testfit.fluctuate(f_paramsTEST);

  m_fitValues.params.clear();
  m_fitValues.sigmas.clear();
  m_fitValues.names.clear();
  TH2D* errMatHist =
      new TH2D("Error Matrix", "ErrorMatrix", pfLH.getNfitSyst(), 0, 1, pfLH.getNfitSyst(), 0, 1);
  TAxis* ax = errMatHist->GetXaxis();
  TAxis* ay = errMatHist->GetYaxis();

  for (int i = 0; i < pfLH.getNfitSyst(); i++) {
    pfLH.getMinosError(i, errP, errM, parab, globCC);
    m_fitValues.params.push_back(pfLH.getFitParamVal(i));
    //   m_fitValues.sigmas.push_back((fabs(errP)+fabs(errM))/2.0);
    m_fitValues.sigmas.push_back(pfLH.getFitParamErr(i));
    m_fitValues.names.push_back(pfLH.getFitSystName(i));

    vector<double> row;
    for (int j = 0; j < pfLH.getNfitSyst(); j++) {
      row.push_back(errMatrix[i][j]);
      errMatHist->SetBinContent(i + 1, j + 1, errMatrix[i][j]);
      ax->SetBinLabel(i + 1, pfLH.getFitSystName(i).c_str());
      ay->SetBinLabel(j + 1, pfLH.getFitSystName(j).c_str());
    }
    m_fitValues.emat.push_back(row);
  }

  TH1D* sigScaleValue = new TH1D("Signal Scale Factor", "Signal Scale Factor", 2, 0, 1);
  sigScaleValue->SetBinContent(1, m_fitValues.sigScale1);
  sigScaleValue->SetBinContent(2, m_fitValues.sigScale2);
  sigScaleValue->GetXaxis()->SetBinLabel(1, testfit.getSigName(0).c_str());
  sigScaleValue->GetXaxis()->SetBinLabel(2, testfit.getSigName(1).c_str());

  //Calculate statistics-only uncertainty
  testfit.setBaselineModel();
  for (int i = 0; i < testfit.getNsyst(); i++)
    testfit.setSystFitFlag(i, false);
  pfLH.setModel(&testfit);
  testfit.fluctuate(f_paramsTEST);
  testfit.setBaselineModel();
  testfit.fluctuate(f_paramsTEST);
  pfLH.setModel(&testfit);
  pfLH.fitProfile();
  testfit.setBaselineModel();

  for (int i = 0; i < testfit.getNsyst(); i++) {
    testfit.setSystFluctValue(i, f_paramsTEST[i]);
    testfit.setSystCentralValue(i, f_paramsTEST[i]);
  }
  pfLH.fitProfile();
  pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy

  //  printf("Stat only chi2: %f\n",pfLH.getFuncMin());
  m_fitValues.sigScaleErrStat1 = pfLH.getSignalSFerr1() * m_fitValues.sigScale1;
  m_fitValues.sigScaleErrStat2 = pfLH.getSignalSFerr2() * m_fitValues.sigScale2;

  //  printf("Stat SigScales: %.3f +/- %.3f; %.3f +/- %.3f\n", pfLH.getSignalSF1(), m_fitValues.sigScaleErrStat1,pfLH.getSignalSF2(), m_fitValues.sigScaleErrStat2);

  double psyst = m_fitValues.sigScaleErrTot1 * m_fitValues.sigScaleErrTot1 -
                 m_fitValues.sigScaleErrStat1 * m_fitValues.sigScaleErrStat1;

  if (psyst > 0)
    m_fitValues.sigScaleErrSyst1 = sqrt(psyst);
  else
    m_fitValues.sigScaleErrSyst1 = -1.0;

  psyst = m_fitValues.sigScaleErrTot2 * m_fitValues.sigScaleErrTot2 -
          m_fitValues.sigScaleErrStat2 * m_fitValues.sigScaleErrStat2;

  if (psyst > 0)
    m_fitValues.sigScaleErrSyst2 = sqrt(psyst);
  else
    m_fitValues.sigScaleErrSyst2 = -1.0;

  testfit.setBaselineModel();
  for (int i = 0; i < testfit.getNsyst(); i++)
    testfit.setSystFitFlag(i, true);
  pfLH.setModel(&testfit);
  testfit.setBaselineModel();
  pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy
  pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy

  if (doPlots) {
    testfit.generateFitHistos(sigScaleValue, testParamHist, errMatHist, "TEST");
    nullfit.generateFitHistos(sigScaleValue, nullParamHist, errMatHist, "NULL");
  }

  return true;
}

bool CrossSectionCalc2D::testFitPE(const SigBkgdDist& dist, double sigVal1, double sigVal2,
                                   int nPE) {

  sigVal1 += 1e-6;
  sigVal2 += 1e-6;

  collie::PoissonRandom poisson;

  SigBkgdDist pfit(dist);
  //Set up the fitter
  ProfileLH_2D pfLH;
  pfLH.sigLLR(1e6);
  pfLH.setLUN(99);
  pfLH.fitSignal(true);
  pfLH.floatSignal(true);
  pfLH.setSignalSF1(sigVal1);
  pfLH.setSignalSF2(sigVal2);
  pfLH.setModel(&pfit);
  pfit.setBaselineModel(); //zeros nuisance param central values...
  pfLH.fitProfile();       //fit twice on the first fit to make sure MINUIT is happy
  pfit.setBaselineModel();
  pfLH.fitProfile();
  pfit.setBaselineModel();

  TH2D* vals = new TH2D("Fitted Xsec", "Fitted Xsec", 500, 0, 5, 500, 0, 5);
  TH1D* sigR1 = new TH1D("Signal Rate 1", "Signal Rate 1", 200, 0, 1000);
  //  TH1D* sigR2 = new TH1D("Signal Rate 2","Signal Rate 2",200,0,1000);
  TH1D* bkgR = new TH1D("Background Rate", "Background Rate", 1000, 0, 20000);
  TH1D* dataR = new TH1D("Data Rate", "Data Rate", 1000, 0, 20000);
  TH1D* errs1 = new TH1D("Fit Error 1", "Fit Error 1", 200, 0, 1.5);
  TH1D* errs2 = new TH1D("Fit Error 2", "Fit Error 2", 200, 0, 1.5);
  double pdata[pfit.nbins()];
  for (int i = 0; i < nPE; i++) {
    if ((i % 1000) == 0)
      printf("processed %d events\n", i);

    //generate pseudodata...
    pfit.setBaselineModel();
    pfit.scaleSignal(pfit.getSigName(0), sigVal1);
    pfit.scaleSignal(pfit.getSigName(1), sigVal2);
    pfit.fluctuate();
    sigR1->Fill(pfit.totSignal());
    bkgR->Fill(pfit.totBkgd());

    //fire off random poisson...
    for (int j = 0; j < pfit.nbins(); j++) {
      pdata[j] = poisson.fire(pfit.bkgd(j) + pfit.signal(j));
    }
    pfit.scaleSignal(pfit.getSigName(0), 1.0 / sigVal1);
    pfit.scaleSignal(pfit.getSigName(1), 1.0 / sigVal2);

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
    //    printf("Xsec: %f, %f\n",pfLH.getSignalSF1(),pfLH.getSignalSF2());
    vals->Fill(pfLH.getSignalSF1(), pfLH.getSignalSF2());
    errs1->Fill(pfLH.getSignalSFerr1());
    errs2->Fill(pfLH.getSignalSFerr2());
  }

  return true;
}

TH2D* CrossSectionCalc2D::get2DContour(const SigBkgdDist& dist, double xMin, double xMax, int binsx,
                                       double yMin, double yMax, int binsy) {

  TH2D* h = new TH2D("2D Cross Section ChiSq Contour", "2D Cross Section ChiSq Contour", binsx,
                     xMin, xMax, binsy, yMin, yMax);

  //Set up the fitter
  SigBkgdDist pfit(dist);
  h->GetXaxis()->SetTitle(pfit.getSigName(0).c_str());
  h->GetYaxis()->SetTitle(pfit.getSigName(1).c_str());

  ProfileLH_2D pfLH;
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

  for (int x = 1; x <= h->GetNbinsX(); x++) {
    for (int y = 1; y <= h->GetNbinsY(); y++) {

      double sigVal1 = h->GetXaxis()->GetBinCenter(x) + 1e-5;
      double sigVal2 = h->GetYaxis()->GetBinCenter(y) + 1e-5;

      pfit.scaleSignal(pfit.getSigName(0), sigVal1);
      pfit.scaleSignal(pfit.getSigName(1), sigVal2);

      pfit.fluctuate();
      pfit.setBaselineModel();
      pfLH.fitProfile(); //fit twice on the first fit to make sure MINUIT is happy
      //      pfit.setBaselineModel();
      pfLH.fitProfile();

      h->SetBinContent(x, y, pfLH.getFuncMin());
      //      printf("value: %.3f, %.3f, %.3f, %.2f, %.2f\n",sigVal1,sigVal2, pfLH.getFuncMin(),pfit.totBkgd(),pfit.totSignal());
      pfit.scaleSignal(pfit.getSigName(0), 1.0 / sigVal1);
      pfit.scaleSignal(pfit.getSigName(1), 1.0 / sigVal2);
    }
  }

  return h;
}

void CrossSectionCalc2D::printEMAT() {
  printf("\n****************************************\n");
  printf("CrossSectionCalc2D error matrix:\n");
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

void CrossSectionCalc2D::print() {

  printf("\n****************************************\n");
  printf("CrossSectionCalc2D results:\n");
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
  printf("> Signal1 Scale Factor: %.3f, Err: %.3f [%.3f (stat), %.3f (syst)]\n",
         m_fitValues.sigScale1, m_fitValues.sigScaleErrTot1, m_fitValues.sigScaleErrStat1,
         m_fitValues.sigScaleErrSyst1);
  printf("> Signal2 Scale Factor: %.3f, Err: %.3f [%.3f (stat), %.3f (syst)]\n",
         m_fitValues.sigScale2, m_fitValues.sigScaleErrTot2, m_fitValues.sigScaleErrStat2,
         m_fitValues.sigScaleErrSyst2);

  /*
  for(uint i=0; i<m_fitValues.names.size(); i++)
    printf("systNames.push_back(\"%s\"); fitVals.push_back(%f); fitConstraints.push_back(%f);\n",m_fitValues.names[i].c_str(), m_fitValues.params[i], m_fitValues.sigmas[i]);
  */
  printf("****************************************\n");

  return;
}
