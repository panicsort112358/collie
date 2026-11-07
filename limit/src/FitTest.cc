#include "FitTest.hh"

#include <cstdio>
#include <cstdlib>

namespace {
void requireFlagUpdate(bool success) {
  if (!success) {
    fprintf(stderr, "FitTest: failed to update a fit flag\n");
    std::abort();
  }
}
} // namespace

FitTest::FitTest() {
  m_iterations = 5000;
  m_fillPE = false;
  m_maxSyst = 0;
  m_maxSystB = 0;
  m_maxSystSB = 0;
}

void FitTest::runNM1Test(ProfileLH& pfLH, SigBkgdDist* sbd, double llrcut) {
  printf("Collie Fit Test: Performing N-1 tests...\n");

  sbd->setBaselineModel();
  int bins = sbd->nbins();
  double min = sbd->min();
  double max = sbd->max();

  vector<string> chanNames = sbd->getChannelNames();
  TH2D* nMinusOneSystChanSB = new TH2D(
      "S+B N Minus One Test, Channels & Systs", "S+B N Minus One Test, Channels & Systs",
      m_maxSyst + 1, 0.5, m_maxSyst + 1.5, chanNames.size() + 1, 0.5, chanNames.size() + 1.5);
  TH2D* nMinusOneSystChanB = new TH2D(
      "B-Only N Minus One Test, Channels & Systs", "B-Only N Minus One Test, Channels & Systs",
      m_maxSyst + 1, 0.5, m_maxSyst + 1.5, chanNames.size() + 1, 0.5, chanNames.size() + 1.5);

  TH1D* nMinusOneIndivB[chanNames.size() + 1];
  TH1D* nMinusOneIndivSB[chanNames.size() + 1];

  TH1D* nMinusOneChannelsB = new TH1D("B-Only N-1 Test, Channels", "B-Only N-1 Test, Channels",
                                      chanNames.size() + 1, 0.5, chanNames.size() + 1.5);
  TH1D* nMinusOneChannelsSB = new TH1D("S+B N-1 Test, Channels", "S+B N-1 Test, Channels",
                                       chanNames.size() + 1, 0.5, chanNames.size() + 1.5);

  TH1D* nMinusOneBinsB = new TH1D("B-Only N-1 Test, Bins", "B-Only N-1 Test, Bins", bins, min, max);
  TH1D* nMinusOneBinsSB = new TH1D("S+B N-1 Test, Bins", "S+B N-1 Test, Bins", bins, min, max);

  TH1D* nMinusOneSystLLR =
      new TH1D("Systematics DeltaChi2 Test - 1", "Systematics DeltaChi2 Test - 1", m_maxSyst + 2,
               0.5, m_maxSyst + 2.5);
  TH2D* nMinusTwoSystLLR =
      new TH2D("Systematics DeltaChi2 Test - 2", "Systematics DeltaChi2 Test - 2", m_maxSyst + 2,
               0.5, m_maxSyst + 2.5, m_maxSyst + 2, 0.5, m_maxSyst + 2.5);

  int maxSources = sbd->getNbkgds() + sbd->getNsignals();

  TH1D* nMinusOneSourceLLR =
      new TH1D("Event Source DeltaChi2 Test - 1", "Event Source DeltaChi2 Test - 1", maxSources + 1,
               0.5, maxSources + 1.5);
  TH2D* nMinusTwoSourceLLR =
      new TH2D("Event Source DeltaChi2 Test - 2", "Event Source DeltaChi2 Test - 2", maxSources + 1,
               0.5, maxSources + 1.5, m_maxSyst + 1, 0.5, m_maxSyst + 1.5);

  nMinusOneIndivB[0] = new TH1D("B-Only N-1 Test, Systematics", "B-Only N-1 Test, Systematics",
                                m_maxSyst + 1, 0.5, m_maxSyst + 1.5);
  nMinusOneIndivSB[0] = new TH1D("S+B N-1 Test, Systematics", "S+B N-1 Test, Systematics",
                                 m_maxSyst + 1, 0.5, m_maxSyst + 1.5);

  char dest[256];
  for (uint c = 0; c < chanNames.size(); ++c) {
    sprintf(dest, "B-Only N-1 Test, Remove %s", chanNames[c].c_str());
    nMinusOneIndivB[c + 1] = new TH1D(dest, dest, m_maxSyst + 1, 0.5, m_maxSyst + 1.5);
    sprintf(dest, "S+B N-1 Test, Remove %s", chanNames[c].c_str());
    nMinusOneIndivSB[c + 1] = new TH1D(dest, dest, m_maxSyst + 1, 0.5, m_maxSyst + 1.5);
  }

  //Perform B-Only N-1 tests
  printf("Collie Fit Test: Performing B-Only N-1 syst/channel tests...\n");
  pfLH.setModel(sbd);
  pfLH.fitSignal(false);
  pfLH.sigLLR(llrcut);
  pfLH.fitProfile();
  sbd->setBaselineModel(); //wf108 add//go back to nominal prediction
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3) {
    printf("-->FitTest, Failed fit on B-Only nominal chi2 test!  Failing...\n");
    return;
  }

  double nomChi2 = pfLH.getFuncMin();
  sbd->setBaselineModel(); //wf108 add//go back to nominal prediction

  //Test deltaChi2 dropping one event source fit at a time...
  for (int source = -1; source < sbd->getNbkgds(); ++source) {
    if (source > -1)
      requireFlagUpdate(sbd->setBkgdFitFlag(source, false));

    pfLH.setModel(sbd);
    sbd->setBaselineModel();
    pfLH.fitProfile();

    sbd->setBaselineModel();
    for (int b = 0; b < sbd->nbins(); ++b)
      sbd->data(b) = sbd->bkgd(b) + sbd->signal(b);
    pfLH.fitProfile();
    pfLH.fitProfile();
    //B-only fit to S+B data
    nMinusOneSourceLLR->Fill(source + 2, pfLH.getFuncMin());

    if (source == -1)
      nMinusOneSourceLLR->GetXaxis()->SetBinLabel(source + 2, "None");
    else
      nMinusOneSourceLLR->GetXaxis()->SetBinLabel(source + 2, sbd->getBkgdName(source).c_str());

    if (source == -1)
      nMinusTwoSourceLLR->GetXaxis()->SetBinLabel(source + 2, "None");
    else
      nMinusTwoSourceLLR->GetXaxis()->SetBinLabel(source + 2, sbd->getBkgdName(source).c_str());

    for (int iSyst = -1; iSyst < sbd->getNsyst(); ++iSyst) {
      if (iSyst > -1)
        requireFlagUpdate(sbd->setSystFitFlag(iSyst, false));

      pfLH.setModel(sbd);
      sbd->setBaselineModel();
      pfLH.fitProfile();

      sbd->setBaselineModel();
      for (int b = 0; b < sbd->nbins(); ++b)
        sbd->data(b) = sbd->bkgd(b) + sbd->signal(b);
      pfLH.fitProfile();
      pfLH.fitProfile();
      //B-only fit to S+B data
      nMinusTwoSourceLLR->Fill(source + 2, iSyst + 2, pfLH.getFuncMin());

      if (iSyst > -1)
        requireFlagUpdate(sbd->setSystFitFlag(iSyst, true));
    }

    if (source > -1)
      requireFlagUpdate(sbd->setBkgdFitFlag(source, true));
  }

  int offset = 2 + sbd->getNbkgds();
  for (int source = 0; source < sbd->getNsignals(); ++source) {
    requireFlagUpdate(sbd->setSigFitFlag(source, false));
    pfLH.setModel(sbd);
    sbd->setBaselineModel();
    pfLH.fitProfile();

    sbd->setBaselineModel();
    for (int b = 0; b < sbd->nbins(); ++b)
      sbd->data(b) = sbd->bkgd(b) + sbd->signal(b);
    pfLH.fitProfile();
    pfLH.fitProfile();
    //B-only fit to S+B data
    nMinusOneSourceLLR->Fill(source + offset, pfLH.getFuncMin());

    nMinusOneSourceLLR->GetXaxis()->SetBinLabel(source + offset, sbd->getSigName(source).c_str());
    nMinusTwoSourceLLR->GetXaxis()->SetBinLabel(source + offset, sbd->getSigName(source).c_str());

    for (int iSyst = -1; iSyst < sbd->getNsyst(); ++iSyst) {
      if (iSyst > -1)
        requireFlagUpdate(sbd->setSystFitFlag(iSyst, false));

      if (iSyst > -1)
        nMinusTwoSourceLLR->GetYaxis()->SetBinLabel(iSyst + 2, sbd->getSystName(iSyst).c_str());
      else
        nMinusTwoSourceLLR->GetYaxis()->SetBinLabel(iSyst + 2, "None");

      pfLH.setModel(sbd);
      sbd->setBaselineModel();
      pfLH.fitProfile();

      sbd->setBaselineModel();
      for (int b = 0; b < sbd->nbins(); ++b)
        sbd->data(b) = sbd->bkgd(b) + sbd->signal(b);
      pfLH.fitProfile();
      pfLH.fitProfile();
      //B-only fit to S+B data
      nMinusTwoSourceLLR->Fill(source + offset, iSyst + 2, pfLH.getFuncMin());

      if (iSyst > -1)
        requireFlagUpdate(sbd->setSystFitFlag(iSyst, true));
    }

    requireFlagUpdate(sbd->setSigFitFlag(source, true));
  }

  int nChan = chanNames.size();
  for (int iSyst = -1; iSyst < sbd->getNsyst(); ++iSyst) {
    for (int iChan = -1; iChan < nChan; ++iChan) {
      if (nChan == 1 && iChan >= 0)
        continue;
      if (nChan > 1)
        sbd->excludeChannel(iChan);

      if (iSyst > -1)
        requireFlagUpdate(sbd->setSystFitFlag(iSyst, false));

      pfLH.setModel(sbd);
      sbd->setBaselineModel();
      pfLH.setModel(sbd);
      sbd->setBaselineModel();
      pfLH.fitProfile();
      pfLH.fitProfile();

      if (pfLH.getStatus() != 3) {
        sbd->setBaselineModel();
        pfLH.fitProfile();
        if (pfLH.getStatus() != 3) {
          printf("-->FitTest, Failed fit on B-Only N-1 channel/syst test!  Continuing...\n");
          if (iSyst >= 0)
            printf("---->Failed while excluding syst=%s\n", sbd->getSystName(iSyst).c_str());
          if (iChan >= 0)
            printf("---->Failed while excluding channel=%s\n", sbd->getChannelName(iChan).c_str());
          continue;
        }
      }

      nMinusOneSystChanB->Fill(iSyst + 2, iChan + 2, pfLH.getFuncMin() - nomChi2);
      nMinusOneIndivB[iChan + 1]->Fill(iSyst + 2, pfLH.getFuncMin() - nomChi2);

      if (iChan == -1) {
        //B-only fit to S+B data
        sbd->setBaselineModel();
        for (int b = 0; b < sbd->nbins(); ++b)
          sbd->data(b) = sbd->bkgd(b) + sbd->signal(b);
        pfLH.fitProfile();
        pfLH.fitProfile();
        nMinusOneSystLLR->Fill(iSyst + 2, pfLH.getFuncMin());

        //Look at 2D systematics correlations and their impact...
        for (int systB = -1; systB < sbd->getNsyst(); systB++) {
          if (systB < iSyst)
            continue;
          if (systB > -1 && iSyst != systB)
            requireFlagUpdate(sbd->setSystFitFlag(systB, false));

          pfLH.setModel(sbd);
          sbd->setBaselineModel();
          pfLH.fitProfile();

          //B-only fit to S+B data
          sbd->setBaselineModel();
          for (int b = 0; b < sbd->nbins(); ++b)
            sbd->data(b) = sbd->bkgd(b) + sbd->signal(b);
          pfLH.fitProfile();
          pfLH.fitProfile();
          nMinusTwoSystLLR->Fill(iSyst + 2, systB + 2, pfLH.getFuncMin());

          if (systB > -1 && iSyst != systB)
            requireFlagUpdate(sbd->setSystFitFlag(systB, true));
        }
      }
      if (iSyst > -1) {
        nMinusOneSystChanB->GetXaxis()->SetBinLabel(iSyst + 2, sbd->getSystName(iSyst).c_str());
        nMinusOneIndivB[iChan + 1]->GetXaxis()->SetBinLabel(iSyst + 2,
                                                            sbd->getSystName(iSyst).c_str());
      } else {
        nMinusOneSystChanB->GetXaxis()->SetBinLabel(iSyst + 2, "None");
        nMinusOneIndivB[iChan + 1]->GetXaxis()->SetBinLabel(iSyst + 2, "None");
        if (iChan > -1)
          nMinusOneChannelsB->GetXaxis()->SetBinLabel(iChan + 2, chanNames[iChan].c_str());
        else
          nMinusOneChannelsB->GetXaxis()->SetBinLabel(iChan + 2, "None");
        nMinusOneChannelsB->Fill(iChan + 2, pfLH.getFuncMin() - nomChi2);
      }

      if (iChan > -1)
        nMinusOneSystChanB->GetYaxis()->SetBinLabel(iChan + 2, chanNames[iChan].c_str());
      else
        nMinusOneSystChanB->GetYaxis()->SetBinLabel(iChan + 2, "All");
      sbd->excludeChannel(-1);
    }
    if (iSyst > -1)
      sbd->setSystFitFlag(iSyst, true);
    sbd->excludeChannel(-1);
    sbd->setBaselineModel(); //wf108 add//go back to nominal prediction
  }
  sbd->excludeChannel(-1);
  for (int iSyst = -1; iSyst < sbd->getNsyst(); ++iSyst)
    sbd->setSystFitFlag(iSyst, true);
  sbd->setBaselineModel(); //wf108 add//go back to nominal prediction

  printf("Collie Fit Test: Performing B-Only N-1 bin tests...\n");
  pfLH.setModel(sbd);
  sbd->setBaselineModel(); //wf108 add//go back to nominal prediction
  pfLH.fitSignal(false);
  pfLH.sigLLR(llrcut);
  pfLH.fitProfile();

  vector<uint> delBins;
  for (int bin = 0; bin < sbd->nbins(); ++bin) {
    delBins.clear();
    delBins.push_back(bin);
    sbd->excludeBins(delBins);
    sbd->setBaselineModel();
    pfLH.setModel(sbd);
    sbd->setBaselineModel();
    pfLH.fitProfile();
    if (pfLH.getStatus() != 3) {
      sbd->setBaselineModel();
      pfLH.fitProfile();
      if (pfLH.getStatus() != 3) {
        printf("-->FitTest, Failed fit on B-Only N-1 bin test!  Continuing...\n");
        continue;
      }
    }
    nMinusOneBinsB->SetBinContent(bin + 1, pfLH.getFuncMin() - nomChi2);
  }

  delBins.clear();
  sbd->clearExludedBins();
  //sbd->excludeBins(delBins);
  sbd->setBaselineModel(); //wf108 add//go back to nominal prediction

  printf("Collie Fit Test: Performing S+B N-1 syst/channel tests...\n");
  //Perform S+B N-1 tests
  pfLH.setModel(sbd);
  sbd->setBaselineModel(); //wf108 add//go back to nominal prediction
  pfLH.fitSignal(true);
  pfLH.fitProfile();
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3) {
    printf("-->FitTest, Failed fit on S+B nominal chi2 test!  Failing...\n");
    return;
  }
  nomChi2 = pfLH.getFuncMin();
  sbd->setBaselineModel(); //wf108 add//go back to nominal prediction

  //Test deltaChi2 dropping one event source fit at a time...
  for (int source = -1; source < sbd->getNbkgds(); ++source) {
    if (source > -1)
      requireFlagUpdate(sbd->setBkgdFitFlag(source, false));

    pfLH.setModel(sbd);
    sbd->setBaselineModel();
    pfLH.fitProfile();

    sbd->setBaselineModel();
    for (int b = 0; b < sbd->nbins(); ++b)
      sbd->data(b) = sbd->bkgd(b);
    pfLH.fitProfile();
    pfLH.fitProfile();
    //S+B fit to B-only data
    nMinusOneSourceLLR->Fill(source + 2, pfLH.getFuncMin());

    for (int iSyst = -1; iSyst < sbd->getNsyst(); ++iSyst) {
      if (iSyst > -1)
        requireFlagUpdate(sbd->setSystFitFlag(iSyst, false));

      pfLH.setModel(sbd);
      sbd->setBaselineModel();
      pfLH.fitProfile();

      sbd->setBaselineModel();
      for (int b = 0; b < sbd->nbins(); ++b)
        sbd->data(b) = sbd->bkgd(b);
      pfLH.fitProfile();
      pfLH.fitProfile();
      //S+B fit to B-only data
      nMinusTwoSourceLLR->Fill(source + 2, iSyst + 2, pfLH.getFuncMin());

      if (iSyst > -1)
        requireFlagUpdate(sbd->setSystFitFlag(iSyst, true));
    }

    if (source > -1)
      requireFlagUpdate(sbd->setBkgdFitFlag(source, true));
  }

  offset = 2 + sbd->getNbkgds();
  for (int source = 0; source < sbd->getNsignals(); ++source) {
    requireFlagUpdate(sbd->setSigFitFlag(source, false));

    pfLH.setModel(sbd);
    sbd->setBaselineModel();
    pfLH.fitProfile();

    sbd->setBaselineModel();
    for (int b = 0; b < sbd->nbins(); ++b)
      sbd->data(b) = sbd->bkgd(b);
    pfLH.fitProfile();
    pfLH.fitProfile();
    //S+B fit to B-only data
    nMinusOneSourceLLR->Fill(source + offset, pfLH.getFuncMin());

    for (int iSyst = -1; iSyst < sbd->getNsyst(); ++iSyst) {
      if (iSyst > -1)
        requireFlagUpdate(sbd->setSystFitFlag(iSyst, false));

      pfLH.setModel(sbd);
      sbd->setBaselineModel();
      pfLH.fitProfile();

      sbd->setBaselineModel();
      for (int b = 0; b < sbd->nbins(); ++b)
        sbd->data(b) = sbd->bkgd(b);
      pfLH.fitProfile();
      pfLH.fitProfile();
      //S+B fit to B-only data
      nMinusTwoSourceLLR->Fill(source + offset, iSyst + 2, pfLH.getFuncMin());

      if (iSyst > -1)
        requireFlagUpdate(sbd->setSystFitFlag(iSyst, true));
    }

    requireFlagUpdate(sbd->setSigFitFlag(source, true));
  }

  for (int syst = -1; syst < sbd->getNsyst(); ++syst) {
    for (int chan = -1; chan < nChan; ++chan) {
      if (nChan == 1 && chan >= 0)
        continue;
      sbd->excludeChannel(chan);

      if (syst > -1)
        sbd->setSystFitFlag(syst, false);

      pfLH.setModel(sbd);
      sbd->setBaselineModel();
      pfLH.fitProfile();
      pfLH.fitProfile();

      if (pfLH.getStatus() != 3) {
        sbd->setBaselineModel();
        pfLH.fitProfile();
        if (pfLH.getStatus() != 3) {
          printf("-->FitTest, Failed fit on S+B N-1 channel/syst test!  Continuing...\n");
          if (syst >= 0)
            printf("---->Failed while excluding syst=%s\n", sbd->getSystName(syst).c_str());
          if (chan >= 0)
            printf("---->Failed while excluding channel=%s\n", sbd->getChannelName(chan).c_str());
          continue;
        }
      }

      nMinusOneSystChanSB->Fill(syst + 2, chan + 2, pfLH.getFuncMin() - nomChi2);
      nMinusOneIndivSB[chan + 1]->Fill(syst + 2, pfLH.getFuncMin() - nomChi2);

      if (chan == -1) {
        sbd->setBaselineModel();
        for (int b = 0; b < sbd->nbins(); ++b)
          sbd->data(b) = sbd->bkgd(b) + sbd->signal(b);
        nMinusOneSystLLR->SetBinContent(nMinusOneSystLLR->GetNbinsX(), -1 * sbd->calculateLLR());

        sbd->setBaselineModel();
        for (int b = 0; b < sbd->nbins(); ++b)
          sbd->data(b) = sbd->bkgd(b);
        nMinusOneSystLLR->AddBinContent(nMinusOneSystLLR->GetNbinsX(), sbd->calculateLLR());
        pfLH.fitProfile();
        pfLH.fitProfile();
        //S+B fit to B-only data
        nMinusOneSystLLR->Fill(syst + 2, pfLH.getFuncMin());

        //Look at 2D correlations and their impact...
        for (int systB = -1; systB < sbd->getNsyst(); ++systB) {
          if (systB < syst)
            continue;
          if (systB > -1 && syst != systB)
            requireFlagUpdate(sbd->setSystFitFlag(systB, false));

          pfLH.setModel(sbd);
          sbd->setBaselineModel();
          pfLH.fitProfile();

          sbd->setBaselineModel();
          for (int b = 0; b < sbd->nbins(); ++b)
            sbd->data(b) = sbd->bkgd(b) + sbd->signal(b);
          double nomL = -1 * sbd->calculateLLR();

          sbd->setBaselineModel();
          for (int b = 0; b < sbd->nbins(); ++b)
            sbd->data(b) = sbd->bkgd(b);
          nomL += sbd->calculateLLR();
          nMinusTwoSystLLR->SetBinContent(nMinusOneSystLLR->GetNbinsX(),
                                          nMinusOneSystLLR->GetNbinsX(), nomL);
          pfLH.fitProfile();
          pfLH.fitProfile();
          //S+B fit to B-only data
          nMinusTwoSystLLR->Fill(syst + 2, systB + 2, pfLH.getFuncMin());

          if (systB > -1 && syst != systB)
            requireFlagUpdate(sbd->setSystFitFlag(systB, true));
        }
      }

      if (syst > -1) {
        nMinusOneSystChanSB->GetXaxis()->SetBinLabel(syst + 2, sbd->getSystName(syst).c_str());
        nMinusOneIndivSB[chan + 1]->GetXaxis()->SetBinLabel(syst + 2,
                                                            sbd->getSystName(syst).c_str());
        nMinusOneSystLLR->GetXaxis()->SetBinLabel(syst + 2, sbd->getSystName(syst).c_str());
        nMinusTwoSystLLR->GetXaxis()->SetBinLabel(syst + 2, sbd->getSystName(syst).c_str());
        nMinusTwoSystLLR->GetYaxis()->SetBinLabel(syst + 2, sbd->getSystName(syst).c_str());
      } else {
        nMinusOneSystLLR->GetXaxis()->SetBinLabel(syst + 2, "None");
        nMinusTwoSystLLR->GetXaxis()->SetBinLabel(syst + 2, "None");
        nMinusTwoSystLLR->GetYaxis()->SetBinLabel(syst + 2, "None");
        nMinusOneSystChanSB->GetXaxis()->SetBinLabel(syst + 2, "None");
        nMinusOneIndivSB[chan + 1]->GetXaxis()->SetBinLabel(syst + 2, "None");
        if (chan > -1)
          nMinusOneChannelsSB->GetXaxis()->SetBinLabel(chan + 2, chanNames[chan].c_str());
        else
          nMinusOneChannelsSB->GetXaxis()->SetBinLabel(chan + 2, "None");
        nMinusOneChannelsSB->Fill(chan + 2, pfLH.getFuncMin() - nomChi2);
      }

      if (chan > -1)
        nMinusOneSystChanSB->GetYaxis()->SetBinLabel(chan + 2, chanNames[chan].c_str());
      else
        nMinusOneSystChanSB->GetYaxis()->SetBinLabel(chan + 2, "All");
      sbd->excludeChannel(-1);
    }
    if (syst > -1)
      sbd->setSystFitFlag(syst, true);
    sbd->excludeChannel(-1);
    sbd->setBaselineModel(); //wf108 add//go back to nominal prediction
  }
  sbd->excludeChannel(-1);
  for (int syst = -1; syst < sbd->getNsyst(); ++syst)
    sbd->setSystFitFlag(syst, true);
  sbd->setBaselineModel(); //wf108 add//go back to nominal prediction

  nMinusOneSystLLR->GetYaxis()->SetTitle("% Reduction in Limit");
  nMinusOneSourceLLR->GetYaxis()->SetTitle("% Reduction in Limit");
  nMinusTwoSourceLLR->GetZaxis()->SetTitle("% Reduction in Limit");

  double nomLLR = nMinusOneSystLLR->GetBinContent(1) + 1e-9;
  for (int i = 1; i <= nMinusOneSystLLR->GetNbinsX(); ++i) {
    //   printf("Bin %d: dLLR: %f\n",i,nMinusOneSystLLR->GetBinContent(i));
    nMinusOneSystLLR->SetBinContent(i,
                                    (sqrt(nMinusOneSystLLR->GetBinContent(i) / nomLLR) - 1) * 100);
    for (int j = 1; j <= nMinusTwoSystLLR->GetNbinsY(); ++j) {
      if (nMinusTwoSystLLR->GetBinContent(i, j) == 0)
        continue;
      nMinusTwoSystLLR->SetBinContent(
          i, j, (sqrt(nMinusTwoSystLLR->GetBinContent(i, j) / nomLLR) - 1) * 100);
    }
  }

  int i = nMinusOneSystLLR->GetNbinsX();
  nMinusOneSystLLR->GetXaxis()->SetBinLabel(i, "All");
  nMinusOneSystLLR->SetBinContent(i, (sqrt(nMinusOneSystLLR->GetBinContent(i) / nomLLR) - 1) * 100);
  nMinusTwoSystLLR->SetBinContent(i, i,
                                  (sqrt(nMinusTwoSystLLR->GetBinContent(i, i) / nomLLR) - 1) * 100);

  for (int i = 1; i <= nMinusOneSourceLLR->GetNbinsX(); ++i) {
    nMinusOneSourceLLR->SetBinContent(i, (sqrt(nMinusOneSourceLLR->GetBinContent(i) / nomLLR) - 1) *
                                             100);
  }

  for (int bx = 1; bx <= nMinusTwoSourceLLR->GetNbinsX(); ++bx) {
    for (int by = 1; by <= nMinusTwoSourceLLR->GetNbinsY(); ++by) {
      nMinusTwoSourceLLR->SetBinContent(
          bx, by, (sqrt(nMinusTwoSourceLLR->GetBinContent(bx, by) / nomLLR) - 1) * 100);
    }
  }

  printf("Collie Fit Test: Performing S+B N-1 bin tests...\n");
  pfLH.setModel(sbd);
  sbd->setBaselineModel(); //wf108 add//go back to nominal prediction
  pfLH.fitSignal(true);
  pfLH.fitProfile();
  pfLH.fitProfile();

  for (int bin = 0; bin < sbd->nbins(); ++bin) {
    delBins.clear();
    delBins.push_back(bin);
    sbd->excludeBins(delBins);
    pfLH.setModel(sbd);
    sbd->setBaselineModel();
    pfLH.fitProfile();

    if (pfLH.getStatus() != 3) {
      sbd->setBaselineModel();
      pfLH.fitProfile();
      if (pfLH.getStatus() != 3) {
        printf("-->FitTest, Failed fit on S+B N-1 bin test!  Continuing...\n");
        continue;
      }
    }
    nMinusOneBinsSB->SetBinContent(bin + 1, pfLH.getFuncMin() - nomChi2);
  }

  sbd->clearExludedBins();
  sbd->setBaselineModel(); //wf108 add//go back to nominal prediction

  return;
}

void FitTest::runPETest(ProfileLH& pfLH, SigBkgdDist* sbd, double llrcut) {
  printf("Collie Fit Test: Performing pseudo-experiment tests...\n");

  int bins = sbd->nbins();
  double min = sbd->min();
  double max = sbd->max();
  double pseudoData[bins];
  SigBkgdDist* refsbd = new SigBkgdDist(*sbd);

  TH1D* sigPull_D[m_maxSyst];
  TH1D* sigPull_PE[m_maxSyst];
  TH1D* chi2_sb[m_maxSystB];
  TH1D* chi2_b[m_maxSystSB];
  TH1D* p1S_bkgd[m_maxSyst];
  TH1D* m1S_bkgd[m_maxSyst];
  TH1D* p1S_sig[m_maxSyst];
  TH1D* m1S_sig[m_maxSyst];

  TH1D* bkgPull_D[m_maxSyst];
  TH1D* bkgPull_PE[m_maxSyst];

  char dest[256];
  char dest2[256];
  for (int i = 0; i < m_maxSyst; ++i) {
    sprintf(dest, "S+B Fit to Data: pull %d", i);
    sprintf(dest2, "S+B Fit to Data: pull %s", sbd->getSystName(i).c_str());
    sigPull_D[i] = new TH1D(dest, dest2, 100, -4.5, 4.5);

    sprintf(dest, "S+B Fit to PEs: pull %d", i);
    sprintf(dest2, "S+B Fit to PEs: pull %s", sbd->getSystName(i).c_str());
    sigPull_PE[i] = new TH1D(dest, dest2, 100, -4.5, 4.5);

    sprintf(dest, "S+B Fit chi2: %d", i);
    sprintf(dest2, "S+B Fit chi2: %s", sbd->getSystName(i).c_str());
    chi2_sb[i] = new TH1D(dest, dest2, 600, -6, 6);

    sprintf(dest, "B-Only Fit chi2: %d", i);
    sprintf(dest2, "B-Only Fit chi2: %s", sbd->getSystName(i).c_str());
    chi2_b[i] = new TH1D(dest, dest2, 600, -6, 6);

    sprintf(dest, "Plus 1 Sigma Bkgd: %d", i);
    sprintf(dest2, "Plus 1 Sigma Bkgd: %s", sbd->getSystName(i).c_str());
    p1S_bkgd[i] = new TH1D(dest, dest2, bins, min, max);

    sprintf(dest, "Minus 1 Sigma Bkgd: %d", i);
    sprintf(dest2, "Minus 1 Sigma Bkgd: %s", sbd->getSystName(i).c_str());
    m1S_bkgd[i] = new TH1D(dest, dest2, bins, min, max);

    sprintf(dest, "Plus 1 Sigma Signal: %d", i);
    sprintf(dest2, "Plus 1 Sigma Signal: %s", sbd->getSystName(i).c_str());
    p1S_sig[i] = new TH1D(dest, dest2, bins, min, max);

    sprintf(dest, "Minus 1 Sigma Signal: %d", i);
    sprintf(dest2, "Minus 1 Sigma Signal: %s", sbd->getSystName(i).c_str());
    m1S_sig[i] = new TH1D(dest, dest2, bins, min, max);

    sprintf(dest, "B-Only Fit to Data: pull %d", i);
    sprintf(dest2, "B-Only Fit to Data: pull %s", sbd->getSystName(i).c_str());
    bkgPull_D[i] = new TH1D(dest, dest2, 100, -5, 5);

    sprintf(dest, "B-Only Fit to PEs: pull %d", i);
    sprintf(dest2, "B-Only Fit to PEs: pull %s", sbd->getSystName(i).c_str());
    bkgPull_PE[i] = new TH1D(dest, dest2, 100, -4.5, 4.5);
  }

  float maxChi = m_chiSB;
  if (m_chiB > m_chiSB)
    maxChi = m_chiB;
  if (maxChi < 30)
    maxChi = 30;
  maxChi *= 2.5;
  TH1D* dataChi2SB = new TH1D("Data Chi-Square S+B", "Data Chi-Square S+B", 500, 0, maxChi);
  TH1D* dataChi2B = new TH1D("Data Chi-Square B-Only", "Data Chi-Square B-Only", 500, 0, maxChi);

  TH1D* dataChi2SB_Ref =
      new TH1D("Data Chi-Square S+B, Ref", "Data Chi-Square S+B, Ref", 500, 0, maxChi);
  TH1D* dataChi2B_Ref =
      new TH1D("Data Chi-Square B-Only, Ref", "Data Chi-Square B-Only, Ref", 500, 0, maxChi);

  TH1D* peChi2SB = new TH1D("PE Chi-Square S+B", "PE Chi-Square S+B", 500, 0, maxChi);
  TH1D* peChi2B = new TH1D("PE Chi-Square B-Only", "PE Chi-Square B-Only", 500, 0, maxChi);

  TH1D* sigSyst = new TH1D("Signal Systematic", "Signal Systematic", 100, -0.5, 2);
  TH1D* bkgSyst = new TH1D("Bkgd Systematic", "Bkgd Systematic", 100, -0.5, 2);

  TH1D* fitsigSB =
      new TH1D("Signal Systematic, S+B Fit", "Signal Systematic, S+B Fit", 100, -0.5, 2);
  TH1D* fitbkgSB = new TH1D("Bkgd Systematic, S+B Fit", "Bkgd Systematic, S+B Fit", 100, -0.5, 2);

  TH1D* fitsigB =
      new TH1D("Signal Systematic, B-Only Fit", "Signal Systematic, B-Only Fit", 100, -0.5, 2);
  TH1D* fitbkgB =
      new TH1D("Bkgd Systematic, B-Only Fit", "Bkgd Systematic, B-Only Fit", 100, -0.5, 2);

  TH1D* sigSystTOT = new TH1D("Signal Systematic TOT", "Signal Systematic TOT", 100, -0.5, 2);
  TH1D* bkgSystTOT = new TH1D("Bkgd Systematic TOT", "Bkgd Systematic TOT", 100, -0.5, 2);

  TH1D* fitsigSBTOT =
      new TH1D("Signal Systematic TOT, S+B Fit", "Signal Systematic TOT, S+B Fit", 100, -0.5, 2);
  TH1D* fitbkgSBTOT =
      new TH1D("Bkgd Systematic TOT, S+B Fit", "Bkgd Systematic TOT, S+B Fit", 100, -0.5, 2);

  TH1D* fitsigBTOT = new TH1D("Signal Systematic TOT, B-Only Fit",
                              "Signal Systematic TOT, B-Only Fit", 100, -0.5, 2);
  TH1D* fitbkgBTOT =
      new TH1D("Bkgd Systematic TOT, B-Only Fit", "Bkgd Systematic TOT, B-Only Fit", 100, -0.5, 2);

  TH1D* fitsigSB_PE =
      new TH1D("Signal Systematic, S+B PE Fit", "Signal Systematic, S+B PE Fit", 100, -0.5, 2);
  TH1D* fitbkgSB_PE =
      new TH1D("Bkgd Systematic, S+B PE Fit", "Bkgd Systematic, S+B PE Fit", 100, -0.5, 2);

  TH1D* fitsigB_PE = new TH1D("Signal Systematic, B-Only PE Fit",
                              "Signal Systematic, B-Only PE Fit", 100, -0.5, 2);
  TH1D* fitbkgB_PE =
      new TH1D("Bkgd Systematic, B-Only PE Fit", "Bkgd Systematic, B-Only PE Fit", 100, -0.5, 2);

  TH1D* minuitIterSB = new TH1D("S+B Fit Iterations", "S+B Fit Iterations", 2500, 0, 5000);
  TH1D* minuitIterB = new TH1D("B-Only Fit Iterations", "S+B Fit Iterations", 2500, 0, 5000);

  TH1D* minuitStatusSB = new TH1D("S+B Fit Status", "S+B Fit Status", 7, -2.5, 4.5);
  TH1D* minuitStatusB = new TH1D("B-Only Fit Status", "S+B Fit Status", 7, -2.5, 4.5);

  TH1D* fitErrsSB = new TH1D("S+B Fit Errors", "S+B Fit Errors", m_maxSyst, 0.5, m_maxSyst + 0.5);
  TH1D* fitErrsB =
      new TH1D("B-Only Fit Errors", "B-Only Fit Errors", m_maxSyst, 0.5, m_maxSyst + 0.5);

  TH1D* totSBsyst[bins];
  TH1D* totBsyst[bins];
  TH1D* totSyst[bins];

  char title[256];
  for (int b = 0; b < bins; ++b) {
    sprintf(title, "non-fit syst hist bin %d", b + 1);
    totSyst[b] = new TH1D(title, title, 200, 0, (sbd->bkgd(b) + sbd->signal(b)) * 3.0);

    sprintf(title, "sig+bkgd syst hist bin %d", b + 1);
    totSBsyst[b] = new TH1D(title, title, 200, 0, (sbd->bkgd(b) + sbd->signal(b)) * 3.0);

    sprintf(title, "bkgd syst hist bin %d", b + 1);
    totBsyst[b] = new TH1D(title, title, 200, 0, (sbd->bkgd(b) + sbd->signal(b)) * 3.0);
  }

  collie::PoissonRandom poisson;

  ///Perform S+B fit iterations
  printf("Collie Fit Test: Performing S+B pseudo-experiment tests...\n");
  sbd->setBaselineModel();
  pfLH.setModel(sbd);
  pfLH.fitProfile();
  sbd->setBaselineModel();
  pfLH.fitSignal(true);
  pfLH.sigLLR(llrcut);
  pfLH.fitProfile();
  pfLH.fitProfile();

  //Record Reference S+B chi2
  dataChi2SB_Ref->Fill(pfLH.getFuncMin());

  int totIter = 0;
  for (int iter = 0; iter < m_iterations; iter++) {
    if ((totIter % 1000) == 0)
      printf("   %d iterations\n", totIter);

    ++totIter;

    sbd->setBaselineModel(); //go back to nominal prediction
    sbd->fluctuate();
    //non-fitted systematics fluctuations
    for (int j = 0; j < bins; ++j) {
      if (refsbd->signal(j) > 0)
        sigSyst->Fill((sbd->signal(j)) / refsbd->signal(j));
      if (refsbd->bkgd(j) > 0)
        bkgSyst->Fill((sbd->bkgd(j)) / refsbd->bkgd(j));
      totSyst[j]->Fill(sbd->bkgd(j) + sbd->signal(j));
    }
    sigSystTOT->Fill(sbd->totSignal() / refsbd->totSignal());
    bkgSystTOT->Fill(sbd->totBkgd() / refsbd->totBkgd());

    ///Collect S+B fluctuations  wf108 add
    if (m_fillPE) {
      for (int j = 0; j < bins; ++j)
        pseudoData[j] = sbd->bkgd(j) + sbd->signal(j);
    }

    //Record minuit statistics
    pfLH.fitProfile();
    minuitIterSB->Fill(pfLH.getNiterations());
    minuitStatusSB->Fill(pfLH.getStatus());
    //S+B-fit systematics fluctuations in best fit to data
    for (int j = 0; j < bins; ++j) {
      if (refsbd->signal(j) > 0)
        fitsigSB->Fill(sbd->signal(j) / refsbd->signal(j));
      if (refsbd->bkgd(j) > 0)
        fitbkgSB->Fill(sbd->bkgd(j) / refsbd->bkgd(j));
      totSBsyst[j]->Fill(sbd->bkgd(j) + sbd->signal(j));
    }
    fitsigSBTOT->Fill(sbd->totSignal() / refsbd->totSignal());
    fitbkgSBTOT->Fill(sbd->totBkgd() / refsbd->totBkgd());

    //Record fitted nuisance parameters
    for (int s = 0; s < m_maxSyst; ++s)
      sigPull_D[s]->Fill(pfLH.getParamVal(s));

    //Record S+B chi2
    dataChi2SB->Fill(pfLH.getFuncMin());

    if (m_fillPE) {
      sbd->setBaselineModel(); //go back to nominal prediction
      //Simulate poisson trials from best fit to data...
      for (int j = 0; j < bins; ++j)
        sbd->data(j) = poisson.fire(pseudoData[j]);
      pfLH.fitProfile();

      //Record best fit parameters
      for (int s = 0; s < m_maxSyst; ++s) {
        sigPull_PE[s]->Fill(pfLH.getParamVal(s));
      }

      //Record systematics fluctuations
      for (int j = 0; j < bins; ++j) {
        if (refsbd->signal(j) > 0)
          fitsigSB_PE->Fill((sbd->signal(j)) / refsbd->signal(j));
        if (refsbd->bkgd(j) > 0)
          fitbkgSB_PE->Fill((sbd->bkgd(j)) / refsbd->bkgd(j));
      }
      peChi2SB->Fill(pfLH.getFuncMin());
    }
  }

  //Now perform B-Only fit iterations
  printf("Collie Fit Test: Performing B-Only pseudo-experiment tests...\n");
  sbd->setBaselineModel();
  pfLH.fitSignal(false);
  pfLH.setModel(sbd);
  pfLH.fitSignal(false);
  pfLH.fitProfile();
  sbd->setBaselineModel(); //wf108 add//go back to nominal prediction
  pfLH.fitProfile();
  pfLH.fitProfile();

  //Record Reference B-Only chi2
  dataChi2B_Ref->Fill(pfLH.getFuncMin());

  totIter = 0;
  for (int iter = 0; iter < m_iterations; iter++) {
    if ((totIter % 1000) == 0)
      printf("   %d iterations\n", totIter);

    ++totIter;

    sbd->setBaselineModel(); //go back to nominal prediction
    sbd->fluctuate();

    //Collect B-only fluctuations wf108 add
    if (m_fillPE) {
      for (int j = 0; j < bins; ++j)
        pseudoData[j] = sbd->bkgd(j);
    }

    //Record minuit statistics
    pfLH.fitProfile();
    minuitIterB->Fill(pfLH.getNiterations());
    minuitStatusB->Fill(pfLH.getStatus());

    for (int j = 0; j < bins; ++j) {
      if (refsbd->signal(j) > 0)
        fitsigB->Fill(sbd->signal(j) / refsbd->signal(j));
      if (refsbd->bkgd(j) > 0)
        fitbkgB->Fill(sbd->bkgd(j) / refsbd->bkgd(j));
      totBsyst[j]->Fill(sbd->bkgd(j));
    }

    fitsigBTOT->Fill(sbd->totSignal() / refsbd->totSignal());
    fitbkgBTOT->Fill(sbd->totBkgd() / refsbd->totBkgd());
    for (int s = 0; s < m_maxSyst; ++s)
      bkgPull_D[s]->Fill(pfLH.getParamVal(s));

    dataChi2B->Fill(pfLH.getFuncMin());

    if (m_fillPE) {
      sbd->setBaselineModel(); //go back to nominal prediction

      for (int j = 0; j < bins; ++j)
        sbd->data(j) = poisson.fire(pseudoData[j]);
      pfLH.fitProfile();

      for (int s = 0; s < m_maxSyst; ++s) {
        bkgPull_PE[s]->Fill(pfLH.getParamVal(s));
      }

      for (int j = 0; j < bins; ++j) {
        if (refsbd->signal(j) > 0)
          fitsigB_PE->Fill((sbd->signal(j)) / refsbd->signal(j));
        if (refsbd->bkgd(j) > 0)
          fitbkgB_PE->Fill((sbd->bkgd(j)) / refsbd->bkgd(j));
      }

      peChi2B->Fill(pfLH.getFuncMin());
    }
  }

  TH1D* systFinal = new TH1D("Total Non-Fit Systematics per Bin",
                             "Total Non-Fit Systematics per Bin", bins, min, max);

  TH1D* systFinalSB =
      new TH1D("Total S+B Systematics per Bin", "Total S+B Systematics per Bin", bins, min, max);

  TH1D* systFinalB = new TH1D("Total B-Only Systematics per Bin",
                              "Total B-Only Systematics per Bin", bins, min, max);

  for (int b = 0; b < bins; ++b) {
    systFinalSB->SetBinContent(b + 1, totSBsyst[b]->GetRMS());
    systFinalB->SetBinContent(b + 1, totBsyst[b]->GetRMS());
    systFinal->SetBinContent(b + 1, totSyst[b]->GetRMS());
  }

  for (int i = 0; i < m_maxSyst; ++i) {
    fitErrsB->SetBinContent(i + 1, bkgPull_D[i]->GetRMS());
    fitErrsSB->SetBinContent(i + 1, sigPull_D[i]->GetRMS());
    fitErrsB->GetXaxis()->SetBinLabel(i + 1, sbd->getSystName(i).c_str());
    fitErrsSB->GetXaxis()->SetBinLabel(i + 1, sbd->getSystName(i).c_str());
  }

  //M Owen - don't use histo stats for these tests
  const bool usehistostats = sbd->usingHistoStats();
  sbd->useHistoStats(false);
  sbd->setBaselineModel();

  //Perform B-Only fit
  sbd->setBaselineModel(); //go back to nominal prediction
  pfLH.fitSignal(false);
  pfLH.sigLLR(llrcut);
  pfLH.fitProfile();
  pfLH.fitProfile();

  double params[m_maxSyst];
  // USED TO BE: double paramsT[m_maxSyst];
  std::vector<double> paramsT(m_maxSyst);

  for (int s = 0; s < m_maxSyst; ++s)
    params[s] = pfLH.getParamVal(s);

  sbd->setBaselineModel();
  for (int s = 0; s < m_maxSyst; ++s) {
    for (int z = 0; z < m_maxSyst; ++z)
      paramsT[z] = params[z];
    for (int i = 0; i < 600; ++i) {
      double fluct = -6.0 + i * 12.0 / 600.0;
      paramsT[s] = fluct;
      //      sbd->fluctuate(&paramsT[0]);
      sbd->fluctuate(paramsT);
      chi2_b[s]->SetBinContent(i + 1, sbd->calculateChi2LLR(false, 1e6)); //float this parameter...
    }
  }

  //Now perform S+B fit
  sbd->setBaselineModel(); //go back to nominal prediction
  pfLH.fitSignal(true);
  pfLH.fitProfile();
  pfLH.fitProfile();

  for (int s = 0; s < m_maxSyst; ++s)
    params[s] = pfLH.getParamVal(s);

  sbd->setBaselineModel();
  for (int s = 0; s < m_maxSyst; ++s) {
    for (int z = 0; z < m_maxSyst; ++z)
      paramsT[z] = params[z];
    for (int i = 0; i < 600; ++i) {
      double fluct = -6.0 + i * 12.0 / 600.0;
      paramsT[s] = fluct;
      sbd->fluctuate(paramsT);
      chi2_sb[s]->SetBinContent(i + 1, sbd->calculateChi2LLR(true, 1e6)); //float this parameter...
    }
  }

  for (int s = 0; s < m_maxSyst; ++s) {
    for (int z = 0; z < m_maxSyst; ++z)
      paramsT[z] = 0;
    sbd->setBaselineModel();
    paramsT[s] = 1;
    sbd->fluctuate(paramsT);
    for (int b = 1; b <= bins; ++b) {
      p1S_bkgd[s]->SetBinContent(b, sbd->bkgd(b - 1));
      p1S_sig[s]->SetBinContent(b, sbd->signal(b - 1));
    }

    sbd->setBaselineModel();
    paramsT[s] = -1;
    sbd->fluctuate(paramsT);
    for (int b = 1; b <= bins; ++b) {
      m1S_bkgd[s]->SetBinContent(b, sbd->bkgd(b - 1));
      m1S_sig[s]->SetBinContent(b, sbd->signal(b - 1));
    }
  }

  delete refsbd;
  refsbd = NULL;

  sbd->useHistoStats(usehistostats);
  sbd->setBaselineModel(); //wf108 add//go back to nominal prediction

  return;
}

void FitTest::runTest(SigBkgdDist* sbd, double llrcut, bool doPE, bool doNM) {

  sbd->setSignalScale(1.0);
  sbd->setBackgroundScale(1.0);
  sbd->setDataScale(1.0);
  sbd->setBaselineModel();

  ProfileLH pfLH;
  pfLH.setFitTest(true);
  pfLH.fitSignal(true);
  pfLH.sigLLR(llrcut);
  pfLH.setModel(sbd);
  sbd->setBaselineModel();
  pfLH.fitProfile();
  sbd->setBaselineModel();
  pfLH.fitProfile();

  //Create some histograms and containers
  m_maxSystSB = pfLH.getNfitSyst();
  int bins = sbd->nbins();
  double min = sbd->min();
  double max = sbd->max();

  pfLH.fitSignal(false);
  pfLH.sigLLR(llrcut);
  sbd->setBaselineModel();
  pfLH.fitProfile();
  sbd->setBaselineModel();

  m_maxSystB = pfLH.getNfitSyst();
  m_maxSyst = m_maxSystSB;
  if (m_maxSystB > m_maxSyst)
    m_maxSyst = m_maxSystB;

  double errMatrixSB[m_maxSystSB][m_maxSystSB];
  double errMatrixB[m_maxSystB][m_maxSystB];

  TH1D* dataDist = new TH1D("Total Data", "Total Data", bins, min, max);

  TH1D* bkgdDist = new TH1D("Total Background, Prefit", "Total Background, prefit", bins, min, max);
  TH1D* bkgdDistSB =
      new TH1D("Total Background, S+B Fit", "Total Background, S+B Fit", bins, min, max);
  TH1D* bkgdDistB =
      new TH1D("Total Background, B-Only Fit", "Total Background, B-Only Fit", bins, min, max);

  TH1D* sigDist = new TH1D("Total Signal, Prefit", "Total Signal, prefit", bins, min, max);
  TH1D* sigDistSB = new TH1D("Total Signal, S+B Fit", "Total Signal, S+B Fit", bins, min, max);
  TH1D* sigDistB = new TH1D("Total Signal, B-Only Fit", "Total Signal, B-Only Fit", bins, min, max);

  TH1D* dataCDF = new TH1D("Data CDF", "Data CDF", bins, min, max);
  TH1D* bkgdCDF = new TH1D("Background CDF, Prefit", "Background CDF, prefit", bins, min, max);
  TH1D* bkgdCDFSB = new TH1D("Background CDF, S+B Fit", "Background CDF, S+B Fit", bins, min, max);
  TH1D* bkgdCDFB =
      new TH1D("Background CDF, B-Only Fit", "Background CDF, B-Only Fit", bins, min, max);

  TH1D* sigCDF = new TH1D("Signal CDF, Prefit", "Signal CDF, prefit", bins, min, max);
  TH1D* sigCDFSB = new TH1D("Signal CDF, S+B Fit", "Signal CDF, S+B Fit", bins, min, max);
  TH1D* sigCDFB = new TH1D("Signal CDF, B-Only Fit", "Signal CDF, B-Only Fit", bins, min, max);

  TH1D* fitParamsSB = new TH1D("S+B Fit Params", "S+B Fit Params", m_maxSyst, 0.5, m_maxSyst + 0.5);
  TH1D* fitParamsB =
      new TH1D("B-Only Fit Params", "B-Only Fit Params", m_maxSyst, 0.5, m_maxSyst + 0.5);

  TH2D* minuitErrorMatrixSB = new TH2D("S+B Fit Error Matrix", "S+B Fit Error Matrix", m_maxSystSB,
                                       0, m_maxSystSB, m_maxSystSB, 0, m_maxSystSB);
  TH2D* minuitErrorMatrixB = new TH2D("B-Only Fit Error Matrix", "B-Only Fit Error Matrix",
                                      m_maxSystB, 0, m_maxSystB, m_maxSystB, 0, m_maxSystB);
  TH1D* fitErrsM =
      new TH1D("MINUIT Fit Errors", "MINUIT Fit Errors", m_maxSyst, 0.5, m_maxSyst + 0.5);

  ///General testing
  ///  Get "best fit" integrals
  printf("\n*********************   Collie Fit Test   *********************\n");
  printf("\nTesting %d fit parameters in %d bins with Log(1+s/b)<%f\n", m_maxSyst, bins, llrcut);

  ///Now run the tests...
  printf("Before fit========>  Sig: %.4f, Bkgd: %.4f, Data: %.1f\n", sbd->totSignal(),
         sbd->totBkgd(), sbd->totData());
  for (int i = 1; i <= bins; ++i) {
    bkgdDist->SetBinContent(i, sbd->bkgd(i - 1));
    sigDist->SetBinContent(i, sbd->signal(i - 1));
    dataDist->SetBinContent(i, sbd->data(i - 1));
  }

  for (int i = 1; i <= bins; ++i) {
    bkgdCDF->SetBinContent(bins - (i - 1), bkgdDist->Integral(i, bins));
    sigCDF->SetBinContent(bins - (i - 1), sigDist->Integral(i, bins));
    dataCDF->SetBinContent(bins - (i - 1), dataDist->Integral(i, bins));
  }

  ///Do B-Only best fit to data
  sbd->zeroModel();        //go back to nominal prediction //wf108sbd->zeroModel();
  sbd->setBaselineModel(); //go back to nominal prediction //wf108sbd->zeroModel();
  pfLH.fitSignal(false);
  pfLH.sigLLR(llrcut);
  pfLH.fitProfile();
  sbd->setBaselineModel(); //go back to nominal prediction //wf108sbd->zeroModel();
  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("Error! B-Only fit did not converge!(s=%d)\n", pfLH.getStatus());

  for (int s = 0; s < m_maxSyst; ++s) {
    fitParamsB->SetBinContent(s + 1, pfLH.getParamVal(s));
    fitParamsB->GetXaxis()->SetBinLabel(s + 1, sbd->getSystName(s).c_str());
  }

  m_chiB = pfLH.getFuncMin();
  pfLH.print();
  printf("After B-Only fit==>  Sig: %.4f, Bkgd: %.4f, Data: %.1f, Chi2: %f (stat=%d)\n",
         sbd->totSignal(), sbd->totBkgd(), sbd->totData(), m_chiB, pfLH.getStatus());
  for (int i = 1; i <= bins; ++i) {
    bkgdDistB->SetBinContent(i, sbd->bkgd(i - 1));
    sigDistB->SetBinContent(i, sbd->signal(i - 1));
  }
  for (int i = 1; i <= bins; ++i) {
    bkgdCDFB->SetBinContent(bins - (i - 1), bkgdDistB->Integral(i, bins));
    sigCDFB->SetBinContent(bins - (i - 1), sigDistB->Integral(i, bins));
  }
  pfLH.getErrorMatrix(&errMatrixB[0][0]);
  for (int i = 0; i < m_maxSystB; ++i) {
    for (int j = 0; j < m_maxSystB; ++j) {
      minuitErrorMatrixB->Fill(i, j, errMatrixB[i][j]);
      minuitErrorMatrixB->GetXaxis()->SetBinLabel(i + 1, pfLH.getFitSystName(i).c_str());
      minuitErrorMatrixB->GetYaxis()->SetBinLabel(j + 1, pfLH.getFitSystName(j).c_str());
    }
  }

  ///Do S+B best fit to data
  sbd->zeroModel();        //go back to nominal prediction //wf108sbd->zeroModel();
  sbd->setBaselineModel(); //go back to nominal prediction //wf108sbd->zeroModel();
  pfLH.fitSignal(true);
  pfLH.fitProfile();
  pfLH.fitProfile();
  sbd->setBaselineModel(); //go back to nominal prediction //wf108sbd->zeroModel();
  pfLH.fitProfile();
  //  pfLH.fitProfile();
  if (pfLH.getStatus() != 3)
    printf("Error! S+B fit did not converge!(s=%d)\n", pfLH.getStatus());

  for (int s = 0; s < m_maxSyst; ++s) {
    fitParamsSB->SetBinContent(s + 1, pfLH.getParamVal(s));
    fitErrsM->SetBinContent(s + 1, pfLH.getParamError(s));
    fitParamsSB->GetXaxis()->SetBinLabel(s + 1, sbd->getSystName(s).c_str());
    fitErrsM->GetXaxis()->SetBinLabel(s + 1, sbd->getSystName(s).c_str());
  }

  m_chiSB = pfLH.getFuncMin();
  printf("\n\n\n\n");
  pfLH.print();

  printf("After S+B fit=====>  Sig: %.4f, Bkgd: %.4f, Data: %.1f, Chi2: %f (stat=%d)\n",
         sbd->totSignal(), sbd->totBkgd(), sbd->totData(), m_chiSB, pfLH.getStatus());
  for (int i = 1; i <= bins; ++i) {
    bkgdDistSB->SetBinContent(i, sbd->bkgd(i - 1));
    sigDistSB->SetBinContent(i, sbd->signal(i - 1));
  }

  for (int i = 1; i <= bins; ++i) {
    bkgdCDFSB->SetBinContent(bins - (i - 1), bkgdDistSB->Integral(i, bins));
    sigCDFSB->SetBinContent(bins - (i - 1), sigDistSB->Integral(i, bins));
  }

  pfLH.getErrorMatrix(&errMatrixSB[0][0]);
  for (int i = 0; i < m_maxSystSB; ++i) {
    for (int j = 0; j < m_maxSystSB; ++j) {
      minuitErrorMatrixSB->Fill(i, j, errMatrixSB[i][j]);
      minuitErrorMatrixSB->GetXaxis()->SetBinLabel(i + 1, pfLH.getFitSystName(i).c_str());
      minuitErrorMatrixSB->GetYaxis()->SetBinLabel(j + 1, pfLH.getFitSystName(j).c_str());
    }
  }

  sbd->setBaselineModel(); //go back to nominal prediction

  //PE testsr
  if (doPE)
    runPETest(pfLH, sbd, llrcut);

  //N Minus one tests
  if (doNM)
    runNM1Test(pfLH, sbd, llrcut);

  pfLH.setFitTest(false);

  return;
}

void FitTest::benchMark(SigBkgdDist* sbd, int navg, bool fitSig) {

  if (sbd == NULL)
    return;
  if (navg <= 0)
    return;

  printf("\n*****************************************************************\n");
  printf("Fit Test Benchmark:\n\n");
  sbd->setBaselineModel(); //zeros nuisance param central values...

  ProfileLH pfLH;
  pfLH.setLUN(99);
  if (fitSig) {
    pfLH.sigLLR(1e10);
    pfLH.fitSignal(true);
  }
  pfLH.setModel(sbd);
  timeval a, b;
  timeval c, d;
  gettimeofday(&a, NULL);
  gettimeofday(&c, NULL);
  int niter = 6 * navg;
  for (int i = 0; i < niter; ++i) {

    if ((i % navg) == 0 && i > 0) {
      gettimeofday(&d, NULL);
      double deltat = c.tv_sec * 1000.0 + c.tv_usec / 1000.0;
      deltat = d.tv_sec * 1000.0 + d.tv_usec / 1000.0 - deltat;
      printf(" %d evt avg = %f s/fit\n", (int)navg, deltat / 1000.0 / navg);
      gettimeofday(&c, NULL);
    }

    sbd->setBaselineModel(); //zeros nuisance param central values...
    sbd->fluctuate();
    pfLH.fitProfile();
    if (pfLH.getStatus() != 3)
      printf("....Fit error status: %d\n", pfLH.getStatus());
  }
  gettimeofday(&b, NULL);
  double deltat = a.tv_sec * 1000.0 + a.tv_usec / 1000.0;
  deltat = b.tv_sec * 1000.0 + b.tv_usec / 1000.0 - deltat;
  printf(" Final avg = %f seconds per fit iteration\n", deltat / 1000 / niter);
  printf("\n*****************************************************************\n");
}
