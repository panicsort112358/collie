void plotLLRdists(const char* fname = "CLtest.root") {

  TFile* f1 = new TFile(fname);
  TH1D* LLRsb = (TH1D*)f1->Get("LLR_SB");
  TH1D* LLRb = (TH1D*)f1->Get("LLR_B");
  TH1D* LLRb_1s = (TH1D*)f1->Get("LLR_B_1sigmas");
  TH1D* LLRb_2s = (TH1D*)f1->Get("LLR_B_2sigmas");
  TH1D* LLRD = (TH1D*)f1->Get("LLR_D");

  LLRsb->SetLineWidth(4);
  LLRsb->SetLineColor(2);

  LLRb->SetLineWidth(4);
  LLRb->SetLineColor(3);
  LLRb->GetXaxis()->SetTitle("LLR");
  LLRb->GetYaxis()->SetTitle("PEs / Bin");
  LLRb->SetTitle("");

  LLRb_1s->SetLineWidth(5);
  LLRb_2s->SetLineWidth(5);
  LLRb_1s->SetLineColor(15);
  LLRb_2s->SetLineColor(12);

  LLRD->SetLineWidth(5);
  LLRD->Scale(1e6);
  LLRb_1s->Scale(1e6);
  LLRb_2s->Scale(1e6);

  TLegend* leg = new TLegend(0.65, 0.65, 0.9, 0.9);
  leg->AddEntry(LLRb, "B-Only LLR", "L");
  leg->AddEntry(LLRb_1s, "B-Only LLR #pm1#sigma", "L");
  leg->AddEntry(LLRb_2s, "B-Only LLR #pm2#sigma", "L");
  leg->AddEntry(LLRsb, "S+B LLR", "L");
  leg->AddEntry(LLRD, "Observed LLR", "L");

  TCanvas* c2 = new TCanvas("LLR Distributions");
  LLRb->Draw("hist");
  LLRsb->Draw("samehist");
  LLRb_1s->Draw("samehist");
  LLRb_2s->Draw("samehist");
  LLRD->Draw("samehist");
  leg->Draw("same");

  int integral = 0;
  float clb = 0;
  float clsb = 0;
  for (int i = 0; i < LLRb->GetNbinsX(); i++) {
    if (LLRb->GetBinLowEdge(i) > LLRD->GetMean()) {
      printf("1-LLRB: %e\n", LLRb->Integral(0, i - 1) / LLRb->Integral());
      clb = 1.0 - LLRb->Integral(0, i - 1) / LLRb->Integral();
      break;
    }
  }
  for (int i = 0; i < LLRb->GetNbinsX(); i++) {
    if (LLRsb->GetBinLowEdge(i) > LLRD->GetMean()) {
      printf("1-LLRSB: %e\n", LLRsb->Integral(0, i - 1) / LLRsb->Integral());
      clsb = 1.0 - LLRsb->Integral(0, i - 1) / LLRsb->Integral();
      break;
    }
  }
  printf("1-CLS: %.4f\n", 1.0 - (clsb / clb));
}
