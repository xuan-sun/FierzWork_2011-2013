DrawBranch()
{
  gStyle->SetOptStat(0);
//  TGaxis::SetMaxDigits(3);

  TCanvas *C = new TCanvas("myC", "myC");
  C->cd();

  TFile f0("Evts_1mill_b_0.root");
  TTree *t0 = (TTree*)f0.Get("Evts");
  TH1D *h0 = new TH1D("fierz0", "b=0", 100, 0, 1000);
  t0->Draw("KE >> fierz0");

  TFile f1("Evts_1mill_b_1.root");
  TTree *t1 = (TTree*)f1.Get("Evts");
  TH1D *h1 = new TH1D("fierz1", "b=1", 100, 0, 1000);
  t1->Draw("KE >> fierz1");

  TFile finf("Evts_1mill_b_inf.root");
  TTree *tinf = (TTree*)finf.Get("Evts");
  TH1D *hinf = new TH1D("fierzinf", "b=inf", 100, 0, 1000);
  tinf->Draw("KE >> fierzinf");

  TFile fn1("Evts_1mill_b_-1.root");
  TTree *tn1 = (TTree*)fn1.Get("Evts");
  TH1D *hn1 = new TH1D("fierzn1", "b=-1", 100, 0, 1000);
  tn1->Draw("KE >> fierzn1");

  TCanvas *C2 = new TCanvas("myC2", "myC2");
  C2->cd();

  hinf->SetFillColorAlpha(38, 0.5);
//  hinf->SetFillStyle(3006);
  hinf->SetTitle("Initial Beta Electron Spectrum with Different b");
  hinf->GetXaxis()->SetTitle("Kinetic Energy (keV)");
  hinf->GetXaxis()->CenterTitle();
  hinf->Draw();
  hn1->SetFillColorAlpha(30, 0.5);
//  hn1->SetFillStyle(3007);
  hn1->Draw("SAME");
  h1->SetFillColorAlpha(41, 0.5);
//  h1->SetFillStyle(3004);
  h1->Draw("SAME");
  h0->SetFillColorAlpha(46, 0.5);
//  h0->SetFillStyle(3005);
  h0->Draw("SAME");

  TLegend *l = new TLegend(0.7, 0.6, 0.9, 0.8);
  l->AddEntry(hinf, "b -> inf", "f");
  l->AddEntry(h0, "b = 0", "f");
  l->AddEntry(h1, "b = 1", "f");
  l->AddEntry(hn1, "b = -1", "f");
  l->Draw();


  C2->Print("output.pdf");



}
