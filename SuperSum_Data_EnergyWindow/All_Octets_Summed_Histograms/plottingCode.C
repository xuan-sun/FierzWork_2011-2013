plottingCode()
{
  TCanvas *c = new TCanvas("c", "c");
  c->cd();

  TFile fData("Octets_0-59_ssDataHist_type0_radialCut_0-49mm_endpointCorrected.root");
  TH1D *hData = (TH1D*)fData.Get("totalData");
  hData->SetLineColor(2);

  TFile fBeta("FullBlind_Feb2019_MC_A_0_b_0_Octets_0-59_ssHist_type0_posCut_0-49mm.root");
  TH1D *hBeta = (TH1D*)fBeta.Get("totalBeta");
  hBeta->SetLineColor(4);

/*
  double normData = 0;
  double normBeta = 0;
  for(int i = 17; i < 65; i++)
  {
    normData = normData + hData->GetBinContent(i);
    normBeta = normBeta + hBeta->GetBinContent(i);
  }

  hBeta->Scale(normData / normBeta);
*/

//  hData->Draw();
//  hBeta->Draw("SAME");




}
