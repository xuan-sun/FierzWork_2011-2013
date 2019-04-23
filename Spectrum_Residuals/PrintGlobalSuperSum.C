PrintGlobalSuperSum()
{
  int octetLow = 80;
  int octetHigh = 122;

  // load the sum of all the data octets, with proper error propagation (if ROOT isn't broken).
  TH1D *hTotalData = new TH1D("totalData", "totalData", 120, 0, 1200);
  hTotalData->Sumw2();

  for(int i = octetLow; i < octetHigh; i++)
  {
    if(i == 9 || i == 59 || i == 91 || i == 93 || i == 101 || i == 107 || i == 121)
    {
      continue;
    }
    TFile f(Form("../PositionCuts/radialCut_0-49/Octet_%i_ssDataHist_type0_radialCut_0-49mm.root", i));
//    TFile f(Form("../PositionCuts/radialCut_0-49/Octet_%i_ssDataHist_type0_radialCut_0-49mm_endpointCorrected.root", i));
    TH1D *hTemp = (TH1D*)f.Get("Super sum");
    hTotalData->Add(hTemp);
    f.Close();
  }

  hTotalData->Draw();

  // create a new canvas and save all the b=0 MC octets
  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->cd();

  TH1D *hTotalBeta = new TH1D("totalBeta", "totalBeta", 120, 0, 1200);
  hTotalBeta->Sumw2();

  for(int i = octetLow; i < octetHigh; i++)
  {
    if(i == 9 || i == 59 || i == 91 || i == 93 || i == 101 || i == 107 || i == 121)
    {
      continue;
    }
    TFile f(Form("../PositionCuts/radialCut_0-49/FullBlind_Feb2019_MC_A_0_b_0_Octet_%i_type0_posCut_0-0.049000m.root", i));
    TH1D *hTemp = (TH1D*)f.Get("Super sum");
    hTotalBeta->Add(hTemp);
    f.Close();
  }

  hTotalBeta->Draw();

  // create a new canvas and save all the b=inf MC octets
  TCanvas *c3 = new TCanvas("c3", "c3");
  c3->cd();

  TH1D *hTotalFierz = new TH1D("totalFierz", "totalFierz", 120, 0, 1200);
  hTotalFierz->Sumw2();

  for(int i = octetLow; i < octetHigh; i++)
  {
    if(i == 9 || i == 59 || i == 91 || i == 93 || i == 101 || i == 107 || i == 121)
    {
      continue;
    }
    TFile f(Form("../PositionCuts/radialCut_0-49/MC_A_0_b_inf_Octet_%i_ssHist_type0_posCut_0-0.049000m.root", i));
    TH1D *hTemp = (TH1D*)f.Get("Super sum");
    hTotalFierz->Add(hTemp);
    f.Close();
  }

  hTotalFierz->Draw();

  // Making files and printing histograms.

  TFile fData(Form("Octets_%i-%i_ssDataHist_type0_radialCut_0-49mm.root", octetLow, octetHigh), "RECREATE");
  hTotalData->Write();
  fData.Close();

/*
  TFile fBeta(Form("FullBlind_Feb2019_MC_A_0_b_0_Octets_80-121_ssHist_type0_posCut_0-49mm.root"), "RECREATE");
  hTotalBeta->Write();
  fBeta.Close();
*/
/*
  TFile fFierz(Form("MC_A_0_b_inf_Octets_80-121_ssHist_type0_posCut_0-49mm.root"), "RECREATE");
  hTotalFierz->Write();
  fFierz.Close();
*/

}
