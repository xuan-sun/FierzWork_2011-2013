Residuals()
{
  // load the sum of all the data octets, with proper error propagation (if ROOT isn't broken).
  TH1D *hTotalData = new TH1D("totalData", "totalData", 120, 0, 1200);
  hTotalData->Sumw2();

  for(int i = 0; i < 60; i++)
  {
    if(i == 9 || i == 59)
    {
      continue;
    }
    TFile f(Form("../PositionCuts/radialCut_0-49/Octet_%i_ssDataHist_type0_radialCut_0-49mm_endpointCorrected.root", i));
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

  for(int i = 0; i < 60; i++)
  {
    if(i == 9 || i == 59)
    {
      continue;
    }
    TFile f(Form("../PositionCuts/radialCut_0-49/MC_A_0_b_0_Octet_%i_ssHist_type0_posCut_0-0.049000m.root", i));
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

  for(int i = 0; i < 60; i++)
  {
    if(i == 9 || i == 59)
    {
      continue;
    }
    TFile f(Form("../PositionCuts/radialCut_0-49/MC_A_0_b_inf_Octet_%i_ssHist_type0_posCut_0-0.049000m.root", i));
    TH1D *hTemp = (TH1D*)f.Get("Super sum");
    hTotalFierz->Add(hTemp);
    f.Close();
  }

  hTotalFierz->Draw();



/*
  int numFilesIndexMin = 0;
  int numFilesIndexMax = 100;

  TH1D* mcTheoryHistBeta = new TH1D("mcTheoryHistBeta", "Base SM", 100, 0, 1000);
  int totalEntries = 0;
  for(int j = numFilesIndexMin; j < numFilesIndexMax; j++)
  {     // note that .c_str() converts a std::string into a useable %s in Form(), must like .Data() for TString
    TFile f(Form("/mnt/data2/xuansun/analyzed_files/2011-2012_geom_twiddles/A_0_b_0_baselineHistograms/Hist_noBlind_SimAnalyzed_2011-2012_Beta_paramSet_100_%i_type0_radialCut_0-49mm.root", j));
    TH1D* hTemp = (TH1D*)f.Get("Erecon blinded hist");
    for(int i = 0; i <= mcTheoryHistBeta->GetNbinsX(); i++)
    {
      mcTheoryHistBeta->SetBinContent(i, mcTheoryHistBeta->GetBinContent(i) + hTemp->GetBinContent(i));
    }
    totalEntries = totalEntries + hTemp->GetEntries();
    f.Close();
  }
  mcTheoryHistBeta->SetEntries(totalEntries);
  cout << "Loaded mcTheoryHistBeta with, after cuts, nEvents = " << mcTheoryHistBeta->GetEntries() << endl;

  // using fierz beta spectrum, histogram files
  totalEntries = 0;
  TH1D* mcTheoryHistFierz = new TH1D("mcTheoryHistFierz", "Fierz", 100, 0, 1000);
  for(int i = numFilesIndexMin; i < numFilesIndexMax; i++)
  {
    TFile f(Form("/mnt/data2/xuansun/analyzed_files/2011-2012_geom_twiddles/A_0_b_inf_baselineHistograms/Hist_noBlind_SimAnalyzed_2011-2012_Beta_paramSet_100_%i_type0_radialCut_0-49mm.root", i));
    TH1D* hTemp = (TH1D*)f.Get("Erecon blinded hist");
    for(int i = 0; i <= mcTheoryHistFierz->GetNbinsX(); i++)
    {
      mcTheoryHistFierz->SetBinContent(i, mcTheoryHistFierz->GetBinContent(i) + hTemp->GetBinContent(i));
    }
    totalEntries = totalEntries + hTemp->GetEntries();
    f.Close();
  }
  mcTheoryHistFierz->SetEntries(totalEntries);

  cout << "Loaded mcTheoryHistFierz with nEvents = " << mcTheoryHistFierz->GetEntries() << endl;
*/

}
