Residuals()
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
    TFile f(Form("../PositionCuts/radialCut_0-49/MC_A_0_b_0_Octet_%i_ssHist_type0_posCut_0-0.049000m.root", i));
    TH1D *hTemp = (TH1D*)f.Get("Super sum");
    hTotalBeta->Add(hTemp);
    f.Close();
  }

  hTotalBeta->Draw();

  double entries = 0;

  int fitBinMin = 65;
  int fitBinMax = 100;

  for(int i = fitBinMin; i <= fitBinMax; i++)
  {
    entries = entries + hTotalBeta->GetBinContent(i);
  }

  cout << "Bins " << fitBinMin << " - " << fitBinMax << ": " << entries << endl;

  // normalize all the histograms using the fit ranges
/*  TCanvas *c2011 = new TCanvas("c2011", "c2011");
  c2011->cd();

  int fitBinMin = 17;
  int fitBinMax = 65;

  double yAxisMin = -0.1;
  double yAxisMax = 0.1;

  double N_data = 0;
  double N_beta = 0;

  for(int i = fitBinMin; i <= fitBinMax; i++)
  {
    N_data = N_data + hTotalData->GetBinContent(i);
    N_beta = N_beta + hTotalBeta->GetBinContent(i);
  }

  hTotalBeta->Scale(N_data / N_beta);

  // draw the 2011-2012 residual histogram
  TH1D *hResidual = new TH1D("residual", "residual", 120, 0, 1200);
  hResidual->Sumw2();
  hResidual->Add(hTotalData, hTotalBeta, 1.0, -1.0);

  hResidual->Divide(hTotalBeta);

  hResidual->SetTitle("2011-2012, data octets, (S_{data}-S_{MC}) / S_{MC}");
  hResidual->GetYaxis()->SetRangeUser(yAxisMin, yAxisMax);
  hResidual->GetXaxis()->SetTitle("Reconstructed Energy (keV)");
  hResidual->GetYaxis()->SetTitle("Fractional residual");

  hResidual->Draw();

  TLine *y0 = new TLine(0, 0, 1200, 0);
  y0->Draw();

  TLine *xMin = new TLine(hResidual->GetBinCenter(fitBinMin), yAxisMin, hResidual->GetBinCenter(fitBinMin), yAxisMax);
  xMin->Draw();
  TLine *xMax = new TLine(hResidual->GetBinCenter(fitBinMax), yAxisMin, hResidual->GetBinCenter(fitBinMax), yAxisMax);
  xMax->Draw();
*/
}
