shapeFactor(int octNb)
{
  TCanvas *C = new TCanvas("canvas", "canvas");

  TFile fData(Form("/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/ExtractedHistograms/Data_Hists/Octet_%i_ssDataHist_noCuts.root", octNb));
  TFile fMCSM(Form("/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/ExtractedHistograms/MC_A_0_b_0/MC_A_0_b_0_Octet_%i_ssHist.root", octNb));

  TH1D *hData = new TH1D("Data", "Data", 100, 0, 1000);
  TH1D *hMCSM = new TH1D("MC", "MC", 100, 0, 1000);

  hData = (TH1D*)fData.Get("Super sum");
  hMCSM = (TH1D*)fMCSM.Get("Super sum");

  double hDataNorm = 0;
  double hMCNorm = 0;

  for(int i = 10; i < 65; i++)
  {
    hDataNorm = hDataNorm + hData->GetBinContent(i);
    hMCNorm = hMCNorm + hMCSM->GetBinContent(i);
  }

  hData->Scale(1.0/hDataNorm);
  hMCSM->Scale(1.0/hMCNorm);
//  hData->Scale(1.0/hData->GetBinContent(22));
//  hMCSM->Scale(1.0/hMCSM->GetBinContent(22));


  cout << "Value of DataNorm " << hDataNorm << endl;
  cout << "Value of MCNorm " << hMCNorm << endl;

  vector <double> energy;
  vector <double> energyErr;
  vector <double> shape;
  vector <double> shapeErr;
  double estimatedErr = 0;
  double shapeValue = 0;
  int numPoints = 0;

  for(int i = 0; i < 100; i++)
  {
    energy.push_back(hData->GetXaxis()->GetBinCenter(i));
    energyErr.push_back(5);	// error bar of 5 KeV aka half the bin width
    if(hMCSM->GetBinContent(i) == 0)
    {
      shape.push_back(0);
      shapeErr.push_back(0);
    }
    else
    {
      shapeValue = (hData->GetBinContent(i) - hMCSM->GetBinContent(i)) / hMCSM->GetBinContent(i);
      shape.push_back(shapeValue);

      cout << "After normalizing, bin contents at " << i << " are now " << hMCSM->GetBinContent(i)*hMCNorm << endl;
//      estimatedErr = sqrt(2.0)*TMath::Abs(shapeValue)*(1.0/(sqrt(hMCSM->GetBinContent(i))*hMCSM->GetBinContent(22)));
      estimatedErr = sqrt(2.0)*TMath::Abs(shapeValue)*(1.0 / sqrt(hMCSM->GetBinContent(i) * hMCNorm));
      shapeErr.push_back(estimatedErr);
    }
    numPoints++;
  }

  TGraphErrors *g = new TGraphErrors(numPoints, &(energy[0]), &(shape[0]), &(energyErr[0]), &(shapeErr[0]));
  g->SetTitle(Form("Shape Factor for Octet %i", octNb));
  g->SetMarkerSize(1);
  g->SetMarkerStyle(21);
  g->SetMarkerColor(38);
  g->GetHistogram()->SetMaximum(0.1);
  g->GetHistogram()->SetMinimum(-0.1);
  g->Draw("AP");

  C->Print(Form("ShapeFactor_%i.pdf", octNb));

}

