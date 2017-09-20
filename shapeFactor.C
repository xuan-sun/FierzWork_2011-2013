shapeFactor(int octNb)
{
  TFile fData(Form("/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/ExtractedHistograms/Data_Hists/Octet_%i_ssDataHist_noCuts.root", octNb));
  TFile fMCSM(Form("/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/ExtractedHistograms/MC_A_0_b_0/MC_A_0_b_0_Octet_%i_ssHist.root", octNb));

  TH1D *hData = new TH1D("Data", "Data", 100, 0, 1000);
  TH1D *hMCSM = new TH1D("MC", "MC", 100, 0, 1000);

  hData = (TH1D*)fData.Get("Super sum");
  hData->Scale(1.0/hData->ComputeIntegral());
  hMCSM = (TH1D*)fMCSM.Get("Super sum");
  hMCSM->Scale(1.0/hMCSM->ComputeIntegral());

  vector <double> energy;
  vector <double> shape;
  int numPoints = 0;

  for(int i = 0; i < 100; i++)
  {
    cout << "At index " << i << " the hData bin contents is " << hData->GetBinContent(i)
	<< " and the hMCSM bin contents is " << hMCSM->GetBinContent(i) << endl;

    energy.push_back(hData->GetXaxis()->GetBinCenter(i));
    if(hMCSM->GetBinContent(i) == 0)
    {
      shape.push_back(0);
    }
    else
    {
      shape.push_back((hData->GetBinContent(i) - hMCSM->GetBinContent(i)) / hMCSM->GetBinContent(i));
    }
    numPoints++;
  }

  TGraph *g = new TGraph(numPoints, &(energy[0]), &(shape[0]));
  g->Draw("AP");


}

