{
//  gROOT->SetStyle("Pub");

  TString fileName = "../Asymmetry_Data_Fitting/AsymmetryDataFit_2011-2012.txt";

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2,1);
  c1->cd(1);

  TH1D* hb = new TH1D("b hist", "b hist", 50, -0.2, -0.1);
  hb->GetXaxis()->SetTitle("b");
  hb->GetYaxis()->SetTitle("N");
//  TH1D* hA = new TH1D("A hist", "A hist", 

  // makes a TGraphAsymmErrors to show the effects of fit windows
  vector <double> EWinLow;
  vector <double> EWinHigh;
  vector <double> EFitValue;

  vector <double> bFitValue;
  vector <double> bFitError1;
  vector <double> bFitError2;

  vector <double> AFitValue;
  vector <double> AFitError1;
  vector <double> AFitError2;

  int nPoints = 0;
  double octNb, avg_mE, chi2, ndf, chi2ndf;
  double E, Elow, Ehigh;
  double b, bErr;
  double A, AErr;

  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> octNb >> avg_mE
		>> Elow >> Ehigh
		>> chi2 >> ndf >> chi2ndf
		>> A >> AErr
		>> b >> bErr;

      E = 400;	// center all data points at 400keV

      EFitValue.push_back(E);
      EWinLow.push_back(E - Elow);
      EWinHigh.push_back(Ehigh - E);

      bFitValue.push_back(b);
      bFitError1.push_back(bErr);
      bFitError2.push_back(bErr);

      AFitValue.push_back(A);
      AFitError1.push_back(AErr);
      AFitError2.push_back(AErr);

      hb->Fill(b);

      nPoints++;
    }
  }
  infile.close();

  for(int i = 0; i < EFitValue.size(); i++)
  {
    cout << "EFitValue[" << i << "] = " << EFitValue[i] << endl;
    cout << "EWinLow[" << i << "] = " << EWinLow[i] << endl;
    cout << "EWinHigh[" << i << "] = " << EWinHigh[i] << endl;
    cout << "bFitValue[" << i << "] = " << bFitValue[i] << endl;
    cout << "bFitError1[" << i << "] = " << bFitError1[i] << endl;
    cout << "bFitError2[" << i << "] = " << bFitError2[i] << endl;
  }

  TGraphAsymmErrors *g = new TGraphAsymmErrors(nPoints, &(EFitValue[0]), &(bFitValue[0]), &(EWinLow[0]), &(EWinHigh[0]), &(bFitError1[0]), &(bFitError2[0]));
  g->SetMarkerStyle(21);
  g->SetMarkerColor(1);
  g->SetTitle("Varying fit windows on asymmetry data fits");
  g->GetXaxis()->SetTitle("E (keV)");
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->SetTitle("b");
  g->GetYaxis()->CenterTitle();
  g->Draw("AP");


  c1->cd(2);
  hb->Draw();

//  c1->Print("33_Ratio_bgfg_variousTimeWindows_allCuts.pdf");
}
