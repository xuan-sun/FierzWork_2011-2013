LoadHistograms()
{
  int numFilesIndexMin = 0;
  int numFilesIndexMax = 0;

  TH1D* mcTheoryHistBeta = new TH1D("mcTheoryHistBeta", "Base SM", 100, 0, 1000);
  int totalEntries = 0;
  for(int j = numFilesIndexMin; j < numFilesIndexMax; j++)
  {     // note that .c_str() converts a std::string into a useable %s in Form(), must like .Data() for TString
    TFile f(Form("/mnt/data2/xuansun/analyzed_files/%s_geom_twiddles/A_0_b_0_baselineHistograms/Hist_noBlind_SimAnalyzed_%s_Beta_paramSet_100_2011-2012_type0_radialCut_0-49mm.root"));
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
}
