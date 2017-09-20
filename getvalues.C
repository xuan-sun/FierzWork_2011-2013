getvalues()
{
  TString treeName = Form("Evts");
  TChain *MCTheoryChainBeta = new TChain(treeName.Data());

  TChain *MCTheoryChainFierz = new TChain(treeName.Data());

  for(int i = 0; i < 100; i++)
  {
    MCTheoryChainBeta->Add(Form("/home/xuansun/Documents/Analysis_Code/ucna_g4_2.1/UCN/UK_EventGen_2016/Evts_Files/b_0_300mill/Evts_%d.root", i));
    MCTheoryChainFierz->Add(Form("/home/xuansun/Documents/Analysis_Code/ucna_g4_2.1/UCN/UK_EventGen_2016/Evts_Files/b_inf_100mill/Evts_%d.root", i));
  }
  TString variableName = Form("KE");

  TH1D* mcTheoryHistBeta = new TH1D("mcBeta", "Test of comparehist code", 1000, 0, 1000);
  TH1D* mcTheoryHistFierz = new TH1D("mcFierz", "Test of comparehist code", 1000, 0, 1000);

  MCTheoryChainBeta -> Draw("KE >> mcBeta");
  MCTheoryChainFierz -> Draw("KE >> mcFierz");

  ofstream outfile_0;
  outfile_0.open("b_0_100mill_histogramContents.txt", ios::app);
  for(int i = 0; i < 1000; i++)
  {
    outfile_0 << mcTheoryHistBeta->GetBinCenter(i) << "\t"
	      << mcTheoryHistBeta->GetBinContent(i) << "\n";
  }
  outfile_0.close();

  ofstream outfile_inf;
  outfile_inf.open("b_inf_100mill_histogramContents.txt", ios::app);
  for(int i = 0; i < 1000; i++)
  {
    outfile_inf << mcTheoryHistFierz->GetBinCenter(i) << "\t"
                << mcTheoryHistFierz->GetBinContent(i) << "\n";
  }
  outfile_inf.close();
}

