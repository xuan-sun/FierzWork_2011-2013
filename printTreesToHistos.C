printTreesToHistos()
{
  for(int i = 0; i < 60; i++)
  {
    if(i == 9 || i == 59 )
    {
      continue;
    }

    TFile f(Form("ExtractedHistograms/Data_Hists/Octet_%i_ssDataHist_type0.root", i));
    TH1D* h = (TH1D*)f.Get("Super sum");
    h->SetTitle(Form("Type 0, Octet %i", i));
    h->Draw();
    gPad->Print(Form("Spectrum_Octet_%i_type0_hist.pdf", i));
  }
}
