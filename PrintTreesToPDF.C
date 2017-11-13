PrintTreesToPDF()
{

  for(int i = 0; i < 60; i++)
  {
    if(i == 9 || i == 59 /* || i == 10 || i == 11 */)
    {
      continue;
    }

    TFile f(Form("ExtractedHistograms/Data_Hists/Octet_%i_ssDataHist.root", i));
	cout << "Opened file.." << endl;

    TH1D *h = (TH1D*)f.Get("Super sum");
	cout << "Got histogram..." << endl;

    h->SetTitle(Form("Octet %i", i));
	cout << "Title is set..." << endl;

    h->Draw();
	cout << "The histogram is drawn..." << endl;

    gPad->Print(Form("SpectraForOctet_%i_allTypes.pdf", i));
	cout << "Done printing. End of program." << endl;
  }

}
