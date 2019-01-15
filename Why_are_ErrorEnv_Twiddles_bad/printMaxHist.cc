// Standard C and C++ libraries
#include         <vector>
#include         <iostream>
#include         <algorithm>
#include         <functional>
#include         <fstream>
#include         <sstream>
#include         <stdlib.h>
#include         <math.h>
#include         <string.h>
#include         <time.h>
#include         <assert.h>

// Pretty much all the ROOT libraries I have ever used.
#include         <TROOT.h>
#include         <TSystem.h>
#include         <TMath.h>
#include         <TF1.h>
#include         <TGaxis.h>
#include         <TGraph.h>
#include         <TGraphErrors.h>
#include         <TCanvas.h>
#include         <TApplication.h>
#include         <TH1.h>
#include         <TProfile.h>
#include         <TObjArray.h>
#include         <TStyle.h>
#include         <TMarker.h>
#include         <TPaveStats.h>
#include         <TPaveText.h>
#include         <TFile.h>
#include         <TLegend.h>
#include         <TLegendEntry.h>
#include         <TH2F.h>
#include         <TRandom.h>
#include         <TTree.h>
#include         <TChain.h>
#include         <TObjArray.h>
#include         <TFractionFitter.h>
#include         <TLatex.h>
#include         <TMatrixD.h>
#include         <TRandom3.h>
#include	 <TMinuit.h>

using            namespace std;

//required later for plot_program
//TApplication plot_program("FADC_readin",0,0,0,0);

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (data file index)" << endl;
    return 0;
  }

  int octNb = atoi(argv[1]);

  TCanvas *C = new TCanvas("canvas", "canvas");

  TFile f(Form("/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/ExtractedHistograms/Data_Hists/Octet_%i_ssDataHist_type0.root", octNb));
  TH1D *octHist = (TH1D*)f.Get("Super sum");

  int maxIndex = -1;
  double maxVal = -1;
  for(int i = 0; i <= octHist->GetNbinsX(); i++)
  {
    if(octHist->GetBinContent(i) > maxVal)
    {
      maxVal = octHist->GetBinContent(i);
      maxIndex = i;
    }
  }


  ofstream outfile;
  outfile.open(Form("FindingSpectrumPeak_DataHists_ssDataHist_%s.txt", "type0"), ios::app);
  outfile << octNb << "\t"
          << maxVal << "\t"
	  << octHist->GetBinCenter(maxIndex) << "\n";
  outfile.close();

  // prints the canvas with a dynamic TString name of the name of the file
//  C -> Print("output_newXuanFitter.png");
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

