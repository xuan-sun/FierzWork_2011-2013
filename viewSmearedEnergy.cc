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
#include	 <TLegend.h>

using            namespace std;

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void FillArrays(TString fileName, TH2D *h);
void FitGausToHistSlice(TH2D *h, int xBin);

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  TCanvas *C = new TCanvas("canvas", "canvas");
  C->cd();
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  TH2D *hNvErecon = new TH2D("NvErecon", "N vs Erecon", 110, 0, 1100, 1000, 0, 20000);

  for(int j = 0; j < 400; j++)
  {
    FillArrays(Form("/mnt/Data/xuansun/analyzed_files/All_Twiddles_Are_Baseline/SimAnalyzed_2011-2012_Beta_paramSet_%i_0.root", j), hNvErecon);
    cout << "Completed filling in counts from SimAnalyzed tree indexed by " << j << endl;
  }

  hNvErecon->Draw("COLZ");

  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->cd();

  for(int i = 1; i < 110; i++)
  {
    FitGausToHistSlice(hNvErecon, i);
  }

  //prints the canvas with a dynamic TString name of the name of the file
  C->Print("EreconSmeared_DistributionOfCounts_2D.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

void FillArrays(TString fileName, TH2D* h)
{
  TFile f(fileName);
  TTree *t = (TTree*)f.Get("SimAnalyzed");
  TH1D *hTemp = new TH1D("tmp", "tmp", 110, 0, 1100);
  t->Draw("EreconSmeared >> tmp", "PID == 1 && type == 0 && side < 2 && EreconSmeared > 0");

  for(int i = 0; i < hTemp->GetNbinsX(); i++)
  {
    h->Fill(hTemp->GetBinCenter(i), hTemp->GetBinContent(i));
  }
}

void FitGausToHistSlice(TH2D *h, int xBin)
{
  TH1D *hXSlice = h->ProjectionY("test", xBin, xBin, "");
  hXSlice->Fit("gaus");

  TF1 *f = hXSlice->GetFunction("gaus");

  cout << "Percentage difference for bin " << xBin << " is: "
	<< abs(sqrt(f->GetParameter(1)) - f->GetParameter(2)) / sqrt(f->GetParameter(1)) << endl;

  ofstream outfile;
  outfile.open("GausFitResultsToSpectra.txt", ios::app);
  outfile << xBin << "\t"
	  << f->GetParameter(0) << "\t"
	  << f->GetParameter(1) << "\t"
	  << f->GetParameter(2) << "\t"
	  << abs(sqrt(f->GetParameter(1)) - f->GetParameter(2)) / sqrt(f->GetParameter(1)) << "\n";
  outfile.close();
}
