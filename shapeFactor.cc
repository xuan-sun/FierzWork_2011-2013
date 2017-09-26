#include         <iostream>
#include         <fstream>
#include         <TGaxis.h>
#include         <sstream>
#include         <TGraph.h>
#include         <TGraphErrors.h>
#include         <TCanvas.h>
#include         <TApplication.h>
#include         <stdlib.h>
#include         <TF1.h>
#include         <TH1.h>
#include         <TProfile.h>
#include         <TObjArray.h>
#include         <TStyle.h>
#include         <TMarker.h>
#include         <math.h>
#include         <TStyle.h>
#include         <TPaveStats.h>
#include         <TPaveText.h>
#include         <vector>
#include         <string.h>
#include         <fstream>
#include         <TROOT.h>
#include         <TFile.h>
#include         <TLegend.h>
#include         <TLegendEntry.h>
#include         <time.h>
#include         <TH2F.h>
#include         <assert.h>
#include         <string>
#include         <TRandom.h>
#include         <TTree.h>
#include         <TChain.h>
#include         <TVector.h>
#include         <vector>
#include         <utility>
#include         <TLeaf.h>
using            namespace std;

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

struct entry
{
  int octNb;
  double avg_mE;
  double chisquared;
  double ndf;
  double bErr;
  double b;
  double GluckErr;
  int entryNum;
  double numberOfOriginalEvents;
};

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (octet #)" << endl;
    return 0;
  }

  int octNb = atoi(argv[1]);

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

//      estimatedErr = sqrt(2.0)*TMath::Abs(shapeValue)*(1.0/(sqrt(hMCSM->GetBinContent(i))*hMCSM->GetBinContent(22)));
//      estimatedErr = sqrt(2.0)*TMath::Abs(shapeValue)*(1.0 / sqrt(hMCSM->GetBinContent(i) * hMCNorm));
      estimatedErr = sqrt(2.0)/sqrt(hMCSM->GetBinContent(i) * hMCNorm);
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

  // extract the fitted b value and plot the shape factor with 'b' on the same graph.






  C->Print(Form("ShapeFactor_%i.pdf", octNb));

  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

