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

#define		GEOM	"2012-2013"
#define		TYPE	"type0"
#define		FITMINBIN	17
#define		FITMAXBIN	65

//required later for plot_program
//TApplication plot_program("FADC_readin",0,0,0,0);

// Fundamental constants that get used
const double m_e = 511.00;                                              ///< electron mass, keV/c^2
double ndf = -1;		// value set in code by fitMax - fitMin - 1 (for the 1 parameter, b).

// Plot things for visualisation
void PlotHist(TCanvas *C, int styleIndex, int canvaxIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors* gPlot, TString title, TString command);

// Perform a few useful, simple calculations
double CalculateAveragemOverE(TH1D* gammaSM, int binMin, int binMax);

vector <double> binContentsMC0;
double avg_mE;

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (data file index)" << endl;
    return 0;
  }

  int numFilesIndexMin = 0;
  int numFilesIndexMax = 100;
  // using unblinded base beta spectrum
  TH1D* mcTheoryHistBeta = new TH1D("mcTheoryHistBeta", "Base SM", 100, 0, 1000);
  TChain* betaChain = new TChain("SimAnalyzed");
  for(int i = numFilesIndexMin; i < numFilesIndexMax; i++)
  {
    betaChain->AddFile(Form("/mnt/Data/xuansun/analyzed_files/%s_geom_twiddledAndBaselineSimulations/A_0_b_0/SimAnalyzed_%s_Beta_paramSet_100_%i.root", GEOM, GEOM, i));
  }
  betaChain->Draw("Erecon >> mcTheoryHistBeta", "PID == 1 && Erecon > 0 && type == 0 && side < 2");

  // the work beyond here is unrelated to which data structure you chose
  for(int i = 0; i < mcTheoryHistBeta->GetNbinsX(); i++)
  {
    binContentsMC0.push_back(mcTheoryHistBeta->GetBinContent(i));
  }

/*
  ofstream outfile;
  outfile.open(Form("gaussianTwiddles_noBlinding_newXuanFitter_bFitsForSyst_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), ios::app);
  outfile << octNb << "\t"
          << avg_mE << "\t"
	  << functionMin << "\t"
	  << ndf << "\t"
	  << functionMin/ndf << "\t"
          << fitVal << "\t"
          << fitErr << "\t"
	  << -1 << "\t"		// these -1's are placeholders so the format is same as combinedAbFitter.cc
	  << -1 << "\t"
	  << covMatrixStatus << "\n";
  outfile.close();
*/


  // plot everything and visualize
/*
  TCanvas *C = new TCanvas("canvas", "canvas");
  C->cd();
  gROOT->SetStyle("Plain");	//on my computer this sets background to white, finally!
  dataHist->Scale(1.0 / totalData);
  mcTheoryHistBeta->Scale(1.0 / totalMC0);
  mcTheoryHistFierz->Scale(1.0 / totalMCinf);
  PlotHist(C, 2, 1, mcTheoryHistFierz, "", "");
  PlotHist(C, 3, 1, dataHist, "", "HISTSAME");
  PlotHist(C, 1, 1, mcTheoryHistBeta, "", "SAME");


  TLegend *l = new TLegend(0.7, 0.6, 0.9, 0.8);
  l->AddEntry(dataHist, "Data", "f");
  l->AddEntry(mcTheoryHistBeta, "SM", "f");
  l->AddEntry(mcTheoryHistFierz, "Fierz", "f");
  l->Draw();
*/
  // update the legend to include valuable variables.
/*  TLatex t2;
  t2.SetTextSize(0.03);
  t2.SetTextAlign(13);
  t2.DrawLatex(600, 0.017, Form("#frac{#Chi^{2}}{NDF} = %f / %f = %f", functionMin, ndf, functionMin/ndf));
  TLatex t3;
  t3.SetTextSize(0.03);
  t3.SetTextAlign(13);
  t3.DrawLatex(700, 0.016, Form("b_{input} = %f", 0.1));
  TLatex t4;
  t4.SetTextSize(0.03);
  t4.SetTextAlign(13);
  t4.DrawLatex(700, 0.014, Form("b_{fit} = %f #pm %f", fitVal, fitErr));
*/
  // prints the canvas with a dynamic TString name of the name of the file
//  C -> Print("output_newXuanFitter.png");
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Energy (keV)");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("Counts (N)");
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(0, 0.000004);

  if(styleIndex == 1)
  {
    hPlot -> SetFillColor(46);
    hPlot -> SetFillStyle(3004);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 2)
  {
    hPlot -> SetFillColor(38);
    hPlot -> SetFillStyle(3005);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 3)
  {
    hPlot -> SetFillColor(29);
//    hPlot -> SetFillStyle(3005);
    hPlot -> SetFillStyle(3001);
  }

  hPlot -> Draw(command);
  C -> Update();
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors* gPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  gPlot -> SetTitle(title);
  gPlot -> GetXaxis() -> SetTitle("Energy (keV)");
  gPlot -> GetXaxis() -> SetRangeUser(0, 1000);
  gPlot -> GetXaxis() -> CenterTitle();
  gPlot -> GetYaxis() -> SetTitle("Hits (N)");
  gPlot -> GetYaxis() -> CenterTitle();

  if(styleIndex == 1)
  {
    gPlot -> SetMarkerColor(46);
  }
  if(styleIndex == 2)
  {
    gPlot -> SetMarkerColor(38);
  }
  if(styleIndex == 3)
  {
    gPlot -> SetMarkerColor(30);
  }
  gPlot -> SetMarkerStyle(7);
  gPlot -> SetMarkerSize(1);

  // add a flat line at y = 0
  TLine *line = new TLine(0, 0, 1000, 0);

  gPlot -> Draw(command);
  line -> Draw("SAME");

}


double CalculateAveragemOverE(TH1D* gammaSM, int binMin, int binMax)
{
  double num = 0;
  double denom = 0;

  for(int i = binMin; i < binMax; i++)
  {
    num = num + (m_e*gammaSM->GetBinContent(i)) / (gammaSM->GetBinCenter(i) + m_e);
    denom = denom + gammaSM->GetBinContent(i);
  }

  return num/denom;
}

