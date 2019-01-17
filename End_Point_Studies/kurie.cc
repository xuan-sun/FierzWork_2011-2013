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

vector <double> binContents;
vector <double> binErrors;
vector <double> energy;
vector <double> energyErr;
double avg_mE;

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (data file index)" << endl;
    return 0;
  }

  int twiddleIndex = atoi(argv[1]);

  // this much longer code loads trees and extracts the histograms that we're interested in for fitting. Done for simulation.
  TH1D* dataHist = new TH1D("dataHist", "Twiddle", 100, 0, 1000);
  TChain* dataChain = new TChain("SimAnalyzed");
  dataChain->AddFile(Form("/mnt/data2/xuansun/analyzed_files/%s_geom_twiddles/TwiddledSimFiles_A_0_b_0_matchingParamSet_15/SimAnalyzed_%s_Beta_paramSet_%i_0.root", GEOM, GEOM, twiddleIndex));
  dataChain->Draw("Erecon >> dataHist", "PID == 1 && Erecon > 0 && type == 0 && side < 2");

  cout << "Loaded dataHist with nEvents = " << dataHist->GetEntries() << ", indexed by " << twiddleIndex << endl;

/*
  int numFilesIndexMin = 0;
  int numFilesIndexMax = 100;
  // using unblinded base beta spectrum i.e. no twiddles, no input b
  TH1D* mcTheoryHistBeta = new TH1D("mcTheoryHistBeta", "Base SM", 100, 0, 1000);
  TChain* betaChain = new TChain("SimAnalyzed");
  for(int i = numFilesIndexMin; i < numFilesIndexMax; i++)
  {
    betaChain->AddFile(Form("/mnt/Data/xuansun/analyzed_files/%s_geom_twiddledAndBaselineSimulations/A_0_b_0/SimAnalyzed_%s_Beta_paramSet_100_%i.root", GEOM, GEOM, i));
  }
  betaChain->Draw("Erecon >> mcTheoryHistBeta", "PID == 1 && Erecon > 0 && type == 0 && side < 2");
  cout << "Completed loading betaChain with events equal to " << mcTheoryHistBeta->GetEntries() << endl;
*/


  // create vectors of the basic, loaded in histogram.
  for(int i = 0; i < dataHist->GetNbinsX(); i++)
  {
    energy.push_back(dataHist->GetBinCenter(i));
    energyErr.push_back(5.0);	// 5keV bin error is 1/2 bin width
    binContents.push_back(dataHist->GetBinContent(i));
    binErrors.push_back(dataHist->GetBinError(i));
  }

  // convert the histogram into a kurie plot.
  vector <double> kurie;
  vector <double> kurieErr;

  double W = 0;
  for(unsigned int i = 0; i < binContents.size(); i++)
  {
    W = energy[i] + 511;

    if(W < 511)
    {
      kurie.push_back(0);
      kurieErr.push_back(0);
      continue;
    }

    kurie.push_back(sqrt( binContents[i] / ( W * sqrt(W*W - 511.0*511.0) ) ));
    kurieErr.push_back(sqrt( binErrors[i] / ( W * sqrt(W*W - 511.0*511.0) ) ));
  }

  TGraphErrors *g2 = new TGraphErrors(kurie.size(), &(energy[0]), &(kurie[0]), &(energyErr[0]), &(kurieErr[0]));

  TF1 *fit1 = new TF1("fit1", "[0] + [1]*x", energy[FITMINBIN], energy[FITMAXBIN]);
  g2->Fit(fit1, "R");

  // plot everything and visualize
  TCanvas *C = new TCanvas("canvas", "canvas");
  C->cd();
  gROOT->SetStyle("Plain");	//on my computer this sets background to white, finally!

  PlotGraph(C, 2, 1, g2, "Kurie Plot for 10^8 b = 0 simulations", "AP");

  // draw fit parameters on the plot
  TLatex t0;
  t0.SetTextSize(0.03);
  t0.SetTextAlign(13);
  t0.DrawLatex(700, 1.75, Form("fit = [0]+[1]*x"));
  TLatex t1;
  t1.SetTextSize(0.03);
  t1.SetTextAlign(13);
  t1.DrawLatex(700, 1.5, Form("Range #in [%f, %f]", energy[FITMINBIN], energy[FITMAXBIN]));
  TLatex t2;
  t2.SetTextSize(0.03);
  t2.SetTextAlign(13);
  t2.DrawLatex(700, 1.25, Form("[0] = %f #pm %f", fit1->GetParameter(0), fit1->GetParError(0)));
  TLatex t3;
  t3.SetTextSize(0.03);
  t3.SetTextAlign(13);
  t3.DrawLatex(700, 1, Form("[1] = %f #pm %f", fit1->GetParameter(1), fit1->GetParError(1)));
  TLatex t4;
  t4.SetTextSize(0.03);
  t4.SetTextAlign(13);
  t4.DrawLatex(700, 0.75, Form("E_{endpoint, fit} = %f", -(fit1->GetParameter(0))/(fit1->GetParameter(1)) ));

  ofstream outfile;
  outfile.open(Form("endPointFits_GaussianTwiddles_%s_%s_Bins_%i-%i_index16.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), ios::app);
  outfile << twiddleIndex << "\t"
          << fit1->GetChisquare() << "\t"
          << fit1->GetNDF() << "\t"
          << fit1->GetChisquare() / fit1->GetNDF() << "\t"
          << fit1->GetParameter(0) << "\t"
          << fit1->GetParError(0) << "\t"
          << fit1->GetParameter(1) << "\t"
          << fit1->GetParError(1) << "\t"         // these -1's are placeholders so the format is same as combinedAbFitter.cc
          << -(fit1->GetParameter(0))/(fit1->GetParameter(1)) << "\t"
          << FITMINBIN << "\t"
	  << FITMAXBIN << "\t"
	  << energy[FITMINBIN] << "\t"
	  << energy[FITMAXBIN] << "\n";
  outfile.close();



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
  gPlot -> GetXaxis() -> SetTitle("KE (keV)"/*"Total Energy W (units of m_{e})"*/);
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

