#include	 "BetaSpectrum.hh"

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

// forward declarations for useful functions
TH1D* CalculateAverageSRAsymmetry(vector <TH1D*> asymms);


// required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

// visualization functions
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString xAxis, TString yAxis, TString command);

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  TCanvas *C = new TCanvas("canvas", "canvas");
  gROOT->SetStyle("Plain");	//on my computer this sets background to white, finally!

  cout << "Begin loading all files and histograms into vectors..." << endl;

  vector <TFile*> octetFiles;
  vector <TH1D*> octetAsymm;
  for(int i = 0; i < 60; i++)
  {
    if(i == 9 || i == 59)
    {
      continue;
    }
    octetAsymm.push_back(new TH1D());
    octetFiles.push_back(new TFile(Form("Xuan_asymmetries/MC_asymm_Octet_%i_type0.root", i)));
    octetAsymm[octetAsymm.size() - 1] = (TH1D*)octetFiles[octetFiles.size() - 1]->Get("AofE");
  }

  cout << "After loading, our total number of octets used is " << octetAsymm.size() << endl;

  cout << "Calculating statistically weighted average of asymmetries..." << endl;

  TH1D* weightedAvgAsymm = CalculateAverageSRAsymmetry(octetAsymm);

  cout << "Plotting final asymmetry..." << endl;

  PlotHist(C, 1, 1, weightedAvgAsymm, "Asymmetry as a function of Energy", "Reconstructed Energy (keV)", "Asymmetry", "");

  //prints the canvas with a dynamic TString name of the name of the file
  C -> Print("statsWeighted_SRAsymm_allOctets_2011-2012.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot->GetXaxis()->SetTitle(xAxis);
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle(yAxis);
  hPlot -> GetYaxis() -> CenterTitle();
  hPlot->GetYaxis()->SetRangeUser(0.02, 0.08);

  if(styleIndex == 1)
  {
    hPlot -> SetFillColor(46);
    hPlot -> SetFillStyle(3004);
  }
  if(styleIndex == 2)
  {
    hPlot -> SetFillColor(38);
    hPlot -> SetFillStyle(3005);
  }
  if(styleIndex == 3)
  {
    hPlot->SetFillColor(30);
    hPlot->SetFillStyle(3003);
  }
  if(styleIndex == 4)
  {
    hPlot->SetFillColor(41);
    hPlot->SetFillStyle(3016);
  }

  hPlot -> Draw(command);

  C->Update();

  TLine *yLow = new TLine(180, gPad->GetUymin(), 180, gPad->GetUymax());
  TLine *yHigh = new TLine(780, gPad->GetUymin(), 780, gPad->GetUymax());
  yLow->SetLineWidth(2);
  yHigh->SetLineWidth(2);
  yLow->SetLineColor(1);
  yHigh->SetLineColor(1);
  yLow->Draw("SAME");
  yHigh->Draw("SAME");

}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command)
{
  C->cd(canvasIndex);
  gPlot->SetTitle(title);
  gPlot->GetXaxis()->SetTitle(xAxis);
  gPlot->GetXaxis()->CenterTitle();
  gPlot->GetYaxis()->SetTitle(yAxis);
  gPlot->GetYaxis()->CenterTitle();

  if(styleIndex == 1)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(2);
  }
  if(styleIndex == 2)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(4);
  }

  gPlot->Draw(command);

  C->Update();

}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString xAxis, TString yAxis, TString command)
{
  C->cd(canvasIndex);
  gPlot->SetTitle(title);
  gPlot->GetXaxis()->SetTitle(xAxis);
  gPlot->GetXaxis()->CenterTitle();
  gPlot->GetYaxis()->SetTitle(yAxis);
  gPlot->GetYaxis()->CenterTitle();

  if(styleIndex == 1)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(2);
  }
  if(styleIndex == 2)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(4);
  }
  if(styleIndex == 3)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(3);
  }

  gPlot->Draw(command);

  C->Update();
}


TH1D* CalculateAverageSRAsymmetry(vector <TH1D*> asymms)
{
  TH1D* avgAsymm = new TH1D("avgAsymm", "Statistically weighted asymmetry", 120, 0, 1200);

  double thisBinContent_numerator = 0;
  double thisBinContent_denominator = 0;

  // this calculation of statistically weighted average uses the input error as sigma
  // formula can be found on Wikipedia: https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
  for(int i = 0; i <= avgAsymm->GetNbinsX(); i++)
  {
    for(unsigned int j = 0; j < asymms.size(); j++)
    {
      if(asymms[j]->GetBinContent(i) != 0 && asymms[j]->GetBinError(i) != 0)
      {
        thisBinContent_numerator = thisBinContent_numerator
		     		 + asymms[j]->GetBinContent(i) * (1.0 / ( asymms[j]->GetBinError(i) * asymms[j]->GetBinError(i) ) );
        thisBinContent_denominator = thisBinContent_denominator
		     		   + (1.0 / ( asymms[j]->GetBinError(i) * asymms[j]->GetBinError(i) ) );
      }
    }

    if(thisBinContent_numerator == 0 || thisBinContent_denominator == 0)
    {
      avgAsymm->SetBinContent(i, 0);
      avgAsymm->SetBinError(i, 0);
    }
    else
    {
      avgAsymm->SetBinContent(i, thisBinContent_numerator / thisBinContent_denominator);
      avgAsymm->SetBinError(i, 1.0 / sqrt(thisBinContent_denominator));
    }

    thisBinContent_numerator = 0;
    thisBinContent_denominator = 0;
  }

  return avgAsymm;
}
