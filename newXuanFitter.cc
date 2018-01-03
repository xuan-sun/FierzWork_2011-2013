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


#define		HIST_IMAGE_PRINTOUT_NAME	"TestPlot_newXuanFitter"
#define		TYPE	"allTypes"

//required later for plot_program
//TApplication plot_program("FADC_readin",0,0,0,0);

// Fundamental constants that get used
const double m_e = 511.00;                                              ///< electron mass, keV/c^2

// Plot things for visualisation
void PlotHist(TCanvas *C, int styleIndex, int canvaxIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors* gPlot, TString title, TString command);


// Perform a few useful, simple calculations
double CalculateAveragemOverE(TH1D* gammaSM, int binMin, int binMax);

// functions needed for TMinuit fitter
void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void CalculateChi2(double b);

double bestChi2 = 1000000000;
double bestFierz = 100000000;

vector <double> binContentsMC0;
vector <double> binContentsMCinf;
vector <double> binContentsData;
vector <double> binErrorsData;
double avg_mE;

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (octet #)" << endl;
    return 0;
  }

  int octNb = atoi(argv[1]);

  TFile fData(TString::Format("ExtractedHistograms/Data_Hists/Octet_%i_ssDataHist_%s.root", octNb, TYPE));
  TFile fMC0(TString::Format("/mnt/Data/xuansun/BLIND_MC_files/2011-2012_geom/BLIND_MC_A_0_b_0_Octet_%i_ssHist_%s.root", octNb, TYPE));
//  TFile fMC0(TString::Format("ExtractedHistograms/MC_A_0_b_0/MC_A_0_b_0_Octet_%i_ssHist%s.root", octNb, ""));
  TFile fMCinf(TString::Format("ExtractedHistograms/MC_A_0_b_inf/MC_A_0_b_inf_Octet_%i_ssHist%s.root", octNb, ""));
//  TFile fMCinf(TString::Format("/mnt/Data/xuansun/BLIND_MC_files/2011-2012_geom/BLIND_MC_A_0_b_inf_Octet_%i_ssHist_%s.root", octNb, TYPE));

  TH1D* dataHist = (TH1D*)fData.Get("Super sum");
  TH1D* mcTheoryHistBeta = (TH1D*)fMC0.Get("Super sum");
  TH1D* mcTheoryHistFierz = (TH1D*)fMCinf.Get("Super sum");
  int fitMin = 10;
  int fitMax = 65;

  for(int i = 0; i < dataHist->GetNbinsX(); i++)
  {
    binContentsMC0.push_back(mcTheoryHistBeta->GetBinContent(i));
    binContentsMCinf.push_back(mcTheoryHistFierz->GetBinContent(i));
    binContentsData.push_back(dataHist->GetBinContent(i));
    binErrorsData.push_back(dataHist->GetBinError(i));
  }

  double totalMC0 = 0;
  double totalMCinf = 0;
  for(int i = fitMin; i < fitMax; i++)
  {
    totalMC0 = totalMC0 + binContentsMC0[i];
    totalMCinf = totalMCinf + binContentsMCinf[i];
  }
  for(int i = fitMin; i < fitMax; i++)
  {
    binContentsMC0[i] = binContentsMC0[i] / totalMC0;
    binContentsMCinf[i] = binContentsMCinf[i] / totalMCinf;
  }

  cout << "Finished loading the bin contents and errors into vectors..." << endl;

  avg_mE = CalculateAveragemOverE(mcTheoryHistBeta, fitMin, fitMax);

  cout << "Set average m/E value between fit range " << dataHist->GetBinCenter(fitMin) << " and " << dataHist->GetBinCenter(fitMax)
	<< " at: " << avg_mE << endl;



  for(double i = -10; i <= 10; i = i + 0.1)
  {
    CalculateChi2(i);
  }

  cout << "Completed testing fit values from -10 to 10..." << endl;
  cout << "Best fit b = " << bestFierz << endl;
  cout << "With corresponding chi2/ndf = " << bestChi2 / 54 << endl;


/*
  TMinuit *gMinuit = new TMinuit(1);
  gMinuit->SetFCN(chi2);

  Double_t arglist[10];
  Int_t iflag = 0;

  arglist[0] = 1;

  gMinuit->mnexcm("SET ERR", arglist, 1, iflag);

  gMinuit->mnparm(0, "fierz", 0, 0.01, -100, 100, iflag);

  arglist[0] = 500;
  arglist[1] = 1;

  gMinuit->mnexcm("CALL FCN", arglist, 1, iflag);
  gMinuit->mnexcm("MIGRAD", arglist, 2, iflag);
*/



  // plot everything and visualize
/*  TCanvas *C = new TCanvas("canvas", "canvas");
  gROOT->SetStyle("Plain");	//on my computer this sets background to white, finally!
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat("en");
  gStyle->SetStatH(0.45);
  gStyle->SetStatW(0.45);
  PlotHist(C, 1, 1, dataHist, "Data histogram", "");
//  PlotHist(C, 2, 1, resultHist, "Fit Result histogram", "SAME");


  // update the legend to include valuable variables.
  TPaveStats *ps = (TPaveStats*)C->GetPrimitive("stats");
  ps->SetName("mystats");
  TList *listOfLines = ps->GetListOfLines();
  TLatex *myText1 = new TLatex(0,0,Form("#Chi^{2} = %f", chisquared));
  listOfLines->Add(myText1);
  TLatex *myText2 = new TLatex(0,0,Form("NDF = %d", ndf));
  listOfLines->Add(myText2);
  TLatex *myText3 = new TLatex(0,0,Form("#frac{#Chi^{2}}{NDF} = %f", chisquared/ndf));
  listOfLines->Add(myText3);
  TLatex *myText4 = new TLatex(0,0,Form("frac0 = %f #pm %f", frac0Val, frac0Err));
  listOfLines->Add(myText4);
  TLatex *myText5 = new TLatex(0,0,Form("frac1 = %f #pm %f", frac1Val, frac1Err));
  listOfLines->Add(myText5);
  // the following line is needed to avoid that the automatic redrawing of stats
  dataHist->SetStats(0);
*/
  // prints the canvas with a dynamic TString name of the name of the file
//  C -> Print(Form("%s.pdf", HIST_IMAGE_PRINTOUT_NAME));
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  double totChi2 = 0;
  double b = par[0];
  double fit = 0;

  for(unsigned int i = 10; i < 65; i++)
  {
    fit = ((binContentsMC0[i] + b*avg_mE*binContentsMCinf[i])*binContentsData[i]) / (1 + b*avg_mE);

    totChi2 = totChi2 + pow((binContentsData[i] - fit) / binErrorsData[i], 2.0);
  }

  f = totChi2;
}


void CalculateChi2(double b)
{
  double totChi2 = 0;
  double fit = 0;

  double totContentsData = 0;
  for(unsigned int i = 10; i < 65; i++)
  {
    totContentsData = totContentsData + binContentsData[i];
  }

  for(unsigned int i = 10; i < 65; i++)
  {
    fit = ((binContentsMC0[i] + b*avg_mE*binContentsMCinf[i])*(binContentsData[i]*totContentsData)) / (1 + b*avg_mE);

    totChi2 = totChi2 + pow((binContentsData[i] - fit) / binErrorsData[i], 2.0);
/*
    if(b > -1.71 && b < -1.69) // there's some rounding error so b == -1.7 isn't working
    {
      cout << "For b value " << b << " we have binContentsData[" << i << "] - fit = " << binContentsData[i] - fit << endl;
      cout << "\t binContentsData[" << i << "] = " << binContentsData[i] << endl;
      cout << "\t fit = " << fit << endl;
      cout << "\t binErrorsData[" << i << "] = " << binErrorsData[i] << endl;
    }
*/
  }

  if(totChi2 < bestChi2)
  {
    bestChi2 = totChi2;
    bestFierz = b;
  }

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

