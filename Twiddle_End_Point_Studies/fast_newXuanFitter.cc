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
#include	 <TSystem.h>
#include	 <TCut.h>

using            namespace std;

#define		GEOM	"2011-2012"
#define		TYPE	"type0"
#define		FITMINBIN	17
#define		FITMAXBIN	65
#define		SAVE_INDEX	"index19"

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

// functions needed for TMinuit fitter
void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

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
    cout << "(executable) (data file index)" << endl;
    return 0;
  }

  int index = atoi(argv[1]);

  TString dataFilePath = Form("/mnt/data2/xuansun/analyzed_files/%s_geom_twiddles/TwiddledSimFiles_A_0_b_0_matchingParamSet_19/SimAnalyzed_%s_Beta_paramSet_%i_0.root", GEOM, GEOM, index);

  // checks if the file exists because we threw out a bunch of twiddles
  if(gSystem->AccessPathName(dataFilePath))
  {
    cout << "File does not exist. Exiting program." << endl;
    return 0;
  }
  else
  {
    cout << "File exists. Continuing..." << endl;
  }

  int radialCutLow_input = 0;
  int radialCutHigh_input = 49;

  double radialCutLow = 0.0 * sqrt(1.0 / 0.6);
  double radialCutHigh = 0.049 * sqrt(1.0 / 0.6);

  TCut positionCut = Form("((MWPCPos.MWPCPosE[0]*MWPCPos.MWPCPosE[0] + MWPCPos.MWPCPosE[1]*MWPCPos.MWPCPosE[1] <= %f && MWPCPos.MWPCPosE[0]MWPCPos.MWPCPosE[0] + MWPCPos.MWPCPosE[1]*MWPCPos.MWPCPosE[1] >= %f && MWPCPos.MWPCPosW[0]*MWPCPos.MWPCPosW[0] + MWPCPos.MWPCPosW[1]*MWPCPos.MWPCPosW[1] == 0 ) || (MWPCPos.MWPCPosE[0]*MWPCPos.MWPCPosE[0] + MWPCPos.MWPCPosE[1]*MWPCPos.MWPCPosE[1] == 0 && MWPCPos.MWPCPosW[0]*MWPCPos.MWPCPosW[0] + MWPCPos.MWPCPosW[1]*MWPCPos.MWPCPosW[1] <= %f && MWPCPos.MWPCPosW[0]*MWPCPos.MWPCPosW[0] + MWPCPos.MWPCPosW[1]*MWPCPos.MWPCPosW[1] >= %f ))",
			radialCutHigh*radialCutHigh, radialCutLow*radialCutLow, radialCutHigh*radialCutHigh, radialCutLow*radialCutLow);

  // this much longer code loads trees and extracts the histograms that we're interested in for fitting. Done for simulation.
  TH1D* dataHist = new TH1D("dataHist", "Twiddle", 100, 0, 1000);
  TChain* dataChain = new TChain("SimAnalyzed");
  dataChain->AddFile(dataFilePath);
  dataChain->Draw("Erecon_corr_r49mm >> dataHist", "PID == 1 && Erecon_corr_r49mm > 0 && type == 0 && side < 2");

  cout << "Loaded dataHist with nEvents = " << dataHist->GetEntries() << ", indexed by " << index << endl;


  int numFilesIndexMin = 0;
  int numFilesIndexMax = 100;

  // using base beta spectrum, histogram files
  TH1D* mcTheoryHistBeta = new TH1D("mcTheoryHistBeta", "Base SM", 100, 0, 1000);
  int totalEntries = 0;
  for(int j = numFilesIndexMin; j < numFilesIndexMax; j++)
  {	// note that .c_str() converts a std::string into a useable %s in Form(), must like .Data() for TString
    TFile f(Form("/mnt/data2/xuansun/analyzed_files/%s_geom_twiddles/A_0_b_0_baselineHistograms/Hist_noBlind_SimAnalyzed_%s_Beta_paramSet_100_%i_type0_radialCut_%i-%imm.root", GEOM, GEOM, j, radialCutLow_input, radialCutHigh_input));
    TH1D* hTemp = (TH1D*)f.Get("Erecon blinded hist");
    for(int i = 0; i <= mcTheoryHistBeta->GetNbinsX(); i++)
    {
      mcTheoryHistBeta->SetBinContent(i, mcTheoryHistBeta->GetBinContent(i) + hTemp->GetBinContent(i));
    }
    totalEntries = totalEntries + hTemp->GetEntries();
    f.Close();
  }
  mcTheoryHistBeta->SetEntries(totalEntries);
  cout << "Loaded mcTheoryHistBeta with, after cuts, nEvents = " << mcTheoryHistBeta->GetEntries() << endl;

  // using fierz beta spectrum, histogram files
  totalEntries = 0;
  TH1D* mcTheoryHistFierz = new TH1D("mcTheoryHistFierz", "Fierz", 100, 0, 1000);
  for(int i = numFilesIndexMin; i < numFilesIndexMax; i++)
  {
    TFile f(Form("/mnt/data2/xuansun/analyzed_files/%s_geom_twiddles/A_0_b_inf_baselineHistograms/Hist_noBlind_SimAnalyzed_%s_Beta_paramSet_100_%i_type0_radialCut_%i-%imm.root", GEOM, GEOM, i, radialCutLow_input, radialCutHigh_input));
    TH1D* hTemp = (TH1D*)f.Get("Erecon blinded hist");
    for(int i = 0; i <= mcTheoryHistFierz->GetNbinsX(); i++)
    {
      mcTheoryHistFierz->SetBinContent(i, mcTheoryHistFierz->GetBinContent(i) + hTemp->GetBinContent(i));
    }
    totalEntries = totalEntries + hTemp->GetEntries();
    f.Close();
  }
  mcTheoryHistFierz->SetEntries(totalEntries);

  cout << "Loaded fierzChain with nEvents = " << mcTheoryHistFierz->GetEntries() << endl;


  // the work beyond here is unrelated to which data structure you chose
  for(int i = 0; i < dataHist->GetNbinsX(); i++)
  {
    binContentsMC0.push_back(mcTheoryHistBeta->GetBinContent(i));
    binContentsMCinf.push_back(mcTheoryHistFierz->GetBinContent(i));
    binContentsData.push_back(dataHist->GetBinContent(i));
    binErrorsData.push_back(dataHist->GetBinError(i));
  }

  double totalMC0 = 0;
  double totalMCinf = 0;
  double totalData = 0;
  double counter_ndf = 0;
  for(int i = FITMINBIN; i < FITMAXBIN; i++)
  {
    totalMC0 = totalMC0 + binContentsMC0[i];
    totalMCinf = totalMCinf + binContentsMCinf[i];
    totalData = totalData + binContentsData[i];
    counter_ndf = counter_ndf + 1;
  }
  ndf = counter_ndf - 1;	// -1 is the single parameter in the fit.
  for(int i = FITMINBIN; i < FITMAXBIN; i++)
  {
    binContentsMC0[i] = binContentsMC0[i] / totalMC0;
    binContentsMCinf[i] = binContentsMCinf[i] / totalMCinf;
  }

  cout << "Finished loading the bin contents and errors into vectors..." << endl;

  avg_mE = CalculateAveragemOverE(mcTheoryHistBeta, FITMINBIN, FITMAXBIN);

  cout << "Set average m/E value between fit range " << dataHist->GetBinCenter(FITMINBIN)
	<< " and " << dataHist->GetBinCenter(FITMAXBIN)
	<< " at: " << avg_mE << endl;

  TMinuit *gMinuit = new TMinuit(1);
  gMinuit->SetFCN(chi2);

  Double_t arglist[10];
  Int_t iflag = 0;

  arglist[0] = 2.0;	// here arglist[0] = 2 means use strategy 2, aka make sure the errors are right
  gMinuit->mnexcm("SET STR", arglist, 1, iflag);

  arglist[0] = 1.0;	// 1.0 = 1 standard deviation errors. If you want 2 standard deviations, use 4.0.
  gMinuit->mnexcm("SET ERR", arglist, 1, iflag);

  gMinuit->mnparm(0, "fierz", 0, 0.01, -10, 10, iflag);

  arglist[0] = 500;
  arglist[1] = 1;

  gMinuit->mnexcm("CALL FCN", arglist, 1, iflag);
  gMinuit->mnexcm("MIGRAD", arglist, 2, iflag);

  cout << "------------ Fit is finished --------------" << endl;

  // can also use GetParameter(...);
  //  double paramValue, paramError;
  //  gMinuit->GetParameter(0, paramValue, paramError);

  Int_t internalInt;
  Double_t fitVal, fitErr, lowLimit, highLimit;
  TString paramName;

  // gets the values of parameter 0
  gMinuit->mnpout(0, paramName, fitVal, fitErr, lowLimit, highLimit, internalInt);

  cout << "paramName = " << paramName << endl;
  cout << "fitVal = " << fitVal << endl;
  cout << "fitErr = " << fitErr << endl;
  cout << "lowLimit = " << lowLimit << endl;
  cout << "highLimit = " << highLimit << endl;
  cout << "internalInt = " << internalInt << endl;

  // get the internal statistics aka chi-squared
  Double_t functionMin, distanceToMin, errorParamDef;
  Int_t nVariableParams, nUserParams, covMatrixStatus;
  gMinuit->mnstat(functionMin, distanceToMin, errorParamDef, nVariableParams, nUserParams, covMatrixStatus);

  cout << "Minimum chi-squared: " << functionMin << endl;
  cout << "Estimated vertical distance to min: " << distanceToMin << endl;
  cout << "Value of UP defining parameter uncertainties: " << errorParamDef << endl;
  cout << "Number of variable parameters: " << nVariableParams << endl;
  cout << "Highest number of parameters defined by user: " << nUserParams << endl;
  cout << "Status of covariance matrix: " << covMatrixStatus << endl;


  ofstream outfile;
  outfile.open(Form("asymmetric_gaussianTwiddles_noBlinding_endpointCorr_fast-newXuanFitter_bFitsForSyst_%s_%s_Bins_%i-%i_%s_radialCut_%i-%imm.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN, SAVE_INDEX, radialCutLow_input, radialCutHigh_input), ios::app);
//  outfile.open(Form("CorrectBlindingOct2018_newXuanFitter_bFitsFromData_endpointCorr_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), ios::app);
//  outfile.open(Form("gaussianTwiddles_CorrectBlindingOct2018_newXuanFitter_bFitsForSyst_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), ios::app);
//  outfile.open(Form("Twiddles_CorrectBlindingOct2018_newXuanFitter_bFitForSystError_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), ios::app);
  outfile << index << "\t"
          << avg_mE << "\t"
	  << functionMin << "\t"
	  << ndf << "\t"
	  << functionMin/ndf << "\t"
          << fitVal << "\t"
          << fitErr << "\t"
	  << -1 << "\t"		// these -1's are placeholders so the format is same as combinedAbFitter.cc
	  << -1 << "\t"
	  << covMatrixStatus << "\n";
  // adding 3 new parameters at the end for only the false b testing!
  // these need to be deleted because A) they don't matter later and B) they mess up plotting code
//	  << numFilesIndexMin << "\t"
//	  << numFilesIndexMax << "\t"
//	  << falseb << "\n";
  outfile.close();



  // plot everything and visualize
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

  // update the legend to include valuable variables.
/*  TLatex t2;
  t2.SetTextSize(0.03);
  t2.SetTextAlign(13);
  t2.DrawLatex(600, 0.017, Form("#frac{#Chi^{2}}{NDF} = %f / %f = %f", functionMin, ndf, functionMin/ndf));
  TLatex t3;
  t3.SetTextSize(0.03);
  t3.SetTextAlign(13);
  t3.DrawLatex(700, 0.016, Form("b_{input} = %f", 0.1));
*/  TLatex t4;
  t4.SetTextSize(0.03);
  t4.SetTextAlign(13);
  t4.DrawLatex(700, 0.014, Form("b_{fit} = %f #pm %f", fitVal, fitErr));

  // prints the canvas with a dynamic TString name of the name of the file
//  C -> Print("output_newXuanFitter.png");
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  double totChi2 = 0;
  double b = par[0];
  double fit = 0;

  double totContentsData = 0;
  for(unsigned int i = FITMINBIN; i < FITMAXBIN; i++)
  {
    totContentsData = totContentsData + binContentsData[i];
  }

  for(unsigned int i = FITMINBIN; i < FITMAXBIN; i++)
  {
    fit = ((binContentsMC0[i] + b*avg_mE*binContentsMCinf[i])*totContentsData) / (1 + b*avg_mE);

    totChi2 = totChi2 + pow((binContentsData[i] - fit) / binErrorsData[i], 2.0);
  }

  f = totChi2;
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

