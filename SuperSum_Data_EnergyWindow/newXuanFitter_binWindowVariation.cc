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

#define		GEOM	"2011-2012"
#define		TYPE	"type0"
#define		RADLOW		0
#define		RADHIGH		49
#define		RADLOWFLOAT	0
#define		RADHIGHFLOAT	0.049000

int global_fitMinBin = 0;
int global_fitMaxBin = 0;

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
  if(argc < 4)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (octet number) (fit bin min) (fit bin max)" << endl;
    return 0;
  }

  int octNb = atoi(argv[1]);
  global_fitMinBin = atoi(argv[2]);
  global_fitMaxBin = atoi(argv[3]);

/*
  // this little bit loads the octets once they have already been separated into super sum histograms
  TFile fData(TString::Format("/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/PositionCuts/radialCut_%i-%i/Octet_%i_ssDataHist_%s_radialCut_%i-%imm.root", RADLOW, RADHIGH, octNb, TYPE, RADLOW, RADHIGH));
  TFile fMC0(TString::Format("/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/PositionCuts/radialCut_%i-%i/FullBlind_Feb2019_MC_A_0_b_0_Octet_%i_%s_posCut_%i-%fm.root", RADLOW, RADHIGH, octNb, TYPE, RADLOWFLOAT, RADHIGHFLOAT));
  TFile fMCinf(TString::Format("/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/PositionCuts/radialCut_%i-%i/MC_A_0_b_inf_Octet_%i_ssHist_%s_posCut_%i-%fm.root", RADLOW, RADHIGH, octNb, TYPE, RADLOWFLOAT, RADHIGHFLOAT));
  TH1D* dataHist = (TH1D*)fData.Get("Super sum");
  TH1D* mcTheoryHistBeta = (TH1D*)fMC0.Get("Super sum");
  TH1D* mcTheoryHistFierz = (TH1D*)fMCinf.Get("Super sum");
*/

  // this loads all the histograms that have all the beta events already summed into one histogram. i.e. the integrated dataset
  TFile fData(TString::Format("All_Octets_Summed_Histograms_readNotes/Octets_0-60_ssDataHist_%s_radialCut_%i-%imm_endpointCorr.root", TYPE, RADLOW, RADHIGH));
  TFile fMC0(TString::Format("All_Octets_Summed_Histograms_readNotes/FullBlind_Feb2019_MC_A_0_b_0_Octets_0-59_ssHist_%s_posCut_%i-%imm.root", TYPE, RADLOW, RADHIGH));
  TFile fMCinf(TString::Format("All_Octets_Summed_Histograms_readNotes/MC_A_0_b_inf_Octets_0-59_ssHist_%s_posCut_%i-%imm.root", TYPE, RADLOW, RADHIGH));
  TH1D* dataHist = (TH1D*)fData.Get("totalData");
  TH1D* mcTheoryHistBeta = (TH1D*)fMC0.Get("totalBeta");
  TH1D* mcTheoryHistFierz = (TH1D*)fMCinf.Get("totalFierz");



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
  for(int i = global_fitMinBin; i < global_fitMaxBin; i++)
  {
    totalMC0 = totalMC0 + binContentsMC0[i];
    totalMCinf = totalMCinf + binContentsMCinf[i];
    totalData = totalData + binContentsData[i];
    counter_ndf = counter_ndf + 1;
  }
  ndf = counter_ndf - 1;	// -1 is the single parameter in the fit.
  for(int i = global_fitMinBin; i < global_fitMaxBin; i++)
  {
    binContentsMC0[i] = binContentsMC0[i] / totalMC0;
    binContentsMCinf[i] = binContentsMCinf[i] / totalMCinf;
  }

  cout << "Finished loading the bin contents and errors into vectors..." << endl;

  avg_mE = CalculateAveragemOverE(mcTheoryHistBeta, global_fitMinBin, global_fitMaxBin);

  cout << "Set average m/E value between fit range " << dataHist->GetBinCenter(global_fitMinBin)
	<< " and " << dataHist->GetBinCenter(global_fitMaxBin)
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
  outfile.open(Form("allOctets_positionCuts_%i-%imm_endpointCorrected_withFullBlind_Feb2019_%s_%s_highBinVari_run3.txt", RADLOW, RADHIGH, TYPE, GEOM), ios::app);
//  outfile.open(Form("positionCuts_%i-%imm_noGainCorrection_withFullBlind_Feb2019_%s_%s_binWindowVariations_individualOctets.txt", RADLOW, RADHIGH, TYPE, GEOM), ios::app);
  outfile << /*octNb*/ "ALL" << "\t"
          << avg_mE << "\t"
	  << functionMin << "\t"
	  << ndf << "\t"
	  << functionMin/ndf << "\t"
	  << TMath::Prob(functionMin, ndf) << "\t"
          << fitVal << "\t"
          << fitErr << "\t"
	  << global_fitMinBin << "\t"
	  << dataHist->GetBinCenter(global_fitMinBin) << "\t"
	  << global_fitMaxBin << "\t"
	  << dataHist->GetBinCenter(global_fitMaxBin) << "\t"
	  << covMatrixStatus << "\n";
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
  for(unsigned int i = global_fitMinBin; i < global_fitMaxBin; i++)
  {
    totContentsData = totContentsData + binContentsData[i];
  }

  for(unsigned int i = global_fitMinBin; i < global_fitMaxBin; i++)
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

