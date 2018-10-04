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

// functions needed for TMinuit fitter
void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

vector <double> binContentsMC0;
vector <double> binContentsMCinf;
vector <double> binContentsData;
vector <double> binErrorsData;
double avg_mE;

vector <double> asymmetriesData;
vector <double> asymmErrorsData;
vector <double> blind_asymmetriesData;
vector <double> blind_asymmErrorsData;
vector <double> energies;

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (index #)" << endl;
    return 0;
  }

  int octNb = atoi(argv[1]);

  // this needs to be here or else it won't compile due to ROOT dictionary issues.
  TCanvas *C = new TCanvas("canvas", "canvas");

  // Loading in energy spectra for fit
  // These energy spectra are revCalSim'd data supersums. Aka they are "data-like" and revCal'd according to octet.
//  TFile fData(TString::Format("ExtractedHistograms/Data_Hists/Octet_%i_ssDataHist_%s.root", octNb, TYPE));
  TFile fData(TString::Format("/mnt/Data/xuansun/analyzed_files/All_Twiddles_Are_Baseline/SimAnalyzed_2011-2012_Beta_paramSet_%i_0.root", octNb));
//  TFile fMC0(TString::Format("/mnt/Data/xuansun/BLIND_MC_files/Blinded_Oct2018/BLIND_MC_A_0_b_0_Octet_%i_%s.root", octNb, TYPE));
//  TFile fMC0(TString::Format("ExtractedHistograms/MC_A_0_b_0/MC_A_0_b_0_Octet_%i_ssHist_%s.root", octNb, TYPE));
//  TFile fMCinf(TString::Format("ExtractedHistograms/MC_A_0_b_inf/MC_A_0_b_inf_Octet_%i_ssHist_%s.root", octNb, TYPE));

//  TH1D* dataHist = (TH1D*)fData.Get("Super sum");

  // Load in a data histogram to get the super sum errors correct
  TFile fErrors(TString::Format("ExtractedHistograms/Data_Hists/Octet_%i_ssDataHist_%s.root", 43, TYPE));
  TH1D* dataErrorsHist = (TH1D*)fErrors.Get("Super sum");

  TH1D* dataHist = new TH1D("dataHist", "dataHist", 120, 0, 1200);
  TTree *t = (TTree*)fData.Get("SimAnalyzed");
  t->Draw("Erecon >> dataHist", "PID == 1 && type == 0 && side < 2 && Erecon >= 0");


//  TH1D* mcTheoryHistBeta = (TH1D*)fMC0.Get("Super sum");
//  TH1D* mcTheoryHistFierz = (TH1D*)fMCinf.Get("Super sum");

  // These 2 versions of mcTheoryHist load SimAnalyzed aka twiddle simulations.
  // I have made several million baseline twiddle simulations for the combined fitting.
  TH1D* mcTheoryHistBeta = new TH1D("mcTheoryHistBeta", "Base SM", 100, 0, 1000);
  int totalEntries = 0;
  for(int j = 0; j < 100; j++)
  {
    TFile f(Form("/mnt/Data/xuansun/analyzed_files/A_0_b_0/BLIND_SimAnalyzed_2011-2012_Beta_paramSet_100_%i_type0.root", j));
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

  TH1D* mcTheoryHistFierz = new TH1D("mcTheoryHistFierz", "Fierz", 100, 0, 1000);
  TChain* fierzChain = new TChain("SimAnalyzed");
  for(int i = 0; i < 100; i++)
  {
    fierzChain->AddFile(Form("/mnt/Data/xuansun/analyzed_files/A_0_b_inf/SimAnalyzed_2011-2012_Beta_paramSet_100_%i.root", i));
  }
  fierzChain->Draw("Erecon >> mcTheoryHistFierz", "PID == 1 && Erecon > 0 && type == 0 && side < 2");
  cout << "Loaded fierzChain with nEvents = " << fierzChain->GetEntries() << endl;



  for(int i = 0; i < dataHist->GetNbinsX(); i++)
  {
    binContentsMC0.push_back(mcTheoryHistBeta->GetBinContent(i));
    binContentsMCinf.push_back(mcTheoryHistFierz->GetBinContent(i));
    binContentsData.push_back(dataHist->GetBinContent(i));
    binErrorsData.push_back(dataErrorsHist->GetBinError(i));
//    binErrorsData.push_back(dataHist->GetBinError(i));
  }

  double totalMC0 = 0;
  double totalMCinf = 0;
  double counter_ndf = 0;
  for(int i = FITMINBIN; i < FITMAXBIN; i++)
  {
    totalMC0 = totalMC0 + binContentsMC0[i];
    totalMCinf = totalMCinf + binContentsMCinf[i];
    counter_ndf = counter_ndf + 1;
  }
  ndf = counter_ndf - 2;	// -2 is the 2 parameters in the fit.
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

  // Loading asymmetry for fit
  TString fileName = Form("MB_asymmetries/AsymmFilesFromMB/AllCorr_OctetAsymmetries_AnaChD_Octets_0-59_BinByBin.txt");

  double energy, asymm, asymmErr;

  // opens the file named above
  string buf1;
  ifstream infile1;
  cout << "The file being opened is: " << fileName << endl;
  infile1.open(fileName);

  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << fileName << endl;

  while(true)
  {
    getline(infile1, buf1);
    istringstream bufstream1(buf1);

    if(!infile1.eof())
    {
      bufstream1 >> energy >> asymm >> asymmErr;
      asymmetriesData.push_back(asymm);
      asymmErrorsData.push_back(asymmErr);
      energies.push_back(energy);
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Done loading in the asymmetry values..." << endl;

  double blindingFactor = 0;
  double b_fromBlinding = -0.2/avg_mE;	// this formula I calculated
  for(int i = 0; i <= asymmetriesData.size(); i++)
  {
    blindingFactor = (1.0 + b_fromBlinding*avg_mE) / (1.0 + b_fromBlinding*m_e/(energies[i] + m_e));
    blind_asymmetriesData.push_back(asymmetriesData[i]*blindingFactor);
    blind_asymmErrorsData.push_back(abs(asymmErrorsData[i]*blindingFactor));
  }


  cout << "About to perform TMinuit fit..." << endl;

  TMinuit *gMinuit = new TMinuit(1);
  gMinuit->SetFCN(chi2);

  Double_t arglist[10];
  Int_t iflag = 0;

  arglist[0] = 2.0;	// here arglist[0] = 2 means use strategy 2, aka make sure the errors are right
  gMinuit->mnexcm("SET STR", arglist, 1, iflag);

  arglist[0] = 1.0;	// 1.0 = 1 standard deviation errors. If you want 2 standard deviations, use 4.0.
  gMinuit->mnexcm("SET ERR", arglist, 1, iflag);

  gMinuit->mnparm(0, "fierz", 0, 0.01, -10, 10, iflag);
  gMinuit->mnparm(1, "Asymmetry", 0, 0.01, -10, 10, iflag);

  arglist[0] = 500;
  arglist[1] = 1;

  gMinuit->mnexcm("CALL FCN", arglist, 1, iflag);
  gMinuit->mnexcm("MIGRAD", arglist, 2, iflag);

  cout << "------------ Fit is finished --------------" << endl;

  // can also use GetParameter(...);
  //  double paramValue, paramError;
  //  gMinuit->GetParameter(0, paramValue, paramError);

  Int_t internalInt0, internalInt1;
  Double_t fitVal0, fitErr0, lowLimit0, highLimit0;
  Double_t fitVal1, fitErr1, lowLimit1, highLimit1;
  TString paramName0, paramName1;

  // gets the values of parameter 0
  gMinuit->mnpout(0, paramName0, fitVal0, fitErr0, lowLimit0, highLimit0, internalInt0);
  gMinuit->mnpout(1, paramName1, fitVal1, fitErr1, lowLimit1, highLimit1, internalInt1);


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

  cout << "------- Results printout --------" << endl;
  cout << "OctNb = " << octNb << endl;
  cout << "avg_mE = " << avg_mE << endl;
  cout << "functionMin = " << functionMin << endl;
  cout << "ndf = " << ndf << endl;
  cout << "chisquared/ndf = " << functionMin/ndf << endl;
  cout << "fitVal0 = " << fitVal0 << endl;
  cout << "fitErr0 = " << fitErr0 << endl;
  cout << "fitVal1 = " << fitVal1 << endl;
  cout << "fitErr1 = " << fitErr1 << endl;
  cout << "covMatrixStatus = " << covMatrixStatus << endl;

  ofstream outfile;
  outfile.open(Form("TestingBlinding_UsingSimAnalyzed_CombinedAbFitter_OneOctetModelErrors_OneOctetbAndA_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), ios::app);
  outfile << octNb << "\t"
          << avg_mE << "\t"
          << functionMin << "\t"
          << ndf << "\t"
          << functionMin/ndf << "\t"
          << fitVal0 << "\t"
          << fitErr0 << "\t"
	  << fitVal1 << "\t"
	  << fitErr1 << "\t"
          << covMatrixStatus << "\n";
  outfile.close();



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
  double A = par[1];
  double fitb = 0;
  double fitA = 0;

  double totContentsData = 0;
  for(unsigned int i = FITMINBIN; i < FITMAXBIN; i++)
  {
    totContentsData = totContentsData + binContentsData[i];
  }

  for(unsigned int i = FITMINBIN; i < FITMAXBIN; i++)
  {
    fitb = ((binContentsMC0[i] + b*avg_mE*binContentsMCinf[i])*totContentsData) / (1 + b*avg_mE);

    fitA = (A*(1.0 + b*avg_mE)) / (1.0 + (b*m_e/(energies[i] + m_e)));

    totChi2 = totChi2 + ( pow((binContentsData[i] - fitb) / (binErrorsData[i]), 2.0) + pow((blind_asymmetriesData[i] - fitA) / blind_asymmErrorsData[i], 2.0) );
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

