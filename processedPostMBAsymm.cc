#include	 "BetaSpectrum.hh"

#include	 <iostream>
#include	 <fstream>
#include	 <TGaxis.h>
#include	 <sstream>
#include	 <TGraph.h>
#include	 <TGraphErrors.h>
#include	 <TCanvas.h>
#include	 <TApplication.h>
#include	 <stdlib.h>
#include	 <TF1.h>
#include	 <TH1.h>
#include	 <TProfile.h>
#include	 <TObjArray.h>
#include	 <TStyle.h>
#include	 <TMarker.h>
#include	 <math.h>
#include	 <TStyle.h>
#include	 <TPaveStats.h>
#include	 <TPaveText.h>
#include	 <vector>
#include	 <string.h>
#include	 <fstream>
#include	 <TROOT.h>
#include	 <TFile.h>
#include	 <TLegend.h>
#include         <TLegendEntry.h>
#include	 <time.h>
#include	 <TH2F.h>
#include         <assert.h>
#include	 <string>
#include	 <TRandom.h>
#include 	 <TTree.h>
#include	 <TChain.h>
#include	 <TVector.h>
#include	 <vector>
#include	 <utility>
#include	 <TLeaf.h>
#include	 <TLine.h>
#include	 <TLatex.h>
using		 namespace std;

struct Event
{
  double primKE;
  double Erecon;
  int side;
  int type;
  int eventNum;
  double time[2];
  double tE;
  double tW;
  int pid;
  double xEastPos;
  double yEastPos;
  double xWestPos;
  double yWestPos;
};

// forward declarations for useful functions
TH1D* CalculateAofE(TH1D* R);
TH1D* DivideByMBAsymm(TH1D* AofE);
double CalculatebFromPercentageMixing(TString fileName);
double CalculateAveragemOverE(TH1D* gammaSM, int binMin, int binMax);

// plotting functions
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xTitle, TString yTitle, TString command);

// these are actual beta run indices
const int index_A2 = 0;
const int index_A5 = 1;
const int index_A7 = 2;
const int index_A10 = 3;
const int index_B2 = 4;
const int index_B5 = 5;
const int index_B7 = 6;
const int index_B10 = 7;
double avg_mE = 0;

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

//-------------------------------------------------//
//------------ Start of Program -------------------//
//-------------------------------------------------//

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  // creating canvas for plotting
  TCanvas *C = new TCanvas("canvas", "canvas", 800, 400);

  int octNb = 43;
  double xMin = 165;
  double xMax = 645;

//  TFile fMC0(TString::Format("/mnt/Data/xuansun/BLIND_MC_files/%s_geom/BLIND_MC_A_0_b_0_Octet_%i_ssHist_%s.root", "2011-2012", octNb, "type0"));
  TFile fMC0(TString::Format("ExtractedHistograms/MC_A_0_b_0/MC_A_0_b_0_Octet_%i_ssHist_%s.root", octNb, "type0"));
  TH1D* mcTheoryHistBeta = (TH1D*)fMC0.Get("Super sum");
  avg_mE = CalculateAveragemOverE(mcTheoryHistBeta, mcTheoryHistBeta->FindBin(xMin), mcTheoryHistBeta->FindBin(xMax));

  cout << "Average m/E for octet " << octNb << " is equal to " << avg_mE
 	<< " over a fit range of " << xMin << " to " << xMax << endl;

  double bMixing = CalculatebFromPercentageMixing("ExtractedHistograms/randomMixingSeeds.txt");

  cout << "For octet " << octNb << ", using random seed s0, we get a b value of " << bMixing << endl;


  // make a file and write the output of the calculations to it
//  TH1D* asymm_Erecon = CalculateAofE(superRatio_Erecon);




/*
  TF1 *fit = new TF1("beta fit", "[0]*sqrt(x*x +2*511*x) / (511 + x)", xMin, xMax);

  fit->SetParName(0, "normalization");
  asymm_Erecon->Fit("beta fit");
  TF1 *fitResults = asymm_Erecon->GetFunction("beta fit");
  cout << "Chi squared value is " << fitResults->GetChisquare()
	<< " with ndf of " << fitResults->GetNDF()
	<< ". For a final chisquared/ndf = " << fitResults->GetChisquare() / fitResults->GetNDF() << endl;
*/


//  PlotHist(C, 1, 1, asymm_Erecon, "", "Primary Energy (keV)", "Asymmetry", "");

  TLine *yLow = new TLine(xMin, 0.02, xMin, 0.08);
  TLine *yHigh = new TLine(xMax, 0.02, xMax, 0.08);
  yLow->SetLineWidth(2);
  yHigh->SetLineWidth(2);
  yLow->SetLineColor(1);
  yHigh->SetLineColor(1);
  yLow->Draw("SAME");
  yHigh->Draw("SAME");


  // Save our plot and print it out as a pdf.
//  C -> Print("100xStats_statsWeighted_SR_asymm.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

TH1D* CalculateAofE(TH1D* R)
{
  TH1D* hAofE = new TH1D("AofE", "AofE", 120, 0, 1200);

  double asymm = 0;
  double asymmErr = 0;

  for(int i = 0; i <= R->GetNbinsX(); i++)
  {
    if(R->GetBinContent(i) == 0)
    {
      asymm = 0;
      asymmErr = 0;
    }
    else
    {
      asymm = (1 - sqrt(R->GetBinContent(i)) ) / ( 1 + sqrt(R->GetBinContent(i)) );
      asymmErr = (R->GetBinError(i))
		 / ( sqrt(R->GetBinContent(i)) * (sqrt(R->GetBinContent(i)) + 1) * (sqrt(R->GetBinContent(i)) + 1) );
    }

    hAofE->SetBinContent(i, asymm);
    hAofE->SetBinError(i, asymmErr);
  }

  return hAofE;
}

TH1D* DivideByMBAsymm(TH1D* AofE)
{
  TH1D* hCorrAofE = new TH1D("Corrected AofE", "Corrected AofE", 120, 0, 1200);

  double energy, asymm, asymmErr;
  double bin, binCounts, binError;

  TString fileName = Form("UnCorr_OctetAsymmetries_AnaChA_180-780_Octets_0-59_BinByBin_withEnergyDependence.txt");
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

      bin = AofE->FindBin(energy);
      binCounts = AofE->GetBinContent(bin);
      binError = AofE->GetBinError(bin);

      if(asymm == 0)
      {
        hCorrAofE->SetBinContent(bin, 0);
        hCorrAofE->SetBinError(bin, binError);
      }
      else
      {
        hCorrAofE->SetBinContent(bin, binCounts / abs(asymm));
        hCorrAofE->SetBinError(bin, sqrt(binError*binError + asymmErr*asymmErr));
      }
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }


  return hCorrAofE;
}


void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xTitle, TString yTitle, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle(xTitle);
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle(yTitle);
  hPlot -> GetYaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetRangeUser(0.02, 0.08);


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

double CalculatebFromPercentageMixing(TString fileName)
{
  double b = 0;

  // read in our random mixing seed so I stay pretty blind.
  double s0, s1, s2, s3, s4, s5, s6, s7, s8, s9;

  string buf1;
  ifstream infile1;
  cout << "The file being opened is: " << fileName.Data() << endl;
  infile1.open(fileName);

  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << fileName.Data() << endl;

  while(true)
  {
    getline(infile1, buf1);
    istringstream bufstream1(buf1);

    if(!infile1.eof())
    {
      bufstream1 >> s0 >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8 >> s9;
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  // load all the histograms of east and west, turn them into rates.
  // ALWAYS USE S3 FOR MIXING.
//  vector < vector < TH1D* > > rates_base = CreateRateHistograms(runFiles_base, 1 - s3);
//  vector < vector < TH1D* > > rates_fierz = CreateRateHistograms(runFiles_fierz, s3);

  b = s4 / ( (1 - s4) * (avg_mE) );

  return b;
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


