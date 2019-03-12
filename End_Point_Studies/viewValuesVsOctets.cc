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

#define		TYPE	"type0"
#define		GEOM	"2011-2012"
#define		FITMINBIN	17
#define		FITMAXBIN	65

using            namespace std;

// Fundamental constants that get used
const double m_e = 511.00;                                              ///< electron mass, keV/c^2
double NDF = -1;

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void FillArrays(TString fileName, TH1D *hist1, int flag);
void ConvertOctetsToDate(int octNb);

struct entry
{
  int octNb;
  double avg_mE;
  double chisquared;
  double ndf;
  double chisquaredperndf;
  double b_minuitFit;
  double bErr_minuitFit;
  double A_minuitFit;
  double AErr_minuitFit;
  int fitMatrixStatus;
};

// global vectors for creating TGraphs.
vector < vector <double> > octets;
vector < vector <double> > octetsErr;
vector < vector <double> > chisquared;
vector < vector <double> > bMinuitValues;
vector < vector <double> > bErrMinuitValues;

/*
vector <double> octets2;
vector <double> octetsErr2;
vector <double> chisquared2;
vector <double> bMinuitValues2;
vector <double> bErrMinuitValues2;
*/

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  NDF = FITMAXBIN - FITMINBIN - 1;

  TCanvas *C = new TCanvas("canvas", "canvas");
//  C -> Divide(2,1);
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  TH1D *h1 = new TH1D("fierz minuit 2011-2012", "fierz 2011-2012", 100, -1, 1);
  TH1D *h2 = new TH1D("position cut fierz", "fierz 2011-2013", 100, -1, 1);
//  h1->SetStats(0);

//  FillArrays(Form("../NewXuanFitter/FullBlindFeb2019_newXuanFitter_dataHists_bFit_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), h1);
  FillArrays(Form("positionCuts_0-150mm_withBlind_andMCCuts_newXuanFitter_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), h2, 0);
//  FillArrays(Form("positionCuts_0-150mm_withBlind_andMCCuts_newXuanFitter_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), h2);
//  FillArrays(Form("positionCuts_0-150mm_withBlind_andMCCuts_newXuanFitter_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), h2);

  vector <double> chisquaredError(chisquared.size(), 0.01);

//  TGraphErrors *g1 = new TGraphErrors(octets.size(), &(octets[0]), &(bMinuitValues[0]), &(octetsErr[0]), &(bErrMinuitValues[0]));
//  TGraphErrors *g2 = new TGraphErrors(octets2.size(), &(octets2[0]), &(bMinuitValues2[0]), &(octetsErr2[0]), &(bErrMinuitValues2[0]));

//  g1->GetYaxis()->SetRangeUser(-0.5, 6);

//  PlotGraph(C, 2, 1, g1, Form("b for %s: 49-150mm radius", GEOM), "Octet Number", "b", "AP");
//  PlotGraph(C, 4, 1, g2, Form("b for %s: 49-150mm radius", GEOM), "Octet Number", "b", "PSAME");

//  PlotHist(C, 1, 2, h1, "b for all octets", "N", "b", "");

  C->cd(1);
  TLegend* leg1 = new TLegend(0.6,0.6,0.9,0.8);
  leg1->AddEntry(h1,"b data","p");
  leg1->AddEntry(h2,"b endpoint corr","p");
//  leg1->Draw();


  double xPrint = 45;
  double yPrint = 0.1;

/*
  TLatex t2;
  t2.SetTextSize(0.03);
  t2.SetTextAlign(13);
  t2.DrawLatex(xPrint, yPrint+0.15, Form("red #mu = %f", h1->GetMean()));
  TLatex t3;
  t3.SetTextSize(0.03);
  t3.SetTextAlign(13);
  t3.DrawLatex(xPrint, yPrint, Form("red RMS = %f", h1->GetRMS()));

  TLatex t4;
  t4.SetTextSize(0.03);
  t4.SetTextAlign(13);
  t4.DrawLatex(xPrint, yPrint-0.15, Form("blue #mu = %f", h2->GetMean()));
  TLatex t5;
  t5.SetTextSize(0.03);
  t5.SetTextAlign(13);
  t5.DrawLatex(xPrint, yPrint-0.30, Form("blue RMS = %f", h2->GetRMS()));
*/



  //prints the canvas with a dynamic TString name of the name of the file
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

  hPlot->Draw(command);

  C->Update();
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command)
{
  C->cd(canvasIndex);
  gPlot->SetTitle(title);
  gPlot->GetXaxis()->SetTitle(xAxis);
  gPlot->GetXaxis()->CenterTitle();
  gPlot->GetYaxis()->SetTitle(yAxis);
  gPlot->GetYaxis()->CenterTitle();

  gPlot->SetMarkerStyle(21);
  gPlot->SetMarkerSize(0.75);
  gPlot->SetMarkerColor(styleIndex);

  gPlot->Draw(command);

  C->Update();
}

void FillArrays(TString fileName, TH1D* hist1, int flag)
{

  entry evt;
  int counter = 0;

  //opens the file that I name in DATA_FILE_IN
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
      bufstream1 >> evt.octNb
		>> evt.avg_mE
		>> evt.chisquared
		>> evt.ndf
		>> evt.chisquaredperndf
		>> evt.b_minuitFit
		>> evt.bErr_minuitFit
		>> evt.A_minuitFit
		>> evt.AErr_minuitFit
		>> evt.fitMatrixStatus;
      counter++;

      octets[flag].push_back(evt.octNb);
      octetsErr[flag].push_back(0.5);
      chisquared[flag].push_back(evt.chisquaredperndf);
      bMinuitValues[flag].push_back(evt.b_minuitFit);
      bErrMinuitValues[flag].push_back(evt.bErr_minuitFit);

    }


    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

