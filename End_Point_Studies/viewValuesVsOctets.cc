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
#define		GEOM	"2012-2013"
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
void FillArrays(TString fileName);
void ConvertOctetsToDate(int octNb);

struct entry
{
  int octNb;
  double endpoint;
  double endpointErr;
  double gainFactor;
};

// global vectors for creating TGraphs.
vector < vector <double> > octets;
vector < vector <double> > octetsErr;
vector < vector <double> > endpoints;
vector < vector <double> > endpointsErr;
vector < vector <double> > gainFactors;


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

  FillArrays(Form("endPointFits_noGainCorrection_ssDataHists_%s_radialCut_%i-%imm.txt", GEOM, 0, 30));
  FillArrays(Form("endPointFits_noGainCorrection_ssDataHists_%s_radialCut_%i-%imm.txt", GEOM, 30, 49));
  FillArrays(Form("endPointFits_noGainCorrection_ssDataHists_%s_radialCut_%i-%imm.txt", GEOM, 0, 49));
  FillArrays(Form("endPointFits_noGainCorrection_ssDataHists_%s_radialCut_%i-%imm.txt", GEOM, 49, 150));
  FillArrays(Form("endPointFits_noGainCorrection_ssDataHists_%s_radialCut_%i-%imm.txt", GEOM, 0, 150));

  TGraphErrors *g1 = new TGraphErrors(octets[0].size(), &(octets[0][0]), &(endpoints[0][0]), &(octetsErr[0][0]), &(endpointsErr[0][0]));
  TGraphErrors *g2 = new TGraphErrors(octets[1].size(), &(octets[1][0]), &(endpoints[1][0]), &(octetsErr[1][0]), &(endpointsErr[1][0]));
  TGraphErrors *g3 = new TGraphErrors(octets[2].size(), &(octets[2][0]), &(endpoints[2][0]), &(octetsErr[2][0]), &(endpointsErr[2][0]));
  TGraphErrors *g4 = new TGraphErrors(octets[3].size(), &(octets[3][0]), &(endpoints[3][0]), &(octetsErr[3][0]), &(endpointsErr[3][0]));
  TGraphErrors *g5 = new TGraphErrors(octets[4].size(), &(octets[4][0]), &(endpoints[4][0]), &(octetsErr[4][0]), &(endpointsErr[4][0]));

  g1->GetYaxis()->SetRangeUser(710, 810);

  PlotGraph(C, 1, 1, g1, Form("endpoints for %s with radial cuts", GEOM), "Octet Number", "endpoints (keV)", "AP");
  PlotGraph(C, 2, 1, g2, "", "", "", "PSAME");
  PlotGraph(C, 3, 1, g3, "", "", "", "PSAME");
  PlotGraph(C, 4, 1, g4, "", "", "", "PSAME");
  PlotGraph(C, 6, 1, g5, "", "", "", "PSAME");

  C->cd(1);
  TLegend* leg1 = new TLegend(0.7,0.1,0.9,0.3);
  leg1->AddEntry(g1,"0-30mm","p");
  leg1->AddEntry(g2,"30-49mm","p");
  leg1->AddEntry(g3,"0-49mm","p");
  leg1->AddEntry(g4,"49-150mm","p");
  leg1->AddEntry(g5,"0-150mm","p");
  leg1->Draw();

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

void FillArrays(TString fileName)
{
  vector <double> octetsTemp;
  vector <double> octetsErrTemp;
  vector <double> endpointsTemp;
  vector <double> endpointsErrTemp;
  vector <double> gainFactorsTemp;

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
		>> evt.endpoint
		>> evt.endpointErr
		>> evt.gainFactor;
      counter++;

      octetsTemp.push_back(evt.octNb);
      octetsErrTemp.push_back(0.5);
      endpointsTemp.push_back(evt.endpoint);
      endpointsErrTemp.push_back(evt.endpointErr);
      gainFactorsTemp.push_back(evt.gainFactor);
    }


    if(infile1.eof() == true)
    {
      break;
    }
  }

  octets.push_back(octetsTemp);
  octetsErr.push_back(octetsErrTemp);
  endpoints.push_back(endpointsTemp);
  endpointsErr.push_back(endpointsErrTemp);
  gainFactors.push_back(gainFactorsTemp);

  octetsTemp.clear();
  octetsErrTemp.clear();
  endpointsTemp.clear();
  endpointsErrTemp.clear();
  gainFactorsTemp.clear();

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

