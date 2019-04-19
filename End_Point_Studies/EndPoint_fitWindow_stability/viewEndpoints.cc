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

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void FillArrays(TString fileName);
void ConvertOctetsToDate(int octNb);

struct entry
{
  string year;
  string octets;
  double bin;
  double endpoint;
  double endpointRMS;
};

// global vectors for creating TGraphs.
vector < vector <double> > x;
vector < vector <double> > xErr;
vector < vector <double> > y;
vector < vector <double> > yErr;


int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  TCanvas *C = new TCanvas("canvas", "canvas");
//  C -> Divide(2,1);
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!


  FillArrays(Form("endPointFits_noCorrection_ssDataHists_type0_2011-2012_radialCut_0-49mm_Bins_summary.txt"));
  FillArrays(Form("endPointFits_noCorrection_ssDataHists_type0_2012-2013_radialCut_0-49mm_Bins_summary.txt"));


  TGraphErrors *g1 = new TGraphErrors(x[0].size(), &(x[0][0]), &(y[0][0]), &(xErr[0][0]), &(yErr[0][0]));
  TGraphErrors *g2 = new TGraphErrors(x[1].size(), &(x[1][0]), &(y[1][0]), &(xErr[1][0]), &(yErr[1][0]));

  g1->GetYaxis()->SetRangeUser(770, 800);


  PlotGraph(C, 1, 1, g1, Form("endpoints vs fit binwdows"), "Low window bin number", "endpoints (keV)", "AP");
  PlotGraph(C, 2, 1, g2, "", "", "", "PSAME");

  C->cd(1);
  TLegend* leg1 = new TLegend(0.7,0.1,0.9,0.3);
  leg1->AddEntry(g1,"2011-2012, 0<r<49mm","p");
  leg1->AddEntry(g2,"2012-2013, 0<r<49mm","p");
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
  gPlot->SetLineColor(styleIndex);

  gPlot->Draw(command);

  C->Update();
}

void FillArrays(TString fileName)
{
  vector <double> xTemp;
  vector <double> xErrTemp;
  vector <double> yTemp;
  vector <double> yErrTemp;

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
      bufstream1 >> evt.year
		>> evt.octets
		>> evt.bin
		>> evt.endpoint
		>> evt.endpointRMS;
      counter++;

      xTemp.push_back(evt.bin);
      xErrTemp.push_back(0.5);
      yTemp.push_back(evt.endpoint);
      yErrTemp.push_back(evt.endpointRMS);
    }


    if(infile1.eof() == true)
    {
      break;
    }
  }

  x.push_back(xTemp);
  xErr.push_back(xErrTemp);
  y.push_back(yTemp);
  yErr.push_back(yErrTemp);

  xTemp.clear();
  xErrTemp.clear();
  yTemp.clear();
  yErrTemp.clear();

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

