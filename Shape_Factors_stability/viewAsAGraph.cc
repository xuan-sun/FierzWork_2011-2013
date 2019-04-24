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

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command, double maxBinContents);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void FillArrays(TString fileName, TH1D *h, int hFillOption);

vector < vector <double> > x;
vector < vector <double> > xErr;
vector < vector <double> > y;
vector < vector <double> > yErr;

struct entry
{
  int octNb;
  double zeroCrossing;
  double crossingError;
};

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  TCanvas *C = new TCanvas("canvas", "canvas");
//  C->Divide(2,1);
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  TH1D* h1 = new TH1D("NoEndpointCorrection", "zero crossing energy", 100, 200, 400);
  FillArrays(Form("zero_crossing_shapeFactors_noEndpointCorrection_2012-2013.txt"), h1, 1);

  TH1D* h2 = new TH1D("endpointCorrection", "zero crossing energy", 100, 200, 400);
  FillArrays(Form("zero_crossing_shapeFactors_endpointCorrection_2012-2013.txt"), h2, 1);

  int max = h1->GetMaximum();


  PlotHist(C, 2, 1, h1, Form("straight-line fit to shape factor, zero-crossing: 2012-2013"), "fitted zero-crossing (keV)", "N", "", 1.2*max);

  TPaveStats *stats = (TPaveStats*)C->GetPrimitive("stats");
  stats->SetName("h1stats");
  stats->SetY1NDC(0.8);
  stats->SetY2NDC(1.0);
  stats->SetTextColor(2);

  PlotHist(C, 4, 1, h2, "", "", "", "SAMES", 1.2*max);
  TPaveStats *stats2 = (TPaveStats*)C->GetPrimitive("stats");
  stats2->SetName("h1stats2");
  stats2->SetY1NDC(0.5);
  stats2->SetY2NDC(0.7);
  stats2->SetTextColor(4);

  C->Print("hold.pdf");

/*
  TGraphErrors *g1 = new TGraphErrors(x[0].size(), &(x[0][0]), &(y[0][0]), &(xErr[0][0]), &(yErr[0][0]));
  TGraphErrors *g2 = new TGraphErrors(x[1].size(), &(x[1][0]), &(y[1][0]), &(xErr[1][0]), &(yErr[1][0]));

  g1->GetYaxis()->SetRangeUser(0, 600);
  PlotGraph(C, 2, 1, g1, Form("straight-line fit to shape factor, zero-crossing: 2011-2012"), "Octet Number", "fitted zero-crossing (keV)", "AP");
  PlotGraph(C, 4, 1, g2, "", "", "", "PSAME");

  TLegend* leg1 = new TLegend(0.7,0.1,0.9,0.25);
  leg1->AddEntry(g1, Form("no endpt corr"),"p");
  leg1->AddEntry(g2, Form("endpt corrected"),"p");
  leg1->Draw();
*/


  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

void FillArrays(TString fileName, TH1D* h, int hFillOption)
{
  vector <double> xTemp;
  vector <double> yTemp;
  vector <double> xErrTemp;
  vector <double> yErrTemp;


  entry evt;

  //opens the file that I name in DATA_FILE_IN
  string buf;
  ifstream infile1;
  cout << "The file being opened is: " << fileName << endl;
  infile1.open(fileName);

  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << fileName << endl;


  while(true)
  {
    getline(infile1, buf);
    istringstream bufstream(buf);

    if(!infile1.eof())
    {
      bufstream >> evt.octNb
		>> evt.zeroCrossing
		>> evt.crossingError;

      xTemp.push_back(evt.octNb);
      xErrTemp.push_back(0.5);
      yTemp.push_back(evt.zeroCrossing);
      yErrTemp.push_back(evt.crossingError);

      h->Fill(evt.zeroCrossing);
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  x.push_back(xTemp);
  y.push_back(yTemp);
  xTemp.clear();
  yTemp.clear();
  xErr.push_back(xErrTemp);
  yErr.push_back(yErrTemp);
  xErrTemp.clear();
  yErrTemp.clear();

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}




void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command, double maxBinContents)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot->GetXaxis()->SetTitle(xAxis);
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle(yAxis);
  hPlot -> GetYaxis() -> CenterTitle();
  hPlot->GetYaxis()->SetRangeUser(0, 1.2*maxBinContents);

  hPlot->SetFillColor(styleIndex);
  hPlot->SetLineColor(styleIndex);

  if(styleIndex == 2)
  {
    hPlot->SetFillStyle(3004);
  }
  else if(styleIndex == 4)
  {
    hPlot->SetFillStyle(3005);
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

