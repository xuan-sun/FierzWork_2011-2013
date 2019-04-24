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
void FillArrays(TString fileName, int flag);

struct entry
{
  string octNb;
  double avg_mE;
  double chisquared;
  double ndf;
  double chisquaredperndf;
  double prob;
  double b_minuitFit;
  double bErr_minuitFit;
  double binMin;
  double EMin;
  double binMax;
  double EMax;
  int fitMatrixStatus;
};

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

  FillArrays("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2011-2012_highBinVari_run3.txt", 1);
  FillArrays("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2012-2013_highBinVari_run2.txt", 2);
//  FillArrays("allOctets_positionCuts_0-49mm_floatingEndpoint_withFullBlind_Feb2019_type0_2011-2012_run2.txt", 1);
//  FillArrays("allOctets_positionCuts_0-49mm_floatingEndpoint_withFullBlind_Feb2019_type0_2012-2013_run2.txt", 2);

  TGraph *g1 = new TGraph(x[0].size(), &(x[0][0]), &(y[0][0]));
  TGraph *g2 = new TGraph(x[1].size(), &(x[1][0]), &(y[1][0]));
//  TGraph *g3 = new TGraph(x[2].size(), &(x[2][0]), &(y[2][0]));
//  TGraph *g4 = new TGraph(x[3].size(), &(x[3][0]), &(y[3][0]));

  g1->GetYaxis()->SetRangeUser(1e-4, 1);

  PlotGraph(C, 2, 1, g1, Form("Probability of integrated dataset super-sum fit."), "Energy Window High (keV)", "Probability", "AP");
  PlotGraph(C, 4, 1, g2, "", "", "", "PSAME");
//  PlotGraph(C, 3, 1, g3, "", "", "", "PSAME");
//  PlotGraph(C, 6, 1, g4, "", "", "", "PSAME");


  C->cd(1);
  TLegend* leg1 = new TLegend(0.6,0.7,0.9,0.9);
  leg1->AddEntry(g1, Form("2011-2012 corr endpt"),"p");
  leg1->AddEntry(g2, Form("2012-2013 corr endpt"),"p");
//  leg1->AddEntry(g3, Form("2011-2012 free endpt"), "p");
//  leg1->AddEntry(g4, Form("2012-2013 free endpt"), "p");
  leg1->Draw();


  TLine *yLine = new TLine(gPad->GetUxmin(), 0.01, gPad->GetUxmax(), 0.01);
  yLine->SetLineStyle(10);
  yLine->Draw();



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

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString xAxis, TString yAxis, TString command)
{
  C->cd(canvasIndex);
  gPlot->SetTitle(title);
  gPlot->GetXaxis()->SetTitle(xAxis);
  gPlot->GetXaxis()->CenterTitle();
  gPlot->GetYaxis()->SetTitle(yAxis);
  gPlot->GetYaxis()->CenterTitle();
  C->SetLogy();

  gPlot->SetMarkerStyle(21);
  gPlot->SetMarkerSize(0.75);
  gPlot->SetMarkerColor(styleIndex);
  gPlot->SetLineColor(styleIndex);

  gPlot->Draw(command);

  C->Update();

}

void FillArrays(TString fileName, int flag)
{
  vector <double> xTemp;
  vector <double> yTemp;
  vector <double> xErrTemp;
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
      bufstream1 >> evt.octNb
		>> evt.avg_mE
		>> evt.chisquared
		>> evt.ndf
		>> evt.chisquaredperndf
		>> evt.prob
		>> evt.b_minuitFit
		>> evt.bErr_minuitFit
		>> evt.binMin
		>> evt.EMin
		>> evt.binMax
		>> evt.EMax
		>> evt.fitMatrixStatus;
      counter++;

      xTemp.push_back(evt.EMax);
      xErrTemp.push_back(5);
//      yTemp.push_back(evt.b_minuitFit);
      yTemp.push_back(evt.prob);
      yErrTemp.push_back(evt.bErr_minuitFit);
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

