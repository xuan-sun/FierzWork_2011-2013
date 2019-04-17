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

// Fundamental constants that get used
const double m_e = 511.00;                                              ///< electron mass, keV/c^2
double NDF = -1;

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

// global vectors for creating TGraphs.
vector <double> octets;
vector <double> octetsErr;
vector <double> chisquared;
vector <double> chi2err;
vector <double> bMinuitValues;
vector <double> bErrMinuitValues;

vector <double> x2011;
vector <double> xErr2011;
vector <double> y2011;
vector <double> yErr2011;
vector <double> x2012;
vector <double> xErr2012;
vector <double> y2012;
vector <double> yErr2012;
vector <double> x2013;
vector <double> xErr2013;
vector <double> y2013;
vector <double> yErr2013;


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

  FillArrays("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2011-2012.txt", 1);
  FillArrays("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2012-2013.txt", 2);
  FillArrays("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2012-2013_Sn113Stitch_try3.txt", 3);

//  TGraphErrors *g1 = new TGraphErrors(x2011.size(), &(x2011[0]), &(y2011[0]), &(xErr2011[0]), &(yErr2011[0]));
//  TGraphErrors *g2 = new TGraphErrors(x2012.size(), &(x2012[0]), &(y2012[0]), &(xErr2012[0]), &(yErr2012[0]));

  TGraph *g1 = new TGraph(x2011.size(), &(x2011[0]), &(y2011[0]));
  TGraph *g2 = new TGraph(x2012.size(), &(x2012[0]), &(y2012[0]));
  TGraph *g3 = new TGraph(x2013.size(), &(x2013[0]), &(y2013[0]));

  g1->GetYaxis()->SetRangeUser(1e-8, 1);

  PlotGraph(C, 2, 1, g1, Form("b fit ratios. Energy window upper end is 645 keV."), "Energy Window Low (keV)", "b fit ratios", "AP");
  PlotGraph(C, 4, 1, g2, "", "", "", "PSAME");
  PlotGraph(C, 3, 1, g3, "", "", "", "PSAME");

//  PlotHist(C, 1, 2, h1, "b for all octets", "N", "b", "");


  C->cd(1);
  TLegend* leg1 = new TLegend(0.15,0.7,0.35,0.9);
  leg1->AddEntry(g1,Form("2011-2012"),"p");
  leg1->AddEntry(g2,Form("2012-2013"),"p");
  leg1->AddEntry(g3, Form("2012-2013, Sn Corr"), "p");
  leg1->Draw();


  double xPrint = 45;
  double yPrint = 0.5;

  TLatex t2;
  t2.SetTextSize(0.03);
  t2.SetTextAlign(13);
//  t2.DrawLatex(xPrint, yPrint+0.1, Form("red: %f #pm %f", h1->GetMean(), h1->GetRMS()));
  TLatex t3;
  t3.SetTextSize(0.03);
  t3.SetTextAlign(13);
//  t3.DrawLatex(xPrint, yPrint, Form("red fit: #chi^{2}/ndf = %f", (fit1->GetChisquare() / fit1->GetNDF())));

  TLatex t4;
  t4.SetTextSize(0.03);
  t4.SetTextAlign(13);
//  t4.DrawLatex(xPrint, yPrint-0.1, Form("blue: %f #pm %f", h2->GetMean(), h2->GetRMS()));
  TLatex t5;
  t5.SetTextSize(0.03);
  t5.SetTextAlign(13);
//  t5.DrawLatex(xPrint, yPrint-0.20, Form("blue fit: #chi^{2}/ndf = %f", (fit2->GetChisquare() / fit2->GetNDF())));




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

      if(evt.EMin < 245)
      {
//        continue;
      }

      if(flag == 1)
      {
/*
        octets.push_back(evt.octNb);
        octetsErr.push_back(0.5);
        chisquared.push_back(evt.chisquaredperndf);
        chi2err.push_back(0.1);
        bMinuitValues.push_back(evt.b_minuitFit);
        bErrMinuitValues.push_back(evt.bErr_minuitFit);
*/
        x2011.push_back(evt.EMin);
	xErr2011.push_back(5);
//        y2011.push_back(evt.b_minuitFit / 0.00279833);
        y2011.push_back(evt.prob);
	yErr2011.push_back(evt.bErr_minuitFit / 0.00279833);
      }
      else if(flag == 2)
      {
        x2012.push_back(evt.EMin);
	xErr2012.push_back(5);
//        y2012.push_back(evt.b_minuitFit / 0.0100813);
        y2012.push_back(evt.prob);
	yErr2012.push_back(evt.bErr_minuitFit / 0.0100813);
      }
      else if(flag == 3)
      {
        x2013.push_back(evt.EMin);
	xErr2013.push_back(5);
//        y2013.push_back(evt.b_minuitFit / 0.0100813);
        y2013.push_back(evt.prob);
	yErr2013.push_back(evt.bErr_minuitFit / 0.0100813);
      }
    }


    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

