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

struct entry1
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

struct entry2
{
  string year;
  string twiddleIndex;
  double bin;
  double bFit;
  double bErr;
};

struct entry3
{
  string year;
  string octets;
  double bin;
  double bMean;
  double bRMS;
  double fitb;
  double fitbErr;
  double fitChi2;
  double fitNDF;
  double fitchi2perndf;
  double prob;
};

struct entry4
{
  int octForRef;
  double avg_mE;
  double Emin;
  double Emax;
  double chi2;
  double ndf;
  double chi2perndf;
  double A;
  double AErr;
  double b;
  double bErr;
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

  FillArrays("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2012-2013.txt", 1);
  FillArrays("twiddle_index19_binVariation_positionCuts_0-49mm_endpointCorrected_noBlind_type0_2012-2013_noStatDependence_summary.txt", 2);
  FillArrays("positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2012-2013_binWindowVariations_individualOctetsSummary.txt", 3);
  FillArrays("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2012-2013.txt", 4);
  FillArrays("positionCuts_0-49mm_noGainCorrection_withFullBlind_Feb2019_type0_2012-2013_binWindowVariations_individualOctetsSummary.txt", 5);

  // make some of the "math operations" error bars
  vector <double> xTemp;
  vector <double> yTemp;
  for(unsigned int i = 0; i < x[0].size(); i++)
  {
    xTemp.push_back(x[0][i]);	// pushes back energy values
    yTemp.push_back(sqrt(y[0][i]*y[0][i] + y[1][i]*y[1][i]));	// sqrt(stat^2 + twiddle^2)
  }
  x.push_back(xTemp);
  y.push_back(yTemp);
  xTemp.clear();
  yTemp.clear();

  for(unsigned int i = 0; i < x[0].size(); i++)
  {
    xTemp.push_back(x[0][i]);   // pushes back energy values
    yTemp.push_back(sqrt(y[3][i])*y[0][i]);	// sqrt(chi2(from all octets))*stat
  }
  x.push_back(xTemp);
  y.push_back(yTemp);
  xTemp.clear();
  yTemp.clear();

  for(unsigned int i = 0; i < x[0].size(); i++)
  {
    xTemp.push_back(x[0][i]);   // pushes back energy values
    yTemp.push_back(y[0][i]*sqrt(37));     // stat*sqrt(58 /*37*/) to represent octet RMS 2011-2012 /*2012-2013*/
  }
  x.push_back(xTemp);
  y.push_back(yTemp);
  xTemp.clear();
  yTemp.clear();

  // fill in super ratio error bars
  FillArrays("AsymmetryDataFit_FullBlind_Feb2019_AbParams_2012-2013_fitWindowSummary.txt", 9);
//  FillArrays("AsymmetryDataFit_FullBlind_Feb2019_justbParam_2012-2013_fitWindowSummary.txt", 9);

  // fill in b values for super ratio and super sum
  FillArrays("AsymmetryDataFit_FullBlind_Feb2019_AbParams_2012-2013_fitWindowSummary.txt", 10);
//  FillArrays("AsymmetryDataFit_FullBlind_Feb2019_justbParam_2012-2013_fitWindowSummary.txt", 10);
  FillArrays("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2012-2013.txt", 11);

  for(unsigned int i = 0; i < x[0].size(); i++)
  {
    xTemp.push_back(x[0][i]);   // pushes back energy values
    yTemp.push_back(abs(y[9][i] - y[10][i]));     // abs(b_SR - b_SS)
  }
  x.push_back(xTemp);
  y.push_back(yTemp);
  xTemp.clear();
  yTemp.clear();


  cout << "We have a total of " << x.size() << " arrays for potential graphs." << endl;

  TGraph *g1 = new TGraph(x[0].size(), &(x[0][0]), &(y[0][0]));	// integrated dataset stat err
  TGraph *g2 = new TGraph(x[1].size(), &(x[1][0]), &(y[1][0]));	// twiddle error
  TGraph *g3 = new TGraph(x[2].size(), &(x[2][0]), &(y[2][0]));	// octet RMS error
  TGraph *g4 = new TGraph(x[5].size(), &(x[5][0]), &(y[5][0]));	// stat^2 + twiddle^2
  TGraph *g5 = new TGraph(x[6].size(), &(x[6][0]), &(y[6][0]));	// sqrt(chi2/ndf) * stat
  TGraph *g6 = new TGraph(x[7].size(), &(x[7][0]), &(y[7][0]));	// stat * sqrt(58) to represent octet RMS
  TGraph *g7 = new TGraph(x[4].size(), &(x[4][0]), &(y[4][0]));	// octet RMS (not endpoint-corrected)
  TGraph *g8 = new TGraph(x[8].size(), &(x[8][0]), &(y[8][0]));	// super ratio fit error
  TGraph *g9 = new TGraph(x[11].size(), &(x[11][0]), &(y[11][0]));	// b value offset |SR - SS|

  g3->GetYaxis()->SetRangeUser(0, 0.25);

  PlotGraph(C, 3, 1, g3, Form("b fit errors. Energy window upper end is 645 keV. Geom is 2012-2013."), "Energy Window Low (keV)", "b fit errors", "AP");
  PlotGraph(C, 4, 1, g2, "", "", "", "PSAME");
  PlotGraph(C, 2, 1, g1, "", "", "", "PSAME");
  PlotGraph(C, 6, 1, g4, "", "", "", "PSAME");
  PlotGraph(C, 46, 1, g5, "", "", "", "PSAME");
  PlotGraph(C, 8, 1, g6, "", "", "", "PSAME");
  PlotGraph(C, 32, 1, g7, "", "", "", "PSAME");
  PlotGraph(C, 7, 1, g8, "", "", "", "PSAME");
  PlotGraph(C, 1, 1, g9, "", "", "", "PSAME");

  C->cd(1);
  TLegend* leg1 = new TLegend(0.1,0.6,0.55,0.9);
  leg1->AddEntry(g1, Form("stat error, integrated dataset"),"p");
  leg1->AddEntry(g2, Form("syst error, twiddles (stat dep)"),"p");
  leg1->AddEntry(g3, Form("octet RMS"), "p");
  leg1->AddEntry(g4, Form("total sqrt(stat^2 + syst^2)"), "p");
  leg1->AddEntry(g5, Form("stat*sqrt(Chi2/ndf)"), "p");
  leg1->AddEntry(g6, Form("stat*sqrt(N_{octets})"), "p");
  leg1->AddEntry(g7, Form("octet RMS (not endpt-corrected)"), "p");
  leg1->AddEntry(g8, Form("super-ratio fit err"), "p");
  leg1->AddEntry(g9, Form("|b_{SR} - b_{SS}|"), "p");
//  leg1->Draw();


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
//  C->SetLogy();

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

  entry1 evt1;
  entry2 evt2;
  entry3 evt3;
  entry4 evt4;

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
      if(flag == 1)
      {
        bufstream1 >> evt1.octNb
		>> evt1.avg_mE
		>> evt1.chisquared
		>> evt1.ndf
		>> evt1.chisquaredperndf
		>> evt1.prob
		>> evt1.b_minuitFit
		>> evt1.bErr_minuitFit
		>> evt1.binMin
		>> evt1.EMin
		>> evt1.binMax
		>> evt1.EMax
		>> evt1.fitMatrixStatus;

        xTemp.push_back(evt1.EMin);
	yTemp.push_back(evt1.bErr_minuitFit);
      }
      else if(flag == 2)
      {
        bufstream1 >> evt2.year
                >> evt2.twiddleIndex
                >> evt2.bin
                >> evt2.bFit
                >> evt2.bErr;

	xTemp.push_back((evt2.bin)*10 - 5);
	yTemp.push_back(evt2.bErr);
      }
      else if(flag == 3)
      {
        bufstream1 >> evt3.year
		>> evt3.octets
		>> evt3.bin
		>> evt3.bMean
		>> evt3.bRMS
		>> evt3.fitb
		>> evt3.fitbErr
		>> evt3.fitChi2
		>> evt3.fitNDF
		>> evt3.fitchi2perndf
		>> evt3.prob;

        xTemp.push_back((evt3.bin)*10 - 5);
	yTemp.push_back(evt3.bRMS);
      }
      else if(flag == 4)
      {
        bufstream1 >> evt1.octNb
                >> evt1.avg_mE
                >> evt1.chisquared
                >> evt1.ndf
                >> evt1.chisquaredperndf
                >> evt1.prob
                >> evt1.b_minuitFit
                >> evt1.bErr_minuitFit
                >> evt1.binMin
                >> evt1.EMin
                >> evt1.binMax
                >> evt1.EMax
                >> evt1.fitMatrixStatus;

        xTemp.push_back(evt1.EMin);
        yTemp.push_back(evt1.chisquaredperndf);
      }
      else if(flag == 5)
      {
        bufstream1 >> evt3.year
                >> evt3.octets
                >> evt3.bin
                >> evt3.bMean
                >> evt3.bRMS
                >> evt3.fitb
                >> evt3.fitbErr
                >> evt3.fitChi2
                >> evt3.fitNDF
                >> evt3.fitchi2perndf
                >> evt3.prob;

        xTemp.push_back((evt3.bin)*10 - 5);
        yTemp.push_back(evt3.bRMS);
      }
      else if(flag == 9)
      {
        bufstream1 >> evt4.octForRef
                >> evt4.avg_mE
                >> evt4.Emin
                >> evt4.Emax
                >> evt4.chi2
                >> evt4.ndf
                >> evt4.chi2perndf
                >> evt4.A
                >> evt4.AErr
                >> evt4.b
                >> evt4.bErr;

        xTemp.push_back(evt4.Emin);
        yTemp.push_back(evt4.bErr);
      }
      else if(flag == 10)
      {
        bufstream1 >> evt4.octForRef
                >> evt4.avg_mE
                >> evt4.Emin
                >> evt4.Emax
                >> evt4.chi2
                >> evt4.ndf
                >> evt4.chi2perndf
                >> evt4.A
                >> evt4.AErr
                >> evt4.b
                >> evt4.bErr;

        xTemp.push_back(evt4.Emin);
        yTemp.push_back(evt4.b);
      }
      else if(flag == 11)
      {
        bufstream1 >> evt1.octNb
                >> evt1.avg_mE
                >> evt1.chisquared
                >> evt1.ndf
                >> evt1.chisquaredperndf
                >> evt1.prob
                >> evt1.b_minuitFit
                >> evt1.bErr_minuitFit
                >> evt1.binMin
                >> evt1.EMin
                >> evt1.binMax
                >> evt1.EMax
                >> evt1.fitMatrixStatus;

        xTemp.push_back(evt1.EMin);
        yTemp.push_back(evt1.b_minuitFit);
      }




      counter++;
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

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

