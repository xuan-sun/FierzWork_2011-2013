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

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH2D *hPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void FillArrays(TString fileName, int flag);

double CalcWeightedAvg(double x1, double w1, double x2, double w2, double x3, double w3, double x4, double w4);
double CalcWeightedError(double x1, double w1, double x2, double w2, double x3, double w3, double x4, double w4);
double CalcChi2(double x1, double w1, double x2, double w2, double x3, double w3, double x4, double w4);

// allOctets aka integrated dataset files
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

// twiddle fit variation files
struct entry2
{
  string year;
  string twiddleIndex;
  double bin;
  double bFit;
  double bErr;
};

// super-ratio dataset fit files
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
  string AErr;
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


  // all our data paths that we use
  // high bin variations, fixed lower bin.
//  TString ssDataPath_2011 = Form("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2011-2012_highBinVari_run2.txt");
//  TString ssDataPath_2012 = Form("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2012-2013_highBinVari_run2.txt");
//  TString ssTwiddlePath_2011 = Form("twiddle_index19_highBinVariation_lowBin-27_positionCuts_0-49mm_endpointCorrected_noBlind_type0_2011-2012_noStatDependence_summary.txt");
//  TString ssTwiddlePath_2012 = Form("twiddle_index19_highBinVariation_lowBin-34_positionCuts_0-49mm_endpointCorrected_noBlind_type0_2012-2013_noStatDependence_summary.txt");

  // low bin variations, fixed higher bin.
  TString ssDataPath_2011 = Form("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2011-2012_lowBinVari_run2.txt");
  TString ssDataPath_2012 = Form("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2012-2013_lowBinVari_run2.txt");
  TString ssTwiddlePath_2011 = Form("twiddle_index19_lowBinVariation_highBin-65_positionCuts_0-49mm_endpointCorrected_noBlind_type0_2011-2012_noStatDependence_summary.txt");
  TString ssTwiddlePath_2012 = Form("twiddle_index19_lowBinVariation_highBin-65_positionCuts_0-49mm_endpointCorrected_noBlind_type0_2012-2013_noStatDependence_summary.txt");

  // super ratio data, not varying fit windows.
  TString srDataPath_2011 = Form("SRDataFit_FullBlind_Feb2019_fixedA_bFitted_2011-2012_190-740keV.txt");
  TString srDataPath_2012 = Form("SRDataFit_FullBlind_Feb2019_fixedA_bFitted_2012-2013_190-740keV.txt");

  // fill in super sum error bars
  FillArrays(ssDataPath_2011, 1);
  FillArrays(ssTwiddlePath_2011, 2);
  // fill in super ratio error bars
  FillArrays(srDataPath_2011, 9);
  // repeat 2012-2013
  FillArrays(ssDataPath_2012, 1);
  FillArrays(ssTwiddlePath_2012, 2);
  FillArrays(srDataPath_2012, 9);


  // fill in super sum fit values
  FillArrays(ssDataPath_2011, 11);
  // fill in super ratio fit values
  FillArrays(srDataPath_2011, 99);
  // repeat 2012-2013
  FillArrays(ssDataPath_2012, 11);
  FillArrays(srDataPath_2012, 99);

  vector <double> SSerr_2011;
  vector <double> SSerr_2012;
  vector <double> SRerr_2011;
  vector <double> SRerr_2012;

  vector <double> SSfit_2011;
  vector <double> SSfit_2012;
  vector <double> SRfit_2011;
  vector <double> SRfit_2012;

  for(unsigned int i = 0; i < x[0].size(); i++)
  {
    SSerr_2011.push_back(sqrt(pow(y[0][i], 2.0) + pow(y[1][i], 2.0)));
  }
  for(unsigned int i = 0; i < x[0].size(); i++)
  {
    SRerr_2011.push_back(y[2][i]);
  }

  for(unsigned int i = 0; i < x[3].size(); i++)
  {
    SSerr_2012.push_back(sqrt(pow(y[3][i], 2.0) + pow(y[4][i], 2.0)));
  }
  for(unsigned int i = 0; i < x[3].size(); i++)
  {
    SRerr_2012.push_back(y[5][i]);
  }

  for(unsigned int i = 0; i < x[6].size(); i++)
  {
    SSfit_2011.push_back(y[6][i]);
  }
  for(unsigned int i = 0; i < x[7].size(); i++)
  {
    SRfit_2011.push_back(y[7][i]);
  }
  for(unsigned int i = 0; i < x[8].size(); i++)
  {
    SSfit_2012.push_back(y[8][i]);
  }
  for(unsigned int i = 0; i < x[9].size(); i++)
  {
    SRfit_2012.push_back(y[9][i]);
  }

  // change the initial energy range of super-ratio from it's currently stored 165-645keV
  // to Michael Brown's 190-740 keV.
/*
  SRfit_2011[0] = -0.0180004;
  SRerr_2011[0] = 0.054527;
  SRfit_2012[0] = -0.0214752;
  SRerr_2012[0] = 0.074534;
*/
  // after running above SR values, the optimum and value of error bar doesn't change. Central value probably does.

  // some sample plotting for visualization
  // low fit window optimization
  TH2D* hist = new TH2D("bOpt", "b optimization using fit windows", 26, 150, 410, 26, 150, 410);
  // high fit window optimization
//  TH2D* hist = new TH2D("bOpt", "b optimization using fit windows", 16, 590, 750, 16, 590, 750);
  hist->SetContour(500);
  double chi2perndf = 0;
  double fillValueWeightedError = 0;

  for(unsigned int a = 0; a < SSfit_2011.size(); a++)
  {
    for(unsigned int c = 0; c < SSfit_2012.size(); c++)
    {
      // take only super ratio fit values for full energy range. Scan over the SS fit values.
      chi2perndf = CalcChi2(SSfit_2011[a], 1.0/pow(SSerr_2011[a], 2.0),
                                SRfit_2011[0], 1.0/pow(SRerr_2011[0], 2.0),
                                SSfit_2012[c], 1.0/pow(SSerr_2012[c], 2.0),
                                SRfit_2012[0], 1.0/pow(SRerr_2012[0], 2.0) );

      fillValueWeightedError = CalcWeightedError( SSfit_2011[a], 1.0/pow(SSerr_2011[a], 2.0),
                                SRfit_2011[0], 1.0/pow(SRerr_2011[0], 2.0),
                                SSfit_2012[c], 1.0/pow(SSerr_2012[c], 2.0),
                                SRfit_2012[0], 1.0/pow(SRerr_2012[0], 2.0) );

      if(hist->GetXaxis()->GetBinCenter(a+2) == 255 && hist->GetYaxis()->GetBinCenter(c+2) == 335)
      {
        cout << "fillValueWeightedError = " << fillValueWeightedError
	     << ", chi2/ndf = " << chi2perndf
	     << ". 2011-2012 energy is " << x[0][a]
	     << ", 2012-2013 energy is " << x[0][c] << endl;
      }

//      if(chi2perndf <= 1.0)
//      {
        hist->SetBinContent(hist->GetXaxis()->FindBin(x[0][a]), hist->GetYaxis()->FindBin(x[0][c]), fillValueWeightedError);
//      }
//      else if(chi2perndf > 1.0)
//      {
//        hist->SetBinContent(hist->GetXaxis()->FindBin(x[0][a]), hist->GetYaxis()->FindBin(x[0][c]), fillValueWeightedError * sqrt(chi2perndf) );
//      }
    }
  }



  C->SetRightMargin(0.20);
  gStyle->SetOptStat(0);
  hist->SetTitle("4-point weighted average b error vs fit windows on SS.");
  hist->GetXaxis()->SetTitle("2011-2012 Super-sum, high fit window cut (keV)");
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->SetTitle("2012-2013 Super-sum, high fit window cut (keV)");
  hist->GetYaxis()->CenterTitle();
  hist->GetZaxis()->SetTitle("Magnitude of 4-point weighted average error");
  TGaxis::SetMaxDigits(3);
  hist->Draw("COLZ");


  //prints the canvas with a dynamic TString name of the name of the file
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH2D *hPlot, TString title, TString xAxis, TString yAxis, TString command)
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
      else if(flag == 99)
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

double CalcWeightedAvg(double x1, double w1, double x2, double w2, double x3, double w3, double x4, double w4)
{
  return (x1*w1 + x2*w2 + x3*w3 + x4*w4) / (w1 + w2 + w3 + w4);
}


double CalcWeightedError(double x1, double w1, double x2, double w2, double x3, double w3, double x4, double w4)
{
  return sqrt(1.0 / (w1 + w2 + w3 + w4));
}

double CalcChi2(double x1, double w1, double x2, double w2, double x3, double w3, double x4, double w4)
{
  double xMean = CalcWeightedAvg(x1, w1, x2, w2, x3, w3, x4, w4);

  return (w1*pow(x1-xMean, 2.0) + w2*pow(x2-xMean, 2.0) + w3*pow(x3-xMean, 2.0) + w4*pow(x4-xMean, 2.0)) / 3.0;
}
