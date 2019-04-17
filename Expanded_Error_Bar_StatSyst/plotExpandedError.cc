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

  FillArrays("allOctets_positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2011-2012.txt", 1);
  FillArrays("twiddle_index19_binVariation_positionCuts_0-49mm_endpointCorrected_noBlind_type0_2011-2012_summary.txt", 2);
  FillArrays("positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_2011-2012_binWindowVariations_individualOctetsSummary.txt", 3);

  vector <double> xTemp;
  vector <double> yTemp;
  for(unsigned int i = 0; i < x[0].size(); i++)
  {
    xTemp.push_back(x[0][i]);	// pushes back energy values
    yTemp.push_back(sqrt(y[0][i]*y[0][i] + y[1][i]*y[1][i]));
  }

  x.push_back(xTemp);
  y.push_back(yTemp);
  xTemp.clear();
  yTemp.clear();

  TGraph *g1 = new TGraph(x[0].size(), &(x[0][0]), &(y[0][0]));
  TGraph *g2 = new TGraph(x[1].size(), &(x[1][0]), &(y[1][0]));
  TGraph *g3 = new TGraph(x[2].size(), &(x[2][0]), &(y[2][0]));
  TGraph *g4 = new TGraph(x[3].size(), &(x[3][0]), &(y[3][0]));


  g1->GetYaxis()->SetRangeUser(0, 0.13);

  PlotGraph(C, 2, 1, g1, Form("b fit errors. Energy window upper end is 645 keV. Geom is 2011-2012."), "Energy Window Low (keV)", "b fit errors", "AP");
  PlotGraph(C, 4, 1, g2, "", "", "", "PSAME");
  PlotGraph(C, 3, 1, g3, "", "", "", "PSAME");
  PlotGraph(C, 6, 1, g4, "", "", "", "PSAME");

  C->cd(1);
  TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
  leg1->AddEntry(g1, Form("stat error, integrated dataset"),"p");
  leg1->AddEntry(g2, Form("syst error, twiddles"),"p");
  leg1->AddEntry(g3, Form("octet RMS"), "p");
  leg1->AddEntry(g4, Form("total sqrt(stat^2 + syst^2)"), "p");
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

