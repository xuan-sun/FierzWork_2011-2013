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

using            namespace std;

// Fundamental constants that get used
const double m_e = 511.00;                                              ///< electron mass, keV/c^2


//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString command);
void FillArrays(TString fileName, TH1D *hist);

struct entry
{
  int octNb;
  double chisquared;
  double ndf;
  double xperndf;
};

// global vectors for creating TGraphs.
vector <double> octets;
vector <double> chisquared;

int main()
{
  TCanvas *C = new TCanvas("canvas", "canvas");
  C -> Divide(2,1);
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  TH1D *h1 = new TH1D("myhist", "myhist", 60, 0, 5);

  FillArrays("chisquared_type0_shapeFactor.txt", h1);

  TGraph *g1 = new TGraph(octets.size(), &(octets[0]), &(chisquared[0]));

  PlotHist(C, 1, 1, h1, "Extracted chi squared per dof values, Type 0", "");
  PlotGraph(C, 1, 2, g1, "chi squared values by octet, Type 0", "AP");

  //prints the canvas with a dynamic TString name of the name of the file
  C -> Print(Form("%s.pdf", "plotChisquared"));
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Extracted #Chi^{2}");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("N");
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(0, 0.000004);

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

  hPlot -> Draw(command);
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString command)
{
  C->cd(canvasIndex);
  gPlot->SetTitle(title);
  gPlot->GetXaxis()->SetTitle("Octet Number");
  gPlot->GetXaxis()->CenterTitle();
  gPlot->GetYaxis()->SetTitle("#frac{#Chi^{2}}{n}");
  gPlot->GetYaxis()->CenterTitle();

  if(styleIndex == 1)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(2);
  }

  gPlot->Draw(command);

  C->Update();

  // all the TLine's needed for 2011-2012 calibration periods
/*  TLine *t1 = new TLine(4.5, gPad->GetUymin(), 4.5, gPad->GetUymax());     // Octet 0-4 inclusive
  TLine *t2 = new TLine(6.5, gPad->GetUymin(), 6.5, gPad->GetUymax());     // Octet 5-6 inclusive
  TLine *t3 = new TLine(9.5, gPad->GetUymin(), 9.5, gPad->GetUymax());     // Octet 7-9 inclusive
  TLine *t4 = new TLine(14.5, gPad->GetUymin(), 14.5, gPad->GetUymax());   // Octet 10-14 inclusive
  TLine *t5 = new TLine(23.5, gPad->GetUymin(), 23.5, gPad->GetUymax());   // Octet 15-23 inclusive
  TLine *t6 = new TLine(31.5, gPad->GetUymin(), 31.5, gPad->GetUymax());   // Octet 24-31 inclusive
  TLine *t7 = new TLine(39.5, gPad->GetUymin(), 39.5, gPad->GetUymax());   // Octet 32-39 inclusive
  TLine *t8 = new TLine(46.5, gPad->GetUymin(), 46.5, gPad->GetUymax());   // Octet 40-46 inclusive
  TLine *t9 = new TLine(50.5, gPad->GetUymin(), 50.5, gPad->GetUymax());   // Octet 47-50 inclusive
  TLine *t11 = new TLine(59.5, gPad->GetUymin(), 59.5, gPad->GetUymax());  // Octet 51-59 inclusive

  t1->SetLineStyle(7);
  t1->Draw("SAME");
  t2->SetLineStyle(7);
  t2->Draw("SAME");
  t3->SetLineStyle(7);
  t3->Draw("SAME");
  t4->SetLineStyle(7);
  t4->Draw("SAME");
  t5->SetLineStyle(7);
  t5->Draw("SAME");
  t6->SetLineStyle(7);
  t6->Draw("SAME");
  t7->SetLineStyle(7);
  t7->Draw("SAME");
  t8->SetLineStyle(7);
  t8->Draw("SAME");
  t9->SetLineStyle(7);
  t9->Draw("SAME");
  t11->SetLineStyle(7);
  t11->Draw("SAME");
*/
}

void FillArrays(TString fileName, TH1D* hist)
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
		>> evt.chisquared
		>> evt.ndf
		>> evt.xperndf;
      {
	counter++;
        hist -> Fill(evt.xperndf);
	octets.push_back(evt.octNb);
	chisquared.push_back(evt.xperndf);
      }
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}
