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
  int octNb;
  double peakValue;
  double peakEnergy;
};

// global vectors for creating TGraphs.
vector <double> octets;
vector <double> octetsErr;
vector <double> peakEnergy;
vector <double> peakEnergyErr;

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  TCanvas *C = new TCanvas("canvas", "canvas");
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  FillArrays("FindingSpectrumPeak_DataHists_ssDataHist_type0.txt");

  TGraphErrors *g1 = new TGraphErrors(octets.size(), &(octets[0]), &(peakEnergy[0]), &(octetsErr[0]), &(peakEnergyErr[0]));

/*
  TF1 *fit1 = new TF1("fit1", "[0]", -1, 59);
  g1->Fit(fit1, "R");
*/

//    g1->GetYaxis()->SetRangeUser(-0.3, 0.3);
  PlotGraph(C, 1, 1, g1, "peak energy for all octets", "Octet Number", "peak energy (keV)", "AP");

/*
  C->cd(1);
  TLegend* leg1 = new TLegend(0.6,0.6,0.9,0.8);
  leg1->AddEntry(h1,"b xuanFitter","f");
  leg1->Draw();
*/

  double xPrint = 45;
  double yPrint = -0.3;

//  double xPrint = 45;
//  double yPrint = -0.3;

/*
  TLatex t2;
  t2.SetTextSize(0.03);
  t2.SetTextAlign(13);
  t2.DrawLatex(xPrint, yPrint+0.1, Form("#frac{#Chi^{2}}{NDF} = #frac{%f}{%i} = %f", fit1->GetChisquare(), fit1->GetNDF(), (fit1->GetChisquare() / fit1->GetNDF())));
  TLatex t3;
  t3.SetTextSize(0.03);
  t3.SetTextAlign(13);
  t3.DrawLatex(xPrint, yPrint, Form("b #pm #sigma_{b} = %f #pm %f", fit1->GetParameter(0), fit1->GetParError(0)));

  TLatex t4;
  t4.SetTextSize(0.03);
  t4.SetTextAlign(13);
  t4.DrawLatex(xPrint, yPrint-0.1, Form("error value = %f", 0.2158));
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
//  gPlot->GetXaxis()->SetTimeDisplay(1);
//  gPlot->GetXaxis()->SetTimeFormat("%d-%m");
//  gPlot->GetXaxis()->SetTimeFormat("%d-%m-%y%F2011-10-23 00:00:01");
  gPlot->GetYaxis()->SetTitle(yAxis);
  gPlot->GetYaxis()->CenterTitle();

  if(styleIndex == 1)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(2);
  }
  if(styleIndex == 2)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(4);
  }

  gPlot->Draw(command);

  C->Update();

//  if(GEOM == "2011-2012")
  {
    // all the TLine's needed for 2011-2012 calibration periods
    TLine *t1 = new TLine(4.5, gPad->GetUymin(), 4.5, gPad->GetUymax());     // Octet 0-4 inclusive
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
  }

//  if(GEOM == "2012-2013")
  {
    // all the TLine's needed for 2012-2013 calibration periods
    TLine *t12 = new TLine(79.5, gPad->GetUymin(), 79.5, gPad->GetUymax());     // Octet 80-85 inclusive
    TLine *t13 = new TLine(85.5, gPad->GetUymin(), 85.5, gPad->GetUymax());     // Octet 86-91 inclusive
    TLine *t14 = new TLine(91.5, gPad->GetUymin(), 91.5, gPad->GetUymax());     // Octet 92-95 inclusive
    TLine *t15 = new TLine(95.5, gPad->GetUymin(), 95.5, gPad->GetUymax());   // Octet 96-105 inclusive
    TLine *t16 = new TLine(105.5, gPad->GetUymin(), 105.5, gPad->GetUymax());   // Octet 105-120 inclusive

    t12->SetLineStyle(7);
    t12->Draw("SAME");
    t13->SetLineStyle(7);
    t13->Draw("SAME");
    t14->SetLineStyle(7);
    t14->Draw("SAME");
    t15->SetLineStyle(7);
    t15->Draw("SAME");
    t16->SetLineStyle(7);
    t16->Draw("SAME");
  }
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString xAxis, TString yAxis, TString command)
{
  C->cd(canvasIndex);
  gPlot->SetTitle(title);
  gPlot->GetXaxis()->SetTitle(xAxis);
  gPlot->GetXaxis()->CenterTitle();
  gPlot->GetYaxis()->SetTitle(yAxis);
  gPlot->GetYaxis()->CenterTitle();

  if(styleIndex == 1)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(2);
  }
  if(styleIndex == 2)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(4);
  }
  if(styleIndex == 3)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(3);
  }

  gPlot->Draw(command);

  C->Update();
}


void FillArrays(TString fileName)
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
		>> evt.peakValue
		>> evt.peakEnergy;
      counter++;

      octets.push_back(evt.octNb);
      octetsErr.push_back(0.5);
      peakEnergy.push_back(evt.peakEnergy);
      peakEnergyErr.push_back(5);
    }


    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

void ConvertOctetsToDate(int octNb)
{
  if(octNb == 0)
  {
    TDatime da(2011, 10, 23, 00, 00, 00);
    octets.push_back(da.Convert());
  }
/*  else if(octNb == 1)
  {
    TDatime da(0, 0, 0, 12, 00, 00);
    octets.push_back(da.Convert());
  }
*/

}
