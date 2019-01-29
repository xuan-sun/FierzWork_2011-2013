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
#define		FITMINBIN	17
#define		FITMAXBIN	65

using            namespace std;

// Fundamental constants that get used
const double m_e = 511.00;                                              ///< electron mass, keV/c^2
double NDF = -1;

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command);

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

  vector <double> calPeriod;
  vector <double> calPeriodError;
  vector <double> CeMean;
  vector <double> CeError;
  vector <double> SnMean;
  vector <double> SnError;
  vector <double> Bi1Mean;
  vector <double> Bi1Error;
  vector <double> Bi2Mean;
  vector <double> Bi2Error;

  calPeriod.push_back(1.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-1.13858);
  CeError.push_back(1.11939);
  SnMean.push_back(1.22008);
  SnError.push_back(3.86199);
  Bi1Mean.push_back(0.968208);
  Bi1Error.push_back(5.91032);
  Bi2Mean.push_back(-0.661375);
  Bi2Error.push_back(9.89761);

  calPeriod.push_back(2.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-1.65087);
  CeError.push_back(2.28167);
  SnMean.push_back(0.0995385);
  SnError.push_back(2.73646);
  Bi1Mean.push_back(-2.11854);
  Bi1Error.push_back(3.44558);
  Bi2Mean.push_back(-1.36738);
  Bi2Error.push_back(6.24466);

  calPeriod.push_back(3.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-0.808786);
  CeError.push_back(2.24667);
  SnMean.push_back(0.0175714);
  SnError.push_back(1.98601);
  Bi1Mean.push_back(-3.45408);
  Bi1Error.push_back(4.21966);
  Bi2Mean.push_back(-3.852);
  Bi2Error.push_back(7.11128);

  calPeriod.push_back(4.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-1.27631);
  CeError.push_back(1.69661);
  SnMean.push_back(0.888833);
  SnError.push_back(2.32956);
  Bi1Mean.push_back(-2.34321);
  Bi1Error.push_back(3.09971);
  Bi2Mean.push_back(-1.96563);
  Bi2Error.push_back(4.52879);

  calPeriod.push_back(5.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-0.346286);
  CeError.push_back(1.69075);
  SnMean.push_back(0.161143);
  SnError.push_back(1.63782);
  Bi1Mean.push_back(-2.69142);
  Bi1Error.push_back(3.19002);
  Bi2Mean.push_back(-2.73467);
  Bi2Error.push_back(5.86021);

  calPeriod.push_back(6.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-2.45114);
  CeError.push_back(1.64703);
  SnMean.push_back(-0.298143);
  SnError.push_back(2.45811);
  Bi1Mean.push_back(-1.92533);
  Bi1Error.push_back(2.2737);
  Bi2Mean.push_back(0.0303333);
  Bi2Error.push_back(3.20514);

  calPeriod.push_back(7.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-1.14017);
  CeError.push_back(1.65866);
  SnMean.push_back(2.68067);
  SnError.push_back(2.01684);
  Bi1Mean.push_back(0.959083);
  Bi1Error.push_back(2.75188);
  Bi2Mean.push_back(-1.40958);
  Bi2Error.push_back(2.76427);

  calPeriod.push_back(8.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-2.23493);
  CeError.push_back(1.59867);
  SnMean.push_back(0.742);
  SnError.push_back(1.6105);
  Bi1Mean.push_back(-0.889417);
  Bi1Error.push_back(1.41915);
  Bi2Mean.push_back(-1.31875);
  Bi2Error.push_back(1.79787);

  calPeriod.push_back(9.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-1.86807);
  CeError.push_back(1.29569);
  SnMean.push_back(0.989);
  SnError.push_back(2.22534);
  Bi1Mean.push_back(-0.365417);
  Bi1Error.push_back(2.31083);
  Bi2Mean.push_back(-0.23025);
  Bi2Error.push_back(3.94056);

  calPeriod.push_back(10.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(0);
  CeError.push_back(0);
  SnMean.push_back(0);
  SnError.push_back(0);
  Bi1Mean.push_back(0);
  Bi1Error.push_back(0);
  Bi2Mean.push_back(0);
  Bi2Error.push_back(0);

  calPeriod.push_back(11.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-1.72943);
  CeError.push_back(1.36);
  SnMean.push_back(1.95307);
  SnError.push_back(1.75864);
  Bi1Mean.push_back(-0.101);
  Bi1Error.push_back(3.03533);
  Bi2Mean.push_back(-1.70543);
  Bi2Error.push_back(3.75728);

  calPeriod.push_back(12.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(0);
  CeError.push_back(0);
  SnMean.push_back(1.64942);
  SnError.push_back(1.79246);
  Bi1Mean.push_back(0);
  Bi1Error.push_back(0);
  Bi2Mean.push_back(0);
  Bi2Error.push_back(0);

  calPeriod.push_back(16.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(1.09621);
  CeError.push_back(2.29213);
  SnMean.push_back(-5.22936);
  SnError.push_back(2.85318);
  Bi1Mean.push_back(-1.29221);
  Bi1Error.push_back(3.33494);
  Bi2Mean.push_back(1.87157);
  Bi2Error.push_back(5.54234);

  calPeriod.push_back(17.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-1.36769);
  CeError.push_back(1.57276);
  SnMean.push_back(-1.98808);
  SnError.push_back(2.39383);
  Bi1Mean.push_back(0.2635);
  Bi1Error.push_back(2.93369);
  Bi2Mean.push_back(0.139269);
  Bi2Error.push_back(4.22588);

  calPeriod.push_back(18.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-0.891286);
  CeError.push_back(1.08107);
  SnMean.push_back(-1.87921);
  SnError.push_back(1.66558);
  Bi1Mean.push_back(-0.2915);
  Bi1Error.push_back(2.79704);
  Bi2Mean.push_back(0.210143);
  Bi2Error.push_back(3.96107);

  calPeriod.push_back(19.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-1.20943);
  CeError.push_back(1.44004);
  SnMean.push_back(-1.76043);
  SnError.push_back(1.66915);
  Bi1Mean.push_back(-1.07779);
  Bi1Error.push_back(2.45911);
  Bi2Mean.push_back(0.0297857);
  Bi2Error.push_back(5.42827);

  calPeriod.push_back(20.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-0.706);
  CeError.push_back(1.19499);
  SnMean.push_back(-2.05807);
  SnError.push_back(1.92076);
  Bi1Mean.push_back(-0.724786);
  Bi1Error.push_back(3.40628);
  Bi2Mean.push_back(-0.621429);
  Bi2Error.push_back(3.00804);

  calPeriod.push_back(21.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(0.570643);
  CeError.push_back(2.63438);
  SnMean.push_back(1.07464);
  SnError.push_back(3.65413);
  Bi1Mean.push_back(1.68879);
  Bi1Error.push_back(4.64891);
  Bi2Mean.push_back(0.390571);
  Bi2Error.push_back(9.79826);

  calPeriod.push_back(22.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-1.80114);
  CeError.push_back(1.49944);
  SnMean.push_back(-1.85214);
  SnError.push_back(2.30521);
  Bi1Mean.push_back(-0.912786);
  Bi1Error.push_back(4.3466);
  Bi2Mean.push_back(-0.606214);
  Bi2Error.push_back(7.6206);

  calPeriod.push_back(23.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(-1.60521);
  CeError.push_back(1.90245);
  SnMean.push_back(-2.27064);
  SnError.push_back(2.94631);
  Bi1Mean.push_back(-1.3515);
  Bi1Error.push_back(5.81087);
  Bi2Mean.push_back(-1.44271);
  Bi2Error.push_back(9.15008);

  calPeriod.push_back(24.0);
  calPeriodError.push_back(0.5);
  CeMean.push_back(0);
  CeError.push_back(0);
  SnMean.push_back(-4.375);
  SnError.push_back(0.958078);
  Bi1Mean.push_back(0);
  Bi1Error.push_back(0);
  Bi2Mean.push_back(0);
  Bi2Error.push_back(0);




//  TGraphErrors *g1 = new TGraphErrors(octets.size(), &(octets[0]), &(bMinuitValues[0]), &(octetsErr[0]), &(bErrMinuitValues[0]));
//  TGraphErrors *g1 = new TGraphErrors(octets.size(), &(octets[0]), &(chisquared[0]), &(octetsErr[0]), &(chisquaredError[0]));

//    g1->GetYaxis()->SetRangeUser(-0.3, 0.3);
//  PlotGraph(C, 1, 1, g1, "b for all octets", "Octet Number", "b", "AP");

  C->cd(1);
//  TLegend* leg1 = new TLegend(0.6,0.6,0.9,0.8);
//  leg1->AddEntry(h1,"b xuanFitter","f");
//  leg1->Draw();


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

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command)
{
  C->cd(canvasIndex);
  gPlot->SetTitle(title);
  gPlot->GetXaxis()->SetTitle(xAxis);
  gPlot->GetXaxis()->CenterTitle();
  gPlot->GetYaxis()->SetTitle(yAxis);
  gPlot->GetYaxis()->CenterTitle();

  gPlot->SetMarkerStyle(21);
  gPlot->SetMarkerSize(0.5);
  gPlot->SetMarkerColor(styleIndex);

  gPlot->Draw(command);

  C->Update();

/*
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
*/
}

