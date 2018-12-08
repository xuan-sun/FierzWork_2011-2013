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

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void FillArrays(TString fileName, TH1D *hist1);
void CreateRealTimeVector();


struct entry
{
  int octNb;
  double avg_mE;
  double chisquared;
  double ndf;
  double chisquaredperndf;
  double b_minuitFit;
  double bErr_minuitFit;
  double A_minuitFit;
  double AErr_minuitFit;
  int fitMatrixStatus;
};

// global vectors for creating TGraphs.
vector <double> octets;
vector <double> TIME;
vector <double> octetsErr;
vector <double> chisquared;
vector <double> bMinuitValues;
vector <double> bErrMinuitValues;

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

  TH1D *h1 = new TH1D("fierz minuit", "fierz", 50, -0.5, 0.5);
  h1->SetStats(0);

  FillArrays(Form("../CorrectBlindingOct2018_newXuanFitter_bFit_%s_%s_Bins_%i-%i.txt", TYPE, "2011-2012", FITMINBIN, FITMAXBIN), h1);
//  FillArrays(Form("../CorrectBlindingOct2018_newXuanFitter_bFit_%s_%s_Bins_%i-%i.txt", TYPE, "2012-2013", FITMINBIN, FITMAXBIN), h1);


  CreateRealTimeVector();

  int endpoint = 40;

  vector <double> test;
  for(int i = 0; i < endpoint; i++)
  {
    test.push_back(i);
  }

  vector <double> zeros;
  for(int i = 0; i < endpoint; i++)
  {
    zeros.push_back(0);
  }


//  TGraphErrors *g1 = new TGraphErrors(test.size(), &(TIME[0]), &(test[0]), &(zeros[0]), &(zeros[0]));

  TGraphErrors *g1 = new TGraphErrors(TIME.size(), &(TIME[0]), &(bMinuitValues[0]), &(octetsErr[0]), &(bErrMinuitValues[0]));

//  PlotHist(C, 1, 1, h1, "b for all octets", "b", "N", "");

//    g1->GetYaxis()->SetRangeUser(-0.3, 0.3);
  PlotGraph(C, 1, 1, g1, "b for all octets", "Dates: day-month (ROOT offset by 24hr)", "b", "AP");

  C->cd(1);
  TLegend* leg1 = new TLegend(0.6,0.6,0.9,0.8);
  leg1->AddEntry(h1,"b xuanFitter","f");
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

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command)
{
  C->cd(canvasIndex);
  gPlot->SetTitle(title);
  gPlot->GetXaxis()->SetTitle(xAxis);
  gPlot->GetXaxis()->CenterTitle();
  gPlot->GetXaxis()->SetTimeDisplay(1);
  gPlot->GetXaxis()->SetTimeFormat("%d-%m");
//  gPlot->GetXaxis()->SetTimeFormat("%d-%m%F2011-10-23 00:00:01");
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
  TDatime da5(2011, 11, 4, 0, 00, 00);
  TDatime da7(2011, 11, 6, 0, 00, 00);
  TDatime da10(2011, 11, 11, 0, 00, 00);
  TDatime da15(2011, 12, 3, 0, 00, 00);
  TDatime da24(2011, 12, 17, 0, 00, 00);
  TDatime da32(2012, 1, 13, 0, 00, 00);
  TDatime da40(2012, 1, 18, 0, 00, 00);
  TDatime da47(2012, 2, 11, 0, 00, 00);
  TDatime da51(2012, 2, 16, 0, 00, 00);


    // all the TLine's needed for 2011-2012 calibration periods
    TLine *t1 = new TLine(da5.Convert(), gPad->GetUymin(), da5.Convert(), gPad->GetUymax());     // Octet 0-4 inclusive
    TLine *t2 = new TLine(da7.Convert(), gPad->GetUymin(), da7.Convert(), gPad->GetUymax());     // Octet 5-6 inclusive
    TLine *t3 = new TLine(da10.Convert(), gPad->GetUymin(), da10.Convert(), gPad->GetUymax());     // Octet 7-9 inclusive
    TLine *t4 = new TLine(da15.Convert(), gPad->GetUymin(), da15.Convert(), gPad->GetUymax());   // Octet 10-14 inclusive
    TLine *t5 = new TLine(da24.Convert(), gPad->GetUymin(), da24.Convert(), gPad->GetUymax());   // Octet 15-23 inclusive
    TLine *t6 = new TLine(da32.Convert(), gPad->GetUymin(), da32.Convert(), gPad->GetUymax());   // Octet 24-31 inclusive
    TLine *t7 = new TLine(da40.Convert(), gPad->GetUymin(), da40.Convert(), gPad->GetUymax());   // Octet 32-39 inclusive
    TLine *t8 = new TLine(da47.Convert(), gPad->GetUymin(), da47.Convert(), gPad->GetUymax());   // Octet 40-46 inclusive
    TLine *t9 = new TLine(da51.Convert(), gPad->GetUymin(), da51.Convert(), gPad->GetUymax());   // Octet 47-50 inclusive
//    TLine *t11 = new TLine(59.5, gPad->GetUymin(), 59.5, gPad->GetUymax());  // Octet 51-59 inclusive

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
//    t11->SetLineStyle(7);
//    t11->Draw("SAME");
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

void FillArrays(TString fileName, TH1D* hist1)
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
		>> evt.b_minuitFit
		>> evt.bErr_minuitFit
		>> evt.A_minuitFit
		>> evt.AErr_minuitFit
		>> evt.fitMatrixStatus;
      counter++;

      hist1->Fill(evt.b_minuitFit);

      octets.push_back(evt.octNb);
      octetsErr.push_back(0.5);
      chisquared.push_back(evt.chisquaredperndf);
      bMinuitValues.push_back(evt.b_minuitFit);
      bErrMinuitValues.push_back(evt.bErr_minuitFit);
    }


    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

void CreateRealTimeVector()
{
//----------------------------------------------//
// 2011-2012 octets that are actually used are below
//----------------------------------------------//

  // octet 0
  TDatime da0(2011, 10, 23, 0, 00, 00);
  TIME.push_back(da0.Convert());

  // octet 1
  TDatime da1(2011, 10, 23, 12, 00, 00);
  TIME.push_back(da1.Convert());

  // octet 2
  TDatime da2(2011, 10, 24, 0, 00, 00);
  TIME.push_back(da2.Convert());

  // octet 3
  TDatime da3(2011, 10, 24, 12, 00, 00);
  TIME.push_back(da3.Convert());

  // octet 4
  TDatime da4(2011, 10, 26, 0, 00, 00);
  TIME.push_back(da4.Convert());

  // octet 5
  TDatime da5(2011, 11, 04, 0, 00, 00);
  TIME.push_back(da5.Convert());

  // octet 6
  TDatime da6(2011, 11, 04, 18, 00, 00);
  TIME.push_back(da6.Convert());

  // octet 7
  TDatime da7(2011, 11, 05, 12, 00, 00);
  TIME.push_back(da7.Convert());

  // octet 8
  TDatime da8(2011, 11, 06, 6, 00, 00);
  TIME.push_back(da8.Convert());

  // octet 9
//  TDatime da9(2011, 11, 07, 0, 00, 00);
//  TIME.push_back(da9.Convert());

  // octet 10
  TDatime da10(2011, 11, 11, 0, 00, 00);
  TIME.push_back(da10.Convert());

  // octet 11
  TDatime da11(2011, 11, 11, 12, 00, 00);
  TIME.push_back(da11.Convert());

  // octet 12
  TDatime da12(2011, 11, 12, 0, 00, 00);
  TIME.push_back(da12.Convert());

  // octet 13
  TDatime da13(2011, 11, 12, 12, 00, 00);
  TIME.push_back(da13.Convert());

  // octet 14
  TDatime da14(2011, 11, 13, 0, 00, 00);
  TIME.push_back(da14.Convert());

  // octet 15
  TDatime da15(2011, 12, 03, 0, 00, 00);
  TIME.push_back(da15.Convert());

  // octet 16
  TDatime da16(2011, 12, 03, 8, 00, 00);
  TIME.push_back(da16.Convert());

  // octet 17
  TDatime da17(2011, 12, 03, 16, 00, 00);
  TIME.push_back(da17.Convert());

  // octet 18
  TDatime da18(2011, 12, 04, 0, 00, 00);
  TIME.push_back(da18.Convert());

  // octet 19
  TDatime da19(2011, 12, 04, 8, 00, 00);
  TIME.push_back(da19.Convert());

  // octet 20
  TDatime da20(2011, 12, 04, 16, 00, 00);
  TIME.push_back(da20.Convert());

  // octet 21
  TDatime da21(2011, 12, 05, 0, 00, 00);
  TIME.push_back(da21.Convert());

  // octet 22
  TDatime da22(2011, 12, 05, 8, 00, 00);
  TIME.push_back(da22.Convert());

  // octet 23
  TDatime da23(2011, 12, 05, 16, 00, 00);
  TIME.push_back(da23.Convert());

  // octet 24
  TDatime da24(2011, 12, 17, 0, 00, 00);
  TIME.push_back(da24.Convert());

  // octet 25
  TDatime da25(2011, 12, 17, 12, 00, 00);
  TIME.push_back(da25.Convert());

  // octet 26
  TDatime da26(2011, 12, 18, 0, 00, 00);
  TIME.push_back(da26.Convert());

  // octet 27
  TDatime da27(2011, 12, 18, 12, 00, 00);
  TIME.push_back(da27.Convert());

  // octet 28
  TDatime da28(2011, 12, 19, 0, 00, 00);
  TIME.push_back(da28.Convert());

  // octet 29
  TDatime da29(2011, 12, 19, 12, 00, 00);
  TIME.push_back(da29.Convert());

  // octet 30
  TDatime da30(2011, 12, 20, 0, 00, 00);
  TIME.push_back(da30.Convert());

  // octet 31
  TDatime da31(2011, 12, 21, 0, 00, 00);
  TIME.push_back(da31.Convert());

  // octet 32
  TDatime da32(2012, 01, 13, 0, 00, 00);
  TIME.push_back(da32.Convert());

  // octet 33
  TDatime da33(2012, 01, 13, 12, 00, 00);
  TIME.push_back(da33.Convert());

  // octet 34
  TDatime da34(2012, 01, 14, 0, 00, 00);
  TIME.push_back(da34.Convert());

  // octet 35
  TDatime da35(2012, 01, 14, 12, 00, 00);
  TIME.push_back(da35.Convert());

  // octet 36
  TDatime da36(2012, 01, 15, 0, 00, 00);
  TIME.push_back(da36.Convert());

  // octet 37
  TDatime da37(2012, 01, 15, 12, 00, 00);
  TIME.push_back(da37.Convert());

  // octet 38
  TDatime da38(2012, 01, 16, 0, 00, 00);
  TIME.push_back(da38.Convert());

  // octet 39
  TDatime da39(2012, 01, 16, 12, 00, 00);
  TIME.push_back(da39.Convert());

  // octet 40
  TDatime da40(2012, 01, 18, 18, 00, 00);
  TIME.push_back(da40.Convert());

  // octet 41
  TDatime da41(2012, 01, 19, 18, 00, 00);
  TIME.push_back(da41.Convert());

  // octet 42
  TDatime da42(2012, 01, 20, 0, 00, 00);
  TIME.push_back(da42.Convert());

  // octet 43
  TDatime da43(2012, 01, 20, 18, 00, 00);
  TIME.push_back(da43.Convert());

  // octet 44
  TDatime da44(2012, 01, 21, 12, 00, 00);
  TIME.push_back(da44.Convert());

  // octet 45
  TDatime da45(2012, 01, 22, 8, 00, 00);
  TIME.push_back(da45.Convert());

  // octet 46
  TDatime da46(2012, 01, 23, 0, 00, 00);
  TIME.push_back(da46.Convert());

  // octet 47
  TDatime da47(2012, 02, 11, 0, 00, 00);
  TIME.push_back(da47.Convert());

  // octet 48
  TDatime da48(2012, 02, 11, 18, 00, 00);
  TIME.push_back(da48.Convert());

  // octet 49
  TDatime da49(2012, 02, 12, 12, 00, 00);
  TIME.push_back(da49.Convert());

  // octet 50
  TDatime da50(2012, 02, 13, 8, 00, 00);
  TIME.push_back(da50.Convert());

  // octet 51
  TDatime da51(2012, 02, 16, 0, 00, 00);
  TIME.push_back(da51.Convert());

  // octet 52
  TDatime da52(2012, 02, 17, 0, 00, 00);
  TIME.push_back(da52.Convert());

  // octet 53
  TDatime da53(2012, 02, 17, 18, 00, 00);
  TIME.push_back(da53.Convert());

  // octet 54
  TDatime da54(2012, 02, 18, 12, 00, 00);
  TIME.push_back(da54.Convert());

  // octet 55
  TDatime da55(2012, 02, 19, 8, 00, 00);
  TIME.push_back(da55.Convert());

  // octet 56
  TDatime da56(2012, 02, 20, 0, 00, 00);
  TIME.push_back(da56.Convert());

  // octet 57
  TDatime da57(2012, 02, 20, 16, 00, 00);
  TIME.push_back(da57.Convert());

  // octet 58
  TDatime da58(2012, 02, 21, 12, 00, 00);
  TIME.push_back(da58.Convert());

  // octet 59
//  TDatime da59(2012, 02, 24, 0, 00, 00);
//  TIME.push_back(da59.Convert());

//----------------------------------------------//
// 2012-2013 octets that are actually used are below
//----------------------------------------------//

}
