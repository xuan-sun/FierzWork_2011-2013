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
void CreateRealTimeVector_2011();
void CreateRealTimeVector_2012();


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
vector <TDatime> TIME_2011;
vector <TDatime> TIME_2012;
vector <double> octetsErr;
vector <double> chisquared;
vector <double> bMinuitValues_2011;
vector <double> bErrMinuitValues_2011;
vector <double> bMinuitValues_2012;
vector <double> bErrMinuitValues_2012;

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
  C -> Divide(2,1);
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  TH1D *h1 = new TH1D("fierz minuit", "fierz", 50, -0.5, 0.5);
  h1->SetStats(0);

  // fills values into arrays and makes a real-time date vector
  FillArrays(Form("../NewXuanFitter/CorrectBlindingOct2018_newXuanFitter_bFit_%s_%s_Bins_%i-%i.txt", TYPE, "2011-2012", FITMINBIN, FITMAXBIN), h1);
  FillArrays(Form("../NewXuanFitter/CorrectBlindingOct2018_newXuanFitter_bFit_%s_%s_Bins_%i-%i.txt", TYPE, "2012-2013", FITMINBIN, FITMAXBIN), h1);

  CreateRealTimeVector_2011();
  CreateRealTimeVector_2012();

  vector <double> err2011(TIME_2011.size(), 0.);
  vector <double> err2012(TIME_2012.size(), 0.0);

  vector <double> converted_TIME_2011;
  vector <double> converted_TIME_2012;

  for(unsigned int i = 0; i < TIME_2011.size(); i++)
  {
    converted_TIME_2011.push_back(TIME_2011[i].Convert());
  }
  for(unsigned int j = 0; j < TIME_2012.size(); j++)
  {
    converted_TIME_2012.push_back(TIME_2012[j].Convert());
  }

  // create night/day TGraphErrors
  vector <double> night_time_2011;
  vector <double> night_timeErr_2011;
  vector <double> night_bVal_2011;
  vector <double> night_bErr_2011;
  vector <double> day_time_2011;
  vector <double> day_timeErr_2011;
  vector <double> day_bVal_2011;
  vector <double> day_bErr_2011;

  vector <double> night_time_2012;
  vector <double> night_timeErr_2012;
  vector <double> night_bVal_2012;
  vector <double> night_bErr_2012;
  vector <double> day_time_2012;
  vector <double> day_timeErr_2012;
  vector <double> day_bVal_2012;
  vector <double> day_bErr_2012;

  for(unsigned int i = 0; i < TIME_2011.size(); i++)
  {
    if(TIME_2011[i].GetHour() > 6 && TIME_2011[i].GetHour() <= 18)
    {
      day_time_2011.push_back(TIME_2011[i].Convert());
      day_timeErr_2011.push_back(2);
      day_bVal_2011.push_back(bMinuitValues_2011[i]);
      day_bErr_2011.push_back(bErrMinuitValues_2011[i]);
    }
    else
    {
      night_time_2011.push_back(TIME_2011[i].Convert());
      night_timeErr_2011.push_back(2);
      night_bVal_2011.push_back(bMinuitValues_2011[i]);
      night_bErr_2011.push_back(bErrMinuitValues_2011[i]);
    }
  }

  for(unsigned int i = 0; i < TIME_2012.size(); i++)
  {
    if(TIME_2012[i].GetHour() > 6 && TIME_2012[i].GetHour() <= 18)
    {
      day_time_2012.push_back(TIME_2012[i].Convert());
      day_timeErr_2012.push_back(2);
      day_bVal_2012.push_back(bMinuitValues_2012[i]);
      day_bErr_2012.push_back(bErrMinuitValues_2012[i]);
    }
    else
    {
      night_time_2012.push_back(TIME_2012[i].Convert());
      night_timeErr_2012.push_back(2);
      night_bVal_2012.push_back(bMinuitValues_2012[i]);
      night_bErr_2012.push_back(bErrMinuitValues_2012[i]);
    }
  }

  cout << "day_time_2011.size() = " << day_time_2011.size() << endl;
  cout << "night_time_2011.size() = " << night_time_2011.size() << endl;
  cout << "day_time_2012.size() = " << day_time_2012.size() << endl;
  cout << "night_time_2012.size() = " << night_time_2012.size() << endl;

  TGraphErrors *g_day_2011 = new TGraphErrors(day_time_2011.size(), &(day_time_2011[0]), &(day_bVal_2011[0]), &(day_timeErr_2011[0]), &(day_bErr_2011[0]));
  TGraphErrors *g_night_2011 = new TGraphErrors(night_time_2011.size(), &(night_time_2011[0]), &(night_bVal_2011[0]), &(night_timeErr_2011[0]), &(night_bErr_2011[0]));
  TGraphErrors *g_day_2012 = new TGraphErrors(day_time_2012.size(), &(day_time_2012[0]), &(day_bVal_2012[0]), &(day_timeErr_2012[0]), &(day_bErr_2012[0]));
  TGraphErrors *g_night_2012 = new TGraphErrors(night_time_2012.size(), &(night_time_2012[0]), &(night_bVal_2012[0]), &(night_timeErr_2012[0]), &(night_bErr_2012[0]));


//  TGraphErrors *g2011 = new TGraphErrors(converted_TIME_2011.size(), &(converted_TIME_2011[0]), &(bMinuitValues_2011[0]), &(err2011[0]), &(bErrMinuitValues_2011[0]));
//  TGraphErrors *g2012 = new TGraphErrors(converted_TIME_2012.size(), &(converted_TIME_2012[0]), &(bMinuitValues_2012[0]), &(err2012[0]), &(bErrMinuitValues_2012[0]));

//    g1->GetYaxis()->SetRangeUser(-0.3, 0.3);
//  PlotGraph(C, 1, 1, g2011, "b for 2011-2012", "Dates: day-month (ROOT offset by 24hr)", "b", "AP");
//  PlotGraph(C, 1, 2, g2012, "b for 2012-2013", "Dates: day-month (ROOT offset by 24hr)", "b", "AP");

  PlotGraph(C, 1, 1, g_day_2011, "b for day time 2011-2012", "Dates: day-month (ROOT offset by 24hr)", "b", "AP");
  PlotGraph(C, 2, 1, g_night_2011, "b for night time 2011-2012", "Dates: day-month (ROOT offset by 24hr)", "b", "PSAME");

  PlotGraph(C, 1, 2, g_day_2012, "b for day time 2012-2013", "Dates: day-month (ROOT offset by 24hr)", "b", "AP");
  PlotGraph(C, 2, 2, g_night_2012, "b for night time 2012-2013", "Dates: day-month (ROOT offset by 24hr)", "b", "PSAME");
/*
  C->cd(1);
  TLegend* leg1 = new TLegend(0.6,0.6,0.9,0.8);
  leg1->AddEntry(g_day_2011,"day (06,18)","p");
  leg1->AddEntry(g_night_2011,"night (00,05) (18,24)","p");
  leg1->Draw();
*/
  C->cd(2);
  TLegend* leg2 = new TLegend(0.1,0.8,0.4,0.9);
  leg2->AddEntry(g_day_2012,"day (06,18)","p");
  leg2->AddEntry(g_night_2012,"night (00,05) (18,24)","p");
  leg2->Draw();


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
    gPlot->SetMarkerSize(0.75);
    gPlot->SetMarkerColor(2);
  }
  if(styleIndex == 2)
  {
    gPlot->SetMarkerStyle(22);
    gPlot->SetMarkerSize(0.75);
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
  TDatime da86(2012, 12, 12, 0, 00, 00);
  TDatime da92(2012, 12, 16, 0, 00, 00);
  TDatime da96(2013, 1, 11, 0, 00, 00);
  TDatime da105(2013, 1, 30, 0, 00, 00);

    // all the TLine's needed for 2012-2013 calibration periods
    TLine *t12 = new TLine(da86.Convert(), gPad->GetUymin(), da86.Convert(), gPad->GetUymax());     // Octet 80-85 inclusive
    TLine *t13 = new TLine(da92.Convert(), gPad->GetUymin(), da92.Convert(), gPad->GetUymax());     // Octet 86-91 inclusive
    TLine *t14 = new TLine(da96.Convert(), gPad->GetUymin(), da96.Convert(), gPad->GetUymax());     // Octet 92-95 inclusive
    TLine *t15 = new TLine(da105.Convert(), gPad->GetUymin(), da105.Convert(), gPad->GetUymax());   // Octet 96-104 inclusive
//    TLine *t16 = new TLine(105.5, gPad->GetUymin(), 105.5, gPad->GetUymax());   // Octet 105-120 inclusive

    t12->SetLineStyle(7);
    t12->Draw("SAME");
    t13->SetLineStyle(7);
    t13->Draw("SAME");
    t14->SetLineStyle(7);
    t14->Draw("SAME");
    t15->SetLineStyle(7);
    t15->Draw("SAME");
//    t16->SetLineStyle(7);
//    t16->Draw("SAME");
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

      if(evt.octNb >= 0 && evt.octNb <= 59)
      {
        bMinuitValues_2011.push_back(evt.b_minuitFit);
        bErrMinuitValues_2011.push_back(evt.bErr_minuitFit);
      }
      else if(evt.octNb >= 60 && evt.octNb <= 121)
      {
        bMinuitValues_2012.push_back(evt.b_minuitFit);
        bErrMinuitValues_2012.push_back(evt.bErr_minuitFit);
      }
    }


    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

void CreateRealTimeVector_2011()
{
//----------------------------------------------//
// 2011-2012 octets that are actually used are below
//----------------------------------------------//

  // octet 0
  TDatime da0(2011, 10, 23, 4, 8, 00);
  TIME_2011.push_back(da0);

  // octet 1
  TDatime da1(2011, 10, 23, 12, 8, 00);
  TIME_2011.push_back(da1);

  // octet 2
  TDatime da2(2011, 10, 23, 19, 47, 00);
  TIME_2011.push_back(da2);

  // octet 3
  TDatime da3(2011, 10, 24, 2, 58, 00);
  TIME_2011.push_back(da3);

  // octet 4
  TDatime da4(2011, 10, 26, 4, 32, 00);
  TIME_2011.push_back(da4);

  // octet 5
  TDatime da5(2011, 11, 05, 1, 45, 00);
  TIME_2011.push_back(da5);

  // octet 6
  TDatime da6(2011, 11, 05, 12, 48, 00);
  TIME_2011.push_back(da6);

  // octet 7
  TDatime da7(2011, 11, 05, 23, 5, 00);
  TIME_2011.push_back(da7);

  // octet 8
  TDatime da8(2011, 11, 06, 12, 50, 00);
  TIME_2011.push_back(da8);

  // octet 9
//  TDatime da9(2011, 11, 07, 0, 00, 00);
//  TIME.push_back(da9);

  // octet 10
  TDatime da10(2011, 11, 11, 5, 8, 00);
  TIME_2011.push_back(da10);

  // octet 11
  TDatime da11(2011, 11, 11, 15, 56, 00);
  TIME_2011.push_back(da11);

  // octet 12
  TDatime da12(2011, 11, 12, 0, 17, 00);
  TIME_2011.push_back(da12);

  // octet 13
  TDatime da13(2011, 11, 12, 16, 54, 00);
  TIME_2011.push_back(da13);

  // octet 14
  TDatime da14(2011, 11, 13, 1, 32, 00);
  TIME_2011.push_back(da14);

  // octet 15
  TDatime da15(2011, 12, 2, 23, 37, 00);
  TIME_2011.push_back(da15);

  // octet 16
  TDatime da16(2011, 12, 03, 8, 32, 00);
  TIME_2011.push_back(da16);

  // octet 17
  TDatime da17(2011, 12, 03, 17, 30, 00);
  TIME_2011.push_back(da17);

  // octet 18
  TDatime da18(2011, 12, 04, 3, 10, 00);
  TIME_2011.push_back(da18);

  // octet 19
  TDatime da19(2011, 12, 04, 12, 50, 00);
  TIME_2011.push_back(da19);

  // octet 20
  TDatime da20(2011, 12, 04, 21, 35, 00);
  TIME_2011.push_back(da20);

  // octet 21
  TDatime da21(2011, 12, 05, 6, 56, 00);
  TIME_2011.push_back(da21);

  // octet 22
  TDatime da22(2011, 12, 05, 16, 22, 00);
  TIME_2011.push_back(da22);

  // octet 23
  TDatime da23(2011, 12, 06, 1, 49, 00);
  TIME_2011.push_back(da23);

  // octet 24
  TDatime da24(2011, 12, 17, 2, 17, 00);
  TIME_2011.push_back(da24);

  // octet 25
  TDatime da25(2011, 12, 17, 14, 12, 00);
  TIME_2011.push_back(da25);

  // octet 26
  TDatime da26(2011, 12, 18, 3, 4, 00);
  TIME_2011.push_back(da26);

  // octet 27
  TDatime da27(2011, 12, 18, 13, 18, 00);
  TIME_2011.push_back(da27);

  // octet 28
  TDatime da28(2011, 12, 19, 1, 32, 00);
  TIME_2011.push_back(da28);

  // octet 29
  TDatime da29(2011, 12, 19, 12, 36, 00);
  TIME_2011.push_back(da29);

  // octet 30
  TDatime da30(2011, 12, 20, 0, 34, 00);
  TIME_2011.push_back(da30);

  // octet 31
  TDatime da31(2011, 12, 21, 23, 16, 00);
  TIME_2011.push_back(da31);

  // octet 32
  TDatime da32(2012, 01, 14, 0, 32, 00);
  TIME_2011.push_back(da32);

  // octet 33
  TDatime da33(2012, 01, 14, 10, 42, 00);
  TIME_2011.push_back(da33);

  // octet 34
  TDatime da34(2012, 01, 14, 20, 50, 00);
  TIME_2011.push_back(da34);

  // octet 35
  TDatime da35(2012, 01, 15, 7, 54, 00);
  TIME_2011.push_back(da35);

  // octet 36
  TDatime da36(2012, 01, 15, 18, 4, 00);
  TIME_2011.push_back(da36);

  // octet 37
  TDatime da37(2012, 01, 16, 4, 33, 00);
  TIME_2011.push_back(da37);

  // octet 38
  TDatime da38(2012, 01, 16, 14, 41, 00);
  TIME_2011.push_back(da38);

  // octet 39
  TDatime da39(2012, 01, 17, 0, 58, 00);
  TIME_2011.push_back(da39);

  // octet 40
  TDatime da40(2012, 01, 18, 21, 59, 00);
  TIME_2011.push_back(da40);

  // octet 41
  TDatime da41(2012, 01, 19, 23, 37, 00);
  TIME_2011.push_back(da41);

  // octet 42
  TDatime da42(2012, 01, 20, 23, 22, 00);
  TIME_2011.push_back(da42);

  // octet 43
  TDatime da43(2012, 01, 21, 15, 48, 00);
  TIME_2011.push_back(da43);

  // octet 44
  TDatime da44(2012, 01, 22, 1, 53, 00);
  TIME_2011.push_back(da44);

  // octet 45
  TDatime da45(2012, 01, 22, 12, 1, 00);
  TIME_2011.push_back(da45);

  // octet 46
  TDatime da46(2012, 01, 22, 22, 48, 00);
  TIME_2011.push_back(da46);

  // octet 47
  TDatime da47(2012, 02, 11, 3, 53, 00);
  TIME_2011.push_back(da47);

  // octet 48
  TDatime da48(2012, 02, 11, 14, 18, 00);
  TIME_2011.push_back(da48);

  // octet 49
  TDatime da49(2012, 02, 12, 6, 33, 00);
  TIME_2011.push_back(da49);

  // octet 50
  TDatime da50(2012, 02, 12, 22, 54, 00);
  TIME_2011.push_back(da50);

  // octet 51
  TDatime da51(2012, 02, 16, 23, 30, 00);
  TIME_2011.push_back(da51);

  // octet 52
  TDatime da52(2012, 02, 18, 1, 28, 00);
  TIME_2011.push_back(da52);

  // octet 53
  TDatime da53(2012, 02, 18, 13, 33, 00);
  TIME_2011.push_back(da53);

  // octet 54
  TDatime da54(2012, 02, 18, 23, 54, 00);
  TIME_2011.push_back(da54);

  // octet 55
  TDatime da55(2012, 02, 19, 10, 00, 00);
  TIME_2011.push_back(da55);

  // octet 56
  TDatime da56(2012, 02, 19, 20, 31, 00);
  TIME_2011.push_back(da56);

  // octet 57
  TDatime da57(2012, 02, 20, 7, 13, 00);
  TIME_2011.push_back(da57);

  // octet 58
  TDatime da58(2012, 02, 20, 17, 45, 00);
  TIME_2011.push_back(da58);

  // octet 59
//  TDatime da59(2012, 02, 24, 0, 00, 00);
//  TIME.push_back(da59);

}

void CreateRealTimeVector_2012()
{
//----------------------------------------------//
// 2012-2013 octets that are actually used are below
//----------------------------------------------//

  // octet 80
  TDatime da80(2012, 12, 7, 20, 29, 00);
  TIME_2012.push_back(da80);

  // octet 81
  TDatime da81(2012, 12, 8, 8, 37, 00);
  TIME_2012.push_back(da81);

  // octet 82
  TDatime da82(2012, 12, 8, 18, 9, 00);
  TIME_2012.push_back(da82);

  // octet 83
  TDatime da83(2012, 12, 9, 3, 45, 00);
  TIME_2012.push_back(da83);

  // octet 84
  TDatime da84(2012, 12, 9, 13, 18, 00);
  TIME_2012.push_back(da84);

  // octet 85
  TDatime da85(2012, 12, 10, 0, 37, 00);
  TIME_2012.push_back(da85);

  // octet 86
  TDatime da86(2012, 12, 12, 0, 3, 00);
  TIME_2012.push_back(da86);

  // octet 87
  TDatime da87(2012, 12, 12, 18, 36, 00);
  TIME_2012.push_back(da87);

  // octet 88
  TDatime da88(2012, 12, 14, 22, 33, 00);
  TIME_2012.push_back(da88);

  // octet 89
  TDatime da89(2012, 12, 15, 8, 38, 00);
  TIME_2012.push_back(da89);

  // octet 90
  TDatime da90(2012, 12, 15, 18, 55, 00);
  TIME_2012.push_back(da90);

  // octet 91
//  TDatime da91(2012, 12, 15, 12, 00, 00);
//  TIME.push_back(da91);

  // octet 92
  TDatime da92(2012, 12, 16, 22, 06, 00);
  TIME_2012.push_back(da92);

  // octet 93
//  TDatime da93(2012, 12, 16, 12, 00, 00);
//  TIME.push_back(da93);

  // octet 94
  TDatime da94(2012, 12, 17, 9, 55, 00);
  TIME_2012.push_back(da94);

  // octet 95
  TDatime da95(2012, 12, 18, 0, 2, 00);
  TIME_2012.push_back(da95);

  // octet 96
  TDatime da96(2013, 1, 12, 0, 10, 00);
  TIME_2012.push_back(da96);

  // octet 97
  TDatime da97(2013, 1, 12, 8, 55, 00);
  TIME_2012.push_back(da97);

  // octet 98
  TDatime da98(2013, 1, 12, 22, 21, 00);
  TIME_2012.push_back(da98);

  // octet 99
  TDatime da99(2013, 1, 13, 11, 59, 00);
  TIME_2012.push_back(da99);

  // octet 100
  TDatime da100(2013, 1, 13, 22, 9, 00);
  TIME_2012.push_back(da100);

  // octet 101
//  TDatime da101(2013, 1, 17, 0, 00, 00);
//  TIME.push_back(da101);

  // octet 102
  TDatime da102(2013, 1, 18, 1, 37, 00);
  TIME_2012.push_back(da102);

  // octet 103
  TDatime da103(2013, 1, 19, 14, 51, 00);
  TIME_2012.push_back(da103);

  // octet 104
  TDatime da104(2013, 1, 20, 1, 35, 00);
  TIME_2012.push_back(da104);

  // octet 105
  TDatime da105(2013, 1, 20, 10, 53, 00);
  TIME_2012.push_back(da105);

  // octet 106
  TDatime da106(2013, 1, 20, 20, 7, 00);
  TIME_2012.push_back(da106);

  // octet 107
//  TDatime da107(2013, 1, 21, 8, 00, 00);
//  TIME.push_back(da107);

  // octet 108
  TDatime da108(2013, 1, 21, 13, 8, 00);
  TIME_2012.push_back(da108);

  // octet 109
  TDatime da109(2013, 1, 21, 23, 29, 00);
  TIME_2012.push_back(da109);

  // octet 110
  TDatime da110(2013, 1, 25, 20, 4, 00);
  TIME_2012.push_back(da110);

  // octet 111
  TDatime da111(2013, 1, 26, 8, 26, 00);
  TIME_2012.push_back(da111);

  // octet 112
  TDatime da112(2013, 1, 27, 1, 42, 00);
  TIME_2012.push_back(da112);

  // octet 113
  TDatime da113(2013, 1, 27, 11, 58, 00);
  TIME_2012.push_back(da113);

  // octet 114
  TDatime da114(2013, 1, 27, 21, 46, 00);
  TIME_2012.push_back(da114);

  // octet 115
  TDatime da115(2013, 1, 31, 0, 14, 00);
  TIME_2012.push_back(da115);

  // octet 116
  TDatime da116(2013, 1, 31, 20, 14, 00);
  TIME_2012.push_back(da116);

  // octet 117
  TDatime da117(2013, 2, 2, 4, 20, 00);
  TIME_2012.push_back(da117);

  // octet 118
  TDatime da118(2013, 2, 2, 17, 46, 00);
  TIME_2012.push_back(da118);

  // octet 119
  TDatime da119(2013, 2, 3, 5, 7, 00);
  TIME_2012.push_back(da119);

  // octet 120
  TDatime da120(2013, 2, 3, 18, 49, 00);
  TIME_2012.push_back(da120);

  // octet 121
//  TDatime da121(2013, 2, 4, 0, 00, 00);
//  TIME.push_back(da121);

}
