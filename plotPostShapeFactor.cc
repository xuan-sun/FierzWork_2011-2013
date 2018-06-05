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
#define		GEOM	"2011-2012"
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
void FillArrays(TString fileName, TH1D *hist1, TH1D* hist2, TH1D* hist3, int codeOption);

struct entry
{
  int octNb;
  int minuitFitStatus;
  double chisquared_minuitShape;
  double ndf_minuitShape;
  double chisquaredperdf_minuitShape;
  double chisquared_hfShape;
  double ndf_hfShape;
  double chisquaredperdf_hfShape;
  double chisquared_minuitFit;
  double ndf_minuitFit;
  double chisquaredperdf_minuitFit;
  double b_minuitFit;
  double bErr_minuitFit;
  double b_hfShape;
  double bErr_hfShape;
};

// global vectors for creating TGraphs.
vector <double> octets;
vector <double> octetsErr;
vector <double> chisquared_minFit;
vector <double> chisquared_minShape;
vector <double> chisquared_hfShape;
vector <double> bHFValues;
vector <double> bErrHFValues;
vector <double> bMinuitValues;
vector <double> bErrMinuitValues;

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (option # = 1 or 2)" << endl;
    return 0;
  }

  int option = atoi(argv[1]);

  NDF = FITMAXBIN - FITMINBIN - 1;

  TCanvas *C = new TCanvas("canvas", "canvas");
  C -> Divide(2,1);
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  if(option == 1)
  {
    TH1D *h1 = new TH1D("fierz minuit", "fierz", 40, -1, 1);
    TH1D *h2 = new TH1D("fierz hf", "fierz", 40, -1, 1);
    TH1D *h6 = new TH1D("filler", "filler", 1, 0, 10);

    FillArrays(Form("FitsUsing_newXuanFitter_ShapeFactors_bAndChisquareds_%s_%s_Bins_%i-%i_shapeFactor.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), h1, h2, h6, option);

    TGraphErrors *g1 = new TGraphErrors(octets.size(), &(octets[0]), &(bMinuitValues[0]), &(octetsErr[0]), &(bErrMinuitValues[0]));
    TGraphErrors *g2 = new TGraphErrors(octets.size(), &(octets[0]), &(bHFValues[0]), &(octetsErr[0]), &(bErrHFValues[0]));

    PlotHist(C, 1, 1, h2, Form("b results for TMinuit fit and TH1::Fit(), %s, %s", TYPE, GEOM), "b", "N", "");
    PlotHist(C, 2, 1, h1, "", "", "", "SAME");

    PlotGraph(C, 1, 2, g1, Form("b results by octet, %s, %s", TYPE, GEOM), "Octet Number", "b", "AP");
    PlotGraph(C, 2, 2, g2, "", "", "", "PSAME");

    C->cd(1);
    TLegend* leg1 = new TLegend(0.6,0.6,0.9,0.8);
    leg1->AddEntry(h1,"b xuanFitter","f");
    leg1->AddEntry(h2,"b Shape Factor fit","f");
    leg1->AddEntry(g1,"b xuanFitter","p");
    leg1->AddEntry(g2,"b Shape Factor fit","p");
    leg1->Draw();

    C -> Print(Form("ReBLINDed_bValues_plotPostShapeFactor_%s_%s_Bins_%i-%i.pdf", TYPE, GEOM, FITMINBIN, FITMAXBIN));
  }
  else if(option == 2)
  {
    TH1D *h3 = new TH1D("chi2 Minuit Fit", "Chi2 Minuit Fit", 80, 0, 4);
    TH1D *h4 = new TH1D("chi2 Minuit Shape", "Chi2 Minuit Shape", 80, 0, 4);
    TH1D *h5 = new TH1D("chi2 HF Shape", "Chi2 HF Shape", 80, 0, 4);

    FillArrays(Form("FitsUsing_newXuanFitter_ShapeFactors_bAndChisquareds_%s_%s_Bins_%i-%i_shapeFactor.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), h3, h4, h5, option);

    TGraph *g3 = new TGraph(octets.size(), &(octets[0]), &(chisquared_minFit[0]));
    TGraph *g4 = new TGraph(octets.size(), &(octets[0]), &(chisquared_minShape[0]));
    TGraph *g5 = new TGraph(octets.size(), &(octets[0]), &(chisquared_hfShape[0]));

    PlotHist(C, 1, 1, h3, Form("#Chi^{2}_{DF} results, %s, %s", TYPE, GEOM), "#frac{#Chi^{2}}{n}", "N", "");
    PlotHist(C, 2, 1, h4, "", "", "", "SAME");
    PlotHist(C, 3, 1, h5, "", "", "", "SAME");

    TF1 *theoryChi = new TF1("theory", Form("-1*(TMath::Prob(x*%f, %f) - TMath::Prob((x-0.1)*%f, %f))", NDF, NDF, NDF, NDF), 0.2, 5);
    TH1D *theoryChiHist = (TH1D*)(theoryChi->GetHistogram());
    double hTot = 0;
    double theoryHTot = 0;
    for(int i = 0; i <= h5->GetNbinsX(); i++)
    {
      hTot = hTot + h5->GetBinContent(i);
      theoryHTot = theoryHTot + theoryChiHist->GetBinContent(i);

    }
    theoryChiHist->Scale(hTot / theoryHTot);
    PlotHist(C, 4, 1, theoryChiHist, "", "", "", "SAME");

    g3->GetYaxis()->SetRangeUser(gPad->GetUymin(), 2.5);
    g4->GetYaxis()->SetRangeUser(gPad->GetUymin(), 2.5);
    g5->GetYaxis()->SetRangeUser(gPad->GetUymin(), 2.5);
    PlotGraph(C, 1, 2, g3, Form("chi squared values by octet, %s, %s", TYPE, GEOM), "Octet Number", "#frac{#Chi^{2}}{n}", "AP");
    PlotGraph(C, 2, 2, g4, "", "", "", "PSAME");
    PlotGraph(C, 3, 2, g5, "", "", "", "PSAME");


    C->cd(1);
    TLegend* leg2 = new TLegend(0.6,0.5,0.9,0.8);
    leg2->AddEntry(h3,"#Chi^{2} xuanFitter","f");
    leg2->AddEntry(h4,"#Chi^{2} b_{xuanFitter} to Shape","f");
    leg2->AddEntry(h5,"#Chi^{2} Shape Function fit","f");
    leg2->AddEntry(theoryChiHist,"Theory #Chi^{2} dist","f");
    leg2->AddEntry(g3,"xuanFitter chi","p");
    leg2->AddEntry(g4,"xuanFitterb to Shape","p");
    leg2->AddEntry(g5,"Shape Function fit","p");
    leg2->Draw();

    C -> Print(Form("ReBLINDed_chisquareds_plotPostShapeFactor_%s_%s_Bins_%i-%i.pdf", TYPE, GEOM, FITMINBIN, FITMAXBIN));
  }

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
    hPlot->GetYaxis()->SetRangeUser(gPad->GetUymin(), hPlot->GetMaximum()*1.2);
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

  hPlot -> Draw(command);

  C->Update();
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command)
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

  gPlot->Draw(command);

  C->Update();

  if(GEOM == "2011-2012")
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

  if(GEOM == "2012-2013")
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


void FillArrays(TString fileName, TH1D* hist1, TH1D* hist2, TH1D* hist3, int codeOption)
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
		>> evt.minuitFitStatus
		>> evt.chisquared_minuitShape
		>> evt.ndf_minuitShape
		>> evt.chisquaredperdf_minuitShape
		>> evt.chisquared_hfShape
		>> evt.ndf_hfShape
		>> evt.chisquaredperdf_hfShape
		>> evt.chisquared_minuitFit
		>> evt.ndf_minuitFit
		>> evt.chisquaredperdf_minuitFit
		>> evt.b_minuitFit
		>> evt.bErr_minuitFit
		>> evt.b_hfShape
		>> evt.bErr_hfShape;
      {
	counter++;

	if(codeOption == 1)
	{
          hist1->Fill(evt.b_minuitFit);
          hist2->Fill(evt.b_hfShape);
	}
	else if(codeOption == 2)
	{
	  hist1->Fill(evt.chisquaredperdf_minuitFit);
	  hist2->Fill(evt.chisquaredperdf_minuitShape);
	  hist3->Fill(evt.chisquaredperdf_hfShape);
	}
	octets.push_back(evt.octNb);
	octetsErr.push_back(0.5);
	chisquared_minFit.push_back(evt.chisquaredperdf_minuitFit);
	chisquared_minShape.push_back(evt.chisquaredperdf_minuitShape);
	chisquared_hfShape.push_back(evt.chisquaredperdf_hfShape);
	bHFValues.push_back(evt.b_hfShape);
	bErrHFValues.push_back(evt.bErr_hfShape);
	bMinuitValues.push_back(evt.b_minuitFit);
	bErrMinuitValues.push_back(evt.bErr_minuitFit);
      }
    }


    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}
