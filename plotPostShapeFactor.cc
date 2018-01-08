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

#define		TYPE	"type0"

using            namespace std;

// Fundamental constants that get used
const double m_e = 511.00;                                              ///< electron mass, keV/c^2
double NDF = 0;

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
  double chisquared_hf;
  double ndf_hf;
  double xperndf_hf;
  double b_hf;
  double bErr_hf;
};

// global vectors for creating TGraphs.
vector <double> octets;
vector <double> octetsErr;
vector <double> chisquared;
vector <double> bValues;
vector <double> bErrValues;

int main()
{
  TCanvas *C = new TCanvas("canvas", "canvas");
  C -> Divide(2,1);
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

//  TH1D *h1 = new TH1D("fierz", "fierz", 80, -2, 2);
  TH1D *h1 = new TH1D("chi2","chi2", 80, 0, 4);

  FillArrays(Form("TFF_HF_chisquared_%s_shapeFactor.txt", TYPE), h1);

  TGraph *g1 = new TGraph(octets.size(), &(octets[0]), &(chisquared[0]));
//  TGraphErrors *g1 = new TGraphErrors(octets.size(), &(octets[0]), &(bValues[0]), &(octetsErr[0]), &(bErrValues[0]));

  PlotHist(C, 1, 1, h1, Form("Extracted chi squared per dof values, %s", TYPE), "");
  PlotGraph(C, 1, 2, g1, Form("chi squared values by octet, %s", TYPE), "AP");
//  PlotHist(C, 1, 1, h1, Form("TH1::Fit() b Results, %s", TYPE), "");
//  PlotGraph(C, 1, 2, g1, Form("TH1::Fit() b results by octet, %s", TYPE), "AP");

/*
  TF1 *theoryChi = new TF1("theory", Form("-1*(TMath::Prob(x*%f, %f) - TMath::Prob((x-0.1)*%f, %f))", NDF, NDF, NDF, NDF), 0.2, 5);
  TH1D *theoryChiHist = (TH1D*)(theoryChi->GetHistogram());
  double h1Tot = 0;
  double theoryHTot = 0;
  for(unsigned int i = 0; i <= 80; i++)
  {
    h1Tot = h1Tot + h1->GetBinContent(i);
    theoryHTot = theoryHTot + theoryChiHist->GetBinContent(i);

  }
  cout << "h1 has content = " << h1Tot << endl;
  cout << "theoryChiHist has content = " << theoryHTot << endl;
  theoryChiHist->Scale(h1Tot / theoryHTot);
  PlotHist(C, 2, 1, theoryChiHist, "", "SAME");
*/
  //prints the canvas with a dynamic TString name of the name of the file
  C -> Print(Form("resultPlotsOfShapeFactorCode_%s.pdf", TYPE));
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot->GetXaxis()->SetTitle("#Chi^{2}_{HF}");
//  hPlot -> GetXaxis() -> SetTitle("b_{HF}");
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
  gPlot->GetYaxis()->SetTitle("#Chi^{2}_{HF}");
//  gPlot->GetYaxis()->SetTitle("b_{HF}");
  gPlot->GetYaxis()->CenterTitle();

  if(styleIndex == 1)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(2);
  }

  gPlot->Draw(command);

  C->Update();
/*
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
*/

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
		>> evt.xperndf
		>> evt.chisquared_hf
		>> evt.ndf_hf
		>> evt.xperndf_hf
		>> evt.b_hf
		>> evt.bErr_hf;
      {
	counter++;
//        hist -> Fill(evt.b_hf);
	hist->Fill(evt.xperndf_hf);
	octets.push_back(evt.octNb);
	octetsErr.push_back(0.5);
	chisquared.push_back(evt.xperndf_hf);
	bValues.push_back(evt.b_hf);
	bErrValues.push_back(evt.bErr_hf);
	NDF = evt.ndf_hf;
      }
    }


    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}
