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

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command, double maxBinContents);
void FillArrays(TString fileName, TH1D *h, int codeOption);

struct entry
{
  int indexNb;
  double avg_mE;
  double chi2;
  double ndf;
  double chi2_ndf;
  double bFitValue;
  double bFitError;
  double holder1;
  double holder2;
  int covMatrixStatus;
};

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
  C->cd();
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  if(option == 1)
  {
    TH1D *hADownbFit = new TH1D("A_-1_b_fit", "A=-1, twiddled b values, 2011-2012", 100, -0.5, 0.5);
    TH1D *hAUpbFit = new TH1D("A_1_b_fit", "A=1, twiddled b values, 2011-2012", 100, -0.5, 0.5);

    FillArrays(Form("TwiddledbValues_NoAsymm100MillBLINDEDBaseline_A_-1_b_0_newXuanFitter_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), hADownbFit, option);
    FillArrays(Form("TwiddledbValues_NoAsymm100MillBLINDEDBaseline_A_1_b_0_newXuanFitter_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), hAUpbFit, option);

    double max = 0;
    if(hADownbFit->GetMaximum() > hAUpbFit->GetMaximum())
    {
      max = hADownbFit->GetMaximum();
    }
    else
    {
      max = hAUpbFit->GetMaximum();
    }

    PlotHist(C, 1, 1, hADownbFit, Form("b fit results, %s, %s", TYPE, GEOM), "b", "N", "", max);
    PlotHist(C, 2, 1, hAUpbFit, Form("b fit results, %s, %s", TYPE, GEOM), "b", "N", "SAME", max);

    TLegend* leg1 = new TLegend(0.1,0.6,0.35,0.8);
    leg1->AddEntry(hADownbFit,"A=-1, b twiddles","f");
    leg1->AddEntry(hAUpbFit,"A=1, b twiddles","f");
    leg1->Draw();

    C -> Print(Form("Twiddledb_bValues_viewTwiddledbAndChi2_%s_%s_Bins_%i-%i.pdf", TYPE, GEOM, FITMINBIN, FITMAXBIN));
  }
  else if(option == 2)
  {
    TH1D *hADownchi2_ndf = new TH1D("A_-1_chi2_ndf", "A=-1, chi2 per ndf values, 2011-2012", 100, 0, 2);
    TH1D *hAUpchi2_ndf = new TH1D("A_1_chi2_ndf", "A=1, chi2 per ndf values, 2011-2012", 100, 0, 2);

    FillArrays(Form("TwiddledbValues_NoAsymm100MillBLINDEDBaseline_A_-1_b_0_newXuanFitter_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), hADownchi2_ndf, option);
    FillArrays(Form("TwiddledbValues_NoAsymm100MillBLINDEDBaseline_A_1_b_0_newXuanFitter_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), hAUpchi2_ndf, option);

    TF1 *theoryChi = new TF1("theory", Form("-1*(TMath::Prob(x*%f, %f) - TMath::Prob((x-0.1)*%f, %f))", NDF, NDF, NDF, NDF), 0.2, 5);
    TH1D *theoryChiHist = (TH1D*)(theoryChi->GetHistogram());
    double hTot = 0;
    double theoryHTot = 0;
    for(int i = 0; i <= hADownchi2_ndf->GetNbinsX(); i++)
    {
      hTot = hTot + (hADownchi2_ndf->GetBinContent(i))*(hADownchi2_ndf->GetBinWidth(i));
      theoryHTot = theoryHTot + (theoryChiHist->GetBinContent(i))*(theoryChiHist->GetBinWidth(i));

    }

    theoryChiHist->Scale(hTot / theoryHTot);

    double max = 0;
    if(hADownchi2_ndf->GetMaximum() > max)
    {
      max = hADownchi2_ndf->GetMaximum();
    }
    if(hAUpchi2_ndf->GetMaximum() > max)
    {
      max = hAUpchi2_ndf->GetMaximum();
    }
    if(theoryChiHist->GetMaximum() > max);
    {
      max = theoryChiHist->GetMaximum();
    }

    PlotHist(C, 1, 1, hADownchi2_ndf, Form("chi2 per ndf results, %s, %s", TYPE, GEOM), "#frac{#Chi^{2}}{n}", "N", "", max);
    PlotHist(C, 2, 1, hAUpchi2_ndf, Form("chi2 per ndf results, %s, %s", TYPE, GEOM), "#frac{#Chi^{2}}{n}", "N", "SAME", max);
    PlotHist(C, 4, 1, theoryChiHist, "", "", "", "SAME", max);

    TLegend* leg2 = new TLegend(0.1,0.5,0.35,0.8);
    leg2->AddEntry(hADownchi2_ndf,"#Chi^{2} for A=-1","f");
    leg2->AddEntry(hAUpchi2_ndf,"#Chi^{2} for A=1","f");
    leg2->AddEntry(theoryChiHist,"Theory #Chi^{2} dist","f");
    leg2->Draw();

    C -> Print(Form("Twiddledb_chi2Values_viewTwiddledbAndChi2_%s_%s_Bins_%i-%i.pdf", TYPE, GEOM, FITMINBIN, FITMAXBIN));
  }

  //prints the canvas with a dynamic TString name of the name of the file
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

void FillArrays(TString fileName, TH1D* h, int codeOption)
{

  entry evt;

  //opens the file that I name in DATA_FILE_IN
  string buf;
  ifstream infile1;
  cout << "The file being opened is: " << fileName << endl;
  infile1.open(fileName);

  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << fileName << endl;


  while(true)
  {
    getline(infile1, buf);
    istringstream bufstream(buf);

    if(!infile1.eof())
    {
      bufstream >> evt.indexNb
		>> evt.avg_mE
		>> evt.chi2
		>> evt.ndf
		>> evt.chi2_ndf
		>> evt.bFitValue
		>> evt.bFitError
		>> evt.holder1
		>> evt.holder2
		>> evt.covMatrixStatus;

      if(codeOption == 1)
      {
	h->Fill(evt.bFitValue);
      }
      else if(codeOption == 2)
      {
	h->Fill(evt.chi2_ndf);
      }
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}




void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command, double maxBinContents)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot->GetXaxis()->SetTitle(xAxis);
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle(yAxis);
  hPlot -> GetYaxis() -> CenterTitle();
  hPlot->GetYaxis()->SetRangeUser(0, 1.2*maxBinContents);

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

