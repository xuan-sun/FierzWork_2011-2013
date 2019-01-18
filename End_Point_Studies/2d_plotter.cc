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
void FillArrays(TString fileName, TH2D *h);

struct entry
{
  int indexNb;
  double chi2;
  double ndf;
  double chi2_ndf;
  double par0;
  double par0Err;
  double par1;
  double par1Err;
  double endpoint;
  double fitbinmin;
  double fitbinmax;
  double lowE;
  double highE;};

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  TCanvas *C = new TCanvas("canvas", "canvas");
  C->cd();
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  TH2D *h2 = new TH2D("bAndChi2", Form("b vs chi2: %s, %s", TYPE, GEOM), 40, -0.5, 1.5, 80, 0, 4);
  h2->GetXaxis()->SetTitle("b");
  h2->GetYaxis()->SetTitle("chi2/ndf");

  FillArrays(Form("endPointFits_GaussianTwiddles_%s_%s_Bins_%i-%i_index15.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN), h2);

  h2->Draw("COLZ");


  //prints the canvas with a dynamic TString name of the name of the file
//  C->Print("2D_bAndChi2_SymmetricTwiddles_newXuanFitter.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

void FillArrays(TString fileName, TH2D* h)
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
                >> evt.chi2
                >> evt.ndf
                >> evt.chi2_ndf
                >> evt.par0
                >> evt.par0Err
                >> evt.par1
                >> evt.par1Err
                >> evt.endpoint
                >> evt.fitbinmin
                >> evt.fitbinmax
                >> evt.lowE
                >> evt.highE;

      h->Fill(evt.endpoint, evt.chi2_ndf);
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

