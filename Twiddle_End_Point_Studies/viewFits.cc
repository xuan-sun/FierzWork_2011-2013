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
#define		RADIALCUTLOW	0
#define		RADIALCUTHIGH	49

using            namespace std;

// Fundamental constants that get used
const double m_e = 511.00;                      ///< electron mass, keV/c^2
double NDF = 47;				// taken from the .txt file. All twiddles fit the same range so same NDF

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command, double maxBinContents);
void FillArrays(TString fileName, TH1D *h, int hFillOption);

struct entry
{
  int indexNb;
  double endpoint;
  double endpointErr;
  double gainFactor;

  int index;
  double avg_mE;
  double functionMin;
  double ndf;
  double reducedchi2;
  double b;
  double bErr;
  int hold1;
  int hold2;
  int status;

};

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

//  int option = atoi(argv[1]);

  TCanvas *C = new TCanvas("canvas", "canvas");
  C->cd();
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  TH1D* h = new TH1D("endpoints", "twiddles", 100, 780, 800);
  FillArrays(Form("endPointFits_noGainCorrection_baselineHistogramsForFitting_2011-2012_radialCut_0-49mm.txt"), h, 1);
//  FillArrays(Form("asymmetric_gaussianTwiddles_noBlinding_endpointCorr_fast-newXuanFitter_bFitsForSyst_type0_2011-2012_Bins_17-65_index19_radialCut_0-49mm.txt"), h, 1);
//  FillArrays(Form("endPointFits_gainAlreadyApplied_asymmTwiddledSpectra_index19_2011-2012_radialCut_0-49mm.txt"), h, 1);
//  FillArrays(Form("endPointFits_noGainCorrection_asymmTwiddledSpectra_index19_2011-2012_radialCut_0-30mm_take2.txt"), h, 1);

  int max = h->GetMaximum();

  PlotHist(C, 2, 1, h, Form("endpoints, twiddles, index 19: %s, %s, %i-%imm", TYPE, GEOM, RADIALCUTLOW, RADIALCUTHIGH), "endpoints (keV)", "N", "", max);

  //prints the canvas with a dynamic TString name of the name of the file
//  C->Print("viewNewXuanFitter_SymmetricTwiddles_finerGrid.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

void FillArrays(TString fileName, TH1D* h, int hFillOption)
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
		>> evt.endpoint
		>> evt.endpointErr;
      h->Fill(evt.endpoint);

/*
      bufstream >> evt.index
		>> evt.avg_mE
		>> evt.functionMin
      		>> evt.ndf
		>> evt.reducedchi2
		>> evt.b
      		>> evt.bErr
		>> evt.hold1
		>> evt.hold2
		>> evt.status;
      h->Fill(evt.b, 1/sqrt(evt.reducedchi2));
*/
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

