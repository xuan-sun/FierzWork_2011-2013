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

#define		GEOM	"2012-2013"

using            namespace std;

// Fundamental constants that get used
const double m_e = 511.00;                                              ///< electron mass, keV/c^2
double NDF = -1;

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void FillArrays(TString fileName, int flag);

struct entry
{
  string octNb;
  double avg_mE;
  double chisquared;
  double ndf;
  double chisquaredperndf;
  double prob;
  double b_minuitFit;
  double bErr_minuitFit;
  double binMin;
  double EMin;
  double binMax;
  double EMax;
  int fitMatrixStatus;
};

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  TCanvas *C = new TCanvas("canvas", "canvas");
//  C -> Divide(2,1);
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  FillArrays(Form("positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_%s_binWindowVariations_individualOctets.txt", GEOM), 1);

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

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString xAxis, TString yAxis, TString command)
{
  C->cd(canvasIndex);
  gPlot->SetTitle(title);
  gPlot->GetXaxis()->SetTitle(xAxis);
  gPlot->GetXaxis()->CenterTitle();
  gPlot->GetYaxis()->SetTitle(yAxis);
  gPlot->GetYaxis()->CenterTitle();
  C->SetLogy();

  gPlot->SetMarkerStyle(21);
  gPlot->SetMarkerSize(0.75);
  gPlot->SetMarkerColor(styleIndex);
  gPlot->SetLineColor(styleIndex);

  gPlot->Draw(command);

  C->Update();

}

void FillArrays(TString fileName, int flag)
{
  ofstream outfile;
  outfile.open(Form("positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_%s_binWindowVariations_individualOctetsSummary.txt", GEOM), ios::app);

  vector <TH1D*> hTemp;

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

  double variableHold = 0;

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
		>> evt.prob
		>> evt.b_minuitFit
		>> evt.bErr_minuitFit
		>> evt.binMin
		>> evt.EMin
		>> evt.binMax
		>> evt.EMax
		>> evt.fitMatrixStatus;
      counter++;

      if(evt.EMin == variableHold)
      {	// extracts histogram at end of vector
        hTemp[hTemp.size() - 1]->Fill(evt.b_minuitFit);
      }
      if(evt.EMin > variableHold)
      {
	// print out the previous bin's RMS and mean.
	if(hTemp.size() != 0)
	{
	  cout << "Bin = " << evt.binMin << ": mu = " << hTemp[hTemp.size() - 1]->GetMean()
	       << ", RMS = " << hTemp[hTemp.size() - 1]->GetRMS() << endl;
	  outfile << "2011-2012" << "\t"
		  << "ALL" << "\t"
		  << evt.binMin << "\t"
	 	  << evt.EMin << "\t"
		  << evt.binMax << "\t"
		  << evt.EMax << "\t"
		  << hTemp[hTemp.size() - 1]->GetMean() << "\t"
		  << hTemp[hTemp.size() - 1]->GetRMS() << "\n";
	}

	// create a new histogram to hold this bin's b fit values
	hTemp.push_back(new TH1D(Form("hist_b_bins_%f", evt.binMin), Form("hist_b_bins_%f", evt.binMin), 200, -1, 1));

        variableHold = evt.EMin;
      }

    }
    if(infile1.eof() == true)
    {
      break;
    }
  }

  outfile.close();

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

