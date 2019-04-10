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

#define		GEOM	"2011-2012"
#define		FITRANGELOW	0	// fit ranges in octet number for diff. geom's
#define		FITRANGEHIGH	60

using            namespace std;

vector <double> vecBinMin;
vector <TH1D*> vecHistOctets;
vector <TGraphErrors*> vecGraphOctets;

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
  double octNb;
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


  ofstream outfile;
  outfile.open(Form("positionCuts_0-49mm_endpointCorrected_withFullBlind_Feb2019_type0_%s_binWindowVariations_individualOctetsSummary.txt", GEOM), ios::app);
  for(unsigned int i = 0; i < vecHistOctets.size(); i++)
  {
    TF1 *fit1 = new TF1(Form("fit_%i", i), "[0]", FITRANGELOW, FITRANGEHIGH);
    vecGraphOctets[i]->Fit(fit1, "R");

    outfile << "2011-2012" << "\t"
            << "ALL" << "\t"
            << vecBinMin[i] << "\t"
            << vecHistOctets[i]->GetMean() << "\t"
            << vecHistOctets[i]->GetRMS() << "\t"
	    << fit1->GetParameter(0) << "\t"
	    << fit1->GetParError(0) << "\t"
	    << fit1->GetChisquare() << "\t"
	    << fit1->GetNDF() << "\t"
	    << fit1->GetChisquare() / fit1->GetNDF() << "\t"
	    << TMath::Prob(fit1->GetChisquare(), fit1->GetNDF()) << "\n";
  }
  outfile.close();


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
  vecHistOctets.push_back(new TH1D(Form("hist_b_bins_%i", 17), Form("hist_b_bins_%i", 17), 200, -1, 1));

  vector <double> octNb;
  vector <double> octNbErr;
  vector <double> bVal;
  vector <double> bErr;

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

  double variableHold = 165;	// note: this has to match the first input to fix a pesky off-by-1 error

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

      if(evt.EMin > variableHold)
      {
	// create a new histogram to hold this bin's b fit values
	vecHistOctets.push_back(new TH1D(Form("hist_b_bins_%f", evt.binMin - 1), Form("hist_b_bins_%f", evt.binMin - 1), 200, -1, 1));

	// create a new TGraph so we can do a straight line fit
        vecGraphOctets.push_back(new TGraphErrors(octNb.size(), &(octNb[0]), &(bVal[0]), &(octNbErr[0]), &(bErr[0])));

	// save the bin values so we can match them in printing later
        vecBinMin.push_back(evt.binMin - 1);

	// reset all the values for the next iteration.
        variableHold = evt.EMin;
	octNb.clear();
	octNbErr.clear();
	bVal.clear();
	bErr.clear();
      }

      // extracts histogram at end of vector and fills it
      vecHistOctets[vecHistOctets.size() - 1]->Fill(evt.b_minuitFit);

      // fill in the vectors that will make the TGraphErrors.
      octNb.push_back(evt.octNb);
      octNbErr.push_back(0.5);
      bVal.push_back(evt.b_minuitFit);
      bErr.push_back(evt.bErr_minuitFit);

    }
    if(infile1.eof() == true)
    {
      break;
    }
  }

  // these are also to address off-by-1 errors cause this method is tricky.
  vecBinMin.push_back(evt.binMin);
  vecGraphOctets.push_back(new TGraphErrors(octNb.size(), &(octNb[0]), &(bVal[0]), &(octNbErr[0]), &(bErr[0])));


  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

