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

using            namespace std;

vector <TH1D*> vecHistTwiddles;
vector <double> vecBinMins;

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void FillArrays(TString fileName, int flag);

struct entry
{
  int index;
  double avg_mE;
  double chisquared;
  double ndf;
  double chisquaredperndf;
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

  for(int fitBinMin = 17; fitBinMin <= 41; fitBinMin++)
  {
    FillArrays(Form("positionCuts_0-49mm_noBlind_newXuanFitter_twiddleHists_bFit_type0_2011-2012_Bins_%i-65_endpointCorrected.txt", fitBinMin), fitBinMin);
  }

  cout << "Number of histograms filled: " << vecHistTwiddles.size() << endl;

  ofstream outfile;
  outfile.open(Form("twiddle_index19_binVariation_positionCuts_0-49mm_endpointCorrected_noBlind_type0_%s_summary.txt", GEOM), ios::app);
  for(unsigned int i = 0; i < vecHistTwiddles.size(); i++)
  {
    outfile << "2011-2012" << "\t"
            << "ALL TWIDDLES" << "\t"
            << vecBinMins[i] << "\t"
            << vecHistTwiddles[i]->GetMean() << "\t"
            << vecHistTwiddles[i]->GetRMS() << "\n";
  }
  outfile.close();


  //prints the canvas with a dynamic TString name of the name of the file
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

void FillArrays(TString fileName, int flag)
{
  vecHistTwiddles.push_back(new TH1D(Form("hist_b_lowBins_%i", flag), Form("hist_b_lowBins_%i", flag), 200, -1, 1));
  vecBinMins.push_back(flag);

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
      bufstream1 >> evt.index
		>> evt.avg_mE
		>> evt.chisquared
		>> evt.ndf
		>> evt.chisquaredperndf
		>> evt.b_minuitFit
		>> evt.bErr_minuitFit
		>> evt.binMin
		>> evt.EMin
		>> evt.binMax
		>> evt.EMax
		>> evt.fitMatrixStatus;
      counter++;

      vecHistTwiddles[vecHistTwiddles.size() - 1]->Fill(evt.b_minuitFit);

    }

    if(infile1.eof() == true)
    {
      break;
    }
  }


  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

