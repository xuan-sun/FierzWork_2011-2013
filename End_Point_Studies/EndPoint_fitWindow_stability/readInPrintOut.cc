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

using            namespace std;

vector <TH1D*> vecHistTwiddles;
vector <double> vecBinMins;

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void FillArrays(TString fileName, int flag);

struct entry
{
  int octNb;
  double endpoint;
  double endpointErr;
  double falseGain;
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
    FillArrays(Form("endPointFits_noCorrection_ssDataHists_type0_2012-2013_radialCut_0-49mm_Bins_%i-65.txt", fitBinMin), fitBinMin);
  }

  cout << "Number of histograms filled: " << vecHistTwiddles.size() << endl;

  ofstream outfile;
  outfile.open(Form("endPointFits_noCorrection_ssDataHists_type0_2012-2013_radialCut_0-49mm_Bins_summary.txt"), ios::app);
  for(unsigned int i = 0; i < vecHistTwiddles.size(); i++)
  {
    outfile << "2011-2012" << "\t"
            << "ALL_OCTETS" << "\t"
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
  vecHistTwiddles.push_back(new TH1D(Form("hist_b_lowBins_%i", flag), Form("hist_b_lowBins_%i", flag), 100, 760, 810));
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
      bufstream1 >> evt.octNb
		>> evt.endpoint
		>> evt.endpointErr
		>> evt.falseGain;
      counter++;

      vecHistTwiddles[vecHistTwiddles.size() - 1]->Fill(evt.endpoint);

    }

    if(infile1.eof() == true)
    {
      break;
    }
  }


  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

