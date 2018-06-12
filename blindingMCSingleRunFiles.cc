#include	 "BetaSpectrum.hh"

#include	 <iostream>
#include	 <fstream>
#include	 <TGaxis.h>
#include	 <sstream>
#include	 <TGraph.h>
#include	 <TGraphErrors.h>
#include	 <TCanvas.h>
#include	 <TApplication.h>
#include	 <stdlib.h>
#include	 <TF1.h>
#include	 <TH1.h>
#include	 <TProfile.h>
#include	 <TObjArray.h>
#include	 <TStyle.h>
#include	 <TMarker.h>
#include	 <math.h>
#include	 <TStyle.h>
#include	 <TPaveStats.h>
#include	 <TPaveText.h>
#include	 <vector>
#include	 <string.h>
#include	 <fstream>
#include	 <TROOT.h>
#include	 <TFile.h>
#include	 <TLegend.h>
#include         <TLegendEntry.h>
#include	 <time.h>
#include	 <TH2F.h>
#include         <assert.h>
#include	 <string>
#include	 <TRandom.h>
#include 	 <TTree.h>
#include	 <TChain.h>
#include	 <TVector.h>
#include	 <vector>
#include	 <utility>
#include	 <TLeaf.h>
#include	 <TRandom3.h>
using		 namespace std;

struct Event
{
  double Erecon;
  int side;
  int type;
  int eventNum;
  double time[2];
  double tE;
  double tW;
  int pid;
  double xEastPos;
  double yEastPos;
  double xWestPos;
  double yWestPos;
};

// forward declarations for useful functions
TH1D* CreateRateHistograms(TChain* runsChains, double percentMix);

// plotting functions
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command);


// these are actual beta run indices
const int index_A2 = 0;
const int index_A5 = 1;
const int index_A7 = 2;
const int index_A10 = 3;
const int index_B2 = 4;
const int index_B5 = 5;
const int index_B7 = 6;
const int index_B10 = 7;

// Used for visualization, keeps the graph on screen.
//TApplication plot_program("FADC_readin",0,0,0,0);

//-------------------------------------------------//
//------------ Start of Program -------------------//
//-------------------------------------------------//

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (index #)" << endl;
    return 0;
  }

  // creating canvas for plotting
  TCanvas *C = new TCanvas("canvas", "canvas", 800, 400);

  // read in arguments.
  Int_t indexNb = atoi(argv[1]);

  // Points TChains at the run files idenified in the octet lists above
  TChain* runFiles_base = new TChain("SimAnalyzed");
  runFiles_base->Add(Form("/mnt/Data/xuansun/analyzed_files/A_0_b_0/SimAnalyzed_2011-2012_Beta_paramSet_100_%i.root", indexNb));

  TChain* runFiles_fierz = new TChain("SimAnalyzed");;
  runFiles_fierz->Add(Form("/mnt/Data/xuansun/analyzed_files/A_0_b_inf/SimAnalyzed_2011-2012_Beta_paramSet_100_%i.root", indexNb));

  // read in our random mixing seed so I stay pretty blind.
  double s0, s1, s2, s3, s4, s5, s6, s7, s8, s9;

  string buf1;
  ifstream infile1;
  cout << "The file being opened is: " << "randomMixingSeeds.txt" << endl;
  infile1.open("ExtractedHistograms/randomMixingSeeds.txt");

  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << "randomMixingSeeds.txt" << endl;

  while(true)
  {
    getline(infile1, buf1);
    istringstream bufstream1(buf1);

    if(!infile1.eof())
    {
      bufstream1 >> s0 >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8 >> s9;
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Loaded mixing seeds" << endl;

  // load all the histograms of east and west, turn them into rates.
  // ALWAYS USE S0 FOR MIXING.
  TH1D* rates_base = CreateRateHistograms(runFiles_base, 1 - s0);
  TH1D* rates_fierz = CreateRateHistograms(runFiles_fierz, s0);

  cout << "Created vectors of rate histograms" << endl;

  // Sum the two files together.
  rates_base->Add(rates_base, rates_fierz, 1, 1);

  cout << "Finished mixing rate histograms." << endl;

  TFile f(TString::Format("/mnt/Data/xuansun/analyzed_files/A_0_b_0/BLIND_SimAnalyzed_2011-2012_Beta_paramSet_100_%i_type0.root", indexNb), "RECREATE");
  // Begin processing the read in data now
  rates_base->Write();

  cout << "About to plot the histogram now..." << endl;

  // YOU NEED THIS PLOT FUNCTION.
  // Or else ROOT doesn't seem to write the histogram at all..
  PlotHist(C, 1, 1, rates_base, "", "");

  // Save our plot and print it out as a pdf.
//  C -> Print("fierz.pdf");
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

TH1D* CreateRateHistograms(TChain* runsChains, double percentMix)
{
  TH1D* rateHistsErecon = new TH1D("Erecon blinded hist", "Erecon blinded hist", 120, 0, 1200);
  TRandom3 *engine = new TRandom3(0);

  // reminder: each evt[i] index corresponds to the indices noted at global scope.
  Event* evt = new Event;

  runsChains->SetBranchAddress("side", &evt->side);
  runsChains->SetBranchAddress("type", &evt->type);
  runsChains->SetBranchAddress("Erecon", &evt->Erecon);
  runsChains->SetBranchAddress("PID", &evt->pid);

  for(unsigned int i = 0; i < runsChains->GetEntries(); i++)
  {
    if(engine -> Rndm() <= percentMix)
    {
      runsChains->GetEntry(i);
      if(evt->pid == 1 && evt->type == 0 && evt->Erecon >= 0 && evt->side < 2)
      {
        rateHistsErecon->Fill(evt->Erecon);
      }
    }
    else
    {
      continue;
    }
  }

  delete engine;
  return rateHistsErecon;
}


void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Energy (KeV)");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("Rate (mHz/KeV)");
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(0, 0.000004);

  if(styleIndex == 1)
  {
    hPlot -> SetFillColor(46);
    hPlot -> SetFillStyle(3004);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 2)
  {
    hPlot -> SetFillColor(38);
    hPlot -> SetFillStyle(3005);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 3)
  {
    hPlot -> SetFillColor(29);
//    hPlot -> SetFillStyle(3005);
    hPlot -> SetFillStyle(3001);
  }

  hPlot -> Draw(command);
  C -> Update();
}
