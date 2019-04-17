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
#include         <TMinuit.h>
#include         <TCut.h>
#include	 <TLeaf.h>


#define         RADIALCUTLOW    0
#define         RADIALCUTHIGH   0.049   //these need to be done in m for simulations
#define         RADIALCUTLOW_NAME	0
#define         RADIALCUTHIGH_NAME	49   //these can be whatever, just using it for file name

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

  double mwpcPosE[3];
  double mwpcPosW[3];
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
  cout << "Input index: " << indexNb << endl;

  TString geom = "2011-2012";

  TString dataFilePath = Form("/mnt/data2/xuansun/analyzed_files/%s_geom_twiddles/TwiddledSimFiles_A_0_b_0_matchingParamSet_19_noStatDependence/SimAnalyzed_%s_Beta_paramSet_%i_0.root", geom.Data(), geom.Data(), indexNb);

  // checks if the file exists because we threw out a bunch of twiddles
  if(gSystem->AccessPathName(dataFilePath))
  {
    cout << "File does not exist. Exiting program." << endl;
    return 0;
  }
  else
  {
    cout << "File exists. Continuing..." << endl;
  }


  TCut positionCut = Form("(MWPCPosE[0]*MWPCPosE[0] + MWPCPosE[1]*MWPCPosE[1] < 0.063*0.063 && MWPCPosW[0]==0 && MWPCPosW[1]==0 || MWPCPosW[0]*MWPCPosW[0] + MWPCPosW[1]*MWPCPosW[1] < 0.063*0.063 && MWPCPosE[0]==0 && MWPCPosE[1]==0)");

  TFile fTwiddle(dataFilePath);
  TTree *t = (TTree*)fTwiddle.Get("SimAnalyzed");

  TH1D *h = new TH1D("Erecon_r49mm_endpointCorr", "Erecon_r49mm_endpointCorr", 120, 0, 1200);

  t->Draw("Erecon_corr_r49mm >> Erecon_r49mm_endpointCorr", "PID == 1 && Erecon > 0 && type == 0 && side < 2" && positionCut);

  TFile f(TString::Format("/mnt/data2/xuansun/analyzed_files/%s_geom_twiddles/TwiddledSimFiles_A_0_b_0_matchingParamSet_19_noStatDependence/Histograms/Hist_SimAnalyzed_%s_Beta_paramSet_%i_0_type0_radialCut_%i-%imm_endpointCorr.root", geom.Data(), geom.Data(), indexNb, RADIALCUTLOW_NAME, RADIALCUTHIGH_NAME), "RECREATE");
  // Begin processing the read in data now
  h->Write();

  cout << "About to plot the histogram now..." << endl;

  // YOU NEED THIS PLOT FUNCTION.
  // Or else ROOT doesn't seem to write the histogram at all..
  PlotHist(C, 1, 1, h, "", "");

  // Save our plot and print it out as a pdf.
//  C -> Print("fierz.pdf");
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

TH1D* CreateRateHistograms(TChain* runsChains, double percentMix)
{
  double radialCutLow = RADIALCUTLOW;   // measured in m
  double radialCutHigh = RADIALCUTHIGH; // measured in m

  radialCutLow = radialCutLow * sqrt(1.0 / 0.6);        // accounts for field expansion region at MWPC
  radialCutHigh = radialCutHigh * sqrt(1.0 / 0.6);      // only applies to simulations because we use MWPCPos

  TH1D* rateHistsErecon = new TH1D("Erecon blinded hist", "Erecon blinded hist", 120, 0, 1200);
  TRandom3 *engine = new TRandom3(0);

  // reminder: each evt[i] index corresponds to the indices noted at global scope.
  Event* evt = new Event;

  runsChains->SetBranchAddress("side", &evt->side);
  runsChains->SetBranchAddress("type", &evt->type);
  runsChains->SetBranchAddress("Erecon", &evt->Erecon);
  runsChains->SetBranchAddress("PID", &evt->pid);

  runsChains->GetBranch("MWPCPos")->GetLeaf("MWPCPosE")->SetAddress(&evt->mwpcPosE);
  runsChains->GetBranch("MWPCPos")->GetLeaf("MWPCPosW")->SetAddress(&evt->mwpcPosW);

  int totalEventNum = 0;

  for(unsigned int i = 0; i < runsChains->GetEntries(); i++)
  {
    if(engine -> Rndm() <= percentMix)
    {
      runsChains->GetEntry(i);
      if(evt->pid == 1 && evt->type == 0 && evt->Erecon >= 0 && evt->side < 2
        && (((pow(evt->mwpcPosE[0], 2.0) + pow(evt->mwpcPosE[1], 2.0) <= pow(radialCutHigh, 2.0))
        && (pow(evt->mwpcPosE[0], 2.0) + pow(evt->mwpcPosE[1], 2.0) >= pow(radialCutLow, 2.0))
        && (pow(evt->mwpcPosW[0], 2.0) + pow(evt->mwpcPosW[1], 2.0) == 0)
        && (pow(evt->mwpcPosW[0], 2.0) + pow(evt->mwpcPosW[1], 2.0) == 0) )
        || ((pow(evt->mwpcPosE[0], 2.0) + pow(evt->mwpcPosE[1], 2.0) == 0)
        && (pow(evt->mwpcPosE[0], 2.0) + pow(evt->mwpcPosE[1], 2.0) == 0)
        && (pow(evt->mwpcPosW[0], 2.0) + pow(evt->mwpcPosW[1], 2.0) <= pow(radialCutHigh, 2.0))
        && (pow(evt->mwpcPosW[0], 2.0) + pow(evt->mwpcPosW[1], 2.0) >= pow(radialCutLow, 2.0)))) )

      {
        rateHistsErecon->Fill(evt->Erecon);
        totalEventNum++;
      }
    }
    else
    {
      continue;
    }
  }

  rateHistsErecon->SetEntries(totalEventNum);

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
