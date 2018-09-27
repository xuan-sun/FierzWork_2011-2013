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
#include	 <math.h>
#include	 <TMath.h>
#include	 <TRandom3.h>

using		 namespace std;

struct Event
{
  int PID;
  double primMomentum[3];
  double primPosition[3];
  int side;
  int type;
  double Erecon;
  double EreconSmeared;
  double Evis[2];
  double Edep[2];
  double EdepQ[2];
  double Eprim;
  double AsymWeight;
  double time[2];
  double MWPCEnergy[2];
  double MWPCPos[2];
  double MWPCPosAdjusted[2];
  double ScintPos[2];
  double ScintPosAdjusted[2];
  double PMT_Evis[16];
};

// forward declarations for useful functions
vector < TChain* > GetChainsOfRuns(vector < pair <string,int> > octetList, TString dataPath);
vector < vector < TH1D* > > CreateRateHistograms(vector <TChain*> runsChains);

// implementing the 2011-2012 error envelope
void ErrorEnvelope_2011();
double converter1top(double *x, double *par);
void LoadEnvelopeHistogram();

// variables related to running the 2011-2012 error envelopes.
TH1D *hEnvelope2011 = new TH1D("2011-2012", "2011-2012", 110, 0, 1100);
TF1* errEnv2011_top_1sigma;


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
  cout << "Executing program on input index " << indexNb << " ..." << endl;

  // create global variables we will use later
  LoadEnvelopeHistogram();
  ErrorEnvelope_2011();


  Event* evt = new Event;

  TFile *f = new TFile(Form("/mnt/Data/xuansun/analyzed_files/All_Twiddles_Are_Baseline/SimAnalyzed_2011-2012_Beta_paramSet_%i_0.root", indexNb), "UPDATE");

  TTree* tree = (TTree*)f->Get("SimAnalyzed");

  tree->SetBranchAddress("PID", &evt->PID);
  tree->SetBranchAddress("primMomentum", &evt->primMomentum);
  tree->SetBranchAddress("primPosition", &evt->primPosition);
  tree->SetBranchAddress("side", &evt->side);
  tree->SetBranchAddress("type", &evt->type);
  tree->SetBranchAddress("Erecon", &evt->Erecon);
  tree->SetBranchAddress("Evis", &evt->Evis);
  tree->SetBranchAddress("Edep", &evt->Edep);
  tree->SetBranchAddress("EdepQ", &evt->EdepQ);
  tree->SetBranchAddress("Eprim", &evt->Eprim);
  tree->SetBranchAddress("AsymWeight", &evt->AsymWeight);
  tree->SetBranchAddress("time", &evt->time);
  tree->SetBranchAddress("MWPCEnergy", &evt->MWPCEnergy);
  tree->SetBranchAddress("MWPCPos", &evt->MWPCPos);
  tree->SetBranchAddress("MWPCPosAdjusted", &evt->MWPCPosAdjusted);
  tree->SetBranchAddress("ScintPos", &evt->ScintPos);
  tree->SetBranchAddress("ScintPosAdjusted", &evt->ScintPosAdjusted);
  tree->SetBranchAddress("PMT_Evis", &evt->PMT_Evis);


  TBranch *branch = tree->Branch("EreconSmeared", &evt->EreconSmeared, "EreconSmeared/D");

  TRandom3 *gRandom = new TRandom3();
  gRandom->SetSeed(0);


  for(int i = 0; i < tree->GetEntries(); i++)
  {
    tree->GetEntry(i);
    evt->EreconSmeared = gRandom->Gaus(evt->Erecon, errEnv2011_top_1sigma->Eval(evt->Erecon));
    branch->Fill();
  }

  tree->Write();

  delete f;


  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}


void ErrorEnvelope_2011()
{
  errEnv2011_top_1sigma = new TF1("2011-2012_error_envelope_top1sigma", converter1top, 0, 1050, 0);
}

double converter1top(double *x, double *par)
{
  double energy = x[0];
  return 1*hEnvelope2011->GetBinContent(hEnvelope2011->FindBin(energy));
}

void LoadEnvelopeHistogram()
{
  TString fileName = "/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/MB_errorEnvelopes_publication/MB_ErrorEnvelope_2011-2012.txt";

  int binNum = 0;
  double energy = 0;
  double error = 0;
  double errorHigh = 0;
  double errorLow = 0;

  // opens the file named above
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
      bufstream1 >> binNum
                 >> energy
                 >> error
                 >> errorHigh
                 >> errorLow;

      hEnvelope2011->SetBinContent(hEnvelope2011->FindBin(energy), errorHigh);
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Done filling in data from " << fileName.Data() << endl;

}


