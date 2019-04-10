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

using		 namespace std;

// forward declarations for useful functions
void AddBranchToTreeInFile(TString fName);
vector < pair <int, int> > ReadInFileOfRunNumbers(TString fileName);

// Used for visualization, keeps the graph on screen.
//TApplication plot_program("FADC_readin",0,0,0,0);

//-------------------------------------------------//
//------------ Start of Program -------------------//
//-------------------------------------------------//

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  // creating canvas for plotting
  TCanvas *C = new TCanvas("canvas", "canvas", 800, 400);


  vector < pair <int, int> > runNumbers = ReadInFileOfRunNumbers(Form("bFits_OctetsUsed_RunNumbers_2012-2013.txt"));

  cout << "runNumbers.size() = " << runNumbers.size() << endl;

  for(unsigned int i = 0; i < runNumbers.size(); i++)
  {
    AddBranchToTreeInFile(Form("/mnt/Data/xuansun/replay_pass3_FINALCAL/replay_pass3_%i.root", runNumbers[i].second));
  }


//  C -> Print("fierz.pdf");
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

void AddBranchToTreeInFile(TString fName)
{
  TFile *f = new TFile(fName.Data(), "UPDATE");
  TTree *t = (TTree*)f->Get("pass3");

  double erecon;
  t->SetBranchAddress("Erecon_corr_r49mm", &erecon);

  double erecon_tinStitch;
  TBranch *b = t->Branch("Erecon_corr_r49mm_Sn113Stitch_try3", &erecon_tinStitch, "erecon_tinStitch/D");

  for(long int i = 0; i < t->GetEntries(); i++)
  {
    t->GetEntry(i);

    if(erecon > 0 && erecon < 368.5)
    {   // scaling pre-Sn energies to the shifted (miscalibrated) Sn energy
      erecon_tinStitch = erecon*(371.6 / 368.5);
    }
    else if(erecon > 368.5 && erecon < 1000)
    {   // scales all values post-Sn to stretch to new-Sn and endpoint (787keV, with resolution effects)
      erecon_tinStitch = erecon*(787.0 / 792.9) + 5.8;
    }

    b->Fill();
  }

  t->Write();
  delete f;

  cout << "Done updating runNumber = " << fName.Data() << endl;
}

vector < pair <int, int> > ReadInFileOfRunNumbers(TString fileName)
{
  vector < pair <int, int> > runs;

  int octNb = -1;
  int runNb = -1;
  string runType;

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
      bufstream >> octNb >> runNb >> runType;

      runs.push_back(make_pair(octNb, runNb));
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  return runs;
}

