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
void AddBranchToTreeInFile(TString fName, double fittedOctetEndpoint, double meanEndpoint);
vector < pair <int, int> > ReadInFileOfRunNumbers(TString fileName);
vector < pair <int, double> > ReadInOctetEndpoints(TString fileName, vector < pair <int, double> > updateThis);

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


  vector < pair <int, int> > runNumbers = ReadInFileOfRunNumbers(Form("bFits_OctetsUsed_RunNumbers.txt"));

  vector < pair <int, double> > fittedEndpointsForOctets;
  fittedEndpointsForOctets = ReadInOctetEndpoints(Form("endPointFits_noGainCorrection_ssDataHists_2011-2012_radialCut_0-30mm.txt"), fittedEndpointsForOctets);
  fittedEndpointsForOctets = ReadInOctetEndpoints(Form("endPointFits_noGainCorrection_ssDataHists_2012-2013_radialCut_0-30mm.txt"), fittedEndpointsForOctets);

  cout << "runNumbers.size() = " << runNumbers.size() << endl;

  double meanEndpoint = -1;
  double octetEndpoint = -1;
  for(unsigned int i = 0; i < runNumbers.size(); i++)
  {
    // use 2011-2012 octet fitted endpoints; 0 < r < 30mm radial cut.
    if(runNumbers[i].first >= 0 && runNumbers[i].first <= 59)
    {
      meanEndpoint = 788.8;
    }
    else if(runNumbers[i].first >= 80 && runNumbers[i].first <= 121)
    {
      meanEndpoint = 785.8;
    }

    for(unsigned int j = 0; j < fittedEndpointsForOctets.size(); j++)
    {
      if(runNumbers[i].first == fittedEndpointsForOctets[j].first)
      {
        octetEndpoint = fittedEndpointsForOctets[j].second;
        break;
      }
    }

    cout << "octetEndpoint = " << octetEndpoint << "; meanEndpoint = " << meanEndpoint << endl;

    AddBranchToTreeInFile(Form("/mnt/Data/xuansun/replay_pass3_FINALCAL/replay_pass3_%i.root", runNumbers[i].second), octetEndpoint, meanEndpoint);
  }


//  C -> Print("fierz.pdf");
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

void AddBranchToTreeInFile(TString fName, double fittedOctetEndpoint, double meanEndpoint)
{
  TFile *f = new TFile(fName.Data(), "UPDATE");
  TTree *t = (TTree*)f->Get("pass3");

  double erecon;
  t->SetBranchAddress("Erecon", &erecon);

  double ereconCorr_r30mm;
  TBranch *b = t->Branch("Erecon_corr_r30mm", &ereconCorr_r30mm, "ereconCorr_r30mm/D");

  for(long int i = 0; i < t->GetEntries(); i++)
  {
    t->GetEntry(i);
    ereconCorr_r30mm = erecon*(meanEndpoint / fittedOctetEndpoint);
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

vector < pair <int, double> > ReadInOctetEndpoints(TString fileName, vector < pair <int, double> > updateThis)
{
  int octNb = -1;
  double endpoint = -1;
  double endpointErr = -1;
  double fakeGain = -1;

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
      bufstream >> octNb >> endpoint >> endpointErr >> fakeGain;

      updateThis.push_back(make_pair(octNb, endpoint));
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  return updateThis;
}
