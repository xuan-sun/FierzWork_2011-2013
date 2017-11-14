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
#include	 <math.h>
using		 namespace std;

struct Event
{
  double Erecon;
  int side;
  int type;
  int eventNum;
  double time;
  double tE;
  double tW;
  int pid;
  int timeFlag;
  double xEastPos;
  double yEastPos;
  double xWestPos;
  double yWestPos;
};

// forward declarations for useful functions
vector < pair <string,int> >  LoadOctetList(TString fileName, vector < pair <string,int> > pairs);
vector < TChain* > GetChainsOfRuns(vector < pair <string,int> > octetList, TString dataPath);

// global event number to be saved to super sum hist
double numberOfEventsInOctet;

// these are actual beta run indices
const int index_A2 = 0;
const int index_A5 = 1;
const int index_A7 = 2;
const int index_A10 = 3;
const int index_B2 = 4;
const int index_B5 = 5;
const int index_B7 = 6;
const int index_B10 = 7;

// these are the background runs
// they correspond to beta run index + 8
const int index_A1 = 8;
const int index_A4 = 9;
const int index_A9 = 10;
const int index_A12 = 11;
const int index_B1 = 12;
const int index_B4 = 13;
const int index_B9 = 14;
const int index_B12 = 15;

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
    cout << "(executable) (run type)" << endl;
    return 0;
  }

  int runType = atoi(argv[1]);

  TCanvas *C = new TCanvas("canvas", "canvas", 800, 400);

  int index = -1;
  TString fileName;
  if(runType == 8)
  {
    index = index_A1;
    fileName = "runType_A1_BGTree_allTypes.root";
  }
  else if(runType == 9)
  {
    index = index_A4;
    fileName = "runType_A4_BGTree_allTypes.root";
  }
  else if(runType == 10)
  {
    index = index_A9;
    fileName = "runType_A9_BGTree_allTypes.root";
  }
  else if(runType == 11)
  {
    index = index_A12;
    fileName = "runType_A12_BGTree_allTypes.root";
  }
  else if(runType == 12)
  {
    index = index_B1;
    fileName = "runType_B1_BGTree_allTypes.root";
  }
  else if(runType == 13)
  {
    index = index_B4;
    fileName = "runType_B4_BGTree_allTypes.root";
  }
  else if(runType == 14)
  {
    index = index_B9;
    fileName = "runType_B9_BGTree_allTypes.root";
  }
  else if(runType == 15)
  {
    index = index_B12;
    fileName = "runType_B12_BGTree_allTypes.root";
  }
  else
  {
    cout << "Value of index doesn't match a background run. Exiting..." << endl;
    return 0;
  }

  TFile f(fileName, "RECREATE");

  cout << "Run type input used to get index... " << index << endl;
  cout << "Creating file name: " << fileName.Data() << endl;

  vector < pair <string,int> > octetIndices;
  vector < TChain* > runsChains;
  for(int octNb = 0; octNb < 60; octNb++)
  {
    if(octNb == 9 || octNb == 59)
    {
      continue;
    }
    octetIndices = LoadOctetList(TString::Format("%s/octet_list_%i.dat", "OctetLists", octNb), octetIndices);
    runsChains = GetChainsOfRuns(octetIndices, "/mnt/Data/xuansun/replay_pass3_FINALCAL/");
  }

  Event* evt = new Event;
  runsChains[index]->SetBranchAddress("EvtN", &evt->eventNum);
  runsChains[index]->SetBranchAddress("Time", &evt->time);
  runsChains[index]->SetBranchAddress("TimeE", &evt->tE);
  runsChains[index]->SetBranchAddress("TimeW", &evt->tW);
  runsChains[index]->SetBranchAddress("Side", &evt->side);
  runsChains[index]->SetBranchAddress("Type", &evt->type);
  runsChains[index]->SetBranchAddress("Erecon", &evt->Erecon);
  runsChains[index]->SetBranchAddress("PID", &evt->pid);
  runsChains[index]->SetBranchAddress("badTimeFlag", &evt->timeFlag);
  runsChains[index]->GetBranch("xE")->GetLeaf("center")->SetAddress(&evt->xEastPos);
  runsChains[index]->GetBranch("yE")->GetLeaf("center")->SetAddress(&evt->yEastPos);
  runsChains[index]->GetBranch("xW")->GetLeaf("center")->SetAddress(&evt->xWestPos);
  runsChains[index]->GetBranch("yW")->GetLeaf("center")->SetAddress(&evt->yWestPos);

  Event* evtWrite = new Event;
  TTree *tFill = new TTree("Background", "Background");
  tFill->Branch("Erecon", &evtWrite->Erecon, "Erecon/D");
  tFill->Branch("EvtN", &evtWrite->eventNum, "EvtN/I");
  tFill->Branch("Time", &evtWrite->time, "Time/D");
  tFill->Branch("TimeE", &evtWrite->tE, "TimeE/D");
  tFill->Branch("TimeW", &evtWrite->tW, "TimeW/D");
  tFill->Branch("Side", &evtWrite->side, "Side/I");
  tFill->Branch("Type", &evtWrite->type, "Type/I");
  tFill->Branch("PID", &evtWrite->pid, "PID/I");
  tFill->Branch("badTimeFlag", &evtWrite->timeFlag, "badTimeFlag/I");

  for(unsigned int i = 0; i < runsChains[index]->GetEntries(); i++)
  {
    if(i % 100000 == 0)
      cout << "Copied entry number " << i << "/" << runsChains[index]->GetEntriesFast() <<  endl;

    runsChains[index]->GetEntry(i);

    if(evt->pid == 1 && evt->type < 4 && evt->Erecon >= 0 && evt->timeFlag == 0)
    {
      evtWrite->Erecon = evt->Erecon;
      evtWrite->eventNum = evt->eventNum;
      evtWrite->time = evt->time;
      evtWrite->tE = evt->tE;
      evtWrite->tW = evt->tW;
      evtWrite->side = evt->side;
      evtWrite->type = evt->type;
      evtWrite->pid = evt->pid;
      evtWrite->timeFlag = evt->timeFlag;

      tFill->Fill();
    }
  }

  f.Write();
  f.Close();



  cout << "-------------- End of Program ---------------" << endl;

  return 0;
}

vector < pair <string,int> >  LoadOctetList(TString fileName, vector <pair <string,int> > pairs)
{
  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  string runType;
  int runIndex;
//  vector <pair <string, int> > pairs;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> runType >> runIndex;
      if(runType == "A2")
      {
        pairs.push_back(make_pair("A2", runIndex));
      }
      else if(runType == "A5")
      {
        pairs.push_back(make_pair("A5", runIndex));
      }
      else if(runType == "A7")
      {
        pairs.push_back(make_pair("A7", runIndex));
      }
      else if(runType == "A10")
      {
        pairs.push_back(make_pair("A10", runIndex));
      }
      else if(runType == "B2")
      {
        pairs.push_back(make_pair("B2", runIndex));
      }
      else if(runType == "B5")
      {
        pairs.push_back(make_pair("B5", runIndex));
      }
      else if(runType == "B7")
      {
        pairs.push_back(make_pair("B7", runIndex));
      }
      else if(runType == "B10")
      {
        pairs.push_back(make_pair("B10", runIndex));
      }
      // Save the background runs as well
      else if(runType == "A1")
      {
	pairs.push_back(make_pair("A1", runIndex));
      }
      else if(runType == "A4")
      {
        pairs.push_back(make_pair("A4", runIndex));
      }
      else if(runType == "A9")
      {
        pairs.push_back(make_pair("A9", runIndex));
      }
      else if(runType == "A12")
      {
        pairs.push_back(make_pair("A12", runIndex));
      }
      else if(runType == "B1")
      {
        pairs.push_back(make_pair("B1", runIndex));
      }
      else if(runType == "B4")
      {
        pairs.push_back(make_pair("B4", runIndex));
      }
      else if(runType == "B9")
      {
        pairs.push_back(make_pair("B9", runIndex));
      }
      else if(runType == "B12")
      {
        pairs.push_back(make_pair("B12", runIndex));
      }
    }
  }

  infile.close();       // ensure you're closing the read-in file.

  return pairs;
}

vector < TChain* > GetChainsOfRuns(vector < pair <string,int> > octetList, TString dataPath)
{
  vector < TChain* > runs;

  for(unsigned int i = 0; i < 16; i++)
  {
    runs.push_back(new TChain("pass3"));
  }

  for(unsigned int l = 0; l < octetList.size(); l++)
  {
    if(octetList[l].first == "A2")
    {
      runs[index_A2] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A5")
    {
      runs[index_A5] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A7")
    {
      runs[index_A7] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A10")
    {
      runs[index_A10] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B2")
    {
      runs[index_B2] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B5")
    {
      runs[index_B5] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B7")
    {
      runs[index_B7] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B10")
    {
      runs[index_B10] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A1")
    {
      runs[index_A1] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A4")
    {
      runs[index_A4] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A9")
    {
      runs[index_A9] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A12")
    {
      runs[index_A12] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B1")
    {
      runs[index_B1] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B4")
    {
      runs[index_B4] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B9")
    {
      runs[index_B9] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B12")
    {
      runs[index_B12] -> Add(TString::Format("%s/replay_pass3_%i.root", dataPath.Data(), octetList[l].second));
    }
  }

  return runs;
}

