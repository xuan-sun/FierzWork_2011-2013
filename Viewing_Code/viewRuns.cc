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
#include 	 <TRandom3.h>
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
vector < pair <string,int> >  LoadOctetList(TString fileName);
vector < TChain* > GetChainsOfRuns(vector < pair <string,int> > octetList, TString dataPath);

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
/*
//This generates the numbers used for b=0 MC "mixing".
  TRandom3* engine = new TRandom3(0);
  ofstream outfile;
  outfile.open("randomMixingSeeds.txt", ios::app);
  outfile << engine->Rndm()/5.0 << "\t"
          << engine->Rndm()/5.0 << "\t"
          << engine->Rndm()/5.0 << "\t"
          << engine->Rndm()/5.0 << "\t"
	  << engine->Rndm()/5.0 << "\t"
          << engine->Rndm()/5.0 << "\t"
          << engine->Rndm()/5.0 << "\t"
          << engine->Rndm()/5.0 << "\t"
          << engine->Rndm()/5.0 << "\t"
          << engine->Rndm()/5.0 << "\n";
  outfile.close();
  return 0;
*/

  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (octet #)" << endl;
    return 0;
  }

  // creating canvas for plotting
  TCanvas *C = new TCanvas("canvas", "canvas", 800, 400);

  // read in arguments.
  Int_t octNb = atoi(argv[1]);

  // Reads in the octet list and saves the run files indices corresponding to an octet number
  vector < pair <string,int> > octetIndices = LoadOctetList(TString::Format("%s/octet_list_%i.dat", "OctetLists", octNb));
  // Points TChains at the run files idenified in the octet lists above
  vector < TChain* > runsChains = GetChainsOfRuns(octetIndices, "/mnt/Data/xuansun/replay_pass3_FINALCAL/");


  // reminder: each evt[i] index corresponds to the indices noted at global scope.
  vector <Event*> evt;

  vector <TH1D*> histos;
  histos.push_back(new TH1D("A2", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("A5", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("A7", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("A10", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("B2", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("B5", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("B7", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("B10", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("A1", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("A4", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("A9", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("A12", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("B1", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("B4", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("B9", "Erecon", 120, 0, 1200));
  histos.push_back(new TH1D("B12", "Erecon", 120, 0, 1200));

  for(unsigned int i = 0; i < runsChains.size(); i++)
  {
    evt.push_back(new Event);

    runsChains[i]->SetBranchAddress("EvtN", &evt[i]->eventNum);
    runsChains[i]->SetBranchAddress("Time", &evt[i]->time);
    runsChains[i]->SetBranchAddress("TimeE", &evt[i]->tE);
    runsChains[i]->SetBranchAddress("TimeW", &evt[i]->tW);
    runsChains[i]->SetBranchAddress("Side", &evt[i]->side);
    runsChains[i]->SetBranchAddress("Type", &evt[i]->type);
    runsChains[i]->SetBranchAddress("Erecon", &evt[i]->Erecon);
    runsChains[i]->SetBranchAddress("PID", &evt[i]->pid);
    runsChains[i]->SetBranchAddress("badTimeFlag", &evt[i]->timeFlag);

    // this additional syntax is needed to get the right leaf inside branch inside tree named "pass3"
    runsChains[i]->GetBranch("xE")->GetLeaf("center")->SetAddress(&evt[i]->xEastPos);
    runsChains[i]->GetBranch("yE")->GetLeaf("center")->SetAddress(&evt[i]->yEastPos);
    runsChains[i]->GetBranch("xW")->GetLeaf("center")->SetAddress(&evt[i]->xWestPos);
    runsChains[i]->GetBranch("yW")->GetLeaf("center")->SetAddress(&evt[i]->yWestPos);
  }


  for(unsigned int j = 0; j < runsChains.size(); j++)
  {
    for(unsigned int i = 0; i < runsChains[j]->GetEntriesFast(); i++)
    {
      runsChains[j]->GetEntry(i);
      if(evt[j]->pid == 1 && evt[j]->type < 4 && evt[j]->Erecon >= 0 && evt[j]->timeFlag == 0)
      {
        histos[j]->Fill(evt[j]->Erecon);
      }
    }
  }



  C -> Divide(4,4);

  // for some mega unclear reason, this for-loop doesn't work...
/*  for(int i = 0; i < 4; i++)
  {
    PlotHist(C, 1, i, histos[i], "", "");
  }
*/
  // ... but individually printing all the histograms does.
  PlotHist(C, 1, 1, histos[0], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 2, histos[1], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 3, histos[2], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 4, histos[3], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 5, histos[4], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 6, histos[5], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 7, histos[6], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 8, histos[7], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 9, histos[8], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 10, histos[9], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 11, histos[10], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 12, histos[11], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 13, histos[12], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 14, histos[13], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 15, histos[14], Form("Octet %i", octNb), "");
  PlotHist(C, 1, 16, histos[15], Form("Octet %i", octNb), "");


  // Save our plot and print it out as a pdf.
  C -> Print(Form("RunsForOctet_%i_allTypes.pdf", octNb));
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

vector < pair <string,int> >  LoadOctetList(TString fileName)
{
  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  string runType;
  int runIndex;
  vector <pair <string, int> > pairs;

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

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Energy (KeV)");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("Counts");
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
