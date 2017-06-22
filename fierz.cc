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
  double xEastPos;
  double yEastPos;
  double xWestPos;
  double yWestPos;
};

// forward declarations for useful functions
//void ExtractA(TString dataPath, Int_t octetNum, Int_t structType, TString analysisChoice, vector < pair <string, int> > octetList);
vector < pair <string,int> >  LoadOctetList(TString fileName);
vector < TChain* > GetChainsOfRuns(vector < pair <string,int> > octetList, TString dataPath);
vector < vector < TH1D* > > CreateRateHistograms(vector <TChain*> runsChains);
TH1D* CreateSuperSum(vector < vector < TH1D* > > sideRates);

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
TApplication plot_program("FADC_readin",0,0,0,0);

//-------------------------------------------------//
//------------ Start of Program -------------------//
//-------------------------------------------------//

int main(int argc, char* argv[])
{
  // creating canvas for plotting
  TCanvas *C = new TCanvas("canvas", "canvas", 800, 400);

  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (octet #)" << endl;
    return 0;
  }

  // read in arguments.
  Int_t octNb = atoi(argv[1]);

  // Reads in the octet list and saves the run files indices corresponding to an octet number
  vector < pair <string,int> > octetIndices = LoadOctetList(TString::Format("%s/octet_list_%i.dat", "OctetLists", octNb));
  // Points TChains at the run files idenified in the octet lists above
  vector < TChain* > runFiles = GetChainsOfRuns(octetIndices, "Data/Octet40/");
  // load all the histograms of east and west, turn them into rates.
  vector < vector < TH1D* > > rates = CreateRateHistograms(runFiles);

  for(int i = 0; i <= 1; i++)
  {
    for(int j = 0; j <= rates[i][index_B2]->GetNbinsX(); j++)
    {
      cout << "Side " << i << " has rates " << rates[i][index_B2]->GetBinContent(j) << " at bin " << j << endl;
    }
  }


  // Begin processing the read in data now
  TH1D* SS_Erecon = CreateSuperSum(rates);
  SS_Erecon->Scale(1000.0/10.0);	// creates mHz/KeV bins

  PlotHist(C, 1, 1, SS_Erecon, "", "");

  // Save our plot and print it out as a pdf.
  C -> Print("fierz.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

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

  for(unsigned int i = 0; i < octetList.size(); i++)
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

vector < vector < TH1D* > > CreateRateHistograms(vector <TChain*> runsChains)
{
  // reminder: each evt[i] index corresponds to the indices noted at global scope.
  vector <Event*> evt;

  vector < vector <TH1D*> > rateHists;

  vector <TH1D*> rateHistsEast;
  vector <TH1D*> rateHistsWest;

  for(unsigned int i = 0; i < runsChains.size(); i++)
  {
    evt.push_back(new Event);
    rateHistsEast.push_back(new TH1D(TString::Format("East Rate %i", i), "East Rate", 120, 0, 1200));
    rateHistsWest.push_back(new TH1D(TString::Format("West Rate %i", i), "West Rate", 120, 0, 1200));

    runsChains[i]->SetBranchAddress("EvtN", &evt[i]->eventNum);
    runsChains[i]->SetBranchAddress("Time", &evt[i]->time);
    runsChains[i]->SetBranchAddress("TimeE", &evt[i]->tE);
    runsChains[i]->SetBranchAddress("TimeW", &evt[i]->tW);
    runsChains[i]->SetBranchAddress("Side", &evt[i]->side);
    runsChains[i]->SetBranchAddress("Type", &evt[i]->type);
    runsChains[i]->SetBranchAddress("Erecon", &evt[i]->Erecon);
    runsChains[i]->SetBranchAddress("PID", &evt[i]->pid);

    // this additional syntax is needed to get the right leaf inside branch inside tree named "pass3"
    runsChains[i]->GetBranch("xE")->GetLeaf("center")->SetAddress(&evt[i]->xEastPos);
    runsChains[i]->GetBranch("yE")->GetLeaf("center")->SetAddress(&evt[i]->yEastPos);
    runsChains[i]->GetBranch("xW")->GetLeaf("center")->SetAddress(&evt[i]->xWestPos);
    runsChains[i]->GetBranch("yW")->GetLeaf("center")->SetAddress(&evt[i]->yWestPos);
  }

  vector <double> liveTimeEast;
  vector <double> liveTimeWest;

  for(unsigned int j = 0; j < runsChains.size(); j++)
  {
    for(unsigned int i = 0; i < runsChains[j]->GetEntriesFast(); i++)
    {
      runsChains[j]->GetEntry(i);
      if(evt[j]->pid == 1 && evt[j]->type < 4 && evt[j]->Erecon >= 0)
      {
        if(evt[j]->side == 0)
        {
          rateHistsEast[j]->Fill(evt[j]->Erecon);
        }
        else if(evt[j]->side == 1)
        {
          rateHistsWest[j]->Fill(evt[j]->Erecon);
        }
      }
      if(i == runsChains[j]->GetEntriesFast() - 1)
      {
        liveTimeEast.push_back(evt[j]->tE);
        liveTimeWest.push_back(evt[j]->tW);
      }
    }
  }


  // here we loop back over the event histograms in east and west
  // and scale them down by time of last event aka live time.
  for(unsigned int i = 0; i < rateHistsEast.size(); i++)
  {
    rateHistsEast[i]->Scale(1.0/liveTimeEast[i]);
  }

  for(unsigned int i = 0; i < rateHistsWest.size(); i++)
  {
    rateHistsWest[i]->Scale(1.0/liveTimeWest[i]);
  }

  rateHists.push_back(rateHistsEast);
  rateHists.push_back(rateHistsWest);

  return rateHists;

}

TH1D* CreateSuperSum(vector < vector < TH1D* > > sideRates)
{
  TH1D* hist = new TH1D("Super sum", "Super sum Erecon spectrum", 120, 0, 1200);

  // Do the background subtraction
  for(unsigned int j = 0; j <= 1; j++)
  {
    for(unsigned int i = 0; i < (sideRates[j].size())/2; i++)
    {
      sideRates[j][i]->Add(sideRates[j][i], sideRates[j][i+8], 1, -1);
    }
  }

  // sum the "like" histograms without any statistical weight
  TH1D* eastPlusRates = new TH1D("East Plus", "East Plus", 120, 0, 1200);
  eastPlusRates->Add(sideRates[0][index_A5]);
  eastPlusRates->Add(sideRates[0][index_A7]);
  eastPlusRates->Add(sideRates[0][index_B2]);
  eastPlusRates->Add(sideRates[0][index_B10]);
  eastPlusRates->Scale(1.0/4.0);

  TH1D* westPlusRates =new TH1D("West Plus", "West Plus", 120, 0, 1200);
  westPlusRates->Add(sideRates[1][index_A5]);
  westPlusRates->Add(sideRates[1][index_A7]);
  westPlusRates->Add(sideRates[1][index_B2]);
  westPlusRates->Add(sideRates[1][index_B10]);
  westPlusRates->Scale(1.0/4.0);

  TH1D* eastMinusRates = new TH1D("East Minus", "East Minus", 120, 0, 1200);
  eastMinusRates->Add(sideRates[0][index_A5]);
  eastMinusRates->Add(sideRates[0][index_A7]);
  eastMinusRates->Add(sideRates[0][index_B2]);
  eastMinusRates->Add(sideRates[0][index_B10]);
  eastMinusRates->Scale(1.0/4.0);

  TH1D* westMinusRates =new TH1D("West Minus", "West Minus", 120, 0, 1200);
  westMinusRates->Add(sideRates[1][index_A5]);
  westMinusRates->Add(sideRates[1][index_A7]);
  westMinusRates->Add(sideRates[1][index_B2]);
  westMinusRates->Add(sideRates[1][index_B10]);
  eastMinusRates->Scale(1.0/4.0);

  // add the histograms together to create a super sum
  for(int i = 0; i <= hist->GetNbinsX(); i++)
  {
    if( (eastPlusRates->GetBinContent(i)*westMinusRates->GetBinContent(i) < 0)
	|| (eastMinusRates->GetBinContent(i)*westPlusRates->GetBinContent(i) < 0))
    {
      hist->SetBinContent(i, 0);
    }
    else
    {
      hist->SetBinContent(i, sqrt(eastPlusRates->GetBinContent(i)*westMinusRates->GetBinContent(i))
			+ sqrt(eastMinusRates->GetBinContent(i)*westPlusRates->GetBinContent(i)));
    }
  }

  return hist;
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



/*
void ExtractA(TString dataPath, Int_t octetNum, Int_t structType, TString analysisChoice, vector < pair <string, int> > octetList)
{
  // create vectors to store event objects, chains to read in files, histograms to store Erecon values
  vector <ProcEvent*> proc_evts;	// need to declare vectors of both types
  vector <PureEvent*> pure_evts;	// but analysis will all be done on one type, indexed by structType
  vector <TChain*> chains;
  vector <TH1D*> East_hists;
  vector <TH1D*> West_hists;
  for(int n = 0; n < 8; n++)
  {
    if(structType == 13)
    {
      pure_evts.push_back(new PureEvent);
    }
    else if(structType == 20)
    {
      proc_evts.push_back(new ProcEvent);
    }
    else
    {
      cout << "Improper type passed to ExtractA(...). StructType needs to be 13 for pure simulation, 20 for processed simulation." << endl;
      return;
    }
    chains.push_back(new TChain("revCalSim"));
    East_hists.push_back(new TH1D(TString::Format("East Erecon %i", n), "East Erecon", 120, 0, 1200));
    West_hists.push_back(new TH1D(TString::Format("West Erecon %i", n), "West Erecon", 120, 0, 1200));
  }

  // set chains
  for(unsigned int l = 0; l < octetList.size(); l++)
  {
    if(octetList[l].first == "A2")
    {
      chains[0] -> Add(TString::Format("%s/revCalSim_%i_Beta.root", dataPath.Data(), octetList[l].second));      // A2 is index 0
    }
    else if(octetList[l].first == "A5")
    {
      chains[1] -> Add(TString::Format("%s/revCalSim_%i_Beta.root", dataPath.Data(), octetList[l].second));      // A5 is index 1
    }
    else if(octetList[l].first == "A7")
    {
      chains[2] -> Add(TString::Format("%s/revCalSim_%i_Beta.root", dataPath.Data(), octetList[l].second));      // A7 is index 2
    }
    else if(octetList[l].first == "A10")
    {
      chains[3] -> Add(TString::Format("%s/revCalSim_%i_Beta.root", dataPath.Data(), octetList[l].second));      // A10 is index 3
    }
    else if(octetList[l].first == "B2")
    {
      chains[4] -> Add(TString::Format("%s/revCalSim_%i_Beta.root", dataPath.Data(), octetList[l].second));      // B2 is index 4
    }
    else if(octetList[l].first == "B5")
    {
      chains[5] -> Add(TString::Format("%s/revCalSim_%i_Beta.root", dataPath.Data(), octetList[l].second));      // B5 is index 5
    }
    else if(octetList[l].first == "B7")
    {
      chains[6] -> Add(TString::Format("%s/revCalSim_%i_Beta.root", dataPath.Data(), octetList[l].second));      // B7 is index 6
    }
    else if(octetList[l].first == "B10")
    {
      chains[7] -> Add(TString::Format("%s/revCalSim_%i_Beta.root", dataPath.Data(), octetList[l].second));      // B10 is index 7
    }
    else
    {
      cout << "Bizarre error. Trying to save a run that isn't a beta." << endl;
    }
  }

  bool checkAnalysisChoice = false;	// whether we satisfy analysis choice or not.

  // perform actions for pure GEANT4 events below
  if(structType == 13)
  {
    for(unsigned int t = 0; t < chains.size(); t++)
    {
      chains[t] -> SetBranchAddress("primKE", &(pure_evts[t]->ke));
      chains[t] -> SetBranchAddress("primTheta", &(pure_evts[t]->theta));
      chains[t] -> SetBranchAddress("time", &(pure_evts[t]->hitTimes));
      chains[t] -> SetBranchAddress("MWPCEnergy", &(pure_evts[t]->edepWC));
      chains[t] -> SetBranchAddress("Edep", &(pure_evts[t]->edepScint));
    }

    for(unsigned int i = 0; i < chains.size(); i++)
    {
      for(int j = 0; j < (chains[i] -> GetEntriesFast()); j++)
      {
        chains[i] -> GetEntry(j);
        pure_evts[i]->eventID = j;
        pure_evts[i]->ptclSpecies = 1;

        // below is the "truth table" to do type identification
        if((pure_evts[i]->edepScint[0] > 0) && (pure_evts[i]->edepWC[0] > 0) && (pure_evts[i]->edepScint[1] == 0) && (pure_evts[i]->edepWC[1] == 0))
        {
          pure_evts[i]->type = 0;      // Type 0 East
        }
        else if((pure_evts[i]->edepScint[0] == 0) && (pure_evts[i]->edepWC[0] == 0) && (pure_evts[i]->edepScint[1] > 0) && (pure_evts[i]->edepWC[1] > 0))
        {
          pure_evts[i]->type = 0;      // Type 0 West
        }
        else if((pure_evts[i]->edepScint[0] > 0) && (pure_evts[i]->edepWC[0] > 0) && (pure_evts[i]->edepScint[1] > 0) && (pure_evts[i]->edepWC[1] > 0))
        {
          pure_evts[i]->type = 1;      	// Type 1 East or West. We need to look at timing structure to differentiate.
     	                           	// Except since we know initial direction, this isn't important for sorting East/West histograms
        	                        // Hence we categorize it into all type 1's.
        }
        else if((pure_evts[i]->edepScint[0] > 0) && (pure_evts[i]->edepWC[0] > 0) && (pure_evts[i]->edepScint[1] == 0) && (pure_evts[i]->edepWC[1] > 0))
        {
          pure_evts[i]->type = 2;      // This is a type 2/3 East trigger
        }
        else if((pure_evts[i]->edepScint[0] == 0) && (pure_evts[i]->edepWC[0] > 0) && (pure_evts[i]->edepScint[1] > 0) && (pure_evts[i]->edepWC[1] > 0))
        {
          pure_evts[i]->type = 2;      // Type 2/3 West trigger
        }
        else
        {
          pure_evts[i]->type = 4;      // Type 4 means unidentified, to match M.Brown's naming convention.
        }

	// check if our analysis choice is satisfied i.e. the types are chosen correctly
	// strcmp( , ) returns 0 if the two strings match. Hence the == 0.
	if(strcmp(analysisChoice.Data(), "Type-0") == 0 && pure_evts[i]->type == 0) checkAnalysisChoice = true;
	else if(strcmp(analysisChoice.Data(), "Type-1") == 0 && pure_evts[i]->type == 1) checkAnalysisChoice = true;
	else if(strcmp(analysisChoice.Data(), "Type-2-3") == 0 && pure_evts[i]->type == 2) checkAnalysisChoice = true;
	else if(strcmp(analysisChoice.Data(), "Misidentified") == 0 && pure_evts[i]->type == 4) checkAnalysisChoice = true;
	else if(strcmp(analysisChoice.Data(), "All") == 0) checkAnalysisChoice = true;
	else if(strcmp(analysisChoice.Data(), "1") == 0) checkAnalysisChoice = true;  // same as "All"
	else if(strcmp(analysisChoice.Data(), "2") == 0 && (pure_evts[i]->type == 0 || pure_evts[i]->type == 1)) checkAnalysisChoice = true;
	else if(strcmp(analysisChoice.Data(), "3") == 0) checkAnalysisChoice = true;  // same as "All" and choice 1
	else if(strcmp(analysisChoice.Data(), "4") == 0 && pure_evts[i]->type == 0) checkAnalysisChoice = true;
	else if(strcmp(analysisChoice.Data(), "5") == 0) checkAnalysisChoice = true;  // same as "All", choice 3, choice 1
	else if(strcmp(analysisChoice.Data(), "6") == 0 && pure_evts[i]->type == 1) checkAnalysisChoice = true;
	else if(strcmp(analysisChoice.Data(), "7") == 0 && pure_evts[i]->type == 2) checkAnalysisChoice = true;
	else if(strcmp(analysisChoice.Data(), "8") == 0 && pure_evts[i]->type == 2) checkAnalysisChoice = true;  // same as 7
	else if(strcmp(analysisChoice.Data(), "9") == 0 && pure_evts[i]->type == 2) checkAnalysisChoice = true;  // same as 8, 7

        // this "if statement" controls which cuts we do.
        // first note that we are assuming all ptcl species are 1. Hence all electrons.
        if( checkAnalysisChoice && ((pure_evts[i]->theta) >= M_PI/2.) && ((pure_evts[i]->theta) <= M_PI) )
        {
          East_hists[i] -> Fill(pure_evts[i] -> ke);
        }
        else if( checkAnalysisChoice && ((pure_evts[i]->theta) > 0) && ((pure_evts[i]->theta) <= M_PI/2.) )
        {
          West_hists[i] -> Fill(pure_evts[i] -> ke);
        }

	checkAnalysisChoice = false;
      }
    }
  }
  // perform actions for processed events below
  else if(structType == 20)
  {
    for(unsigned int t = 0; t < chains.size(); t++)
    {
      chains[t] -> SetBranchAddress("PID", &(proc_evts[t]->pid));
      chains[t] -> SetBranchAddress("side", &(proc_evts[t]->side));
      chains[t] -> SetBranchAddress("type", &(proc_evts[t]->type));
      chains[t] -> SetBranchAddress("Erecon", &(proc_evts[t]->Er));
    }

    for(unsigned int i = 0; i < chains.size(); i++)
    {
      for(int j = 0; j < (chains[i] -> GetEntriesFast()); j++)
      {
        chains[i] -> GetEntry(j);

	// check if our analysis choice is satisfied i.e. the types are chosen correctly
        // strcmp( , ) returns 0 if the two strings match. Hence the == 0.
        if(strcmp(analysisChoice.Data(), "Type-0") == 0 && proc_evts[i]->type == 0) checkAnalysisChoice = true;
        else if(strcmp(analysisChoice.Data(), "Type-1") == 0 && proc_evts[i]->type == 1) checkAnalysisChoice = true;
        else if(strcmp(analysisChoice.Data(), "Type-2-3") == 0 && proc_evts[i]->type == 2) checkAnalysisChoice = true;
        else if(strcmp(analysisChoice.Data(), "Misidentified") == 0 && proc_evts[i]->type == 4) checkAnalysisChoice = true;
        else if(strcmp(analysisChoice.Data(), "All") == 0) checkAnalysisChoice = true;
	else if(strcmp(analysisChoice.Data(), "1") == 0) checkAnalysisChoice = true;  // same as "All"                                         
	else if(strcmp(analysisChoice.Data(), "2") == 0 && (proc_evts[i]->type == 0 || proc_evts[i]->type == 1)) checkAnalysisChoice = true;
	else if(strcmp(analysisChoice.Data(), "3") == 0) checkAnalysisChoice = true;  // same as "All" and choice 1                            
        else if(strcmp(analysisChoice.Data(), "4") == 0 && proc_evts[i]->type == 0) checkAnalysisChoice = true;
        else if(strcmp(analysisChoice.Data(), "5") == 0) checkAnalysisChoice = true;  // same as "All", choice 3, choice 1                     
	else if(strcmp(analysisChoice.Data(), "6") == 0 && proc_evts[i]->type == 1) checkAnalysisChoice= true;
	else if(strcmp(analysisChoice.Data(), "7") == 0 && proc_evts[i]->type == 2) checkAnalysisChoice= true;
	else if(strcmp(analysisChoice.Data(), "8") == 0 && proc_evts[i]->type == 2) checkAnalysisChoice= true; // same as 7                   
	else if(strcmp(analysisChoice.Data(), "9") == 0 && proc_evts[i]->type == 2) checkAnalysisChoice= true; // same as 8, 7                


        // this "if statement" controls which cuts we do.
        if( checkAnalysisChoice && ((proc_evts[i]->pid) == 1) && ((proc_evts[i]->side) == 0) )
        {
          East_hists[i] -> Fill(proc_evts[i] -> Er);
        }
        else if( checkAnalysisChoice && ((proc_evts[i]->pid) == 1) && ((proc_evts[i]->side) == 1) )
        {
          West_hists[i] -> Fill(proc_evts[i] -> Er);
        }

	checkAnalysisChoice = false;
      }
    }
  }

  vector < vector <double> > RatesE_X, RatesW_X, RatesE_Y, RatesW_Y;
  vector <double> tmpEX, tmpWX, tmpEY, tmpWY;
  for(unsigned int m = 0; m < chains.size(); m++)
  {
    // get the bin contents in a vector of doubles
    for(int n = 0; n < (East_hists[m]->GetNbinsX()); n++)
    {
      tmpEY.push_back(East_hists[m]->GetBinContent(n));
      tmpWY.push_back(West_hists[m]->GetBinContent(n));
      tmpEX.push_back(East_hists[m]->GetXaxis()->GetBinCenter(n));
      tmpWX.push_back(West_hists[m]->GetXaxis()->GetBinCenter(n));
    }
    // store the vectors in a vector of size equal to number of files we are indexing over.
    RatesE_Y.push_back(tmpEY);
    RatesW_Y.push_back(tmpWY);
    RatesE_X.push_back(tmpEX);
    RatesW_X.push_back(tmpWX);

    // clear the tmp vectors otherwise we'll keep appending at the end
    tmpEY.clear();
    tmpWY.clear();
    tmpEX.clear();
    tmpWX.clear();
  }

  // super ratio calculation for quartet. 1 = East, 2 = West in Brad Plaster paper.
  vector <double> A;
  vector <double> A_error;
  vector <double> E_error;
  double N_E_minus = 0;
  double N_W_plus = 0;
  double N_E_plus = 0;
  double N_W_minus = 0;
  double w_E_minus = 0;
  double w_W_plus = 0;
  double w_E_plus = 0;
  double w_W_minus = 0;
  double ratioValue = 0;
  for(unsigned int l = 0; l < RatesE_Y[0].size(); l++)
  {
    // first index (from left) is the file or chain name. Second index is the energy window.
    if(RatesE_Y[index_A2][l] == 0 || RatesE_Y[index_A10][l] == 0 || RatesE_Y[index_B5][l] == 0 || RatesE_Y[index_B7][l] == 0)
    {
      N_E_minus = 0;
    }
    else
    {
      w_E_minus = 1./RatesE_Y[index_A2][l] + 1./RatesE_Y[index_A10][l] + 1./RatesE_Y[index_B5][l] + 1./RatesE_Y[index_B7][l];
      N_E_minus = 4./w_E_minus;
    }
    if(RatesW_Y[index_A5][l] == 0 || RatesW_Y[index_A7][l] == 0 || RatesW_Y[index_B2][l] == 0 || RatesW_Y[index_B10][l] == 0)
    {
      N_W_plus = 0;
    }
    else
    {
      w_W_plus = 1./RatesW_Y[index_A5][l] + 1./RatesW_Y[index_A7][l] + 1./RatesW_Y[index_B2][l] + 1./RatesW_Y[index_B10][l];
      N_W_plus = 4./w_W_plus;
    }
    if(RatesE_Y[index_A5][l] == 0 || RatesE_Y[index_A7][l] == 0 || RatesE_Y[index_B2][l] == 0 || RatesE_Y[index_B10][l] == 0)
    {
      N_E_plus = 0;
    }
    else
    {
      w_E_plus = 1./RatesE_Y[index_A5][l] + 1./RatesE_Y[index_A7][l] + 1./RatesE_Y[index_B2][l] + 1./RatesE_Y[index_B10][l];
      N_E_plus = 4./w_E_plus;
    }
    if(RatesW_Y[index_A2][l] == 0 || RatesW_Y[index_A10][l] == 0 || RatesW_Y[index_B5][l] == 0 || RatesW_Y[index_B7][l] == 0)
    {
      N_W_minus = 0;
    }
    else
    {
      w_W_minus = 1./RatesW_Y[index_A2][l] + 1./RatesW_Y[index_A10][l] + 1./RatesW_Y[index_B5][l] + 1./RatesW_Y[index_B7][l];
      N_W_minus = 4./w_W_minus;
    }

    // calculate super ratio
    if(N_E_minus == 0 || N_W_plus == 0 || N_E_plus == 0 || N_W_minus == 0)
    {
      ratioValue = 1;
    }
    else
    {
      ratioValue = (N_E_minus * N_W_plus) / (N_E_plus * N_W_minus);
    }

    A.push_back( (1 - sqrt(ratioValue)) / (1 + sqrt(ratioValue)) );

    if(ratioValue == 1)
    {
      A_error.push_back(0);
    }
    else
    {
      A_error.push_back( (sqrt(ratioValue) / (4*(sqrt(ratioValue) + 1)*(sqrt(ratioValue) + 1)))
                        *sqrt(w_E_minus + w_W_plus + w_W_minus + w_E_plus) );
    }
    E_error.push_back(5);	// 10keV bins -> +/-5 keV error?
//    cout << "Asym value: " << A[l] << endl;
  }

  ofstream outfile;

  char tmp[250];
  sprintf(tmp, "Asymm_Octet-%i_Analysis-%s_%i.txt", octetNum, analysisChoice.Data(), structType);
  outfile.open(tmp, ios::app);
  outfile << "Energy midpoint \t Asymm value \t E error \t A error" << endl;
  for(unsigned h = 0; h < A.size(); h++)
  {
    // any index from 0 to 7 will do, since they all are the same energy bins
    outfile << RatesE_X[0][h] << "\t"
            << A[h] << "\t"
            << E_error[h] << "\t"
            << A_error[h] << "\n";
  }
  outfile.close();

}

*/
