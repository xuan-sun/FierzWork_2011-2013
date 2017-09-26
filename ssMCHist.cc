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
  vector < TChain* > runFiles = GetChainsOfRuns(octetIndices, "fromSept2017Onwards/A_0_b_0");
  // load all the histograms of east and west, turn them into rates.
  vector < vector < TH1D* > > rates = CreateRateHistograms(runFiles);

  TFile f(TString::Format("MC_A_0_b_0_Octet_%i_ssHist.root", octNb), "RECREATE");
  // Begin processing the read in data now
  TH1D* SS_Erecon = CreateSuperSum(rates);
  SS_Erecon->Write();

  cout << "About to plot the histogram now..." << endl;

  PlotHist(C, 1, 1, SS_Erecon, "", "");

  // Save our plot and print it out as a pdf.
//  C -> Print("fierz.pdf");
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
    runs.push_back(new TChain("revCalSim"));
  }

  for(unsigned int l = 0; l < octetList.size(); l++)
  {
    if(octetList[l].first == "A2")
    {
      runs[index_A2] -> Add(TString::Format("%s/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A5")
    {
      runs[index_A5] -> Add(TString::Format("%s/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A7")
    {
      runs[index_A7] -> Add(TString::Format("%s/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A10")
    {
      runs[index_A10] -> Add(TString::Format("%s/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B2")
    {
      runs[index_B2] -> Add(TString::Format("%s/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B5")
    {
      runs[index_B5] -> Add(TString::Format("%s/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B7")
    {
      runs[index_B7] -> Add(TString::Format("%s/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B10")
    {
      runs[index_B10] -> Add(TString::Format("%s/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
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

    runsChains[i]->SetBranchAddress("PrimaryID", &evt[i]->eventNum);
    runsChains[i]->SetBranchAddress("side", &evt[i]->side);
    runsChains[i]->SetBranchAddress("type", &evt[i]->type);
    runsChains[i]->SetBranchAddress("Erecon", &evt[i]->Erecon);
    runsChains[i]->SetBranchAddress("PID", &evt[i]->pid);

    runsChains[i]->SetBranchAddress("time", &evt[i]->time);

//    runsChains[i]->GetBranch("time")->GetLeaf("timeE")->SetAddress(&(evt[i])->tE);
//    runsChains[i]->GetBranch("time")->GetLeaf("timeW")->SetAddress(&(evt[i])->tW);

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
        liveTimeEast.push_back(evt[j]->time[0]);
        liveTimeWest.push_back(evt[j]->time[1]);
      }
    }
  }


  // here we loop back over the event histograms in east and west
  // and scale them down by time of last event aka live time.
/*  for(unsigned int i = 0; i < rateHistsEast.size(); i++)
  {
    rateHistsEast[i]->Scale(1.0/liveTimeEast[i]);
  }
  for(unsigned int i = 0; i < rateHistsWest.size(); i++)
  {
    rateHistsWest[i]->Scale(1.0/liveTimeWest[i]);
  }
*/
  rateHists.push_back(rateHistsEast);
  rateHists.push_back(rateHistsWest);

  return rateHists;

}

TH1D* CreateSuperSum(vector < vector < TH1D* > > sideRates)
{
  TH1D* hist = new TH1D("Super sum", "Super sum Erecon spectrum", 120, 0, 1200);

  // sum the "like" histograms without any statistical weight
  TH1D* eastPlusRates = new TH1D("East Plus", "East Plus", 120, 0, 1200);
  eastPlusRates->Add(sideRates[0][index_A5]);
  eastPlusRates->Add(sideRates[0][index_A7]);
  eastPlusRates->Add(sideRates[0][index_B2]);
  eastPlusRates->Add(sideRates[0][index_B10]);
//  eastPlusRates->Scale(1.0/4.0);

  TH1D* westPlusRates =new TH1D("West Plus", "West Plus", 120, 0, 1200);
  westPlusRates->Add(sideRates[1][index_A5]);
  westPlusRates->Add(sideRates[1][index_A7]);
  westPlusRates->Add(sideRates[1][index_B2]);
  westPlusRates->Add(sideRates[1][index_B10]);
//  westPlusRates->Scale(1.0/4.0);

  TH1D* eastMinusRates = new TH1D("East Minus", "East Minus", 120, 0, 1200);
  eastMinusRates->Add(sideRates[0][index_A2]);
  eastMinusRates->Add(sideRates[0][index_A10]);
  eastMinusRates->Add(sideRates[0][index_B5]);
  eastMinusRates->Add(sideRates[0][index_B7]);
//  eastMinusRates->Scale(1.0/4.0);

  TH1D* westMinusRates =new TH1D("West Minus", "West Minus", 120, 0, 1200);
  westMinusRates->Add(sideRates[1][index_A2]);
  westMinusRates->Add(sideRates[1][index_A10]);
  westMinusRates->Add(sideRates[1][index_B5]);
  westMinusRates->Add(sideRates[1][index_B7]);
//  eastMinusRates->Scale(1.0/4.0);

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
      if(i == 22)
      {
        cout << "Bin contents are " << eastPlusRates->GetBinContent(i) << endl;
      }
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
