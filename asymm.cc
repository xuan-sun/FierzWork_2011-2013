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
vector < vector < TH1D* > > CreateMCCountsHistograms(vector <TChain*> runsChains);
TH1D* CreateSuperRatio(vector < vector < TH1D* > > sideCounts);
TH1D* CalculateAofE(TH1D* R);

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

  cout << "Running executable asymm for octet " << octNb << " ..." << endl;

  // Reads in the octet list and saves the run files indices corresponding to an octet number
  vector < pair <string,int> > octetIndices = LoadOctetList(TString::Format("%s/octet_list_%i.dat", "OctetLists", octNb));
  // Points TChains at the run files idenified in the octet lists above
  vector < TChain* > runFiles = GetChainsOfRuns(octetIndices, "/mnt/Data/xuansun/G4Sims_AbFit/AbFit_Fierz_2011-2013/2011-2012_geom/");
  // load all the histograms of east and west given the cuts of interest
  vector < vector < TH1D* > > counts = CreateMCCountsHistograms(runFiles);

  cout << "Completed loading all counts from simulations..." << endl;

  // make a file and write the output of the calculations to it
  TFile f(TString::Format("MC_asymm_Octet_%i_type0.root", octNb), "RECREATE");
  TH1D* superRatio_Erecon = CreateSuperRatio(counts);
  TH1D* asymm_Erecon = CalculateAofE(superRatio_Erecon);
  asymm_Erecon->Write();



//  PlotHist(C, 1, 1, SS_Erecon, "", "");

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

  for(unsigned int i = 0; i < 8; i++)
  {
    runs.push_back(new TChain("revCalSim"));
  }

  for(unsigned int l = 0; l < octetList.size(); l++)
  {
    if(octetList[l].first == "A2")
    {
      runs[index_A2] -> Add(TString::Format("%s/A_-1_b_0/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A5")
    {
      runs[index_A5] -> Add(TString::Format("%s/A_1_b_0/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A7")
    {
      runs[index_A7] -> Add(TString::Format("%s/A_1_b_0/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A10")
    {
      runs[index_A10] -> Add(TString::Format("%s/A_-1_b_0/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B2")
    {
      runs[index_B2] -> Add(TString::Format("%s/A_1_b_0/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B5")
    {
      runs[index_B5] -> Add(TString::Format("%s/A_-1_b_0/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B7")
    {
      runs[index_B7] -> Add(TString::Format("%s/A_-1_b_0/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B10")
    {
      runs[index_B10] -> Add(TString::Format("%s/A_1_b_0/revCalSim_%i.root", dataPath.Data(), octetList[l].second));
    }
  }

  return runs;
}

vector < vector < TH1D* > > CreateMCCountsHistograms(vector <TChain*> runsChains)
{
  // reminder: each evt[i] index corresponds to the indices noted at global scope.
  vector <Event*> evt;

  vector < vector <TH1D*> > countHists;

  vector <TH1D*> countHistsEast;
  vector <TH1D*> countHistsWest;

  for(unsigned int i = 0; i < runsChains.size(); i++)
  {
    evt.push_back(new Event);
    countHistsEast.push_back(new TH1D(TString::Format("East Counts %i", i), "East Counts", 120, 0, 1200));
    countHistsWest.push_back(new TH1D(TString::Format("West Counts %i", i), "West Counts", 120, 0, 1200));

    runsChains[i]->SetBranchAddress("PrimaryID", &evt[i]->eventNum);
    runsChains[i]->SetBranchAddress("side", &evt[i]->side);
    runsChains[i]->SetBranchAddress("type", &evt[i]->type);
    runsChains[i]->SetBranchAddress("Erecon", &evt[i]->Erecon);
    runsChains[i]->SetBranchAddress("PID", &evt[i]->pid);

    runsChains[i]->SetBranchAddress("time", &evt[i]->time);
  }

  for(unsigned int j = 0; j < runsChains.size(); j++)
  {
    for(unsigned int i = 0; i < runsChains[j]->GetEntriesFast(); i++)
    {
      runsChains[j]->GetEntry(i);
      if(evt[j]->pid == 1 && evt[j]->type == 0 && evt[j]->Erecon >= 0)
      {
        if(evt[j]->side == 0)
        {
          countHistsEast[j]->Fill(evt[j]->Erecon);
        }
        else if(evt[j]->side == 1)
        {
          countHistsWest[j]->Fill(evt[j]->Erecon);
        }
      }
    }
  }


  countHists.push_back(countHistsEast);
  countHists.push_back(countHistsWest);

  return countHists;
}

TH1D* CreateSuperRatio(vector < vector < TH1D* > > sideCounts)
{
  // here I want to calculate the super ratio using the counts in East, West sides
  // as a function of energy, for the 8 runs that for an octet.
  // we assume background subtraction is already done (makes sense, it's a sim)
  // and that rates and counts are effectively the same thing (makes sense, no blinded live times in a sim).


  TH1D* hist = new TH1D("Super sum", "Super sum Erecon spectrum", 120, 0, 1200);

/*
  vector <double> mySetErrors;

  // sum the "like" histograms without any statistical weight
  TH1D* eastPlusRates = new TH1D("East Plus", "East Plus", 120, 0, 1200);
  sideRates[0][index_A5]->Sumw2();
  eastPlusRates->Add(sideRates[0][index_A5]);
  sideRates[0][index_A7]->Sumw2();
  eastPlusRates->Add(sideRates[0][index_A7]);
  sideRates[0][index_B2]->Sumw2();
  eastPlusRates->Add(sideRates[0][index_B2]);
  sideRates[0][index_B10]->Sumw2();
  eastPlusRates->Add(sideRates[0][index_B10]);

  TH1D* westPlusRates =new TH1D("West Plus", "West Plus", 120, 0, 1200);
  sideRates[1][index_A5]->Sumw2();
  westPlusRates->Add(sideRates[1][index_A5]);
  sideRates[1][index_A7]->Sumw2();
  westPlusRates->Add(sideRates[1][index_A7]);
  sideRates[1][index_B2]->Sumw2();
  westPlusRates->Add(sideRates[1][index_B2]);
  sideRates[1][index_B10]->Sumw2();
  westPlusRates->Add(sideRates[1][index_B10]);

  TH1D* eastMinusRates = new TH1D("East Minus", "East Minus", 120, 0, 1200);
  sideRates[0][index_A2]->Sumw2();
  eastMinusRates->Add(sideRates[0][index_A2]);
  sideRates[0][index_A10]->Sumw2();
  eastMinusRates->Add(sideRates[0][index_A10]);
  sideRates[0][index_B5]->Sumw2();
  eastMinusRates->Add(sideRates[0][index_B5]);
  sideRates[0][index_B7]->Sumw2();
  eastMinusRates->Add(sideRates[0][index_B7]);

  TH1D* westMinusRates =new TH1D("West Minus", "West Minus", 120, 0, 1200);
  sideRates[1][index_A2]->Sumw2();
  westMinusRates->Add(sideRates[1][index_A2]);
  sideRates[1][index_A10]->Sumw2();
  westMinusRates->Add(sideRates[1][index_A10]);
  sideRates[1][index_B5]->Sumw2();
  westMinusRates->Add(sideRates[1][index_B5]);
  sideRates[1][index_B7]->Sumw2();
  westMinusRates->Add(sideRates[1][index_B7]);

  double errorValueEachBin = 0;
  // add the histograms together to create a super sum
  for(int i = 0; i <= hist->GetNbinsX(); i++)
  {
    if( (eastPlusRates->GetBinContent(i)*westMinusRates->GetBinContent(i) < 0)
	|| (eastMinusRates->GetBinContent(i)*westPlusRates->GetBinContent(i) < 0))
    {
      hist->SetBinContent(i, 0);
      mySetErrors.push_back(0);
    }
    else
    {
      hist->SetBinContent(i, sqrt(eastPlusRates->GetBinContent(i)*westMinusRates->GetBinContent(i))
			+ sqrt(eastMinusRates->GetBinContent(i)*westPlusRates->GetBinContent(i)));

      if(eastPlusRates->GetBinContent(i) == 0 || westMinusRates->GetBinContent(i) == 0
	|| eastMinusRates->GetBinContent(i) == 0 || westPlusRates->GetBinContent(i) == 0)
      {
        errorValueEachBin = 0;
      }
      else
      {
	errorValueEachBin = sqrt( 0.25*westMinusRates->GetBinContent(i)*eastPlusRates->GetBinError(i)
				*eastPlusRates->GetBinError(i)/eastPlusRates->GetBinContent(i)
			  + 0.25*eastPlusRates->GetBinContent(i)*westMinusRates->GetBinError(i)
                                *westMinusRates->GetBinError(i)/westMinusRates->GetBinContent(i)
                          + 0.25*westPlusRates->GetBinContent(i)*eastMinusRates->GetBinError(i)
                                *eastMinusRates->GetBinError(i)/eastMinusRates->GetBinContent(i)
                          + 0.25*eastMinusRates->GetBinContent(i)*westPlusRates->GetBinError(i)
                                *westPlusRates->GetBinError(i)/westPlusRates->GetBinContent(i) );
      }
      mySetErrors.push_back(errorValueEachBin);
    }
  }

  hist->SetError(&(mySetErrors[0]));
*/
  return hist;
}

TH1D* CalculateAofE(TH1D* R)
{
  TH1D* hAofE = new TH1D("AofE", "AofE", 120, 0, 1200);
 



  return hAofE;
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
