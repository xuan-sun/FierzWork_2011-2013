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
#include	 <TLine.h>
#include	 <TLatex.h>
using		 namespace std;

struct Event
{
  int num;
  int PID;
  double KE;
  double vertex[3];
  double direction[3];
  double time;
  double weight;
};

// forward declarations for useful functions
vector < pair <string,int> >  LoadOctetList(TString fileName);
vector < TChain* > GetChainsOfRuns(vector < pair <string,int> > octetList, TString dataPath);
vector < vector < TH1D* > > CreateMCCountsHistograms(vector <TChain*> runsChains);
TH1D* CreateSuperRatio(vector < vector < TH1D* > > sideCounts);
TH1D* CalculateAofE(TH1D* R);
TH1D* DivideMPMEffects(TH1D* AofE);

// debugging functions
TH1D* SingleRunAsymm(vector < vector < TH1D* > > sideCounts);

// plotting functions
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xTitle, TString yTitle, TString command);

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
TApplication plot_program("FADC_readin",0,0,0,0);

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
  vector < TChain* > runFiles = GetChainsOfRuns(octetIndices, "/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/Evts_Files/testing_asymm_2");
  // load all the histograms of east and west given the cuts of interest
  vector < vector < TH1D* > > counts = CreateMCCountsHistograms(runFiles);

  cout << "Completed loading all counts from simulations..." << endl;

  // make a file and write the output of the calculations to it
//  TFile f(TString::Format("MC_asymm_Octet_%i_type0.root", octNb), "RECREATE");
  TH1D* superRatio_Erecon = CreateSuperRatio(counts);
  TH1D* asymm_Erecon = CalculateAofE(superRatio_Erecon);
  TH1D* flattenedAsymm = DivideMPMEffects(asymm_Erecon);

  double xMin = 180;
  double xMax = 780;

  TF1 *fit = new TF1("asymm fit", "[0] + [1]*x", xMin, xMax);

  fit->SetParName(0, "A0");
  fit->SetParName(1, "weak magnetism");
  flattenedAsymm->Fit("asymm fit");
  TF1 *fitResults = flattenedAsymm->GetFunction("asymm fit");
  cout << "Chi squared value is " << fitResults->GetChisquare()
	<< " with ndf of " << fitResults->GetNDF()
	<< ". For a final chisquared/ndf = " << fitResults->GetChisquare() / fitResults->GetNDF() << endl;


//  asymm_Erecon->Write();

  PlotHist(C, 1, 1, flattenedAsymm, "", "Primary Energy (keV)", "Asymmetry", "");

  // Save our plot and print it out as a pdf.
  C -> Print("preMichael_flattenedAsymm_noMPMEffects.pdf");
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
    }
  }

  infile.close();       // ensure you're closing the read-in file.

  return pairs;
}

vector < TChain* > GetChainsOfRuns(vector < pair <string,int> > octetList, TString dataPath)
{
  vector < TChain* > runs;

  TString asymmDown = "A_-1_b_0";
  TString asymmUp = "A_1_b_0";

  for(unsigned int i = 0; i < 8; i++)
  {
    runs.push_back(new TChain("Evts"));
  }

  for(unsigned int l = 0; l < octetList.size(); l++)
  {
    if(octetList[l].first == "A2")
    {
      runs[index_A2] -> Add(TString::Format("%s/%s/Evts_%i.root", dataPath.Data(), asymmDown.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A5")
    {
      runs[index_A5] -> Add(TString::Format("%s/%s/Evts_%i.root", dataPath.Data(), asymmUp.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A7")
    {
      runs[index_A7] -> Add(TString::Format("%s/%s/Evts_%i.root", dataPath.Data(), asymmUp.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "A10")
    {
      runs[index_A10] -> Add(TString::Format("%s/%s/Evts_%i.root", dataPath.Data(), asymmDown.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B2")
    {
      runs[index_B2] -> Add(TString::Format("%s/%s/Evts_%i.root", dataPath.Data(), asymmUp.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B5")
    {
      runs[index_B5] -> Add(TString::Format("%s/%s/Evts_%i.root", dataPath.Data(), asymmDown.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B7")
    {
      runs[index_B7] -> Add(TString::Format("%s/%s/Evts_%i.root", dataPath.Data(), asymmDown.Data(), octetList[l].second));
    }
    else if(octetList[l].first == "B10")
    {
      runs[index_B10] -> Add(TString::Format("%s/%s/Evts_%i.root", dataPath.Data(), asymmUp.Data(), octetList[l].second));
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

    runsChains[i]->SetBranchAddress("num", &evt[i]->num);
    runsChains[i]->SetBranchAddress("PID", &evt[i]->PID);
    runsChains[i]->SetBranchAddress("KE", &evt[i]->KE);
    runsChains[i]->SetBranchAddress("vertex", &evt[i]->vertex);
    runsChains[i]->SetBranchAddress("direction", &evt[i]->direction);
    runsChains[i]->SetBranchAddress("time", &evt[i]->time);
    runsChains[i]->SetBranchAddress("weight", &evt[i]->weight);
  }

  for(unsigned int j = 0; j < runsChains.size(); j++)
  {
    cout << "For chain j = " << j << " we have number of events = " << runsChains[j]->GetEntries() << endl;

    for(unsigned int i = 0; i < runsChains[j]->GetEntries(); i++)
    {
      runsChains[j]->GetEntry(i);

      if(evt[j]->PID == 11)
      {
        if(evt[j]->direction[2] < 0)
        {
	    countHistsEast[j]->Fill(evt[j]->KE);
        }
        else if(evt[j]->direction[2] > 0)
        {
	    countHistsWest[j]->Fill(evt[j]->KE);
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

  for(unsigned int i = 0; i < 2; i++)
  {
    for(unsigned int j = 0; j < sideCounts[i].size(); j++)
    {
      // ensure that all the errors are gonna get propagated correctly later.
      // This is a mandatory ROOT trick.
      sideCounts[i][j]->Sumw2();
    }
  }


  TH1D* hist = new TH1D("Super ratio", "Super ratio", 120, 0, 1200);
  TH1D* eastPlusCounts = new TH1D("East Plus", "East Plus", 120, 0, 1200);
  TH1D* westPlusCounts =new TH1D("West Plus", "West Plus", 120, 0, 1200);
  TH1D* eastMinusCounts = new TH1D("East Minus", "East Minus", 120, 0, 1200);
  TH1D* westMinusCounts =new TH1D("West Minus", "West Minus", 120, 0, 1200);

  vector <double> mySetErrors;

  double nEP, nWP, nEM, nWM;
  double sumEP2, sumWP2, sumEM2, sumWM2;

  double nEA5, nEA7, nEB2, nEB10, nWA5, nWA7, nWB2, nWB10;
  double nEA2, nEA10, nEB5, nEB7, nWA2, nWA10, nWB5, nWB7;

  for(int i = 0; i < sideCounts[0][index_A2]->GetNbinsX(); i++)
  {
    nEA5 = sideCounts[0][index_A5]->GetBinContent(i);
    nEA7 = sideCounts[0][index_A7]->GetBinContent(i);
    nEB2 = sideCounts[0][index_B2]->GetBinContent(i);
    nEB10 = sideCounts[0][index_B10]->GetBinContent(i);

    nWA5 = sideCounts[1][index_A5]->GetBinContent(i);
    nWA7 = sideCounts[1][index_A7]->GetBinContent(i);
    nWB2 = sideCounts[1][index_B2]->GetBinContent(i);
    nWB10 = sideCounts[1][index_B10]->GetBinContent(i);

    nEA2 = sideCounts[0][index_A2]->GetBinContent(i);
    nEA10 = sideCounts[0][index_A10]->GetBinContent(i);
    nEB5 = sideCounts[0][index_B5]->GetBinContent(i);
    nEB7 = sideCounts[0][index_B7]->GetBinContent(i);

    nWA2 = sideCounts[1][index_A2]->GetBinContent(i);
    nWA10 = sideCounts[1][index_A10]->GetBinContent(i);
    nWB5 = sideCounts[1][index_B5]->GetBinContent(i);
    nWB7 = sideCounts[1][index_B7]->GetBinContent(i);

    sumEP2 = ( 1.0/nEA5 + 1.0/nEA7 + 1.0/nEB2 + 1.0/nEB10 );
    sumWP2 = ( 1.0/nWA5 + 1.0/nWA7 + 1.0/nWB2 + 1.0/nWB10 );
    sumEM2 = ( 1.0/nEA2 + 1.0/nEA10 + 1.0/nEB5 + 1.0/nEB7 );
    sumWM2 = ( 1.0/nWA2 + 1.0/nWA10 + 1.0/nWB5 + 1.0/nWB7 );

    if(nEA5 == 0 || nEA7 == 0 || nEB2 == 0 || nEB10 == 0 || nWA5 == 0 || nWA7 == 0 || nWB2 == 0 || nWB10 == 0
       || nEA2 == 0 || nEA10 == 0 || nEB5 == 0 || nEB7 == 0 || nWA2 == 0 || nWA10 == 0 || nWB5 == 0 || nWB7 == 0)
    {
      nEP = 0;
      nWP = 0;
      nEM = 0;
      nWM = 0;
    }
    else
    {
      // using weight formula N_avg = (sum) sigma^-2 * N / (sum) sigma^-2
      nEP = 4.0 / sumEP2;
      nWP = 4.0 / sumWP2;
      nEM = 4.0 / sumEM2;
      nWM = 4.0 / sumWM2;
    }

    eastPlusCounts->SetBinContent(i, nEP);
    eastPlusCounts->SetBinError(i, sqrt(1.0 / sumEP2));
    westPlusCounts->SetBinContent(i, nWP);
    westPlusCounts->SetBinError(i, sqrt(1.0 / sumWP2));
    eastMinusCounts->SetBinContent(i, nEM);
    eastMinusCounts->SetBinError(i, sqrt(1.0 / sumEM2));
    westMinusCounts->SetBinContent(i, nWM);
    westMinusCounts->SetBinError(i, sqrt(1.0 / sumWM2));
  }

  double r1up, r1down, r2up, r2down;
  double er1up, er1down, er2up, er2down;
  for(int i = 0; i <= hist->GetNbinsX(); i++)
  {
    r1up = eastPlusCounts->GetBinContent(i);
    r1down = eastMinusCounts->GetBinContent(i);
    r2up = westPlusCounts->GetBinContent(i);
    r2down = westMinusCounts->GetBinContent(i);
    er1up = eastPlusCounts->GetBinError(i);
    er1down = eastMinusCounts->GetBinError(i);
    er2up = westPlusCounts->GetBinError(i);
    er2down = westMinusCounts->GetBinError(i);


    if(r1up <= 0 || r1down <= 0 || r2up <= 0 || r2down <= 0 )
    {
      hist->SetBinContent(i, 0);
      mySetErrors.push_back(0);
    }
    else
    {
      hist->SetBinContent(i, (r1down*r2up) / (r1up*r2down));

      mySetErrors.push_back( ( (r1down*r2up) / (r1up*r2down) )
                           * sqrt( (er1down/r1down)*(er1down/r1down)
				 + (er2up/r2up)*(er2up/r2up)
				 + (er1up/r1up)*(er1up/r1up)
				 + (er2down/r2down)*(er2down/r2down) ) );
    }
  }

  hist->SetError(&(mySetErrors[0]));

  return hist;
}

TH1D* CalculateAofE(TH1D* R)
{
  TH1D* hAofE = new TH1D("AofE", "AofE", 120, 0, 1200);

  double asymm = 0;
  double asymmErr = 0;

  for(int i = 0; i <= R->GetNbinsX(); i++)
  {
    if(R->GetBinContent(i) == 0)
    {
      asymm = 0;
      asymmErr = 0;
    }
    else
    {
      asymm = (1 - sqrt(R->GetBinContent(i)) ) / ( 1 + sqrt(R->GetBinContent(i)) );
      asymmErr = (R->GetBinError(i))
		 / ( sqrt(R->GetBinContent(i)) * (sqrt(R->GetBinContent(i)) + 1) * (sqrt(R->GetBinContent(i)) + 1) );
    }

    hAofE->SetBinContent(i, asymm);
    hAofE->SetBinError(i, asymmErr);
  }

  return hAofE;
}

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xTitle, TString yTitle, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle(xTitle);
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle(yTitle);
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(0.02, 0.08);


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

  TLine *yLow = new TLine(180, gPad->GetUymin(), 180, gPad->GetUymax());
  TLine *yHigh = new TLine(780, gPad->GetUymin(), 780, gPad->GetUymax());
  yLow->SetLineWidth(2);
  yHigh->SetLineWidth(2);
  yLow->SetLineColor(1);
  yHigh->SetLineColor(1);
  yLow->Draw("SAME");
  yHigh->Draw("SAME");

}

TH1D* DivideMPMEffects(TH1D* AofE)
{
  cout << "Dividing out the corrections from BetaSpectrum.hh ... " << endl;

  TH1D* hCorrectionDivided = new TH1D("corrected", "MPM corrected", 120, 0, 1200);

  double energyEffects = 0;

  for(int i = 0; i <= AofE->GetNbinsX(); i++)
  {
    energyEffects = (beta(AofE->GetBinCenter(i)) / 2.0);
//		  * (1.0 + WilkinsonACorrection((AofE->GetBinCenter(i) + m_e) / m_e) + shann_h_minus_g_a2pi((AofE->GetBinCenter(i) + m_e) / m_e));
    hCorrectionDivided->SetBinContent(i, AofE->GetBinContent(i) / energyEffects);
    hCorrectionDivided->SetBinError(i, AofE->GetBinError(i) / energyEffects);
  }

  cout << "Created new histogram with energy effects divided out..." << endl;

  return hCorrectionDivided;

}



TH1D* SingleRunAsymm(vector < vector < TH1D* > > sideCounts)
{
  for(unsigned int i = 0; i < 2; i++)
  {
    for(unsigned int j = 0; j < sideCounts[i].size(); j++)
    {
      // ensure that all the errors are gonna get propagated correctly later.
      // This is a mandatory ROOT trick.
      sideCounts[i][j]->Sumw2();
    }
  }

  TH1D *hAofE = new TH1D("AofE", "AofE", 120, 0, 1200);
  hAofE->Sumw2();

  vector <double> mySetErrors;
  int index = index_A2;

  double asymm = 0;
  for(int i = 0; i <= sideCounts[0][index]->GetNbinsX(); i++)
  {
    cout << "sideCounts[0][index_A2]->GetBinContent(" << i << ") = " << sideCounts[0][index]->GetBinContent(i) << endl;

    if(sideCounts[0][index]->GetBinContent(i) == 0 || sideCounts[1][index]->GetBinContent(i) == 0)
    {
      asymm = 0;
    }
    else
    {
      asymm = ( sideCounts[0][index]->GetBinContent(i) - sideCounts[1][index]->GetBinContent(i) )
            / ( sideCounts[0][index]->GetBinContent(i) + sideCounts[1][index]->GetBinContent(i) );


    }

    hAofE->SetBinContent(i, asymm);
  }

  return hAofE;

}

