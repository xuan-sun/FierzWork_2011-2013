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
void AddBranchToTreeInFile(TString fName, double fittedEndpoint, double meanEndpoint);
vector < pair <int, double> > ReadInEndpoints(TString fileName, vector < pair <int, double> > updateThis);

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

  vector < pair <int, double> > fittedEndpoints;
  fittedEndpoints = ReadInEndpoints(Form("endPointFits_noGainCorrection_asymmTwiddledSpectra_index19_2011-2012_radialCut_0-49mm.txt"), fittedEndpoints);

  cout << "Done reading in endpoints... " << endl;

  double meanEndpoint = 786.7;	// keV. From 2011-2012 baseline fitting histograms, all with b_input = 0. radial cut 49mm.
//  double meanEndpoint = 787.8;	// keV. From 2011-2012 baseline fitting histograms, all with b_input = 0. radial cut 30mm.
//  double meanEndpoint = 789.6;	// keV. From 2011-2012 asymmetric error envelope, twiddle index 19, radial cut 30mm. WRONG.
//  double meanEndpoint = 782.6;	// keV. From 2011-2012 asymmetric error envelope, twiddle index 19. WRONG.
  double twiddleEndpoint = -1;
  for(unsigned int i = 0; i < fittedEndpoints.size(); i++)
  {
    twiddleEndpoint = fittedEndpoints[i].second;
    AddBranchToTreeInFile(Form("/mnt/data2/xuansun/analyzed_files/2011-2012_geom_twiddles/TwiddledSimFiles_A_0_b_0_matchingParamSet_19/SimAnalyzed_2011-2012_Beta_paramSet_%i_0.root", fittedEndpoints[i].first), twiddleEndpoint, meanEndpoint);
  }


//  C -> Print("fierz.pdf");
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

void AddBranchToTreeInFile(TString fName, double fittedEndpoint, double meanEndpoint)
{
  TFile *f = new TFile(fName.Data(), "UPDATE");
  TTree *t = (TTree*)f->Get("SimAnalyzed");

  double erecon;
  t->SetBranchAddress("Erecon", &erecon);

  double ereconCorr_r49mm;
  TBranch *b = t->Branch("Erecon_corr_r49mm", &ereconCorr_r49mm, "ereconCorr_r49mm/D");

  for(long int i = 0; i < t->GetEntries(); i++)
  {
    t->GetEntry(i);
    ereconCorr_r49mm = erecon*(meanEndpoint / fittedEndpoint);
    b->Fill();
  }

  t->Write();
  delete f;

  cout << "Done updating runNumber = " << fName.Data() << endl;
}


vector < pair <int, double> > ReadInEndpoints(TString fileName, vector < pair <int, double> > updateThis)
{
  int indexNb = -1;
  double endpoint = -1;
  double endpointErr = -1;

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
      bufstream >> indexNb >> endpoint >> endpointErr;

      updateThis.push_back(make_pair(indexNb, endpoint));
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  return updateThis;
}
