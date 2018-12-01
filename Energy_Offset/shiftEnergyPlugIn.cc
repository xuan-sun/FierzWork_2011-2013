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
vector <int>  LoadOctetList(TString fileName);
void AddBranchToTreeInFile(int fileNumber);


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

  // create a vector with all the run numbers in this octet
  vector <int> runIndices = LoadOctetList(TString::Format("../OctetLists/octet_list_%i.dat", octNb));

  for(unsigned int i =0; i < runIndices.size(); i++)
  {
    AddBranchToTreeInFile(runIndices[i]);
    cout << "Done updating runNumber = " << runIndices[i] << endl;
  }

//  C -> Print("fierz.pdf");
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

vector <int>  LoadOctetList(TString fileName)
{
  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  string runType;
  int runIndex;
  vector <int> runNumbers;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> runType >> runIndex;

      // take only runs that actually exist (beta or background, we don't have depol runs).
      if(runType == "A1" || runType == "A2" || runType == "A4" || runType == "A5" ||
         runType == "A7" || runType == "A9" || runType == "A10" || runType == "A12" ||
         runType == "B1" || runType == "B2" || runType == "B4" || runType == "B5" ||
         runType == "B7" || runType == "B9" || runType == "B10" || runType == "B12")
      {
        runNumbers.push_back(runIndex);
      }
    }
  }

  infile.close();       // ensure you're closing the read-in file.

  return runNumbers;
}

void AddBranchToTreeInFile(int fileNumber)
{
  TFile *f = new TFile(Form("/mnt/Data/xuansun/replay_pass3_FINALCAL/replay_pass3_%i.root", fileNumber), "UPDATE");
  TTree *t = (TTree*)f->Get("pass3");

  double Erecon;
  t->SetBranchAddress("Erecon", &Erecon);

  double shiftedUp_2012_Erecon;
  TBranch *b = t->Branch("shiftedUp_2012_Erecon", &shiftedUp_2012_Erecon, "shiftedUp_2012_Erecon/D");

  for(long int i = 0; i < t->GetEntries(); i++)
  {
    t->GetEntry(i);
//    shiftedUp_2011_Erecon = Erecon + 3.438;	// this is the residual keV from the 4 calibration source points.
    shiftedUp_2012_Erecon = Erecon + (3.438/4.0) ;	// incredibly, the 2012-2013 summed error residual is the same as 2011-2012.

    b->Fill();
  }

  t->Write();
  delete f;
}
