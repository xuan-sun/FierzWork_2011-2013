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

  // read in arguments.
//  Int_t octNb = atoi(argv[1]);

  TString fileName = "Evts_A_0_b_0_half.root";
  AddBranchToTreeInFile(fileName);
  cout << "Done updating runNumber = " << fileName.Data() << endl;

//  C -> Print("fierz.pdf");
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

void AddBranchToTreeInFile(TString fName)
{
  TFile *f = new TFile(fName.Data(), "UPDATE");
  TTree *t = (TTree*)f->Get("Evts");

  double ke;
  t->SetBranchAddress("KE", &ke);

  double ke_99;
  TBranch *b = t->Branch("99KE", &ke_99, "ke_99/D");

  for(long int i = 0; i < t->GetEntries(); i++)
  {
    t->GetEntry(i);
    ke_99 = 0.99*ke;
    b->Fill();
  }

  t->Write();
  delete f;
}
