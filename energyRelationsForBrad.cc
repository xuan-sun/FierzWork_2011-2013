#include         <iostream>
#include         <fstream>
#include         <TGaxis.h>
#include         <sstream>
#include         <TGraph.h>
#include         <TGraphErrors.h>
#include         <TCanvas.h>
#include         <TApplication.h>
#include         <stdlib.h>
#include         <TF1.h>
#include         <TH1.h>
#include         <TProfile.h>
#include         <TObjArray.h>
#include         <TStyle.h>
#include         <TMarker.h>
#include         <math.h>
#include         <TStyle.h>
#include         <TPaveStats.h>
#include         <TPaveText.h>
#include         <vector>
#include         <string.h>
#include         <fstream>
#include         <TROOT.h>
#include         <TFile.h>
#include         <TLegend.h>
#include         <TLegendEntry.h>
#include         <time.h>
#include         <TH2F.h>
#include         <assert.h>
#include         <string>
#include         <TRandom.h>
#include         <TTree.h>
#include         <TChain.h>
#include         <TVector.h>
#include         <vector>
#include         <utility>
#include         <TLeaf.h>
#include         <math.h>
using            namespace std;

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

int main(int argc, char* argv[])
{

  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (run number)" << endl;
    return 0;
  }

  int runNb = atoi(argv[1]);

  TCanvas *C = new TCanvas("canvas", "canvas");

  TChain *chain = new TChain("revCalSim");
  chain->Add(Form("fromSept2017Onwards/A_0_b_0/revCalSim_%i.root", runNb));
  int side = -1;
  int type = -1;
  double erecon = -1;
  double pke = -1;
  int pid = -1;

  vector <double> ereconVec;
  vector <double> pkeVec;

  chain->SetBranchAddress("PID", &pid);
  chain->SetBranchAddress("side", &side);
  chain->SetBranchAddress("type", &type);
  chain->SetBranchAddress("primaryKE", &pke);
  chain->SetBranchAddress("Erecon", &erecon);

  for(int i = 0; i < chain->GetEntries(); i++)
  {
    chain->GetEntry(i);
    if(pid == 1 && type < 4 && side < 2 && erecon >= 0)
    {
      ereconVec.push_back(erecon);
      pkeVec.push_back(pke);
    }
  }

  if(ereconVec.size() == 0)
  {
    cout << "No events passed the cut. Exiting.." << endl;
    return 0;
  }

  TGraph *graph = new TGraph(ereconVec.size(), &(ereconVec[0]), &(pkeVec[0]));

  graph->GetXaxis()->SetTitle("Erecon (keV)");
  graph->GetYaxis()->SetTitle("Primary KE (keV)");

  graph->SetTitle(Form("revCalSim_%i, N_{entries} = %i", runNb, ereconVec.size()));

  graph->Draw("AP");

  // prints the canvas with a dynamic TString name of the name of the file
  C -> Print(Form("energy_revCalSim_%i.pdf", runNb));
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;

}
