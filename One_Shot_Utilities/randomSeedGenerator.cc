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
#include	 <TLatex.h>
#include	 <TRandom3.h>
using            namespace std;


int main(int argc, char* argv[])
{
  TRandom3 *engine = new TRandom3(0);

  ofstream outfile;
  TString fileName = Form("randomMixingSeeds.txt");
  outfile.open(fileName.Data(), ios::app);

  double flagCheck = engine->Rndm();

  cout << flagCheck << endl;

  if(flagCheck < 0.5)
  {
    outfile << 1000 << "\t";
    outfile << 0.06 << "\n";	// 0.06 -> b = -0.1
//    outfile << 0.06*engine->Rndm() << "\n";
  }
  else if(flagCheck >= 0.5)
  {
    outfile << -1 << "\t";
    outfile << 0.04 << "\n";	// 0.04 -> b = 0.1
//    outfile << 0.04*engine->Rndm() << "\n";
  }

  outfile.close();


  return 0;
}

