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
  TString fileName = Form("randomMixingSeeds_newFormatThatWorks.txt");
  outfile.open(fileName.Data(), ios::app);

  for(int i = 0; i < 9; i++)
  {
    outfile << 0.1 + (0.06*engine->Rndm() - 0.03) << "\t";
  }
  outfile << 0.1 + (0.06*engine->Rndm() - 0.03) << "\n";

  for(int j = 0; j < 10; j++)
  {
    cout << "For verification, our random seed generator is giving " << j << ": "
	 << 0.1 + (0.06*engine->Rndm() - 0.03) << endl;
  }


  outfile.close();


  return 0;
}

