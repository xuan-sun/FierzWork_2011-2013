#include	 "BetaSpectrum.hh"

#include	 <complex.h>
#include	 <stdio.h>
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
#include	 <TRandom3.h>
#include 	 <TTree.h>
#include	 <TMath.h>
using		 namespace std;

//required later for plot_program
//TApplication plot_program("FADC_readin",0,0,0,0);

int main(int argc, char* argv[])
{
  if(argc != 1)
  {
    cout << "Incorrect format. Execute with: \n";
    cout << "(executable)" << endl;
    return 0;
  }

  TRandom3* engine = new TRandom3(0);


  double tau = 880.0;	// neutron lifetime in seconds

  double integrationEndpoint = 782;
  double integrationStartpoint = 200;
  double numberOfSteps = 582;	// approximately corresponds to 1keV steps
  double stepWidth = (integrationEndpoint - integrationStartpoint) / numberOfSteps;

  double integratedArea = 0;

  // endpoint is chosen such that we integrate over the entire beta decay spectrum. DO NOT go over 782 or else divide by 0 error.
  for(double ke = integrationStartpoint; ke <= integrationEndpoint; ke = ke + stepWidth)
  {
    integratedArea = integratedArea + ( (1.0/(m_e+ke))*(1.0/(m_e+ke))*neutronCorrectedBetaSpectrum(ke) ) * stepWidth;

    cout << "At ke = " << ke << ", function value is: " << integratedArea << endl;
  }

  integratedArea = integratedArea*tau*m_e*m_e;

  cout << "Integrating in steps of " << stepWidth << "keV, final area is: " << integratedArea << endl;



  cout << "-------------- End of Program ---------------" << endl;
  return 0;
}

