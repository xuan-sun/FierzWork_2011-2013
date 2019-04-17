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

void CreateEvts(TRandom3* factor, TString outFile, double pol, double b, int n_events);

int main(int argc, char* argv[])
{
  if(argc != 4)
  {
    cout << "Incorrect format. Execute with: \n";
    cout << "(executable) (-1 for East, +1 for West, 0 for none) (value of b) (nb events)" << endl;
    return 0;
  }

  int nbEvts = atoi(argv[3]);		// converts argument to int. C lib command. Getting number of events to gen
  double polArg = atof(argv[1]);		// polarization value
  double fierz_b = atof(argv[2]);	// value of fierz term

  TRandom3* engine = new TRandom3(0);

  for(int i = 0; i < 10; i++)
  {
    CreateEvts(engine, TString::Format("Evts_A_0_b_0_statTest_%i.root", i), polArg, fierz_b, nbEvts);
  }


  cout << "-------------- End of Program ---------------" << endl;
  return 0;
}

void CreateEvts(TRandom3* factor, TString outFile, double pol, double b, int n_events)
{

  if(pol == 0 || pol == -1 || pol == +1) { cout << "Value of polarization: " << pol << endl; }
  else { cout << "Polarization flag is incorrect." << endl; }

  cout << "Value of b: " << b << endl;

  cout << "Saving initial events kinematics in file " << outFile << endl;
  TFile fTree(outFile, "RECREATE");
  TTree* Evts = new TTree("Evts", "initial events kinematics");

  Int_t event_id = -1;  	// event ID will be incremented
  Int_t event_ptclID = 11;      // PDG flag 11 means electron

  // randomly throw the kinetic energy
  Double_t test_prob, pdf_value, Te_test, /*theta_test,*/ phi_test, cosTheta_test;
  Double_t event_KE = -1;       // keV
  Double_t event_theta = 0;
  Double_t event_phi = 0;

  Double_t event_pos[3];        // position in m
  Double_t event_dir[3];        // momentum direction, unit vector
  Double_t event_time;          // time in ns or s, it's all zero anyway
  Double_t event_weight;

  Evts -> Branch("num", &event_id, "event_id/I");
  Evts -> Branch("PID", &event_ptclID, "event_ptclID/I");
  Evts -> Branch("KE", &event_KE, "event_KE/D");
  Evts -> Branch("vertex", event_pos, "event_pos[3]/D");
  Evts -> Branch("direction", event_dir, "event_dir[3]/D");
  Evts -> Branch("time", &event_time, "event_time/D");
  Evts -> Branch("weight", &event_weight, "event_weight/D");

  Double_t normalizer = -1;     // find max value of prob distribution, to normalize beta PDF
  for(int i = 0; i < 9999; i++)
  {
    // Below is the max value, in 10000 steps, of the 1 + asymm spectra for normalization later
//    Double_t value = neutronCorrectedBetaSpectrum((neutronBetaEp*i)/10000)*(1 + correctedAsymmetry((neutronBetaEp*i)/10000, -1));

    // Below is max value of 1 + (spectral index on Fierz spectra)
    // si = 0 is fierz, si = 1 is SM.
    Double_t value =  1*neutronCorrectedSpectralBetaSpectrum(((neutronBetaEp*i)/10000), 1)
			+ pol*neutronCorrectedBetaSpectrum((neutronBetaEp*i)/10000)*correctedAsymmetry((neutronBetaEp*i)/10000, -1)
			+ b*neutronCorrectedSpectralBetaSpectrum(((neutronBetaEp*i)/10000), 0);

    if(value > normalizer)
    {
      normalizer = value;
    }
  }

  Double_t maxPDF = 0;

  // fetch number of events equal to max t in for loop
  for(int t = 0; t < n_events; t++)
  {
    test_prob = 1;
    pdf_value = 0;
    // rejection-acceptance sampling using BetaSpectrum.cc neutronCorrectedBetaSpectrum(...), correctedAsymmetry(...)
    while(test_prob > pdf_value)
    {
      cosTheta_test = 2*(factor -> Rndm()) - 1;      // uniformly sample cos(theta) from -1 to 1

      phi_test = 2*M_PI*(factor -> Rndm());
      Te_test = neutronBetaEp*(factor -> Rndm());    // sample random kinetic energy

      // normalizer is calculated to within 1/10000 precision and used to set the max value of the energy distribution.
      // Done for Standard Model spectra + asymmetric term
//      pdf_value = (neutronCorrectedBetaSpectrum(Te_test)*(1 + pol*correctedAsymmetry(Te_test, cosTheta_test))) / normalizer;

      // Now done for Standard Model spectra + Fierz term (depending on value of SI)
      pdf_value = ( 1*neutronCorrectedSpectralBetaSpectrum(Te_test, 1)
			+ pol*correctedAsymmetry(Te_test, cosTheta_test)*neutronCorrectedBetaSpectrum(Te_test)
			+ b*neutronCorrectedSpectralBetaSpectrum(Te_test, 0) ) / normalizer;

      test_prob = (factor -> Rndm());
    }

    if(pdf_value > maxPDF)
    {
      maxPDF = pdf_value;
    }

    event_theta = TMath::ACos(cosTheta_test);
    event_phi = phi_test;
    // set momentum direction. Sufficient to use spherical angles since it is unit vector
    event_dir[0] = TMath::Sin(event_theta)*TMath::Cos(event_phi);
    event_dir[1] = TMath::Sin(event_theta)*TMath::Sin(event_phi);
    event_dir[2] = TMath::Cos(event_theta);

    // set the values of each event for the TTree storage.
    event_pos[0] = 1;
    event_pos[1] = 1;
    while(event_pos[0]*event_pos[0] + event_pos[1]*event_pos[1] >= 0.0034129)
    {
      event_pos[0] = 0.05842*(2*(factor -> Rndm())-1);       // randomly choose x, y between -0.058 and 0.058m
      event_pos[1] = 0.05842*(2*(factor -> Rndm())-1);       // once it passes the radial cut, exit loop
    }
    event_pos[2] = 1.5*(2*(factor -> Rndm())-1);     // randomly choose z between -1.5 to 1.5m

    event_KE = Te_test;
    event_ptclID = 11;  // stays as electron
    event_id = t;       // one unique event ID # for each ptcl
    event_time = 0;     // ns or s, doesn't matter since it's zero
    event_weight = 1;   // all equally weighted

    Evts -> Fill();
  }

  if(maxPDF > 1)
  {
    cout << "Badness. Max value of pdf is " << maxPDF << endl;
  }


  fTree.Write();
  fTree.Close();

}
