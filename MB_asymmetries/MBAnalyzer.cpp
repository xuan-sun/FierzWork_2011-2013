/*
Compilable code which houses the heart of the analysis. 
Here are classes which handle calculating asymmetries of either data 
or simulation. There will be flags to be passed dictating what type of 
asymmetry, what type of data, and what octet/runs to use. 

The code should be able to do BG subtraction, construct asymmetries, 
or do both depending on the flags passed.

*/

#include "EvtRateHandler.hh"
#include "Asymmetries.hh"
#include "SystematicCorrections.hh"
#include "MBUtils.hh"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>

#include "BetaSpectrum.hh"

TString anaChoices[10] = {"A","B","C","D","E","F","G","H","J","K"};


//Types of Corrections to apply
std::string corr ("UnCorr");//{"UnCorr","DeltaExpOnly","DeltaTheoryOnly","AllCorr"};
                             

Double_t POL_minus = 0.9981;
Double_t POL_plus = 0.9937;
Double_t delta_POL = POL_plus-POL_minus;
Double_t POL_ave = (POL_plus+POL_minus) / 2.;

bool withPOL = false; //Set this to true to correct DATA for the polarimetry measurement


std::vector <Int_t> badOct {7,60,61,62,63,64,65,66}; // 9,59 used to be part of this, but not totally sure why... they are just not complete octets
                                  // Octet 7 had W anode dead for part of run
                                  // Either need to discard, or apply the 
                                  // Charge cloud method to determine a 
                                  // coincidence
                                  // 60,61,62,63,64,65,66 have bad TDC spectra
                                  // 70 and 92 just have low statistics so they had bad corrections. Put them back in when using high statistics corrections

// The Process functions will calculate all the raw asymmetries on a bin-by-bin basis and write them to file
void ProcessOctets(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool UNBLIND=false);
void ProcessQuartets(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool UNBLIND=false);
void ProcessPairs(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool UNBLIND=false); 

// makes plots of 2*A/Beta for each octet (pair, quartet, and octet as a whole) and fits over the range specified, for whatever grouping provided ("Octet", "Quartet", "Pair")
void PlotAsymmetriesByGrouping(std::string groupType, Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t Elow=220., Double_t Ehigh=680., Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool UNBLIND=false, int key=0);

//Returns the values of the systematic corrections and statistical error on this for now
std::vector < std::vector <Double_t> >  LoadOctetSystematics(Int_t octet, std::string anaChoice, std::vector <Double_t> enBinMidpoint);

//Returns a vector containing all the theory corrections to A0 for a particular bin
std::vector <Double_t> LoadTheoryCorrections(std::vector <Double_t> enBinMidpoint);

// Collects all the asymmetries, bin-by-bin, and produces a final asymmetry plot, both A_SR and 2*A/Beta, for whatever grouping provided ("Octet", "Quartet", "Pair"). Also makes a plot of the 
// integrated asymmetry vs octet/quartet/pair number
void PlotFinalAsymmetries(std::string groupType, Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t Elow=220., Double_t Ehigh=680., Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool UNBLIND=false, int key=0);

// notes on key value...
// key==0 : Normal asymmetries
// key>=10 : Quadrant asymmtries, where 10 is quadrant0, 11 is quadrant 1, etc
// key>=20 : Radial Asymmetries, where 20 is radial cut 0, etc (in 10 mm incremements)



// Function to return beta when given the kinetic energy of an electron
Double_t returnBeta(Double_t En) { 
  Double_t me = 510.998928; //rest mass energy of electron in keV
  return sqrt(En*En+2.*En*me)/(En+me);
};



int main(int argc, char* argv[])
{
  
  //if ( corr.size()>7 && corr.compare(corr.size()-8,8,"_withPOL") == 0 ) withPOL = true; 

  /*  std::cout << "Is this thing broken???\n";
  std::vector<double> vec{1.,3.,2.,7.,8.,4.,5.,6.,5.,10.};

    for ( auto i : vec ) {
      std::cout << i << " ";
    }
    std::cout << std::endl;
    
    std::vector<double> vecSorted = sortVecDouble(vec,true);

    for ( auto i : vecSorted ) {
      std::cout << i << " ";
    }

    std::cout << std::endl;
  */	

  if (argc==1) {
    std::cout << "USAGE: ./MBAnalyzer.exe [anaChoice] [octmin] [octmax] [analysis window min] [analysis window max] [corrections] [key]\n";
    exit(0);
  }
    
  std::string analysisChoice = argc>1 ? std::string(argv[1]) : "A";
  Int_t octBegin = argc>2 ? atoi(argv[2]) : 0;
  Int_t octEnd = argc>2 ? atoi(argv[3]) : 1;
  Double_t enBinWidth = 10.;
  Double_t Elow = argc>4 ? atoi(argv[4]) : 220.;//220
  Double_t Ehigh = argc>4 ? atoi(argv[5]) : 680.;//680
  if ( argc>6 ) corr = std::string(argv[6]);
  int key = 0;
  if ( argc==8 ) key = atoi(argv[7]);
  bool UKdata = true;//true;
  bool simulation = false;
  bool applyAsymm = false;

  if (simulation) withPOL=false;

  //I should keep track of the raw asymmetry, Experimental systematic corrected asymmetry, and theoretical 
  // Systematic asymmetry in a file each time I run to see the percent correction of each...

  //NEED TO ADD IN ABILITY TO READ IN SYSTEMATICS FOR QUARTET AND PAIR

  //****************************************************************
  //****************************************************************
  // ONLY TURN THIS ON WHEN READY TO UNBLIND. IT WILL USE THE TRUE 
  // CLOCK TIMES.
  //****************************************************************
  //****************************************************************
  bool UNBLIND = false;


  if (UNBLIND) {
    std::string decision;
    std::cout << "YOU ARE ABOUT TO DO SOME SORT OF UNBLINDING!!!! \n\n";
    std::cout << "To continue, type YES: ";
    std::cin >> decision;

    if (decision!=std::string("YES")) exit(0);

  }

  
  try {
    
    /*std::vector < Double_t > enBinMedian; //Holds the center of the energy bins
    for (Int_t i=0; i<120; i++) {
      Double_t En = i*enBinWidth+enBinWidth/2.;
      enBinMedian.push_back(En);
    }    
    std::vector <Double_t> theoryCorr = LoadTheoryCorrections(enBinMedian);    
    for (UInt_t i=0; i<theoryCorr.size(); i++) std::cout << enBinMedian[i] << " " << theoryCorr[i] << "\n";*/
    
    
    /*TString aCh[5] = {"C","D","F","J","K"};//{"A","B","G","H"};//{"C","J","K","H"};//"A","D"
    for (auto ach : aCh) {
      //ProcessOctets(octBegin, octEnd, std::string(ach.Data()), enBinWidth, UKdata, simulation, UNBLIND);
      ProcessPairs(octBegin, octEnd, std::string(ach.Data()), enBinWidth, UKdata, simulation, UNBLIND);
      //PlotAsymmetriesByGrouping("Octet",octBegin, octEnd, std::string(ach.Data()), Elow, Ehigh, enBinWidth, UKdata, simulation, UNBLIND);
      //PlotFinalAsymmetries("Octet",octBegin, octEnd, std::string(ach.Data()), Elow, Ehigh, enBinWidth, UKdata, simulation, UNBLIND);
      }*/

    // Loop over keys
    /*int keys[] {0,10,11,12,13,20,21,22,23,24,25};

    for ( auto& k : keys ) {
      PlotAsymmetriesByGrouping("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, UNBLIND, k);
      PlotFinalAsymmetries("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, UNBLIND, k);
      }*/

   
    ProcessOctets(octBegin, octEnd, analysisChoice, enBinWidth, UKdata, simulation, UNBLIND);
    //PlotAsymmetriesByGrouping("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, UNBLIND, key);
    //PlotFinalAsymmetries("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, UNBLIND, key);
    //ProcessPairs(octBegin, octEnd, analysisChoice, enBinWidth, UKdata, simulation, UNBLIND);
    //*************** DONT FORGET TO CHANGE FIDUCIAL CUT TO 50!

    //ProcessQuartets(octBegin, octEnd, analysisChoice, enBinWidth, UKdata, simulation, applyAsymm, UNBLIND);
    //PlotAsymmetriesByGrouping("Quartet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, UNBLIND);
    //PlotFinalAsymmetries("Quartet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, UNBLIND);
     
    //    ProcessPairs(octBegin, octEnd, analysisChoice, enBinWidth, UKdata, simulation, UNBLIND);
    //PlotAsymmetriesByGrouping("Pair",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, UNBLIND);
    //PlotFinalAsymmetries("Pair",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, UNBLIND);
    
  }
  catch(const char* ex){
    std::cerr << "Error: " << ex << std::endl;
  }
  
  
  return 0;
}


void ProcessOctets(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t enBinWidth, bool UKdata, bool simulation, bool UNBLIND) {

  //std::ofstream octval("octvalUK.dat");
  //octval << "oct" << "\t" 
  //	 << "Asymm" << "\t" 
  //	 << "Error" << "\t" 
  //	 << "Pull" << "\n";
  

  for ( Int_t octet=octBegin; octet<=octEnd; octet++ ) {
    if ( std::find(badOct.begin(), badOct.end(),octet) != badOct.end() ) continue;  //Checking if octet should be ignored for data quality reasons
    try {
      OctetAsymmetry oct(octet,anaChoice,enBinWidth, 50., UKdata, simulation, UNBLIND);
      //oct.calcTotalAsymmetry(180.,780.);
      oct.calcAsymmetryBinByBin(); 
      //oct.calcNCSUSumAsymmetryBinByBin(); 
      //oct.calcSuperSum();
      //oct.calcSuperSumNCSUstyle();
      oct.writeAsymToFile();
      //oct.writeSuperSumToFile(); 
      
      //  octval << octet << "\t" 
      //     << oct.returnTotalAsymmetry() << "\t" 
      //     << oct.returnTotalAsymmetryError() << "\t" 
      //     << 0. << "\n";
      
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }

  }

  //octval.close();

};

void ProcessQuartets(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t enBinWidth, bool UKdata, bool simulation, bool UNBLIND) {
  
  for (Int_t octet=octBegin; octet<=octEnd; octet++) {
    if ( std::find(badOct.begin(), badOct.end(),octet) != badOct.end() ) continue;  //Checking if octet should be ignored for data quality reasons
    try {
      QuartetAsymmetry quart(octet,anaChoice,enBinWidth, 50., UKdata, simulation, UNBLIND);
      quart.calcAsymmetryBinByBin();
      //quart.calcSuperSum();
      quart.writeAsymToFile();
      //quart.writeSuperSumToFile();
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }
  
};

void ProcessPairs(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t enBinWidth, bool UKdata, bool simulation, bool UNBLIND) {
  
  for (Int_t octet=octBegin; octet<=octEnd; octet++) {
    if ( std::find(badOct.begin(), badOct.end(),octet) != badOct.end() ) continue;  //Checking if octet should be ignored for data quality reasons
    try {
      PairAsymmetry pair(octet, anaChoice, enBinWidth, 50., UKdata, simulation, UNBLIND);
      //pair.calcAsymmetryBinByBin();
      //pair.calcSuperSum();
      //pair.writeAsymToFile();
      //pair.writeSuperSumToFile();
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }
  
};


void PlotFinalAsymmetries(std::string groupType, Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t Elow, Double_t Ehigh, Double_t enBinWidth, bool UKdata, bool simulation, bool UNBLIND, int key) {

  if (groupType!="Quartet" && groupType!="Octet" && groupType!="Pair") throw "Bad group type given to PlotFinalAsymmetries. Options are \"Octet\", \"Quartet\", or \"Pair\""; 

  key = ( groupType==std::string("Octet") ? key : 0 );

  Int_t numBins = 1200/enBinWidth;

  std::vector < Double_t > enBinMedian; //Holds the center of the energy bins
  for (Int_t i=0; i<numBins; i++) {
    Double_t En = i*enBinWidth+enBinWidth/2.;
    enBinMedian.push_back(En);
  }

  //Loading theory systematics... These aren't dependent on grouping, only energy of bin
  std::vector <Double_t> theoryCorr = LoadTheoryCorrections(enBinMedian);

  std::vector < std::vector <Double_t> > rawAsymByGroup(3,std::vector <Double_t> (0));
  std::vector < std::vector <Double_t> > AsymByGroup(3,std::vector <Double_t> (0));
  
  std::vector < std::vector <Double_t> > groupRawAsymByBin(3,std::vector <Double_t> (numBins,0.));
  std::vector < std::vector <Double_t> > groupAsymByBin(3,std::vector <Double_t> (numBins,0.));

  
  //bool isOctet = groupType=="Octet" ? true : false;
  bool isQuartet = groupType=="Quartet" ? true : false;
  bool isPair = groupType=="Pair" ? true : false;
  std::map < std::string, Int_t > numQuartetLoops = {{"Octet",1},{"Quartet",2},{"Pair",2}};
  std::map < std::string, Int_t > numPairLoops = {{"Octet",1},{"Quartet",1},{"Pair",2}};
  std::string quartetName[2] = {"A","B"};
  
  std::string basePath = simulation ? getenv("SIM_ANALYSIS_RESULTS") : UKdata ? getenv("ANALYSIS_RESULTS") : getenv("MPM_ANALYSIS_RESULTS");
  
  std::string txt;
  Double_t binEdge, Asym, AsymError;

  Int_t pair = 0, quartet = 0;

  for (Int_t octet=octBegin; octet<=octEnd; octet++) {

    if (std::find(badOct.begin(), badOct.end(),octet) != badOct.end()) { pair+=4; quartet+=2; continue; } //Checking if octet should be ignored for data quality reasons
    
    for (Int_t q = 0; q<numQuartetLoops[groupType]; q++) {

      for (Int_t p=0; p<numPairLoops[groupType]; p++) {

	std::string path = ( basePath+"Octet_"+itos(octet)+"/"+groupType+"Asymmetry/"+(UNBLIND?"UNBLINDED_":"") + 
			     corr + "_" + (withPOL?"withPOL_":"") +"FittedAsymmetry_Octet"+itos(octet)+"_AnaCh"+anaChoice+"_"+
			     (isQuartet?("Quartet_"+quartetName[q]+"_"):isPair?("Pair_"+quartetName[q]+itos(p)+"_"):"")+itos((int)Elow)+"-"+itos((int)Ehigh) + 
			     ( key/10==1 ? "_Quadrant"+itos(key-10) : ( key/10==2 ? "_RadialRing"+itos(key-20) :"" ) ) + ".dat");
	//std::cout << path << std::endl;

	std::ifstream infile(path.c_str());  
   
	if (infile.is_open()) {
	  infile >> txt >> Asym >> AsymError;
	  rawAsymByGroup[0].push_back(isQuartet?quartet:isPair?pair:octet);
	  rawAsymByGroup[1].push_back(-Asym); // Note the negative sign here to turn the raw asymmetry into purely a positive number (by definition, it is negative)
	  rawAsymByGroup[2].push_back(AsymError);

	  //std::cout << octet << " " << Asym << " " << AsymError << std::endl;
      
	  infile >> txt >> Asym >> AsymError;
	  AsymByGroup[0].push_back(isQuartet?quartet:isPair?pair:octet);
	  AsymByGroup[1].push_back(Asym);
	  AsymByGroup[2].push_back(AsymError);
	  infile.close();
	
	  path = basePath + "Octet_" + itos(octet) + "/" + groupType + "Asymmetry/" + (UNBLIND?"UNBLINDED_":"") + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice  +
	    (isQuartet?("_Quartet_"+quartetName[q]):isPair?("_Pair_"+quartetName[q]+itos(p)):"")+( key/10==1 ? "_Quadrant"+itos(key-10) : ( key/10==2 ? "_RadialRing"+itos(key-20) :"" ) ) + ".dat"; 
	  infile.open(path.c_str());
	  //std::cout << path << std::endl;
	  
	  //Apply systematics if run is an octet
	  //TODO: put in systematics for other types of runs...
	  std::vector < std::vector <Double_t> > deltaSys(enBinMedian.size(),std::vector<Double_t>(2,1.));
	  if (groupType=="Octet") deltaSys = LoadOctetSystematics(octet,anaChoice,enBinMedian);
	  
	  
	  Int_t i = 0;
	  while (infile >> binEdge >> Asym >> AsymError) {
	    groupRawAsymByBin[0][i] = binEdge;
	    groupRawAsymByBin[1][i] += AsymError>0. ? 1./power(AsymError*deltaSys[i][0]/theoryCorr[i],2)*Asym*deltaSys[i][0]/theoryCorr[i] : 0.; //Applying Delta_exp here.. Should this affect AsymError???
	    groupRawAsymByBin[2][i] += AsymError>0. ? 1/power(AsymError*deltaSys[i][0]/theoryCorr[i],2) : 0.;         // Since Asym is multiplied by delta, AsymError would be as well...
	    //std::cout << binEdge << " " << groupRawAsymByBin[1][i] << " " << groupRawAsymByBin[2][i] << std::endl;
	    i++;
	  }
	  infile.close();
	}
	pair++;
      }
      quartet++;
    }
  }
  //Do final calculations of the rates in each bin and their associated errors
  for (unsigned int i=0; i<groupRawAsymByBin[1].size(); i++) {
    groupRawAsymByBin[1][i] = groupRawAsymByBin[2][i]>0. ? groupRawAsymByBin[1][i] / groupRawAsymByBin[2][i] : 0.; // This is sum of weights*Asym / sum of weights 
    groupRawAsymByBin[2][i] = groupRawAsymByBin[2][i]>0. ? 1./sqrt(groupRawAsymByBin[2][i]) : 0.; // Sqrt of sum of weights
    
    groupAsymByBin[0][i] = groupRawAsymByBin[0][i];
    groupAsymByBin[1][i] = 2*groupRawAsymByBin[1][i]/returnBeta(enBinMedian[i]); //Divide out the energy dependence...
    groupAsymByBin[2][i] = 2*groupRawAsymByBin[2][i]/returnBeta(enBinMedian[i]);

    // Apply polarimetry correction if necessary
    if (withPOL) {
      groupRawAsymByBin[1][i] = groupRawAsymByBin[1][i] / ( POL_ave ); 
      groupRawAsymByBin[2][i] = groupRawAsymByBin[2][i] / ( POL_ave );
      groupAsymByBin[1][i] = groupAsymByBin[1][i] / ( POL_ave ); 
      groupAsymByBin[2][i] = groupAsymByBin[2][i] / ( POL_ave );
    }
    
    std::cout << enBinMedian[i] << " " << groupAsymByBin[1][i] << " " << groupAsymByBin[2][i] << std::endl;
  }
  
  // Plotting stuff
  std::string outFile = basePath + "Asymmetries/"+(UNBLIND?"UNBLINDED_":"") + corr + "_" + (withPOL?"withPOL_":"") + groupType +"Asymmetries_AnaCh" + anaChoice 
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) + "_Octets_" +itos(octBegin)+"-"+itos(octEnd)+( key/10==1 ? "_Quadrant"+itos(key-10) : ( key/10==2 ? "_RadialRing"+itos(key-20) :"" ) );
  
  std::string pdfFile = outFile+std::string(".pdf");
  std::string txtFile = outFile+std::string(".txt");

  std::ofstream asymFile(txtFile.c_str());
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1000., 1400.);
  c1->Divide(1,4);
  
  gStyle->SetOptFit(1111);
  gStyle->SetTitleX(0.25);
  gStyle->SetStatX(0.75);
  gStyle->SetStatY(0.85);
  
  std::vector <double> errorX(1000,0.);
  std::string title;
  
  c1->cd(1);
  
  TGraphErrors *g = new TGraphErrors(rawAsymByGroup[0].size(), &rawAsymByGroup[0][0],&rawAsymByGroup[1][0],&errorX[0], &rawAsymByGroup[2][0]);
  title = groupType+"-By-"+groupType+" Raw Measured Asymmetry " + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  g->SetTitle(title.c_str());
  g->SetMarkerStyle(20);
  g->SetLineWidth(2);
  g->GetXaxis()->SetLimits(rawAsymByGroup[0][0]-2., rawAsymByGroup[0][rawAsymByGroup[0].size()-1]+2.);
  g->GetXaxis()->SetTitle("Number");
  g->GetYaxis()->SetTitle("Raw Asymmetry");
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();
  
  TF1 *fit = new TF1("fit","[0]",rawAsymByGroup[0][0]-1., rawAsymByGroup[0][rawAsymByGroup[0].size()-1]+1.);
  fit->SetLineColor(kRed);
  fit->SetLineWidth(3);
  fit->SetParameter(0,0.05);
  
  g->Fit("fit","R");

  asymFile << "RawA_oct_by_oct\t" << fit->GetParameter(0) << "\t" << fit->GetParError(0) << std::endl;

  //Writing to file the raw asymmetries of each octet
  /*std::ofstream octval("octvalUK.dat");
  octval << "oct" << "\t" 
	 << "Asymm" << "\t" 
	 << "Error" << "\t" 
	 << "Pull" << "\n";
  
  for (UInt_t i=0; i<rawAsymByGroup[0].size(); i++) {
    octval << rawAsymByGroup[0][i] << "\t" 
	   << rawAsymByGroup[1][i] << "\t" 
	   << rawAsymByGroup[2][i] << "\t" 
	   << ( rawAsymByGroup[1][i] - fit->GetParameter(0) ) / rawAsymByGroup[2][i] << "\n";
  } 
  octval.close();*/
  
  g->Draw("AP");
  g->SetMinimum(fit->GetParameter(0)-0.02);//((simulation && !AsymmOn) ? -0.05 : 0.03);
  g->SetMaximum(fit->GetParameter(0)+0.02);//(simulation && !AsymmOn) ? 0.05 : 0.07);
  c1->Update();
  
  c1->cd(2);
  
  std::string corrText = corr=="UnCorr" ? " 2A_{SR}/#beta " : ( corr=="DeltaExpOnly" ? " 2A_{SR}#upoint(1+#Delta_{exp})/#beta " : 
								( corr=="DeltaTheoryOnly" ? " 2A_{SR}#upoint(1+#Delta_{Th})/#beta " : 
								  ( " 2A_{SR}#upoint(1+#Delta_{exp})/(#beta#upoint(1+#Delta_{Th})) ")));
  
  TGraphErrors *gBeta = new TGraphErrors(AsymByGroup[0].size(), &AsymByGroup[0][0],&AsymByGroup[1][0],&errorX[0], &AsymByGroup[2][0]);
  title = groupType+"-By-"+groupType+corrText + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gBeta->SetTitle(title.c_str());
  gBeta->SetMarkerStyle(20);
  gBeta->SetLineWidth(2);
  gBeta->GetXaxis()->SetLimits(rawAsymByGroup[0][0]-2., rawAsymByGroup[0][rawAsymByGroup[0].size()-1]+2.);
  gBeta->GetXaxis()->SetTitle("Number");
  gBeta->GetYaxis()->SetTitle("Asymmetry");
  gBeta->GetXaxis()->CenterTitle();
  gBeta->GetYaxis()->CenterTitle();
  
  TF1 *fitBeta = new TF1("fitBeta","[0]",rawAsymByGroup[0][0]-1., rawAsymByGroup[0][rawAsymByGroup[0].size()-1]+1.);
  fitBeta->SetLineColor(kRed);
  fitBeta->SetLineWidth(3);
  fitBeta->SetParameter(0,0.05);
  
  gBeta->Fit("fitBeta","R");
  
  asymFile << "BetaCorrA_oct_by_oct\t" << fitBeta->GetParameter(0) << "\t" << fitBeta->GetParError(0) << std::endl;

  gBeta->Draw("AP");
  gBeta->SetMinimum(fitBeta->GetParameter(0)-0.02);//(simulation && !AsymmOn) ? -0.5 : -0.14);
  gBeta->SetMaximum(fitBeta->GetParameter(0)+0.03);//(simulation && !AsymmOn) ? 0.5 : -0.09);
  c1->Update();
  
  
  //c1->Print(pdfFile.c_str());

  std::string corrText2 = corr=="UnCorr" ? " A_{SR} " : ( corr=="DeltaExpOnly" ? " A_{SR}#upoint(1+#Delta_{exp}) " : 
								( corr=="DeltaTheoryOnly" ? " A_{SR}#upoint(1+#Delta_{Th})/ " : 
								  ( " A_{SR}#upoint(1+#Delta_{exp})/(#upoint(1+#Delta_{Th})) ")));

  c1->cd(3);
  Int_t offset = 0;
  TGraphErrors *g2 = new TGraphErrors(enBinMedian.size()-offset, &enBinMedian[offset], &groupRawAsymByBin[1][offset], 0, &groupRawAsymByBin[2][offset]);
  title = groupType + corrText2 + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  g2->SetTitle(title.c_str());
  g2->SetMarkerStyle(20);
  g2->SetLineWidth(2);
  g2->GetXaxis()->SetLimits(0., 800.);
  g2->GetXaxis()->SetTitle("Energy (keV)");
  g2->GetYaxis()->SetTitle("Uncorrected Asymmetry");
  g2->GetXaxis()->CenterTitle();
  g2->GetYaxis()->CenterTitle();
  
  g2->Draw("AP");
  g2->SetMinimum(-0.07);//(simulation && !AsymmOn) ? -0.05 : -0.08);
  g2->SetMaximum(-0.02);//(simulation && !AsymmOn) ? 0.05 : -0.01);
  c1->Update();
	
  c1->cd(4);
	
  TGraphErrors *gBeta2 = new TGraphErrors(enBinMedian.size()-offset, &enBinMedian[offset], &groupAsymByBin[1][offset], 0, &groupAsymByBin[2][offset]);
  title = groupType + corrText + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gBeta2->SetTitle(title.c_str());
  gBeta2->SetMarkerStyle(20);
  gBeta2->SetLineWidth(2);
  gBeta2->GetXaxis()->SetLimits(0., 800.);
  gBeta2->GetXaxis()->SetTitle("Energy (keV)");
  gBeta2->GetYaxis()->SetTitle("Asymmetry");
  gBeta2->GetXaxis()->CenterTitle();
  gBeta2->GetYaxis()->CenterTitle();
  
  TF1 *fitBeta2 = new TF1("fitBeta2","[0]",Elow, Ehigh);
  fitBeta2->SetLineColor(kRed);
  fitBeta2->SetLineWidth(3);
  fitBeta2->SetParameter(0,-0.12);
  
  gBeta2->Fit("fitBeta2","R");

  asymFile << "BetaCorrA_binSummed\t" << fitBeta2->GetParameter(0) << "\t" << fitBeta2->GetParError(0) << std::endl;

  
  gBeta2->Draw("AP");
  gBeta2->SetMinimum(fitBeta2->GetParameter(0)-0.04);//(simulation && !AsymmOn) ? -0.5 : -0.16);
  gBeta2->SetMaximum(fitBeta2->GetParameter(0)+0.05);//(simulation && !AsymmOn) ? 0.5 : -0.07);
  c1->Update();
  c1->Print(pdfFile.c_str());

  asymFile.close();

  //Write out corrected (but energy dependent) bin-by-bin asymmetry for use in calculating effect from doing corrections
  outFile = basePath + "Asymmetries/"+(UNBLIND?"UNBLINDED_":"") + corr + "_" + (withPOL?"withPOL_":"") + groupType +"Asymmetries_AnaCh" + anaChoice 
    + std::string("_")  + "Octets_" +itos(octBegin)+"-"+itos(octEnd)+( key/10==1 ? "_Quadrant"+itos(key-10) : ( key/10==2 ? "_RadialRing"+itos(key-20) :"" ) ) ;
  txtFile = outFile+std::string("_BinByBin_withEnergyDependence.txt");
  asymFile.open(txtFile.c_str());

  for (UInt_t n=0; n<enBinMedian.size(); n++) {
    asymFile << enBinMedian[n] << "\t" << groupRawAsymByBin[1][n] << "\t" << groupRawAsymByBin[2][n] << "\n";
  }
  asymFile.close();

  //Write out corrected bin-by-bin asymmetry for use in calculating effect from doing corrections
  txtFile = outFile+std::string("_BinByBin.txt");
  asymFile.open(txtFile.c_str());

  for (UInt_t n=0; n<enBinMedian.size(); n++) {
    asymFile << enBinMedian[n] << "\t" << groupAsymByBin[1][n] << "\t" << groupAsymByBin[2][n] << "\n";
  }
  asymFile.close();
  
  
  delete c1; delete g; delete fit; delete gBeta; delete fitBeta; 
  delete g2; delete gBeta2; delete fitBeta2;

};


//This will create asymmetry plots for each pair, quartet, and octet, and also write the raw integrated asymmetry and Beta corrected integrated asymmetry
// to file for each grouping to be used by PlotFinalAsymmetries

void PlotAsymmetriesByGrouping(std::string groupType, Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t Elow, Double_t Ehigh, Double_t enBinWidth, bool UKdata, bool simulation, bool UNBLIND,int key) {
  if (groupType!="Quartet" && groupType!="Octet" && groupType!="Pair") throw "Bad group type given to PlotAsymmetriesByGrouping. Options are \"Octet\", \"Quartet\", or \"Pair\""; 
  
  key = ( groupType==std::string("Octet") ? key : 0 );

  std::string basePath = simulation ? getenv("SIM_ANALYSIS_RESULTS") : UKdata ? getenv("ANALYSIS_RESULTS") : getenv("MPM_ANALYSIS_RESULTS");

  for (Int_t octet=octBegin; octet<=octEnd; octet++) {

    if (std::find(badOct.begin(), badOct.end(),octet) != badOct.end()) {  continue; } //Checking if octet should be ignored for data quality reasons

    std::ifstream infile;

    if (groupType=="Octet")
    {
      std::string basePath2 =basePath+ "Octet_"+itos(octet)+"/" + groupType + "Asymmetry/"; 
      std::string infilePath = basePath2 + (UNBLIND?"UNBLINDED_":"") + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice +( key/10==1 ? "_Quadrant"+itos(key-10) : ( key/10==2 ? "_RadialRing"+itos(key-20) :"" ) ) +  ".dat";
      std::string outfilePath = basePath2 + (UNBLIND?"UNBLINDED_":"")+ corr +"_" + (withPOL?"withPOL_":"") +"FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ( key/10==1 ? "_Quadrant"+itos(key-10) : ( key/10==2 ? "_RadialRing"+itos(key-20) :"" ) ) + ".dat";
      
      //Remove old output files in case they were created on accident and aren't filled with good values
      std::string command = "rm " + outfilePath; 
      system(command.c_str());

      //First check that Octet was good
      std::string checkStatus;
      infile.open(infilePath.c_str());

      if ( infile.is_open() ) {
	infile >> checkStatus; 
	infile.close();
	if (checkStatus=="BAD") continue;
      }
      else {
	std::cout << "Could not open binned Asymmetries for Octet " << octet << " so skipping...\n";
	continue; 
      }

      std::vector < std::vector <Double_t > > AsymAndError;
      std::vector < std::vector <Double_t > > RawAsymAndError;
      std::vector < Double_t > enBinMedian;
      
      AsymAndError.resize(2,std::vector <Double_t> (0));
      RawAsymAndError.resize(2,std::vector <Double_t> (0));
	
      Double_t eBinLow, Asym, AsymError;

      infile.open(infilePath.c_str());

      Int_t bin = 0;
      
      while (infile >> eBinLow >> Asym >> AsymError) {
	Double_t En = eBinLow+enBinWidth/2.;
	enBinMedian.push_back(En);
	Double_t Beta = returnBeta(En);
	AsymAndError[0].push_back(Asym*2./Beta); 
	AsymAndError[1].push_back(AsymError*2./Beta);
	RawAsymAndError[0].push_back(Asym);
	RawAsymAndError[1].push_back(AsymError);
	bin++;
      }
      
      infile.close();

      //Read in and apply the systematic corrections for this octet
      std::vector < std::vector <Double_t> > deltaSys = LoadOctetSystematics(octet,anaChoice,enBinMedian);

      //Loading theory systematics... 
      std::vector <Double_t> theoryCorr = LoadTheoryCorrections(enBinMedian);

      for (UInt_t i=0 ; i<AsymAndError[0].size() ; i++) {
	AsymAndError[0][i] *= (deltaSys[i][0] / theoryCorr[i]); //Here is where the corrections to Ameas are made.. Need to add in delta theory
	AsymAndError[1][i] *= (deltaSys[i][0] / theoryCorr[i]);

	if ( withPOL ) {

	  AsymAndError[0][i] /= ( POL_ave ); 
	  AsymAndError[1][i] /= ( POL_ave );

	}	  
      }


      std::string pdfPath = basePath2 + (UNBLIND?"UNBLINDED_":"") + corr + "_" +  "Asymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ( key/10==1 ? "_Quadrant"+itos(key-10) : ( key/10==2 ? "_RadialRing"+itos(key-20) :"" ) ) + ".pdf";
      TCanvas *c1 = new TCanvas("c1", "c1",800, 400);
      gStyle->SetOptFit(1111);
      gStyle->SetTitleX(0.25);
	
      TGraphErrors *gOct = new TGraphErrors(enBinMedian.size(), &enBinMedian[0], &RawAsymAndError[0][0], 0, &RawAsymAndError[1][0]);
      std::string title = std::string("Raw Asymmetry A_{SR} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
      gOct->SetTitle(title.c_str());
      gOct->SetMarkerStyle(20);
      gOct->SetLineWidth(2);
      gOct->GetXaxis()->SetLimits(0., 800.);
      gOct->GetXaxis()->SetTitle("Energy (keV)");
      gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
      gOct->GetXaxis()->CenterTitle();
      gOct->GetYaxis()->CenterTitle();
	
      TF1 *fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
      fitOct->SetLineColor(kRed);
      fitOct->SetLineWidth(3);
      fitOct->SetParameter(0,-0.05);
	
      gOct->Fit("fitOct","R");
	
      std::ofstream ofile(outfilePath.c_str());
      ofile << "RawA_SR " << (fitOct->GetParameter(0)) << " " << fitOct->GetParError(0) << std::endl;
	
      gOct->Draw("AP");
      gOct->SetMinimum(-0.1);
      gOct->SetMaximum(0.);
      c1->Update();
      c1->Print(pdfPath.c_str());
	
	
      delete c1; delete gOct; delete fitOct;
	
      pdfPath = basePath2 + (UNBLIND?"UNBLINDED_":"") + corr + "_"  + (withPOL?"withPOL_":"") + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ( key/10==1 ? "_Quadrant"+itos(key-10) : ( key/10==2 ? "_RadialRing"+itos(key-20) :"" ) ) + ".pdf";
      c1 = new TCanvas("c1", "c1",800, 400);
	
      gOct = new TGraphErrors(enBinMedian.size(), &enBinMedian[0], &AsymAndError[0][0], 0, &AsymAndError[1][0]);
      title = std::string("#frac{2A_{SR}}{#beta}#upoint(1+#Delta_{exp}) ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
      gOct->SetTitle(title.c_str());
      gOct->SetMarkerStyle(20);
      gOct->SetLineWidth(2);
      gOct->GetXaxis()->SetLimits(0., 800.);
      gOct->GetXaxis()->SetTitle("Energy (keV)");
      gOct->GetYaxis()->SetTitle("Corrected Asymmetry");
      gOct->GetXaxis()->CenterTitle();
      gOct->GetYaxis()->CenterTitle();
	
      fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
      fitOct->SetLineColor(kRed);
      fitOct->SetLineWidth(3);
      fitOct->SetParameter(0,-0.12);
	
      gOct->Fit("fitOct","R");
	
      gOct->Draw("AP");
      gOct->SetMinimum(-0.6);
      gOct->SetMaximum(0.4);
      c1->Update();
      c1->Print(pdfPath.c_str());
	
	
      ofile << "2*A/Beta " << fitOct->GetParameter(0) << " " << fitOct->GetParError(0);
      ofile.close();
	
      delete c1; delete gOct; delete fitOct;
      
    }

    //Now the Quartets
    else if (groupType=="Quartet")
    {	  

      std::string basePath2 =basePath+ "Octet_"+itos(octet)+"/" + groupType + "Asymmetry/"; 

      std::string infilePath[2];
      infilePath[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_A" + ".dat";
      infilePath[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_B" + ".dat";

      std::string outfilePath[2];
      outfilePath[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_A" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) +".dat";
      outfilePath[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_B" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) +".dat";

      std::string pdfPathRaw[2];
      pdfPathRaw[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_A" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
      pdfPathRaw[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_B" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";

      std::string pdfPathCorr[2];
      pdfPathCorr[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_A" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
      pdfPathCorr[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_B" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
      
      for (int quart=0; quart<2; quart++) {

	std::string currentQuart = quart==0?"A":"B";

	//Remove old output files in case they were created on accident and aren't filled with good values
	std::string command = "rm " + outfilePath[quart]; 
	system(command.c_str());

	//First check that Quartet was good
	std::string checkStatus;
	infile.open(infilePath[quart].c_str());
	
	if ( infile.is_open() ) {
	  infile >> checkStatus; 
	  infile.close();
	  if (checkStatus=="BAD") continue;
	}
	else {
	  std::cout << "Could not open binned Asymmetries for Quartet " << currentQuart << " in Octet " << octet << " so skipping...\n";
	  continue; 
	}
	    
	infile.open(infilePath[quart].c_str());
      
	std::vector < std::vector <Double_t > > AsymAndError;
	std::vector < std::vector <Double_t > > RawAsymAndError;
	std::vector < Double_t > enBinMedian;
	
	AsymAndError.resize(2,std::vector <Double_t> (0));
	RawAsymAndError.resize(2,std::vector <Double_t> (0));
	
	Double_t eBinLow, Asym, AsymError;
	
	while (infile >> eBinLow >> Asym >> AsymError) {
	  Double_t En = eBinLow+enBinWidth/2.;
	  enBinMedian.push_back(En);
	  Double_t Beta = returnBeta(En);
	  AsymAndError[0].push_back(Asym*2./Beta);
	  AsymAndError[1].push_back(AsymError*2./Beta);
	  RawAsymAndError[0].push_back(Asym);
	  RawAsymAndError[1].push_back(AsymError);
	}

	infile.close();
	
	TCanvas *c1 = new TCanvas("c1", "c1",800, 400);
	gStyle->SetOptFit(1111);
	gStyle->SetTitleX(0.25);
      
	TGraphErrors *gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &RawAsymAndError[0][5], 0, &RawAsymAndError[1][5]);
	std::string title = std::string("Raw Asymmetry A_{SR} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
	gOct->SetTitle(title.c_str());
	gOct->SetMarkerStyle(20);
	gOct->SetLineWidth(2);
	gOct->GetXaxis()->SetLimits(0., 800.);
	gOct->GetXaxis()->SetTitle("Energy (keV)");
	gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
	gOct->GetXaxis()->CenterTitle();
	gOct->GetYaxis()->CenterTitle();

	TF1 *fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
	fitOct->SetLineColor(kRed);
	fitOct->SetLineWidth(3);
	fitOct->SetParameter(0,-0.05);
	
	gOct->Fit("fitOct","R");
	
	std::ofstream ofile(outfilePath[quart].c_str());
	ofile << "RawA_SR " << -(fitOct->GetParameter(0)) << " " << fitOct->GetParError(0) << std::endl;
      
	gOct->Draw("AP");
	gOct->SetMinimum(-0.1);
	gOct->SetMaximum(0.);
	c1->Update();
	c1->Print(pdfPathRaw[quart].c_str());
	
	delete c1; delete gOct; delete fitOct;
	
	c1 = new TCanvas("c1", "c1",800, 400);
	
	gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &AsymAndError[0][5], 0, &AsymAndError[1][5]);
	title = std::string("#frac{2A_{SR}}{#beta} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
	gOct->SetTitle(title.c_str());
	gOct->SetMarkerStyle(20);
	gOct->SetLineWidth(2);
	gOct->GetXaxis()->SetLimits(0., 800.);
	gOct->GetXaxis()->SetTitle("Energy (keV)");
	gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
	gOct->GetXaxis()->CenterTitle();
	gOct->GetYaxis()->CenterTitle();
	
        fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
	fitOct->SetLineColor(kRed);
	fitOct->SetLineWidth(3);
	fitOct->SetParameter(0,-0.12);
	
	gOct->Fit("fitOct","R");
	
	gOct->Draw("AP");
	gOct->SetMinimum(-0.6);
	gOct->SetMaximum(0.4);
	c1->Update();
	c1->Print(pdfPathCorr[quart].c_str());

	ofile << "2*A/Beta " << fitOct->GetParameter(0) << " " << fitOct->GetParError(0);
	ofile.close();

	delete c1; delete gOct; delete fitOct;
      }
    }

    //And last of all the Pairs
    else if (groupType=="Pair")
    {
      std::string basePath2 =basePath+ "Octet_"+itos(octet)+"/" + groupType + "Asymmetry/"; 

      for (int pair=0; pair<2; pair++) {
	std::string infilePath[2];
	infilePath[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") +"rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_A" + itos(pair) + ".dat";
	infilePath[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_B" + itos(pair) + ".dat";

	std::string outfilePath[2];
	outfilePath[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_A" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".dat";
	outfilePath[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_B" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".dat";

	std::string pdfPathRaw[2];
	pdfPathRaw[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_A" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
	pdfPathRaw[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_B" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
	
	std::string pdfPathCorr[2];
	pdfPathCorr[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_A" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
	pdfPathCorr[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_B" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
	

	//Remove old output files in case they were created on accident and aren't filled with good values
	std::string command = "rm " + outfilePath[0]; 
	system(command.c_str());
	command = "rm " + outfilePath[1]; 
	system(command.c_str());
	
	for (int quart=0; quart<2; quart++) {
	  
	  std::string currentPair = (quart==0?"A":"B") + itos(pair);	  

	  //First check that pair was good
	  std::string checkStatus;
	  infile.open(infilePath[quart].c_str());
	  
	  if ( infile.is_open() ) {
	    infile >> checkStatus; 
	    infile.close();
	    if (checkStatus=="BAD") continue;
	  }
	  else {
	    std::cout << "Could not open binned Asymmetries for Pair " << currentPair << " in Octet " << octet << " so skipping...\n";
	    continue; 
	  }


	 
	  
	  infile.open(infilePath[quart].c_str());
	  
	  std::vector < std::vector <Double_t > > AsymAndError;
	  std::vector < std::vector <Double_t > > RawAsymAndError;
	  std::vector < Double_t > enBinMedian;
	  
	  AsymAndError.resize(2,std::vector <Double_t> (0));
	  RawAsymAndError.resize(2,std::vector <Double_t> (0));
	  
	  Double_t eBinLow, Asym, AsymError;
	
	  while (infile >> eBinLow >> Asym >> AsymError) {
	    Double_t En = eBinLow+enBinWidth/2.;
	    enBinMedian.push_back(En);
	    Double_t Beta = returnBeta(En);
	    AsymAndError[0].push_back(Asym*2./Beta);
	    AsymAndError[1].push_back(AsymError*2./Beta);
	    RawAsymAndError[0].push_back(Asym);
	    RawAsymAndError[1].push_back(AsymError);
	  }

	  infile.close();
	  
	  
	  TCanvas *c1 = new TCanvas("c1", "c1",800, 400);
	  gStyle->SetOptFit(1111);
	  gStyle->SetTitleX(0.25);
      
	  TGraphErrors *gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &RawAsymAndError[0][5], 0, &RawAsymAndError[1][5]);
	  std::string title = std::string("Raw Asymmetry A_{SR} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
	  gOct->SetTitle(title.c_str());
	  gOct->SetMarkerStyle(20);
	  gOct->SetLineWidth(2);
	  gOct->GetXaxis()->SetLimits(0., 800.);
	  gOct->GetXaxis()->SetTitle("Energy (keV)");
	  gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
	  gOct->GetXaxis()->CenterTitle();
	  gOct->GetYaxis()->CenterTitle();

	  TF1 *fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
	  fitOct->SetLineColor(kRed);
	  fitOct->SetLineWidth(3);
	  fitOct->SetParameter(0,-0.05);
	  
	  gOct->Fit("fitOct","R");
	  
	  std::ofstream ofile(outfilePath[quart].c_str());
	  ofile << "RawA_SR " << -(fitOct->GetParameter(0)) << " " << fitOct->GetParError(0) << std::endl;
	  
	  gOct->Draw("AP");
	  gOct->SetMinimum(-0.1);
	  gOct->SetMaximum(0.);
	  c1->Update();
	  c1->Print(pdfPathRaw[quart].c_str());
	  
	  delete c1; delete gOct; delete fitOct;
	  
	  c1 = new TCanvas("c1", "c1",800, 400);
	  
	  gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &AsymAndError[0][5], 0, &AsymAndError[1][5]);
	  title = std::string("#frac{2A_{SR}}{#beta} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
	  gOct->SetTitle(title.c_str());
	  gOct->SetMarkerStyle(20);
	  gOct->SetLineWidth(2);
	  gOct->GetXaxis()->SetLimits(0., 800.);
	  gOct->GetXaxis()->SetTitle("Energy (keV)");
	  gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
	  gOct->GetXaxis()->CenterTitle();
	  gOct->GetYaxis()->CenterTitle();
	  
	  fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
	  fitOct->SetLineColor(kRed);
	  fitOct->SetLineWidth(3);
	  fitOct->SetParameter(0,-0.12);
	  
	  gOct->Fit("fitOct","R");
	  
	  gOct->Draw("AP");
	  gOct->SetMinimum(-0.6);
	  gOct->SetMaximum(0.4);
	  c1->Update();
	  c1->Print(pdfPathCorr[quart].c_str());

	  ofile << "2*A/Beta " << fitOct->GetParameter(0) << " " << fitOct->GetParError(0);
	  ofile.close();
	  
	  delete c1; delete gOct; delete fitOct;
	}
      }
    }
  }
};
      
  

//These are summed over the energy range given as Elow-Ehigh

void ProduceRawAsymmetries(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t Elow, Double_t Ehigh, Double_t enBinWidth, bool UKdata, bool simulation) {
 
    //Looking at octet evts and asymmetries
  std::string basePath = simulation ? getenv("SIM_ANALYSIS_RESULTS") : UKdata ? getenv("ANALYSIS_RESULTS") : getenv("MPM_ANALYSIS_RESULTS");
  
  std::string pairFile = basePath + std::string("Asymmetries/Pair_RawAsymmetries_AnaCh") + anaChoice 
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".dat");
  std::string quartetFile = basePath + std::string("Asymmetries/Quartet_RawAsymmetries_AnaCh") + anaChoice
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".dat");
  std::string octetFile = basePath + std::string("Asymmetries/Octet_RawAsymmetries_AnaCh") + anaChoice 
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".dat");

    
    //unsigned int octetNum = 1;
  std::ofstream octAsym(octetFile.c_str());
  std::ofstream quartAsym(quartetFile.c_str());
  std::ofstream pairAsym(pairFile.c_str());
 
  Int_t numPair=0, numQuart=0;
  std::vector < std::vector <Double_t > > octetRawAsymAndError;
  std::vector < std::vector <Double_t > > quartetRawAsymAndError;
  std::vector < std::vector <Double_t > > pairRawAsymAndError;
  octetRawAsymAndError.resize(3,std::vector <Double_t> (0));
  quartetRawAsymAndError.resize(3,std::vector <Double_t> (0));
  pairRawAsymAndError.resize(3,std::vector <Double_t> (0)); 
  
  Double_t Asym=0., AsymError=0.;

  for (Int_t octet=octBegin; octet<=octEnd; octet++) {

    if (std::find(badOct.begin(), badOct.end(),octet) != badOct.end()) { numPair+=4; numQuart+=2; continue; } //Checking if octet should be ignored for data quality reasons

    try {
      //OCTETS
      OctetAsymmetry oct(octet,anaChoice,enBinWidth, 50., UKdata, simulation);     
      oct.calcTotalAsymmetry(Elow,Ehigh);
      //oct.calcAsymmetryBinByBin(1);
      //oct.calcTotalAsymmetry(170.,630.,1);
      
      Asym = oct.returnTotalAsymmetry();
      AsymError = oct.returnTotalAsymmetryError();
      octetRawAsymAndError[0].push_back(octet);
      octetRawAsymAndError[1].push_back(Asym);
      octetRawAsymAndError[2].push_back(AsymError);

      std::cout << "Asymmetry for Octet " << octet << ":\n";
      std::cout << Asym << " +- " << AsymError << std::endl;
      octAsym << octet << " " << Asym << " " << AsymError << "\n";
     
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
    
    try {
      // QUARTETS
      QuartetAsymmetry quart(octet,anaChoice,enBinWidth, 50., UKdata, simulation); 
      quart.calcTotalAsymmetry(Elow,Ehigh);
	//quart.calcAsymmetryBinByBin(1);
	//quart.calcTotalAsymmetry(170.,630.,1);
      
      if (quart.boolGoodQuartet(0)) {
	
	Asym = quart.returnTotalAsymmetry_QuartetA();
	AsymError = quart.returnTotalAsymmetryError_QuartetA();
	quartetRawAsymAndError[0].push_back(numQuart);
	quartetRawAsymAndError[1].push_back(Asym);
	quartetRawAsymAndError[2].push_back(AsymError);
  
	std::cout << "Asymmetry for Quartet " << numQuart << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	quartAsym << numQuart << " " << Asym << " " << AsymError << "\n";
      }
      numQuart++;

      if (quart.boolGoodQuartet(1)) {
	Asym = quart.returnTotalAsymmetry_QuartetB();
	AsymError = quart.returnTotalAsymmetryError_QuartetB();
	quartetRawAsymAndError[0].push_back(numQuart);
	quartetRawAsymAndError[1].push_back(Asym);
	quartetRawAsymAndError[2].push_back(AsymError);

	std::cout << "Asymmetry for Quartet " << numQuart << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	quartAsym << numQuart << " " << Asym << " " << AsymError << "\n";
      }
      numQuart++;
    }
    catch(const char* ex){
      numQuart+=2;
      std::cerr << "Error: " << ex << std::endl;
    }

    try {
      // PAIRS
      PairAsymmetry pair(octet,anaChoice,enBinWidth, 50., UKdata, simulation); 
      pair.calcTotalAsymmetry(Elow,Ehigh);
	//pair.calcAsymmetryBinByBin(1);
	//pair.calcTotalAsymmetry(170.,630.,1);
      
      if (pair.boolGoodPair(0,0)) {
	
	Asym = pair.returnTotalAsymmetry_PairA0();
	AsymError = pair.returnTotalAsymmetryError_PairA0();
	pairRawAsymAndError[0].push_back(numPair);
	pairRawAsymAndError[1].push_back(Asym);
	pairRawAsymAndError[2].push_back(AsymError);
	
	std::cout << "Asymmetry for Pair " << numPair << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	pairAsym << numPair << " " << Asym << " " << AsymError << "\n";
      }
      numPair++;

      if (pair.boolGoodPair(0,1)) {
	Asym = pair.returnTotalAsymmetry_PairA1();
	AsymError = pair.returnTotalAsymmetryError_PairA1();
	pairRawAsymAndError[0].push_back(numPair);
	pairRawAsymAndError[1].push_back(Asym);
	pairRawAsymAndError[2].push_back(AsymError);

	std::cout << "Asymmetry for Pair " << numPair << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	pairAsym << numPair << " " << Asym << " " << AsymError << "\n";
      }
      numPair++;

      if (pair.boolGoodPair(1,0)) {
	
	Asym = pair.returnTotalAsymmetry_PairB0();
	AsymError = pair.returnTotalAsymmetryError_PairB0();
	pairRawAsymAndError[0].push_back(numPair);
	pairRawAsymAndError[1].push_back(Asym);
	pairRawAsymAndError[2].push_back(AsymError);

	std::cout << "Asymmetry for Pair " << numPair << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	pairAsym << numPair << " " << Asym << " " << AsymError << "\n";
      }
      numPair++;

      if (pair.boolGoodPair(1,1)) {
	Asym = pair.returnTotalAsymmetry_PairB1();
	AsymError = pair.returnTotalAsymmetryError_PairB1();
	pairRawAsymAndError[0].push_back(numPair);
	pairRawAsymAndError[1].push_back(Asym);
	pairRawAsymAndError[2].push_back(AsymError);

	std::cout << "Asymmetry for Pair " << numPair << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	pairAsym << numPair << " " << Asym << " " << AsymError << "\n";
      }
      numPair++;
    }
    catch(const char* ex){
      numPair+=4;
      std::cerr << "Error: " << ex << std::endl;
    }
  }

  octAsym.close();
  quartAsym.close();
  pairAsym.close();

  std::string pdfFile = basePath + std::string("Asymmetries/RawAsymmetries_AnaCh") + anaChoice 
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".pdf");

  TCanvas *c1 = new TCanvas("c1", "c1", 1000., 1000.);
  c1->Divide(1,3);

  gStyle->SetOptFit(1111);
  gStyle->SetTitleX(0.25);
 
  std::vector <double> errorX(pairRawAsymAndError[0].size(),0.);
  std::string title;

  c1->cd(1);

  TGraphErrors *gOct = new TGraphErrors(octetRawAsymAndError[0].size(), &octetRawAsymAndError[0][0],&octetRawAsymAndError[1][0],&errorX[0], &octetRawAsymAndError[2][0]);
  title = std::string("Octet-By-Octet Raw Measured Asymmetry ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gOct->SetTitle(title.c_str());
  gOct->SetMarkerStyle(20);
  gOct->SetLineWidth(2);
  gOct->GetXaxis()->SetLimits(-2., octetRawAsymAndError[0][octetRawAsymAndError[0].size()-1]+2.);
  gOct->GetXaxis()->SetTitle("Octet Number");
  gOct->GetYaxis()->SetTitle("Raw Asymmetry");
  gOct->GetXaxis()->CenterTitle();
  gOct->GetYaxis()->CenterTitle();
  
  TF1 *fitOct = new TF1("fitOct","[0]",octetRawAsymAndError[0][0]-1., octetRawAsymAndError[0][octetRawAsymAndError[0].size()-1]+1.);
  fitOct->SetLineColor(kRed);
  fitOct->SetLineWidth(3);
  fitOct->SetParameter(0,0.05);

  gOct->Fit("fitOct","R");
  
  gOct->Draw("AP");
  gOct->SetMinimum(0.03);
  gOct->SetMaximum(0.07);
  c1->Update();

  c1->cd(2);

  TGraphErrors *gQuart = new TGraphErrors(quartetRawAsymAndError[0].size(), &quartetRawAsymAndError[0][0],&quartetRawAsymAndError[1][0],&errorX[0], &quartetRawAsymAndError[2][0]);
  title = std::string("Quartet-By-Quartet Raw Measured Asymmetry ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gQuart->SetTitle(title.c_str());
  gQuart->SetMarkerStyle(20);
  gQuart->SetLineWidth(2);
  gQuart->GetXaxis()->SetLimits(-2., quartetRawAsymAndError[0][quartetRawAsymAndError[0].size()-1]+2.);
  gQuart->GetXaxis()->SetTitle("Quartet Number");
  gQuart->GetYaxis()->SetTitle("Raw Asymmetry");
  gQuart->GetXaxis()->CenterTitle();
  gQuart->GetYaxis()->CenterTitle();
  
  TF1 *fitQuart = new TF1("fitQuart","[0]",quartetRawAsymAndError[0][0]-1., quartetRawAsymAndError[0][quartetRawAsymAndError[0].size()-1]+1.);
  fitQuart->SetLineColor(kRed);
  fitQuart->SetLineWidth(3);
  fitQuart->SetParameter(0,0.05);


  gQuart->Fit("fitQuart","R");
  
  gQuart->Draw("AP");
  gQuart->SetMinimum(0.03);
  gQuart->SetMaximum(0.07);
  c1->Update();

  c1->cd(3);

  TGraphErrors *gPair = new TGraphErrors(pairRawAsymAndError[0].size(), &pairRawAsymAndError[0][0],&pairRawAsymAndError[1][0],&errorX[0], &pairRawAsymAndError[2][0]);
  title = std::string("Pair-By-Pair Raw Measured Asymmetry ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gPair->SetTitle(title.c_str());
  gPair->SetMarkerStyle(20);
  gPair->SetLineWidth(2);
  gPair->GetXaxis()->SetLimits(-2., pairRawAsymAndError[0][pairRawAsymAndError[0].size()-1]+2.);
  gPair->GetXaxis()->SetTitle("Pair Number");
  gPair->GetYaxis()->SetTitle("Raw Asymmetry");
  gPair->GetXaxis()->CenterTitle();
  gPair->GetYaxis()->CenterTitle();
  
  TF1 *fitPair = new TF1("fitPair","[0]",pairRawAsymAndError[0][0]-1., pairRawAsymAndError[0][pairRawAsymAndError[0].size()-1]+1.);
  fitPair->SetLineColor(kRed);
  fitPair->SetLineWidth(3);
  fitPair->SetParameter(0,0.05);


  gPair->Fit("fitPair","R");
  
  gPair->Draw("AP");
  gPair->SetMinimum(0.03);
  gPair->SetMaximum(0.07);
  c1->Update();
  
  c1->Print(pdfFile.c_str());

  delete gPair; delete gQuart; delete gOct; delete fitPair; delete fitQuart; delete fitOct; delete c1;
																     

};
  
  
std::vector < std::vector <Double_t> > LoadOctetSystematics(Int_t octet, std::string anaChoice, std::vector <Double_t> enBinMidpoint) {

  Int_t iAnaChoice;

  for (UInt_t i = 0; i<10; i++) {
    
    if (anaChoices[i]==anaChoice) iAnaChoice = i+1;

  }
  
  //TString filename = TString::Format("%s/Octet_%i/OctetAsymmetry/Systematics/ThOverProc_Octet-%i_Analysis-%i.txt",getenv("ANALYSIS_RESULTS"),octet,octet,iAnaChoice);
  TString filename = TString::Format("%s/systematics/MC_Corrections/DeltaExp_OctetByOctetCorrections/ThOverProc_Octet-%i_Analysis-%i.txt",getenv("ANALYSIS_CODE"),octet,iAnaChoice);
  //std::cout << filename.Data() << std::endl;
  std::vector < std::vector <Double_t> > syst(enBinMidpoint.size(), std::vector<Double_t>(2,1.));

  if ( corr!=std::string("DeltaExpOnly") && corr!=std::string("AllCorr") ) return syst;
  
  std::ifstream infile(filename.Data());

  if (!infile.is_open()) throw "Couldn't open file in LoadOctetSystematics!";

  //Read in the header crap
  //  std::string firstline[8];
  //for (int i=0; i<8; i++) {
  //  infile >> firstline[i];
  //  std::cout << firstline[i] << " " ;
  // }
  std::cout << "\n";

  //Read in the systematics
  Double_t mid; //midpoint of bin
  Double_t ratio;  //systematics correction (Apure/Aproc)
  Double_t midErr; //BinWidth/2
  Double_t ratioErr; //Systematic error on ratio

  Int_t it = 0;

  while (infile >> mid >> ratio >> ratioErr) {
    if (mid==enBinMidpoint[it]) {
      syst[it][0] = ratio!=0. ? ratio : 1.;
      syst[it][1] = ratioErr;
      std::cout << mid << " " << ratio << " " << ratioErr << "\n";
      it++;
    }
  }

  std::cout << "Loaded Systematic Errors for Octet " << octet << " ...\n";

  return syst;
};

std::vector <Double_t> LoadTheoryCorrections(std::vector <Double_t> enBinMidpoint) {

  std::vector <Double_t> syst(enBinMidpoint.size(), 1.);

  if ( corr!=std::string("DeltaTheoryOnly") && corr!=std::string("AllCorr") ) return syst;

  for (UInt_t i=0; i<syst.size(); i++) {
   
    syst[i] = asymmetryCorrectionFactor(enBinMidpoint[i]); //As defined in BetaSpectrum.hh by MPM

  }

  return syst;

};
