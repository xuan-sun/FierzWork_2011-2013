/*
Defines different types of Asymmetries to be calculated within an Octet.
ppair, quartet, octet
*/

#include "EvtRateHandler.hh"
#include "Asymmetries.hh"
#include "MBUtils.hh"

#include <fstream>
#include <cstdlib>
#include <cmath>
#include <TMath.h>

void calcBinByBinSuperRatio(std::vector< std::vector<double> > &A2, 
			    std::vector< std::vector<double> > &A5,
			    std::vector< std::vector<double> > &A7, 
			    std::vector< std::vector<double> > &A10,
			    std::vector< std::vector<double> > &B2, 
			    std::vector< std::vector<double> > &B5,
			    std::vector< std::vector<double> > &B7,
			    std::vector< std::vector<double> > &B10,
			    std::vector< std::vector<double> > &A2err, 
			    std::vector< std::vector<double> > &A5err,
			    std::vector< std::vector<double> > &A7err, 
			    std::vector< std::vector<double> > &A10err,
			    std::vector< std::vector<double> > &B2err, 
			    std::vector< std::vector<double> > &B5err,
			    std::vector< std::vector<double> > &B7err, 
			    std::vector< std::vector<double> > &B10err,
			    std::vector<double> &asymmetry, std::vector<double> &asymmetryError);

std::map<int,int> BGrunReplace = {{23017,23006},{21175,21164},{21176,21187}};
// 21175, 21176 beam out in runlog
// 23017 very low statistics

AsymmetryBase::AsymmetryBase(int oct, std::string anaCh, double enBinWidth, double fidCut, bool ukdata, bool simulation, bool unblind) : UKdata(ukdata), Simulation(simulation),  UNBLIND(unblind), octet(oct), energyBinWidth(enBinWidth), numEnBins(0), fiducialCut(fidCut), runsInOctet(0), analysisChoice(anaCh) {
  
  numEnBins = (unsigned int)(1200./energyBinWidth);

  if (UNBLIND) {
    std::cout << "***************************************************************\n"
	      << "***************************************************************\n\n"
	      << "                UNBLINDING OCTET " << " " << octet << "\n\n"
	      << "***************************************************************\n"
	      << "***************************************************************\n\n";
  }


  A2.resize(2,std::vector <double> (numEnBins,0.));
  B2.resize(2,std::vector <double> (numEnBins,0.));
  A5.resize(2,std::vector <double> (numEnBins,0.));
  B5.resize(2,std::vector <double> (numEnBins,0.));
  A7.resize(2,std::vector <double> (numEnBins,0.));
  B7.resize(2,std::vector <double> (numEnBins,0.));
  A10.resize(2,std::vector <double> (numEnBins,0.));
  B10.resize(2,std::vector <double> (numEnBins,0.));
  A2err.resize(2,std::vector <double> (numEnBins,0.));
  B2err.resize(2,std::vector <double> (numEnBins,0.));
  A5err.resize(2,std::vector <double> (numEnBins,0.));
  B5err.resize(2,std::vector <double> (numEnBins,0.));
  A7err.resize(2,std::vector <double> (numEnBins,0.));
  B7err.resize(2,std::vector <double> (numEnBins,0.));
  A10err.resize(2,std::vector <double> (numEnBins,0.));
  B10err.resize(2,std::vector <double> (numEnBins,0.));

  A2Quad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B2Quad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A5Quad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B5Quad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A7Quad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B7Quad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A10Quad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B10Quad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A2errQuad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B2errQuad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A5errQuad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B5errQuad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A7errQuad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B7errQuad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A10errQuad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B10errQuad.resize(4,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  
  A2Rad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B2Rad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A5Rad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B5Rad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A7Rad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B7Rad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A10Rad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B10Rad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A2errRad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B2errRad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A5errRad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B5errRad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A7errRad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B7errRad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  A10errRad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  B10errRad.resize(6,std::vector< std::vector<double> > (2,std::vector <double> (numEnBins,0.)));
  
  A2len.resize(2,std::vector <double> (2,0.));
  A5len.resize(2,std::vector <double> (2,0.));
  A7len.resize(2,std::vector <double> (2,0.));
  A10len.resize(2,std::vector <double> (2,0.));
  B2len.resize(2,std::vector <double> (2,0.));
  B5len.resize(2,std::vector <double> (2,0.));
  B7len.resize(2,std::vector <double> (2,0.));
  B10len.resize(2,std::vector <double> (2,0.));

  binLowerEdge.resize(numEnBins,0.);
  binUpperEdge.resize(numEnBins,0.);
  
  for (unsigned int i=0; i<numEnBins; i++) {
    binLowerEdge[i] = (double)(i)*energyBinWidth;
    binUpperEdge[i] = (double)(i+1)*energyBinWidth;
  }  

  readOctetFile();
 
};

void AsymmetryBase::readOctetFile() {
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  int numRuns = 0;
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {
    runType[runNumberHold] = runTypeHold;
    numRuns++;
  }
  infile.close();
  runsInOctet = numRuns;
  std::cout << "Read in octet file for octet " << octet << " with " << runsInOctet << " runs\n";
};

int AsymmetryBase::getBGrunNumber(std::string type) {
 
  std::string bgType = type=="A2"?"A1" : ( type=="A5"?"A4" : ( type=="A7"?"A9" : ( type=="A10"?"A12" : ( type=="B2"?"B1" : ( type=="B5"?"B4" : ( type=="B7"?"B9" : ( type=="B10"?"B12" : "BAD" )))))));
  
  if (bgType=="BAD") { throw "PASSED RUN WITH BAD TYPE TO getBGrun() IN ASYMMETRY BASE!"; }
  //std::map <int, std::string>::iterator it;
  for ( auto& it: runType ) { //= runType.begin();it!=runType.end(); it++) {

    if ( it.second==bgType ) {
      if ( BGrunReplace.find(it.first)==BGrunReplace.end() ) return it.first;
      else return BGrunReplace[it.first];
    }
  }

  std::cout<< "NO BACKGROUND RUN IN BG SUBTRACTED RATE FOR RUN " << type << "!\n";
  throw "FAILED IN getBGrunNumber()";
};


int AsymmetryBase::getBGrunNumber(int run) {
  std::string type = runType.at(run);
  std::string bgType = type=="A2"?"A1" : ( type=="A5"?"A4" : ( type=="A7"?"A9" : ( type=="A10"?"A12" : ( type=="B2"?"B1" : ( type=="B5"?"B4" : ( type=="B7"?"B9" : ( type=="B10"?"B12" : "BAD" )))))));
  
  if (bgType=="BAD") { throw "PASSED RUN WITH BAD TYPE TO getBGrun() IN ASYMMETRY BASE!"; }
  //std::map <int, std::string>::iterator it;
  for ( auto& it: runType ) { //= runType.begin();it!=runType.end(); it++) {

    if ( it.second==bgType ) {
      if ( BGrunReplace.find(it.first)==BGrunReplace.end() ) return it.first;
      else return BGrunReplace[it.first];
    }
  }

  std::cout<< "NO BACKGROUND RUN IN BG SUBTRACTED RATE FOR RUN " << run << "!\n";
  throw "FAILED IN getBGrun(int run)";
};


std::string AsymmetryBase::getBGrunLabel(std::string type) {
 
  std::string bgType = type=="A2"?"A1" : ( type=="A5"?"A4" : ( type=="A7"?"A9" : ( type=="A10"?"A12" : ( type=="B2"?"B1" : ( type=="B5"?"B4" : ( type=="B7"?"B9" : ( type=="B10"?"B12" : "BAD" )))))));
  
  if (bgType=="BAD") { throw "PASSED RUN WITH BAD TYPE TO getBGrunLabel() IN ASYMMETRY BASE!"; }
  else return bgType;
 
};


bool AsymmetryBase::isFullOctet() {
  bool A2=false, B2=false, A5=false, B5=false, A7=false, B7=false, A10=false, B10=false;
  bool A2bg=false, B2bg=false, A5bg=false, B5bg=false, A7bg=false, B7bg=false, A10bg=false, B10bg=false;
  std::map <int, std::string>::iterator it;
  for ( it = runType.begin();it!=runType.end(); it++) {
    if (it->second=="A2") {A2=true; continue;}
    if (it->second=="B2") {B2=true; continue;}
    if (it->second=="A5") {A5=true; continue;}
    if (it->second=="B5") {B5=true; continue;}
    if (it->second=="A7") {A7=true; continue;}
    if (it->second=="B7") {B7=true; continue;}
    if (it->second=="A10") {A10=true; continue;}
    if (it->second=="B10") {B10=true; continue;}
    if (it->second=="A1") {A2bg=true; continue;}
    if (it->second=="B1") {B2bg=true; continue;}
    if (it->second=="A4") {A5bg=true; continue;}
    if (it->second=="B4") {B5bg=true; continue;}
    if (it->second=="A9") {A7bg=true; continue;}
    if (it->second=="B9") {B7bg=true; continue;}
    if (it->second=="A12") {A10bg=true; continue;}
    if (it->second=="B12") {B10bg=true; continue;}
  }
  return (A2 && B2 && A5 && B5 && A7 && B7 && A10 && B10 && A2bg && B2bg && A5bg && B5bg && A7bg && B7bg && A10bg && B10bg);
};

bool AsymmetryBase::isPair(int quartNum, int pairNum) {
  bool A2=false, B2=false, A5=false, B5=false, A7=false, B7=false, A10=false, B10=false;
  bool A2bg=false, B2bg=false, A5bg=false, B5bg=false, A7bg=false, B7bg=false, A10bg=false, B10bg=false;
  std::map<int,std::string>::iterator it;
  for ( it = runType.begin();it!=runType.end(); it++) {
    if (it->second=="A2") {A2=true; continue;}
    if (it->second=="B2") {B2=true; continue;}
    if (it->second=="A5") {A5=true; continue;}
    if (it->second=="B5") {B5=true; continue;}
    if (it->second=="A7") {A7=true; continue;}
    if (it->second=="B7") {B7=true; continue;}
    if (it->second=="A10") {A10=true; continue;}
    if (it->second=="B10") {B10=true; continue;}
    if (it->second=="A1") {A2bg=true; continue;}
    if (it->second=="B1") {B2bg=true; continue;}
    if (it->second=="A4") {A5bg=true; continue;}
    if (it->second=="B4") {B5bg=true; continue;}
    if (it->second=="A9") {A7bg=true; continue;}
    if (it->second=="B9") {B7bg=true; continue;}
    if (it->second=="A12") {A10bg=true; continue;}
    if (it->second=="B12") {B10bg=true; continue;}
  }
  if (quartNum==0) {
    if (pairNum==0) return (A2 && A5 && A2bg && A5bg);
    if (pairNum==1) return (A7 && A10 && A7bg && A10bg);
    else throw "BAD PAIR NUMBER GIVEN in AsymmetryBase::isPair(int QuartNum=0,1, int pairNum=0,1)";
  }
  else if (quartNum==1) {
    if (pairNum==0) return (B2 && B5 && B2bg && B5bg);
    if (pairNum==1) return (B7 && B10 && B7bg && B10bg);
    else throw "BAD PAIR NUMBER GIVEN in AsymmetryBase::isPair(int QuartNum=0,1, int pairNum=0,1)";
  }
  else throw "BAD QUARTET NUMBER GIVEN in AsymmetryBase::isPair(int QuartNum=0,1, int pairNum=0,1)";
  
};

bool AsymmetryBase::isFullQuartet(int quartNum) {
  bool A2=false, B2=false, A5=false, B5=false, A7=false, B7=false, A10=false, B10=false;
  bool A2bg=false, B2bg=false, A5bg=false, B5bg=false, A7bg=false, B7bg=false, A10bg=false, B10bg=false;
  std::map<int,std::string>::iterator it;
  for ( it = runType.begin();it!=runType.end(); it++) {
    if (it->second=="A2") {A2=true; continue;}
    if (it->second=="B2") {B2=true; continue;}
    if (it->second=="A5") {A5=true; continue;}
    if (it->second=="B5") {B5=true; continue;}
    if (it->second=="A7") {A7=true; continue;}
    if (it->second=="B7") {B7=true; continue;}
    if (it->second=="A10") {A10=true; continue;}
    if (it->second=="B10") {B10=true; continue;}
    if (it->second=="A1") {A2bg=true; continue;}
    if (it->second=="B1") {B2bg=true; continue;}
    if (it->second=="A4") {A5bg=true; continue;}
    if (it->second=="B4") {B5bg=true; continue;}
    if (it->second=="A9") {A7bg=true; continue;}
    if (it->second=="B9") {B7bg=true; continue;}
    if (it->second=="A12") {A10bg=true; continue;}
    if (it->second=="B12") {B10bg=true; continue;}
  }
  if (quartNum==0) return (A2 && A5 && A7 && A10 && A2bg && A5bg && A7bg && A10bg);
  else if (quartNum==1) return (B2 && B5 && B7 && B10 && B2bg && B5bg && B7bg && B10bg);
  else throw "BAD QUARTET NUMBER GIVEN in AsymemtryBase::isQuartet(int QuartNum)";
  
};


std::vector <int> AsymmetryBase::makeRunVec(std::string rnType) {

  std::vector <int> rns;
   
  for ( auto& it: runType ) {

    if ( it.second==rnType ) {
      if ( BGrunReplace.find(it.first)==BGrunReplace.end() ) rns.push_back(it.first);
      else rns.push_back(BGrunReplace[it.first]);
    }
  }

  if (rns.size()>0) return rns;

  else {

    std::cout<< "NO BACKGROUND RUN IN BG SUBTRACTED RATE FOR RUN " << rnType << "!\n";
    throw "FAILED IN makeRunVec()";
  }
};

void AsymmetryBase::loadRates() {
  std::map<int,std::string>::iterator it = runType.begin();
  BGSubtractedRate *bgSubtr;

  std::vector <std::string> typePrevRun; // holds the runtypes that have already been run so that it doesnt try to run another one.

  while (it!=runType.end()) {
    
    bool skip = false;
    for (auto ch : typePrevRun) {
      if (ch==it->second) skip = true;
    }
        
    if (checkIfBetaRun(it->second) && !skip) {

      typePrevRun.push_back(it->second);

      std::vector <int> fg_runs = makeRunVec(it->second);
      std::vector <int> bg_runs = makeRunVec(getBGrunLabel(it->second));

      bgSubtr = new BGSubtractedRate(fg_runs,bg_runs,analysisChoice,energyBinWidth,fiducialCut,UKdata,Simulation, UNBLIND); //UNBLIND defaults to false

      std::cout << "initialized BGStubtractedRate for run " << it->first << " (BG run " 
		<< getBGrunNumber(it->first) << ") \n";

      bgSubtr->calcBGSubtRates();
	
      if ( it->second=="A2" ) {
	A2[0] = bgSubtr->returnBGSubtRate(0);
	A2err[0] = bgSubtr->returnBGSubtRateError(0);
	A2[1] = bgSubtr->returnBGSubtRate(1);
	A2err[1] = bgSubtr->returnBGSubtRateError(1);

	for (int i=0; i<4; ++i) {
	  A2Quad[i][0] = bgSubtr->returnBGSubtRateQuad(0,i);
	  A2errQuad[i][0] = bgSubtr->returnBGSubtRateErrorQuad(0,i);
	  A2Quad[i][1] = bgSubtr->returnBGSubtRateQuad(1,i);
	  A2errQuad[i][1] = bgSubtr->returnBGSubtRateErrorQuad(1,i);
	}

	for (int i=0; i<6; ++i) {
	  A2Rad[i][0] = bgSubtr->returnBGSubtRateRad(0,i);
	  A2errRad[i][0] = bgSubtr->returnBGSubtRateErrorRad(0,i);
	  A2Rad[i][1] = bgSubtr->returnBGSubtRateRad(1,i);
	  A2errRad[i][1] = bgSubtr->returnBGSubtRateErrorRad(1,i);
	}
	  
	A2len[0] = bgSubtr->returnRunLengths(true);
	A2len[1] = bgSubtr->returnRunLengths(false);
      }
      else if ( it->second=="A5" ) {
	A5[0] = bgSubtr->returnBGSubtRate(0);
	A5err[0] = bgSubtr->returnBGSubtRateError(0);
	A5[1] = bgSubtr->returnBGSubtRate(1);
	A5err[1] = bgSubtr->returnBGSubtRateError(1);

	for (int i=0; i<4; ++i) {
	  A5Quad[i][0] = bgSubtr->returnBGSubtRateQuad(0,i);
	  A5errQuad[i][0] = bgSubtr->returnBGSubtRateErrorQuad(0,i);
	  A5Quad[i][1] = bgSubtr->returnBGSubtRateQuad(1,i);
	  A5errQuad[i][1] = bgSubtr->returnBGSubtRateErrorQuad(1,i);
	}

	for (int i=0; i<6; ++i) {
	  A5Rad[i][0] = bgSubtr->returnBGSubtRateRad(0,i);
	  A5errRad[i][0] = bgSubtr->returnBGSubtRateErrorRad(0,i);
	  A5Rad[i][1] = bgSubtr->returnBGSubtRateRad(1,i);
	  A5errRad[i][1] = bgSubtr->returnBGSubtRateErrorRad(1,i);
	}
	  
	A5len[0] = bgSubtr->returnRunLengths(true);
	A5len[1] = bgSubtr->returnRunLengths(false);
      }
      else if ( it->second=="A7" ) {
	A7[0] = bgSubtr->returnBGSubtRate(0);
	A7err[0] = bgSubtr->returnBGSubtRateError(0);
	A7[1] = bgSubtr->returnBGSubtRate(1);
	A7err[1] = bgSubtr->returnBGSubtRateError(1);

	for (int i=0; i<4; ++i) {
	  A7Quad[i][0] = bgSubtr->returnBGSubtRateQuad(0,i);
	  A7errQuad[i][0] = bgSubtr->returnBGSubtRateErrorQuad(0,i);
	  A7Quad[i][1] = bgSubtr->returnBGSubtRateQuad(1,i);
	  A7errQuad[i][1] = bgSubtr->returnBGSubtRateErrorQuad(1,i);
	}

	for (int i=0; i<6; ++i) {
	  A7Rad[i][0] = bgSubtr->returnBGSubtRateRad(0,i);
	  A7errRad[i][0] = bgSubtr->returnBGSubtRateErrorRad(0,i);
	  A7Rad[i][1] = bgSubtr->returnBGSubtRateRad(1,i);
	  A7errRad[i][1] = bgSubtr->returnBGSubtRateErrorRad(1,i);
	}
	  
	A7len[0] = bgSubtr->returnRunLengths(true);
	A7len[1] = bgSubtr->returnRunLengths(false);
      }
      else if ( it->second=="A10" ) {
	A10[0] = bgSubtr->returnBGSubtRate(0);
	A10err[0] = bgSubtr->returnBGSubtRateError(0);
	A10[1] = bgSubtr->returnBGSubtRate(1);
	A10err[1] = bgSubtr->returnBGSubtRateError(1);

	for (int i=0; i<4; ++i) {
	  A10Quad[i][0] = bgSubtr->returnBGSubtRateQuad(0,i);
	  A10errQuad[i][0] = bgSubtr->returnBGSubtRateErrorQuad(0,i);
	  A10Quad[i][1] = bgSubtr->returnBGSubtRateQuad(1,i);
	  A10errQuad[i][1] = bgSubtr->returnBGSubtRateErrorQuad(1,i);
	}

	for (int i=0; i<6; ++i) {
	  A10Rad[i][0] = bgSubtr->returnBGSubtRateRad(0,i);
	  A10errRad[i][0] = bgSubtr->returnBGSubtRateErrorRad(0,i);
	  A10Rad[i][1] = bgSubtr->returnBGSubtRateRad(1,i);
	  A10errRad[i][1] = bgSubtr->returnBGSubtRateErrorRad(1,i);
	}
	  
	A10len[0] = bgSubtr->returnRunLengths(true);
	A10len[1] = bgSubtr->returnRunLengths(false);
      }
      else if ( it->second=="B2" ) {
	B2[0] = bgSubtr->returnBGSubtRate(0);
	B2err[0] = bgSubtr->returnBGSubtRateError(0);
	B2[1] = bgSubtr->returnBGSubtRate(1);
	B2err[1] = bgSubtr->returnBGSubtRateError(1);
	
	for (int i=0; i<4; ++i) {
	  B2Quad[i][0] = bgSubtr->returnBGSubtRateQuad(0,i);
	  B2errQuad[i][0] = bgSubtr->returnBGSubtRateErrorQuad(0,i);
	  B2Quad[i][1] = bgSubtr->returnBGSubtRateQuad(1,i);
	  B2errQuad[i][1] = bgSubtr->returnBGSubtRateErrorQuad(1,i);
	}

	for (int i=0; i<6; ++i) {
	  B2Rad[i][0] = bgSubtr->returnBGSubtRateRad(0,i);
	  B2errRad[i][0] = bgSubtr->returnBGSubtRateErrorRad(0,i);
	  B2Rad[i][1] = bgSubtr->returnBGSubtRateRad(1,i);
	  B2errRad[i][1] = bgSubtr->returnBGSubtRateErrorRad(1,i);
	}
	  
	B2len[0] = bgSubtr->returnRunLengths(true);
	B2len[1] = bgSubtr->returnRunLengths(false);
      }
      else if ( it->second=="B5" ) {
	B5[0] = bgSubtr->returnBGSubtRate(0);
	B5err[0] = bgSubtr->returnBGSubtRateError(0);
	B5[1] = bgSubtr->returnBGSubtRate(1);
	B5err[1] = bgSubtr->returnBGSubtRateError(1);

	for (int i=0; i<4; ++i) {
	  B5Quad[i][0] = bgSubtr->returnBGSubtRateQuad(0,i);
	  B5errQuad[i][0] = bgSubtr->returnBGSubtRateErrorQuad(0,i);
	  B5Quad[i][1] = bgSubtr->returnBGSubtRateQuad(1,i);
	  B5errQuad[i][1] = bgSubtr->returnBGSubtRateErrorQuad(1,i);
	}

	for (int i=0; i<6; ++i) {
	  B5Rad[i][0] = bgSubtr->returnBGSubtRateRad(0,i);
	  B5errRad[i][0] = bgSubtr->returnBGSubtRateErrorRad(0,i);
	  B5Rad[i][1] = bgSubtr->returnBGSubtRateRad(1,i);
	  B5errRad[i][1] = bgSubtr->returnBGSubtRateErrorRad(1,i);
	}
	  
	B5len[0] = bgSubtr->returnRunLengths(true);
	B5len[1] = bgSubtr->returnRunLengths(false);
      }
      else if ( it->second=="B7" ) {
	B7[0] = bgSubtr->returnBGSubtRate(0);
	B7err[0] = bgSubtr->returnBGSubtRateError(0);
	B7[1] = bgSubtr->returnBGSubtRate(1);
	B7err[1] = bgSubtr->returnBGSubtRateError(1);

	for (int i=0; i<4; ++i) {
	  B7Quad[i][0] = bgSubtr->returnBGSubtRateQuad(0,i);
	  B7errQuad[i][0] = bgSubtr->returnBGSubtRateErrorQuad(0,i);
	  B7Quad[i][1] = bgSubtr->returnBGSubtRateQuad(1,i);
	  B7errQuad[i][1] = bgSubtr->returnBGSubtRateErrorQuad(1,i);
	}

	for (int i=0; i<6; ++i) {
	  B7Rad[i][0] = bgSubtr->returnBGSubtRateRad(0,i);
	  B7errRad[i][0] = bgSubtr->returnBGSubtRateErrorRad(0,i);
	  B7Rad[i][1] = bgSubtr->returnBGSubtRateRad(1,i);
	  B7errRad[i][1] = bgSubtr->returnBGSubtRateErrorRad(1,i);
	}
	  
	B7len[0] = bgSubtr->returnRunLengths(true);
	B7len[1] = bgSubtr->returnRunLengths(false);
      }
      else if ( it->second=="B10" ) {
	B10[0] = bgSubtr->returnBGSubtRate(0);
	B10err[0] = bgSubtr->returnBGSubtRateError(0);
	B10[1] = bgSubtr->returnBGSubtRate(1);
	B10err[1] = bgSubtr->returnBGSubtRateError(1);

	for (int i=0; i<4; ++i) {
	  B10Quad[i][0] = bgSubtr->returnBGSubtRateQuad(0,i);
	  B10errQuad[i][0] = bgSubtr->returnBGSubtRateErrorQuad(0,i);
	  B10Quad[i][1] = bgSubtr->returnBGSubtRateQuad(1,i);
	  B10errQuad[i][1] = bgSubtr->returnBGSubtRateErrorQuad(1,i);
	}

	for (int i=0; i<6; ++i) {
	  B10Rad[i][0] = bgSubtr->returnBGSubtRateRad(0,i);
	  B10errRad[i][0] = bgSubtr->returnBGSubtRateErrorRad(0,i);
	  B10Rad[i][1] = bgSubtr->returnBGSubtRateRad(1,i);
	  B10errRad[i][1] = bgSubtr->returnBGSubtRateErrorRad(1,i);
	}
	  
	B10len[0] = bgSubtr->returnRunLengths(true);
	B10len[1] = bgSubtr->returnRunLengths(false);
      }
	
      else throw "Run misidentified in loadRates"; 
	  
      delete bgSubtr;
    }
    it++;
  }
  
};

bool AsymmetryBase::checkIfBetaRun(std::string type) {
  if (type=="A2" || type=="A5" || type=="A7" || type=="A10" || type=="B2" || type=="B5" || type=="B7" || type=="B10") return true;
  else return false;
};


std::vector < std::vector<double> > AsymmetryBase::returnBGsubtractedRate(std::string runType) {
  if (runType=="A2") return A2;
  else if (runType=="A5") return A5;
  else if (runType=="A7") return A7;
  else if (runType=="A10") return A10;
  else if (runType=="B2") return B2;
  else if (runType=="B5") return B5;
  else if (runType=="B7") return B7;
  else if (runType=="B10") return B10;
  else throw "Bad Run Type given to returnBGsubtractedRate(string runType)";

};

std::vector < std::vector<double> > AsymmetryBase::returnBGsubtractedRateError(std::string runType) {
  if (runType=="A2") return A2err;
  else if (runType=="A5") return A5err;
  else if (runType=="A7") return A7err;
  else if (runType=="A10") return A10err;
  else if (runType=="B2") return B2err;
  else if (runType=="B5") return B5err;
  else if (runType=="B7") return B7err;
  else if (runType=="B10") return B10err;
  else throw "Bad Run Type given to returnBGsubtractedRateError(string runType)";

};

std::vector < std::vector<double> > AsymmetryBase::returnBGsubtractedRateQuad(std::string runType,int quad) {
  if (runType=="A2") return A2Quad[quad];
  else if (runType=="A5") return A5Quad[quad];
  else if (runType=="A7") return A7Quad[quad];
  else if (runType=="A10") return A10Quad[quad];
  else if (runType=="B2") return B2Quad[quad];
  else if (runType=="B5") return B5Quad[quad];
  else if (runType=="B7") return B7Quad[quad];
  else if (runType=="B10") return B10Quad[quad];
  else throw "Bad Run Type given to returnBGsubtractedRate(string runType)";

};

std::vector < std::vector<double> > AsymmetryBase::returnBGsubtractedRateErrorQuad(std::string runType, int quad) {
  if (runType=="A2") return A2errQuad[quad];
  else if (runType=="A5") return A5errQuad[quad];
  else if (runType=="A7") return A7errQuad[quad];
  else if (runType=="A10") return A10errQuad[quad];
  else if (runType=="B2") return B2errQuad[quad];
  else if (runType=="B5") return B5errQuad[quad];
  else if (runType=="B7") return B7errQuad[quad];
  else if (runType=="B10") return B10errQuad[quad];
  else throw "Bad Run Type given to returnBGsubtractedRateError(string runType)";

};

std::vector < std::vector<double> > AsymmetryBase::returnBGsubtractedRateRad(std::string runType,int rad) {
  if (runType=="A2") return A2Rad[rad];
  else if (runType=="A5") return A5Rad[rad];
  else if (runType=="A7") return A7Rad[rad];
  else if (runType=="A10") return A10Rad[rad];
  else if (runType=="B2") return B2Rad[rad];
  else if (runType=="B5") return B5Rad[rad];
  else if (runType=="B7") return B7Rad[rad];
  else if (runType=="B10") return B10Rad[rad];
  else throw "Bad Run Type given to returnBGsubtractedRate(string runType)";

};

std::vector < std::vector<double> > AsymmetryBase::returnBGsubtractedRateErrorRad(std::string runType, int rad) {
  if (runType=="A2") return A2errRad[rad];
  else if (runType=="A5") return A5errRad[rad];
  else if (runType=="A7") return A7errRad[rad];
  else if (runType=="A10") return A10errRad[rad];
  else if (runType=="B2") return B2errRad[rad];
  else if (runType=="B5") return B5errRad[rad];
  else if (runType=="B7") return B7errRad[rad];
  else if (runType=="B10") return B10errRad[rad];
  else throw "Bad Run Type given to returnBGsubtractedRateError(string runType)";

};

///////////////////////////////////////////////////////////////////////////////////////////

OctetAsymmetry::OctetAsymmetry(int oct, std::string anaCh, double enBinWidth, double fidCut, bool ukdata, bool simulation, bool unblind) : AsymmetryBase(oct,anaCh,enBinWidth,fidCut,ukdata,simulation,unblind), totalAsymmetry(0.), totalAsymmetryError(0.), boolSuperSum(false), boolAsymmetry(false) {
  if (isFullOctet()) {
    //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
    asymmetry.resize(numEnBins,0.);
    asymmetryError.resize(numEnBins,0.);

    asymmetryQuad.resize(4,std::vector<double> (numEnBins,0.) );
    asymmetryErrorQuad.resize(4,std::vector<double> (numEnBins,0.) );
    
    asymmetryRad.resize(6,std::vector<double> (numEnBins,0.) );
    asymmetryErrorRad.resize(6,std::vector<double> (numEnBins,0.) );

    superSum.resize(numEnBins,0.);
    superSumError.resize(numEnBins,0.);
    loadRates(); // load the rates in the rate vectors for each run
    std::cout <<"//////////////////////////////////////////////////////////////////\n"
	      <<"OctetAsymmetry for octet " << octet << std::endl;
  }
  else { 
    boolAsymmetry=boolSuperSum=true; //Writing "bad octet" to all files so that it isn't used in the future
    writeAsymToFile();
    writeSuperSumToFile();
  
    throw "OCTET NOT A COMPLETE OCTET";
  }
}; 

void OctetAsymmetry::calcAsymmetryBinByBin() {
  
  calcBinByBinSuperRatio(A2,A5,A7,A10,B2,B5,B7,B10,
			 A2err,A5err,A7err,A10err,B2err,B5err,B7err,B10err,
			 asymmetry,asymmetryError);
  for (int i=0; i<4; ++i) {
    calcBinByBinSuperRatio(A2Quad[i],A5Quad[i],A7Quad[i],A10Quad[i],
			   B2Quad[i],B5Quad[i],B7Quad[i],B10Quad[i],
			   A2errQuad[i],A5errQuad[i],A7errQuad[i],A10errQuad[i],
			   B2errQuad[i],B5errQuad[i],B7errQuad[i],B10errQuad[i],
			   asymmetryQuad[i],asymmetryErrorQuad[i]);
  }
  for (int i=0; i<6; ++i) {
    calcBinByBinSuperRatio(A2Rad[i],A5Rad[i],A7Rad[i],A10Rad[i],
			   B2Rad[i],B5Rad[i],B7Rad[i],B10Rad[i],
			   A2errRad[i],A5errRad[i],A7errRad[i],A10errRad[i],
			   B2errRad[i],B5errRad[i],B7errRad[i],B10errRad[i],
			   asymmetryRad[i],asymmetryErrorRad[i]);
  }
  boolAsymmetry = true;
};



void OctetAsymmetry::calcNCSUSumAsymmetryBinByBin() {
  
  //These are the rates and their errors that go into the Super-ratio
  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};


  for (unsigned int bin=0; bin<asymmetry.size(); bin++) {

    double R = 0.;
    double deltaR = 0.;

    for (unsigned int side=0; side<2; side++) {

      double numer = 0.;
      double denom = 0.;
      double numer_Err = 0.;

      if ( A2err[side][bin] > 0. && A2[side][bin] > 0. ) {
	numer += A2[side][bin]*A2len[0][side];
	numer_Err += power(A2err[side][bin]*A2len[0][side],2);
	denom += A2len[0][side];
      }
      if ( A10err[side][bin] > 0. && A10[side][bin] > 0. ) {
	numer += A10[side][bin]*A10len[0][side];
	numer_Err += power(A10err[side][bin]*A10len[0][side],2);
	denom += A10len[0][side];
      }
      if ( B5err[side][bin] > 0. && B5[side][bin] > 0. ) {
	numer += B5[side][bin]*B5len[0][side];
	numer_Err += power(B5err[side][bin]*B5len[0][side],2);
	denom += B5len[0][side];
      }
      if ( B7err[side][bin] > 0. && B7[side][bin] > 0. ) {
	numer += B7[side][bin]*B7len[0][side];
	numer_Err += power(B7err[side][bin]*B7len[0][side],2);
	denom += B7len[0][side];
      }

      sfOFF[side] = denom > 0. ? numer/denom : 0.;
      sfOFF_err[side] = denom > 0. ? sqrt(numer_Err)/denom : 0.;

      // Flipper On
      numer = denom = numer_Err = 0.;

      if ( B2err[side][bin] > 0. && B2[side][bin] > 0. ) {
	numer += B2[side][bin]*B2len[0][side];
	numer_Err += power(B2err[side][bin]*B2len[0][side],2);
	denom += B2len[0][side];
      }
      if ( B10err[side][bin] > 0. && B10[side][bin] > 0. ) {
	numer += B10[side][bin]*B10len[0][side];
	numer_Err += power(B10err[side][bin]*B10len[0][side],2);
	denom += B10len[0][side];
      }
      if ( A5err[side][bin] > 0. && A5[side][bin] > 0. ) {
	numer += A5[side][bin]*A5len[0][side];
	numer_Err += power(A5err[side][bin]*A5len[0][side],2);
	denom += A5len[0][side];
      }
      if ( A7err[side][bin] > 0. && A7[side][bin] > 0. ) {
	numer += A7[side][bin]*A7len[0][side];
	numer_Err += power(A7err[side][bin]*A7len[0][side],2);
	denom += A7len[0][side];
      }

      sfON[side] = denom > 0. ? numer/denom : 0.;
      sfON_err[side] = denom > 0. ? sqrt(numer_Err)/denom : 0.;

      
    }

    R = sfOFF[0]*sfON[1]/(sfON[0]*sfOFF[1]);

    if ( R!=0. && sfON[0]>0. && sfOFF[1]>0. ) { 
     
      deltaR = sqrt( power( sfOFF_err[0]*sfON[1]/(sfOFF[1]*sfON[0]), 2 ) + 
		     power( sfON_err[1]*sfOFF[0]/(sfOFF[1]*sfON[0]), 2 ) +
		     power( sfOFF_err[1]*sfON[1]*sfOFF[0]/(sfOFF[1]*sfOFF[1]*sfON[0]), 2 ) +
		     power( sfON_err[0]*sfON[1]*sfOFF[0]/(sfOFF[1]*sfON[0]*sfON[0]), 2 ) );
      //deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      asymmetry[bin] = (1.-sqrt(std::abs(R)))/(1+sqrt(std::abs(R)));
      asymmetryError[bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      asymmetry[bin] = 0.;
      asymmetryError[bin] = 1.;
    }

  } 
 
  boolAsymmetry = true;
};

void OctetAsymmetry::calcTotalAsymmetry(double enWinLow, double enWinHigh) {

  unsigned int binLow = (unsigned int)(enWinLow/energyBinWidth);
  unsigned int binHigh = (unsigned int)(enWinHigh/energyBinWidth)-1;

  double sfON[2]={0.}, sfOFF[2]={0.};
  double sfON_err[2]={0.}, sfOFF_err[2]={0.};
  double sumA2[2]={0.}, sumA5[2]={0.}, sumA7[2]={0.}, sumA10[2]={0.}, sumB2[2]={0.}, sumB5[2]={0.}, sumB7[2]={0.}, sumB10[2]={0.};
  double sumA2_err[2]={0.}, sumA5_err[2]={0.}, sumA7_err[2]={0.}, sumA10_err[2]={0.}, sumB2_err[2]={0.}, sumB5_err[2]={0.}, sumB7_err[2]={0.}, sumB10_err[2]={0.};

  double R = 0.;
  double deltaR = 0.;

  for (unsigned int side=0; side<2; side++) {
    for (unsigned int bin=binLow; bin<=binHigh; bin++) {
      sumA2[side]+=A2[side][bin];
      sumA2_err[side]+=power(A2err[side][bin],2);
      sumA5[side]+=A5[side][bin];
      sumA5_err[side]+=power(A5err[side][bin],2);
      sumA7[side]+=A7[side][bin];
      sumA7_err[side]+=power(A7err[side][bin],2);
      sumA10[side]+=A10[side][bin];
      sumA10_err[side]+=power(A10err[side][bin],2);
      sumB2[side]+=B2[side][bin];
      sumB2_err[side]+=power(B2err[side][bin],2);
      sumB5[side]+=B5[side][bin];
      sumB5_err[side]+=power(B5err[side][bin],2);
      sumB7[side]+=B7[side][bin];
      sumB7_err[side]+=power(B7err[side][bin],2);
      sumB10[side]+=B10[side][bin];
      sumB10_err[side]+=power(B10err[side][bin],2);

      /*if (side==1) {
	std::cout << A10[side][bin] << " " << A10err[side][bin] << std::endl;
	}*/ 
    }
    sumA2_err[side] = sqrt(sumA2_err[side]);
    sumA5_err[side] = sqrt(sumA5_err[side]);
    sumA7_err[side] = sqrt(sumA7_err[side]);
    sumA10_err[side] = sqrt(sumA10_err[side]);
    sumB2_err[side] = sqrt(sumB2_err[side]);
    sumB5_err[side] = sqrt(sumB5_err[side]);
    sumB7_err[side] = sqrt(sumB7_err[side]);
    sumB10_err[side] = sqrt(sumB10_err[side]);
    
  }
    		   
  for (unsigned int side=0; side<2; side++) {
  
   double numer = 0.;
   double denom = 0.;
   double numer_Err = 0.;
   
   if ( sumA2_err[side] > 0. && sumA2[side] > 0. ) {
     numer += sumA2[side]*A2len[0][side];
     numer_Err += power(sumA2_err[side]*A2len[0][side],2);
     denom += A2len[0][side];
   }
   if ( sumA10_err[side] > 0. && sumA10[side] > 0. ) {
     numer += sumA10[side]*A10len[0][side];
     numer_Err += power(sumA10_err[side]*A10len[0][side],2);
     denom += A10len[0][side];
   }
   if ( sumB5_err[side] > 0. && sumB5[side] > 0. ) {
     numer += sumB5[side]*B5len[0][side];
     numer_Err += power(sumB5_err[side]*B5len[0][side],2);
     denom += B5len[0][side];
   }
   if ( sumB7_err[side] > 0. && sumB7[side] > 0. ) {
     numer += sumB7[side]*B7len[0][side];
     numer_Err += power(sumB7_err[side]*B7len[0][side],2);
     denom += B7len[0][side];
   }
   
   sfOFF[side] = denom > 0. ? numer/denom : 0.;
   sfOFF_err[side] = denom > 0. ? sqrt(numer_Err)/denom : 0.;
   
   // Flipper On
   numer = denom = numer_Err = 0.;
   
   if ( sumB2_err[side] > 0. && sumB2[side] > 0. ) {
     numer += sumB2[side]*B2len[0][side];
     numer_Err += power(sumB2_err[side]*B2len[0][side],2);
     denom += B2len[0][side];
   }
   if ( sumB10_err[side] > 0. && sumB10[side] > 0. ) {
     numer += sumB10[side]*B10len[0][side];
     numer_Err += power(sumB10_err[side]*B10len[0][side],2);
     denom += B10len[0][side];
   }
   if ( sumA5_err[side] > 0. && sumA5[side] > 0. ) {
     numer += sumA5[side]*A5len[0][side];
     numer_Err += power(sumA5_err[side]*A5len[0][side],2);
     denom += A5len[0][side];
   }
   if ( sumA7_err[side] > 0. && sumA7[side] > 0. ) {
     numer += sumA7[side]*A7len[0][side];
     numer_Err += power(sumA7_err[side]*A7len[0][side],2);
     denom += A7len[0][side];
   }
   
   sfON[side] = denom > 0. ? numer/denom : 0.;
   sfON_err[side] = denom > 0. ? sqrt(numer_Err)/denom : 0.;
   
   
  }

  R = sfOFF[0]*sfON[1]/(sfON[0]*sfOFF[1]);

  if ( R!=0. && sfON[0]>0. && sfOFF[1]>0. ) { 
     
    deltaR = sqrt( power( sfOFF_err[0]*sfON[1]/(sfOFF[1]*sfON[0]), 2 ) + 
		   power( sfON_err[1]*sfOFF[0]/(sfOFF[1]*sfON[0]), 2 ) +
		   power( sfOFF_err[1]*sfON[1]*sfOFF[0]/(sfOFF[1]*sfOFF[1]*sfON[0]), 2 ) +
		   power( sfON_err[0]*sfON[1]*sfOFF[0]/(sfOFF[1]*sfON[0]*sfON[0]), 2 ) );
    //deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
    totalAsymmetry = (1.-sqrt(std::abs(R)))/(1+sqrt(std::abs(R)));
    totalAsymmetryError = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
  }
  else {
    totalAsymmetry = 0.;
    totalAsymmetryError = 1.;
  }  
  
  std::ofstream ofile(TString::Format("OctetAveRateFiles/UKOctetAvg_%i.txt",octet).Data());
  
  ofile << "RunType\tEast_bgSubt\tWest_bgSubt\tEast_time\tWest_time\n";
  ofile << "A2\t" << sumA2[0] << "\t" << sumA2[1] << "\t" << A2len[0][0] << "\t" << A2len[0][1] << std::endl;
  ofile << "A5\t" << sumA5[0] << "\t" << sumA5[1] << "\t" << A5len[0][0] << "\t" << A5len[0][1] << std::endl;
  ofile << "A7\t" << sumA7[0] << "\t" << sumA7[1] << "\t" << A7len[0][0] << "\t" << A7len[0][1] << std::endl;
  ofile << "A10\t" << sumA10[0] << "\t" << sumA10[1] << "\t" << A10len[0][0] << "\t" << A10len[0][1] << std::endl;
  ofile << "B2\t" << sumB2[0] << "\t" << sumB2[1] << "\t" << B2len[0][0] << "\t" << B2len[0][1] << std::endl;
  ofile << "B5\t" << sumB5[0] << "\t" << sumB5[1] << "\t" << B5len[0][0] << "\t" << B5len[0][1] << std::endl;
  ofile << "B7\t" << sumB7[0] << "\t" << sumB7[1] << "\t" << B7len[0][0] << "\t" << B7len[0][1] << std::endl;
  ofile << "B10\t" << sumB10[0] << "\t" << sumB10[1] << "\t" << B10len[0][0] << "\t" << B10len[0][1] << "\n\n";
  ofile << "sfOFF_E\t" << sfOFF[0] <<"\n";
  ofile << "sfON_W\t" << sfON[1] <<"\n";
  ofile << "sfOFF_W\t" << sfOFF[1] <<"\n";
  ofile << "sfON_E\t" << sfON[0] <<"\n\n";
  ofile << "R\t" << R << "\n\n";
  ofile << "A_ave\t" << totalAsymmetry << "\n";

  ofile.close();
  //std::cout << std::endl << R << " " << deltaR << "\n";
  
};


void OctetAsymmetry::calcSuperSum() {

  //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
  superSum.resize(numEnBins,0.);
  superSumError.resize(numEnBins,0.);

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  for (unsigned int bin=0; bin<superSum.size(); bin++) {
 
    for (unsigned int side=0; side<2; side++) {
      double weightsum=0.;

      sfOFF[side] = ( ( A2err[side][bin]>0. ? power(1./A2err[side][bin],2)*A2[side][bin]  :0. ) + 
		      ( A10err[side][bin]>0.? power(1./A10err[side][bin],2)*A10[side][bin]:0. ) + 
		      ( B5err[side][bin]>0. ? power(1./B5err[side][bin],2)*B5[side][bin]  :0. ) + 
		      ( B7err[side][bin]>0. ? power(1./B7err[side][bin],2)*B7[side][bin]  :0. ) );
      weightsum = ( ( A2err[side][bin]>0. ? power(1./A2err[side][bin],2):0. ) +
		    ( A10err[side][bin]>0.? power(1./A10err[side][bin],2):0.) + 
		    ( B5err[side][bin]>0. ? power(1./B5err[side][bin],2):0. ) + 
		    ( B7err[side][bin]>0. ? power(1./B7err[side][bin],2):0. ) );

      sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
      sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;

      weightsum=0.;

      sfON[side] = ( ( B2err[side][bin] > 0. ? B2[side][bin]/power(B2err[side][bin],2)   : 0. ) +
		      ( B10err[side][bin]> 0. ? B10[side][bin]/power(B10err[side][bin],2) : 0. ) +
		      (	A5err[side][bin] > 0. ? A5[side][bin]/power(A5err[side][bin],2)   : 0. ) +
		      ( A7err[side][bin] > 0. ? A7[side][bin]/power(A7err[side][bin],2)   : 0. ) );

      weightsum = ( ( B2err[side][bin] > 0. ? 1./power(B2err[side][bin],2)  : 0. ) + 
		    ( B10err[side][bin]> 0. ? 1./power(B10err[side][bin],2) : 0. ) +
		    ( A5err[side][bin] > 0. ? 1./power(A5err[side][bin],2)  : 0. ) +
		    ( A7err[side][bin] > 0. ? 1./power(A7err[side][bin],2)  : 0. ) );
      
      sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
      sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;

    }

    //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
    double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
    double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
    double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
      0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
    double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
      0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));

    superSum[bin] = R1 + R2;
    //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
    superSumError[bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
    
  } 
  boolSuperSum = true;
};

void OctetAsymmetry::calcSuperSumNCSUstyle() {

  //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
  superSum.resize(numEnBins,0.);
  superSumError.resize(numEnBins,0.);

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  for (unsigned int bin=0; bin<superSum.size(); bin++) {
 
    for (unsigned int side=0; side<2; side++) {

      double numer = 0.;
      double denom = 0.;
      double numer_Err = 0.;

      if ( A2err[side][bin] > 0. ) {
	numer += A2[side][bin]*A2len[0][side];
	numer_Err += power(A2err[side][bin]*A2len[0][side],2);
	denom += A2len[0][side];
      }
      if ( A10err[side][bin] > 0. ) {
	numer += A10[side][bin]*A10len[0][side];
	numer_Err += power(A10err[side][bin]*A10len[0][side],2);
	denom += A10len[0][side];
      }
      if ( B5err[side][bin] > 0. ) {
	numer += B5[side][bin]*B5len[0][side];
	numer_Err += power(B5err[side][bin]*B5len[0][side],2);
	denom += B5len[0][side];
      }
      if ( B7err[side][bin] > 0. ) {
	numer += B7[side][bin]*B7len[0][side];
	numer_Err += power(B7err[side][bin]*B7len[0][side],2);
	denom += B7len[0][side];
      }

      sfOFF[side] = denom > 0. ? numer/denom : 0.;
      sfOFF_err[side] = denom > 0. ? sqrt(numer_Err)/denom : 1.;

      // Flipper On
      numer = denom = numer_Err = 0.;

      if ( B2err[side][bin] > 0. ) {
	numer += B2[side][bin]*B2len[0][side];
	numer_Err += power(B2err[side][bin]*B2len[0][side],2);
	denom += B2len[0][side];
      }
      if ( B10err[side][bin] > 0. ) {
	numer += B10[side][bin]*B10len[0][side];
	numer_Err += power(B10err[side][bin]*B10len[0][side],2);
	denom += B10len[0][side];
      }
      if ( A5err[side][bin] > 0. ) {
	numer += A5[side][bin]*A5len[0][side];
	numer_Err += power(A5err[side][bin]*A5len[0][side],2);
	denom += A5len[0][side];
      }
      if ( A7err[side][bin] > 0. ) {
	numer += A7[side][bin]*A7len[0][side];
	numer_Err += power(A7err[side][bin]*A7len[0][side],2);
	denom += A7len[0][side];
      }

      sfON[side] = denom > 0. ? numer/denom : 0.;
      sfON_err[side] = denom > 0. ? sqrt(numer_Err)/denom : 1.;


    }

    //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
    double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
    double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
    double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
      0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
    double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
      0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));

    superSum[bin] = R1 + R2;
    //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
    superSumError[bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
    
  } 
  boolSuperSum = true;
};

void OctetAsymmetry::writeAsymToFile() {
  if (!isAsymmetry()) {
    std::cout << "No Asymmetry has been calculated. Not writing to file!\n\n";
    return;
  }

  std::string outpath;  
  
  if (Simulation) outpath = std::string(getenv("SIM_ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+".dat";
  else if (UKdata) outpath = std::string(getenv("ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+".dat";
  else outpath = std::string(getenv("MPM_ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+".dat";
  std::ofstream outfile(outpath.c_str());
  
  if (isFullOctet()) {
    for (unsigned int i=0; i<asymmetry.size(); i++) {
      outfile << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
    }
  }
  else {
    outfile << "BAD OCTET";
  }
  outfile.close();

  // Now for quadrants
  for (int j=0; j<4; ++j) {
    
    if (Simulation) outpath = std::string(getenv("SIM_ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Quadrant"+itos(j)+".dat";
    else if (UKdata) outpath = std::string(getenv("ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Quadrant"+itos(j)+".dat";
    else outpath = std::string(getenv("MPM_ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Quadrant"+itos(j)+".dat";
    outfile.open(outpath.c_str());
    
    if (isFullOctet()) {
      for (unsigned int i=0; i<asymmetryQuad[j].size(); i++) {
	outfile << binLowerEdge[i] << " " << asymmetryQuad[j][i] << " " << asymmetryErrorQuad[j][i] << std::endl;
	//std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
      }
    }
    else {
      outfile << "BAD OCTET";
    }
    outfile.close();
  }

  for (int j=0; j<6; ++j) {
    
    if (Simulation) outpath = std::string(getenv("SIM_ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_RadialRing"+itos(j)+".dat";
    else if (UKdata) outpath = std::string(getenv("ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_RadialRing"+itos(j)+".dat";
    else outpath = std::string(getenv("MPM_ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_RadialRing"+itos(j)+".dat";
    outfile.open(outpath.c_str());
    
    if (isFullOctet()) {
      for (unsigned int i=0; i<asymmetryRad[j].size(); i++) {
	outfile << binLowerEdge[i] << " " << asymmetryRad[j][i] << " " << asymmetryErrorRad[j][i] << std::endl;
	//std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
      }
    }
    else {
      outfile << "BAD OCTET";
    }
    outfile.close();
  }
  
  std::cout << "Wrote Asymmetry to file for " << analysisChoice << " in " << outpath << "\n";
};

void OctetAsymmetry::writeSuperSumToFile() {
  if (!isSuperSum()) {
    std::cout << "No Super-Sum has been calculated. Not writing to file!\n\n";
    return;
  }

  std::string outpath;

  if (Simulation) outpath = std::string(getenv("SIM_ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+analysisChoice+".dat";
  else if (UKdata) outpath = std::string(getenv("ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+analysisChoice+".dat";
  else outpath = std::string(getenv("MPM_ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/superSum_Octet"+ (UNBLIND?"UNBLINDED_":"")+"" + itos(octet)+"_AnaCh"+analysisChoice+".dat";

  
  std::ofstream outfile(outpath.c_str());
  
  if (isFullOctet()) {
    for (unsigned int i=0; i<superSum.size(); i++) {
      outfile << binLowerEdge[i] << " " << superSum[i] << " " << superSumError[i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << superSum[i] << " " << superSumError[i] << std::endl;
    }
  }
  else {
    outfile << "BAD OCTET";
  }

  outfile.close(); 
};

    

///////////////////////////////////////////////////////////////////////////////////////////

QuartetAsymmetry::QuartetAsymmetry(int oct, std::string anaCh, double enBinWidth, double fidCut, bool ukdata, bool simulation, bool unblind) : AsymmetryBase(oct,anaCh,enBinWidth,fidCut,ukdata,simulation,unblind), totalAsymmetryA(0.), totalAsymmetryErrorA(0.), totalAsymmetryB(0.), totalAsymmetryErrorB(0.), boolSuperSum(false), boolAsymmetry(false) {

  isGoodQuartet.push_back(isFullQuartet(0));
  isGoodQuartet.push_back(isFullQuartet(1));

  if (isGoodQuartet[0] || isGoodQuartet[1]) {
    //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
    asymmetry.resize(2, std::vector<double> (numEnBins,0.));
    asymmetryError.resize(2, std::vector<double> (numEnBins,0.));
    superSum.resize(2, std::vector<double> (numEnBins,0.));
    superSumError.resize(2, std::vector<double> (numEnBins,0.));
    loadRates(); // load the rates in the rate vectors for each run
    std::cout <<"//////////////////////////////////////////////////////////////////\n"
	      <<"QuartetAsymmetry for octet " << octet << std::endl
	      <<"\nQuartet Status (0=bad, 1=good): First: " << isGoodQuartet[0] << " " 
	      <<"Second: " << isGoodQuartet[1] << std::endl;
  }
  else {
    boolAsymmetry=boolSuperSum=true; //Writing "bad quartet" to all files so that it isn't used in the future
    writeAsymToFile();
    writeSuperSumToFile();
    
    throw "NO QUARTETS IN THIS OCTET";
  }
};

void QuartetAsymmetry::calcAsymmetryBinByBin() {

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  // A type runs
  if (isGoodQuartet[0]) {
    for (unsigned int bin=0; bin<asymmetry[0].size(); bin++) {
      double R = 0.;
      double deltaR = 0.;
      
      for (unsigned int side=0; side<2; side++) {
	double weightsum=0.;

	//Start with flipper off... Only use rate if it has an error to avoid dividing by 0
	sfOFF[side] = ( ( A2err[side][bin] > 0. ? A2[side][bin]/power(A2err[side][bin],2)   : 0. ) +
			( A10err[side][bin]> 0. ? A10[side][bin]/power(A10err[side][bin],2) : 0. ) );
	
	weightsum = ( ( A2err[side][bin] > 0. ? 1./power(A2err[side][bin],2)  : 0. ) + 
		      ( A10err[side][bin]> 0. ? 1./power(A10err[side][bin],2) : 0. ) ) ;
		      
	sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
		      
		      
	weightsum=0.;

	// Now for flipper on
	sfON[side] = (  ( A5err[side][bin] > 0. ? A5[side][bin]/power(A5err[side][bin],2)   : 0. ) +
			( A7err[side][bin] > 0. ? A7[side][bin]/power(A7err[side][bin],2)   : 0. ) );
	
	weightsum = ( ( A5err[side][bin] > 0. ? 1./power(A5err[side][bin],2)  : 0. ) +
		      ( A7err[side][bin] > 0. ? 1./power(A7err[side][bin],2)  : 0. ) );
	
	sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	
	
      }
      
      //Calculate super-ratio, R. Note that any of the rates that go into this could be negative 
      R = (sfOFF[1]*sfON[0]) / (sfON[1]*sfOFF[0]) ;
      
      //Only use rates which are real and not zero to avoid nan (couldn't get TMath::IsNaN() to work properly :( ...)
      if ( R!=0. && sfON[1]!=0. && sfOFF[0]!=0.) { 
	
	deltaR = sqrt( power( sfOFF_err[1]*sfON[0]/(sfOFF[0]*sfON[1]), 2 ) + 
		       power( sfON_err[0]*sfOFF[1]/(sfOFF[0]*sfON[1]), 2 ) +
		       power( sfOFF_err[0]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfOFF[0]*sfON[1]), 2 ) +
		       power( sfON_err[1]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfON[1]*sfON[1]), 2 ) );
	//deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	asymmetry[0][bin] = (1.-sqrt(std::abs(R)))/(1+sqrt(std::abs(R)));
	asymmetryError[0][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	asymmetry[0][bin] = 0.;
	asymmetryError[0][bin] = 1.;
      } 
    }
  }
	
  
	
  // B type runs
  if (isGoodQuartet[1]) {
    for (unsigned int bin=0; bin<asymmetry[1].size(); bin++) {
      double R = 0.;
      double deltaR = 0.;
      
      for (unsigned int side=0; side<2; side++) {
	double weightsum=0.;
	
	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;


	//Start with flipper off... Only use rate if it has an error to avoid dividing by 0
	sfON[side] = ( ( B2err[side][bin] > 0. ? B2[side][bin]/power(B2err[side][bin],2)   : 0. ) +
			( B10err[side][bin]> 0. ? B10[side][bin]/power(B10err[side][bin],2) : 0. ) );
	
	weightsum = ( ( B2err[side][bin] > 0. ? 1./power(B2err[side][bin],2)  : 0. ) + 
		      ( B10err[side][bin]> 0. ? 1./power(B10err[side][bin],2) : 0. ) ) ;
		      
	sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
		      
		      
	weightsum=0.;

	// Now for flipper on
	sfOFF[side] = (  ( B5err[side][bin] > 0. ? B5[side][bin]/power(B5err[side][bin],2)   : 0. ) +
			 ( B7err[side][bin] > 0. ? B7[side][bin]/power(B7err[side][bin],2)   : 0. ) );
	
	weightsum = ( ( B5err[side][bin] > 0. ? 1./power(B5err[side][bin],2)  : 0. ) +
		      ( B7err[side][bin] > 0. ? 1./power(B7err[side][bin],2)  : 0. ) );
	
	sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	
	
      }
      
      //Calculate super-ratio, R. Note that any of the rates that go into this could be negative 
      R = (sfOFF[1]*sfON[0]) / (sfON[1]*sfOFF[0]) ;
      
      //Only use rates which are real and not zero to avoid nan (couldn't get TMath::IsNaN() to work properly :( ...)
      if ( R!=0. && sfON[1]!=0. && sfOFF[0]!=0.) { 
	
	deltaR = sqrt( power( sfOFF_err[1]*sfON[0]/(sfOFF[0]*sfON[1]), 2 ) + 
		       power( sfON_err[0]*sfOFF[1]/(sfOFF[0]*sfON[1]), 2 ) +
		       power( sfOFF_err[0]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfOFF[0]*sfON[1]), 2 ) +
		       power( sfON_err[1]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfON[1]*sfON[1]), 2 ) );
	//deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	asymmetry[1][bin] = (1.-sqrt(std::abs(R)))/(1+sqrt(std::abs(R)));
	asymmetryError[1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	asymmetry[1][bin] = 0.;
	asymmetryError[1][bin] = 1.;
      } 
    }
  }
  
  boolAsymmetry = true;
};

void QuartetAsymmetry::calcTotalAsymmetry(double enWinLow, double enWinHigh) {

  unsigned int binLow = (unsigned int)(enWinLow/energyBinWidth);
  unsigned int binHigh = (unsigned int)(enWinHigh/energyBinWidth)-1;

  double sfON[2]={0.}, sfOFF[2]={0.};
  double sfON_err[2]={0.}, sfOFF_err[2]={0.};
  double sumA2[2]={0.}, sumA5[2]={0.}, sumA7[2]={0.}, sumA10[2]={0.}, sumB2[2]={0.}, sumB5[2]={0.}, sumB7[2]={0.}, sumB10[2]={0.};
  double sumA2_err[2]={0.}, sumA5_err[2]={0.}, sumA7_err[2]={0.}, sumA10_err[2]={0.}, sumB2_err[2]={0.}, sumB5_err[2]={0.}, sumB7_err[2]={0.}, sumB10_err[2]={0.};

  for (unsigned int side=0; side<2; side++) {
    for (unsigned int bin=binLow; bin<=binHigh; bin++) {
      sumA2[side]+=A2[side][bin];
      sumA2_err[side]+=power(A2err[side][bin],2);
      sumA5[side]+=A5[side][bin];
      sumA5_err[side]+=power(A5err[side][bin],2);
      sumA7[side]+=A7[side][bin];
      sumA7_err[side]+=power(A7err[side][bin],2);
      sumA10[side]+=A10[side][bin];
      sumA10_err[side]+=power(A10err[side][bin],2);
      sumB2[side]+=B2[side][bin];
      sumB2_err[side]+=power(B2err[side][bin],2);
      sumB5[side]+=B5[side][bin];
      sumB5_err[side]+=power(B5err[side][bin],2);
      sumB7[side]+=B7[side][bin];
      sumB7_err[side]+=power(B7err[side][bin],2);
      sumB10[side]+=B10[side][bin];
      sumB10_err[side]+=power(B10err[side][bin],2);

    }
    sumA2_err[side] = sqrt(sumA2_err[side]);
    sumA5_err[side] = sqrt(sumA5_err[side]);
    sumA7_err[side] = sqrt(sumA7_err[side]);
    sumA10_err[side] = sqrt(sumA10_err[side]);
    sumB2_err[side] = sqrt(sumB2_err[side]);
    sumB5_err[side] = sqrt(sumB5_err[side]);
    sumB7_err[side] = sqrt(sumB7_err[side]);
    sumB10_err[side] = sqrt(sumB10_err[side]);
    
  }
    
  double R = 0.;
  double deltaR = 0.; 

  if (isGoodQuartet[0]) {
    for (unsigned int side=0; side<2; side++) {
      
      double weightsum=0.;
      sfOFF[side] = (power(1./sumA2_err[side],2)*sumA2[side]+power(1./sumA10_err[side],2)*sumA10[side]);
      weightsum = power(1./sumA2_err[side],2)+power(1./sumA10_err[side],2);
      
      sfOFF[side] = sfOFF[side]/weightsum;
      sfOFF_err[side] = 1./sqrt(weightsum);
      
      weightsum=0.;
      sfON[side] = (power(1./sumA5_err[side],2)*sumA5[side]+power(1./sumA7_err[side],2)*sumA7[side]);
      weightsum = power(1./sumA5_err[side],2)+power(1./sumA7_err[side],2);
      
      sfON[side] = sfON[side]/weightsum;
      sfON_err[side] = 1./sqrt(weightsum);
      //std::cout << sfOFF[side] << " " << sfON[side] << std::endl;
    }
    
    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      totalAsymmetryA = (1.-sqrt(R))/(1+sqrt(R));
      totalAsymmetryErrorA = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      totalAsymmetryA = 0.;
      totalAsymmetryErrorA= 0.;
    }
    
    //std::cout << std::endl << R << " " << deltaR << "\n";
  }

  R=0.;
  deltaR=0.;

  if (isGoodQuartet[1]) {
    for (unsigned int side=0; side<2; side++) {

      sfOFF[side]=0.;
      sfON[side]=0.;
      sfOFF_err[side]=0.;
      sfON_err[side]=0.;
      
      double weightsum=0.;
      sfOFF[side] = (power(1./sumB5_err[side],2)*sumB5[side]+power(1./sumB7_err[side],2)*sumB7[side]);
      weightsum = power(1./sumB5_err[side],2)+power(1./sumB7_err[side],2);
      
      sfOFF[side] = sfOFF[side]/weightsum;
      sfOFF_err[side] = 1./sqrt(weightsum);
      
      weightsum=0.;
      sfON[side] = (power(1./sumB2_err[side],2)*sumB2[side]+power(1./sumB10_err[side],2)*sumB10[side]);
      weightsum = power(1./sumB2_err[side],2)+power(1./sumB10_err[side],2);
      
      sfON[side] = sfON[side]/weightsum;
      sfON_err[side] = 1./sqrt(weightsum);
      //std::cout << sfOFF[side] << " " << sfON[side] << std::endl;
    }
    
    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      totalAsymmetryB= (1.-sqrt(R))/(1+sqrt(R));
      totalAsymmetryErrorB = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      totalAsymmetryB = 0.;
      totalAsymmetryErrorB= 0.;
    }
    //std::cout << std::endl << R << " " << deltaR << "\n";
  }

  
};


void QuartetAsymmetry::calcSuperSum() {

  //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
  superSum.resize(2, std::vector<double> (numEnBins,0.));
  superSumError.resize(2, std::vector<double> (numEnBins,0.));

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  if (isGoodQuartet[0]) {
    for (unsigned int bin=0; bin<superSum[0].size(); bin++) {
      
      for (unsigned int side=0; side<2; side++) {
	double weightsum=0.;

	// AFP Off
	sfOFF[side] = (A2err[side][bin]>0.?power(1./A2err[side][bin],2)*A2[side][bin]:0.) + (A10err[side][bin]>0.?power(1./A10err[side][bin],2)*A10[side][bin]:0.);
	weightsum = (A2err[side][bin]>0.?power(1./A2err[side][bin],2):0.) + (A10err[side][bin]>0.?power(1./A10err[side][bin],2):0.);
	
	sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = (A5err[side][bin]>0.?power(1./A5err[side][bin],2)*A5[side][bin]:0.) + (A7err[side][bin]>0.?power(1./A7err[side][bin],2)*A7[side][bin]:0.);
	weightsum = (A5err[side][bin]>0.?power(1./A5err[side][bin],2):0.) + (A7err[side][bin]>0.?power(1./A7err[side][bin],2):0.);
	
	sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }
      
      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
	0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      superSum[0][bin] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      superSumError[0][bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
    } 
  }

  if (isGoodQuartet[1]) {
    for (unsigned int bin=0; bin<superSum[1].size(); bin++) {
      
      for (unsigned int side=0; side<2; side++) {

	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;
	double weightsum=0.;

	// AFP Off
	sfOFF[side] = (B5err[side][bin]>0.?power(1./B5err[side][bin],2)*B5[side][bin]:0.) + (B7err[side][bin]>0.?power(1./B7err[side][bin],2)*B7[side][bin]:0.);
	weightsum = (B5err[side][bin]>0.?power(1./B5err[side][bin],2):0.) + (B7err[side][bin]>0.?power(1./B7err[side][bin],2):0.);
	
	sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = (B2err[side][bin]>0.?power(1./B2err[side][bin],2)*B2[side][bin]:0.) + (B10err[side][bin]>0.?power(1./B10err[side][bin],2)*B10[side][bin]:0.);

	weightsum = (B2err[side][bin]>0.?power(1./B2err[side][bin],2):0.) + (B10err[side][bin]>0.?power(1./B10err[side][bin],2):0.);
	
	sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }

      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
	0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      superSum[1][bin] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      superSumError[1][bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
      
    } 
  }
 
  boolSuperSum = true;
};

void QuartetAsymmetry::writeAsymToFile() {
  if (!isAsymmetry()) {
    std::cout << "No Asymmetry has been calculated. Not writing to file!\n\n";
    return;
  }

 
  //Setting paths to output files
  std::string outpathA = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) + 
    "Octet_" + itos(octet) + "/QuartetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Quartet_A.dat";
  std::string outpathB = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) 
    + "Octet_" + itos(octet) + "/QuartetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Quartet_B.dat";
  

  //Open and fill file for quartet A
  std::ofstream outfileA(outpathA.c_str()); 
  
  if (isGoodQuartet[0]) {
    for (unsigned int i=0; i<asymmetry[0].size(); i++) {
      outfileA << binLowerEdge[i] << " " << asymmetry[0][i] << " " << asymmetryError[0][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
    }
  }
  else outfileA << "BAD QUARTET";
  outfileA.close();

  //Open and fill file for quartet B
  std::ofstream outfileB(outpathB.c_str());

  if (isGoodQuartet[1]) {
    for (unsigned int i=0; i<asymmetry[1].size(); i++) {
      outfileB << binLowerEdge[i] << " " << asymmetry[1][i] << " " << asymmetryError[1][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
    }
  }
  else outfileB << "BAD QUARTET";
  outfileB.close();

  std::cout << "Wrote Asymmetry to file for anaChoice " << analysisChoice << " in:\n" << outpathA << std::endl << outpathB << "\n";
};

void QuartetAsymmetry::writeSuperSumToFile() {
  if (!isSuperSum()) {
    std::cout << "No Super-sum has been calculated. Not writing to file!\n\n";
    return;
  }

  //Setting paths to output files
  std::string outpathA = Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS")) + 
    "Octet_" + itos(octet) + "/QuartetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Quartet_A.dat";
  std::string outpathB = Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS")) 
    + "Octet_" + itos(octet) + "/QuartetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Quartet_B.dat";

  //Open and fill file for quartet A
  std::ofstream outfileA(outpathA.c_str()); 
  
  if (isGoodQuartet[0]) {
    for (unsigned int i=0; i<superSum[0].size(); i++) {
      outfileA << binLowerEdge[i] << " " << superSum[0][i] << " " << superSumError[0][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << superSum[i] << " " << superSumError[i] << std::endl;
    }
  }
  else outfileA << "BAD QUARTET";
  outfileA.close();

  //Open and fill file for quartet B
  std::ofstream outfileB(outpathB.c_str());

  if (isGoodQuartet[1]) {
    for (unsigned int i=0; i<superSum[1].size(); i++) {
      outfileB << binLowerEdge[i] << " " << superSum[1][i] << " " << superSumError[1][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << superSum[i] << " " << superSumError[i] << std::endl;
    }
  }
  else outfileB << "BAD QUARTET";
  outfileB.close();

};


///////////////////////////////////////////////////////////////////////////////////////////

PairAsymmetry::PairAsymmetry(int oct, std::string anaCh, double enBinWidth, double fidCut, bool ukdata, bool simulation, bool unblind) : AsymmetryBase(oct,anaCh,enBinWidth,fidCut,ukdata,simulation,unblind), totalAsymmetryA0(0.), totalAsymmetryErrorA0(0.), totalAsymmetryB0(0.), totalAsymmetryErrorB0(0.), totalAsymmetryA1(0.), totalAsymmetryErrorA1(0.), totalAsymmetryB1(0.), totalAsymmetryErrorB1(0.), boolSuperSum(false), boolAsymmetry(false) {

  isGoodPair.resize(2,std::vector <bool> (2));

  isGoodPair[0][0] = isPair(0,0);
  isGoodPair[0][1] = isPair(0,1);
  isGoodPair[1][0] = isPair(1,0);
  isGoodPair[1][1] = isPair(1,1);


  if (isGoodPair[0][0] || isGoodPair[1][0] || isGoodPair[0][1] || isGoodPair[1][1]) {
    //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
    asymmetry.resize(2, std::vector < std::vector <double> > (2, std::vector <double>(numEnBins,0.)));
    asymmetryError.resize(2, std::vector < std::vector <double> > (2, std::vector <double>(numEnBins,0.)));
    superSum.resize(2, std::vector < std::vector <double> > (2, std::vector <double>(numEnBins,0.)));
    superSumError.resize(2, std::vector < std::vector <double> > (2, std::vector <double>(numEnBins,0.)));
    loadRates(); // load the rates in the rate vectors for each run
    std::cout <<"//////////////////////////////////////////////////////////////////\n"
	      <<"PairAsymmetry for octet " << octet << std::endl
	      <<"\nPair Status (0=bad, 1=good): First: " << isGoodPair[0][0] << " " 
	      <<"Second: " << isGoodPair[0][1]  << " " <<"Third: " << isGoodPair[1][0] << " "
	      <<"Fourth: " << isGoodPair[1][1] << std::endl;
  }
  else {
    boolAsymmetry=boolSuperSum=true; //Writing "bad pair" to all files so that it isn't used in the future
    writeAsymToFile();
    writeSuperSumToFile();
    
    throw "NO PAIRS IN THIS OCTET";
  }
};

void PairAsymmetry::calcAsymmetryBinByBin() {

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};
      
  // A type runs, pair 0
  if (isGoodPair[0][0]) {
    for (unsigned int bin=0; bin<asymmetry[0][0].size(); bin++) {
      double R = 0.;
      double deltaR = 0.;
      
      for (unsigned int side=0; side<2; side++) {

	// AFP Off
	sfOFF[side] = A2[side][bin];
	sfOFF_err[side] = A2err[side][bin];
	
	// AFP ON
	sfON[side] = A5[side][bin];
	sfON_err[side] = A5err[side][bin];
	
      }

      //Calculate super-ratio, R. Note that any of the rates that go into this could be negative 
      R = (sfOFF[1]*sfON[0]) / (sfON[1]*sfOFF[0]) ;
      
      //Only use rates which are real and not zero to avoid nan (couldn't get TMath::IsNaN() to work properly :( ...)
      if ( R!=0. && sfON[1]!=0. && sfOFF[0]!=0.) { 
	
	deltaR = sqrt( power( sfOFF_err[1]*sfON[0]/(sfOFF[0]*sfON[1]), 2 ) + 
		       power( sfON_err[0]*sfOFF[1]/(sfOFF[0]*sfON[1]), 2 ) +
		       power( sfOFF_err[0]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfOFF[0]*sfON[1]), 2 ) +
		       power( sfON_err[1]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfON[1]*sfON[1]), 2 ) );
	//deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	asymmetry[0][0][bin] = (1.-sqrt(std::abs(R)))/(1+sqrt(std::abs(R)));
	asymmetryError[0][0][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	asymmetry[0][0][bin] = 0.;
	asymmetryError[0][0][bin] = 1.;
      } 

    }
  }

  if (isGoodPair[0][1]) {
    for (unsigned int bin=0; bin<asymmetry[0][1].size(); bin++) {
      double R = 0.;
      double deltaR = 0.;
      
      for (unsigned int side=0; side<2; side++) {

	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;

	// AFP Off
	sfOFF[side] = A10[side][bin];
	sfOFF_err[side] = A10err[side][bin];
	
	// AFP ON
	sfON[side] = A7[side][bin];
	sfON_err[side] = A7err[side][bin];
	
      }

      //Calculate super-ratio, R. Note that any of the rates that go into this could be negative 
      R = (sfOFF[1]*sfON[0]) / (sfON[1]*sfOFF[0]) ;
      
      //Only use rates which are real and not zero to avoid nan (couldn't get TMath::IsNaN() to work properly :( ...)
      if ( R!=0. && sfON[1]!=0. && sfOFF[0]!=0.) { 
	
	deltaR = sqrt( power( sfOFF_err[1]*sfON[0]/(sfOFF[0]*sfON[1]), 2 ) + 
		       power( sfON_err[0]*sfOFF[1]/(sfOFF[0]*sfON[1]), 2 ) +
		       power( sfOFF_err[0]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfOFF[0]*sfON[1]), 2 ) +
		       power( sfON_err[1]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfON[1]*sfON[1]), 2 ) );
	//deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	asymmetry[0][1][bin] = (1.-sqrt(std::abs(R)))/(1+sqrt(std::abs(R)));
	asymmetryError[0][1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	asymmetry[0][1][bin] = 0.;
	asymmetryError[0][1][bin] = 1.;
      } 
    }
  }

  if (isGoodPair[1][0]) {
    for (unsigned int bin=0; bin<asymmetry[1][0].size(); bin++) {
      double R = 0.;
      double deltaR = 0.;
      
      for (unsigned int side=0; side<2; side++) {
	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;

	// AFP Off
	sfOFF[side] = B5[side][bin];
	sfOFF_err[side] = B5err[side][bin];
	
	// AFP ON
	sfON[side] = B2[side][bin];
	sfON_err[side] = B2err[side][bin];
	
      }

      //Calculate super-ratio, R. Note that any of the rates that go into this could be negative 
      R = (sfOFF[1]*sfON[0]) / (sfON[1]*sfOFF[0]) ;
      
      //Only use rates which are real and not zero to avoid nan (couldn't get TMath::IsNaN() to work properly :( ...)
      if ( R!=0. && sfON[1]!=0. && sfOFF[0]!=0.) { 
	
	deltaR = sqrt( power( sfOFF_err[1]*sfON[0]/(sfOFF[0]*sfON[1]), 2 ) + 
		       power( sfON_err[0]*sfOFF[1]/(sfOFF[0]*sfON[1]), 2 ) +
		       power( sfOFF_err[0]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfOFF[0]*sfON[1]), 2 ) +
		       power( sfON_err[1]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfON[1]*sfON[1]), 2 ) );
	//deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	asymmetry[1][0][bin] = (1.-sqrt(std::abs(R)))/(1+sqrt(std::abs(R)));
	asymmetryError[1][0][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	asymmetry[1][0][bin] = 0.;
	asymmetryError[1][0][bin] = 1.;
      } 
    }
  }

  if (isGoodPair[1][1]) {
    for (unsigned int bin=0; bin<asymmetry[1][1].size(); bin++) {
      double R = 0.;
      double deltaR = 0.;
      
      for (unsigned int side=0; side<2; side++) {
	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;

	// AFP Off
	sfOFF[side] = B7[side][bin];
	sfOFF_err[side] = B7err[side][bin];
	
	// AFP ON
	sfON[side] = B10[side][bin];
	sfON_err[side] = B10err[side][bin];
	
      }

      //Calculate super-ratio, R. Note that any of the rates that go into this could be negative 
      R = (sfOFF[1]*sfON[0]) / (sfON[1]*sfOFF[0]) ;
      
      //Only use rates which are real and not zero to avoid nan (couldn't get TMath::IsNaN() to work properly :( ...)
      if ( R!=0. && sfON[1]!=0. && sfOFF[0]!=0.) { 
	
	deltaR = sqrt( power( sfOFF_err[1]*sfON[0]/(sfOFF[0]*sfON[1]), 2 ) + 
		       power( sfON_err[0]*sfOFF[1]/(sfOFF[0]*sfON[1]), 2 ) +
		       power( sfOFF_err[0]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfOFF[0]*sfON[1]), 2 ) +
		       power( sfON_err[1]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfON[1]*sfON[1]), 2 ) );
	//deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	asymmetry[1][1][bin] = (1.-sqrt(std::abs(R)))/(1+sqrt(std::abs(R)));
	asymmetryError[1][1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	asymmetry[1][1][bin] = 0.;
	asymmetryError[1][1][bin] = 1.;
      } 
    }
  }
  
  boolAsymmetry = true;
};

void PairAsymmetry::calcTotalAsymmetry(double enWinLow, double enWinHigh) {
  
  unsigned int binLow = (unsigned int)(enWinLow/energyBinWidth);
  unsigned int binHigh = (unsigned int)(enWinHigh/energyBinWidth)-1;

  double sfON[2]={0.}, sfOFF[2]={0.};
  double sfON_err[2]={0.}, sfOFF_err[2]={0.};
  double sumA2[2]={0.}, sumA5[2]={0.}, sumA7[2]={0.}, sumA10[2]={0.}, sumB2[2]={0.}, sumB5[2]={0.}, sumB7[2]={0.}, sumB10[2]={0.};
  double sumA2_err[2]={0.}, sumA5_err[2]={0.}, sumA7_err[2]={0.}, sumA10_err[2]={0.}, sumB2_err[2]={0.}, sumB5_err[2]={0.}, sumB7_err[2]={0.}, sumB10_err[2]={0.};

  for (unsigned int side=0; side<2; side++) {
    for (unsigned int bin=binLow; bin<=binHigh; bin++) {
      sumA2[side]+=A2[side][bin];
      sumA2_err[side]+=power(A2err[side][bin],2);
      sumA5[side]+=A5[side][bin];
      sumA5_err[side]+=power(A5err[side][bin],2);
      sumA7[side]+=A7[side][bin];
      sumA7_err[side]+=power(A7err[side][bin],2);
      sumA10[side]+=A10[side][bin];
      sumA10_err[side]+=power(A10err[side][bin],2);
      sumB2[side]+=B2[side][bin];
      sumB2_err[side]+=power(B2err[side][bin],2);
      sumB5[side]+=B5[side][bin];
      sumB5_err[side]+=power(B5err[side][bin],2);
      sumB7[side]+=B7[side][bin];
      sumB7_err[side]+=power(B7err[side][bin],2);
      sumB10[side]+=B10[side][bin];
      sumB10_err[side]+=power(B10err[side][bin],2);

    }
    sumA2_err[side] = sqrt(sumA2_err[side]);
    sumA5_err[side] = sqrt(sumA5_err[side]);
    sumA7_err[side] = sqrt(sumA7_err[side]);
    sumA10_err[side] = sqrt(sumA10_err[side]);
    sumB2_err[side] = sqrt(sumB2_err[side]);
    sumB5_err[side] = sqrt(sumB5_err[side]);
    sumB7_err[side] = sqrt(sumB7_err[side]);
    sumB10_err[side] = sqrt(sumB10_err[side]);
    
  }
    
  
  // A type runs, pair 0
  if (isGoodPair[0][0]) {
    
    double R = 0.;
    double deltaR = 0.;
    
    for (unsigned int side=0; side<2; side++) {
      double weightsum=0.;
      
      // AFP Off
      sfOFF[side] = sumA2[side];
      weightsum = (sumA2_err[side]>0.?power(1./sumA2_err[side],2):0.);
      
      sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      
      weightsum=0.;
      
      // AFP ON
      sfON[side] = sumA5[side];
	weightsum = (sumA5_err[side]>0.?power(1./sumA5_err[side],2):0.);
	
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
    }
    
    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      totalAsymmetryA0= (1.-sqrt(R))/(1+sqrt(R));
      totalAsymmetryErrorA0 = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      totalAsymmetryA0 = 0.;
      totalAsymmetryErrorA0= 0.;
    }
  }
  
  
  if (isGoodPair[0][1]) {
    double R = 0.;
    double deltaR = 0.;
    
    for (unsigned int side=0; side<2; side++) {
      sfOFF[side]=0.;
      sfON[side]=0.;
      sfOFF_err[side]=0.;
      sfON_err[side]=0.;
      double weightsum=0.;
      
      // AFP Off
      sfOFF[side] = sumA10[side];
	weightsum = (sumA10_err[side]>0.?power(1./sumA10_err[side],2):0.);
	
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = sumA7[side];
	weightsum = (sumA7_err[side]>0.?power(1./sumA7_err[side],2):0.);
	
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
    }
    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      totalAsymmetryA1= (1.-sqrt(R))/(1+sqrt(R));
      totalAsymmetryErrorA1 = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      totalAsymmetryA1 = 0.;
      totalAsymmetryErrorA1= 0.;
    }
  }

  
  if (isGoodPair[1][0]) {
    double R = 0.;
    double deltaR = 0.;
    
    for (unsigned int side=0; side<2; side++) {
      sfOFF[side]=0.;
      sfON[side]=0.;
      sfOFF_err[side]=0.;
      sfON_err[side]=0.;
      double weightsum=0.;
      
      // AFP Off
      sfOFF[side] = sumB5[side];
      weightsum = (sumB5_err[side]>0.?power(1./sumB5_err[side],2):0.);
      
      sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      
      weightsum=0.;
      
      // AFP ON
      sfON[side] = sumB2[side];
      weightsum = (sumB2_err[side]>0.?power(1./sumB2_err[side],2):0.);
      
      sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      
    }

    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      totalAsymmetryB0= (1.-sqrt(R))/(1+sqrt(R));
      totalAsymmetryErrorB0 = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      totalAsymmetryB0 = 0.;
      totalAsymmetryErrorB0= 0.;
    }

  }


  if (isGoodPair[1][1]) {
    double R = 0.;
    double deltaR = 0.;
    
    for (unsigned int side=0; side<2; side++) {
      sfOFF[side]=0.;
      sfON[side]=0.;
      sfOFF_err[side]=0.;
      sfON_err[side]=0.;
      double weightsum=0.;
      
      // AFP Off
      sfOFF[side] = sumB7[side];
      weightsum = (sumB7_err[side]>0.?power(1./sumB7_err[side],2):0.);
      
      sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      
      weightsum=0.;
      
      // AFP ON
      sfON[side] = sumB10[side];
      weightsum = (sumB10_err[side]>0.?power(1./sumB10_err[side],2):0.);
      
      sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      
    }

    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      totalAsymmetryB1= (1.-sqrt(R))/(1+sqrt(R));
      totalAsymmetryErrorB1 = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      totalAsymmetryB1 = 0.;
      totalAsymmetryErrorB1= 0.;
    }
    

  }
};


void PairAsymmetry::calcSuperSum() {

  //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
  superSum.resize(2, std::vector < std::vector <double> > (2,std::vector <double> (numEnBins, 0.)));
  superSumError.resize(2, std::vector < std::vector <double> > (2,std::vector <double> (numEnBins, 0.)));;

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  

  // A type runs, pair 0
  if (isGoodPair[0][0]) {
    for (unsigned int bin=0; bin<asymmetry[0][0].size(); bin++) {     
      for (unsigned int side=0; side<2; side++) {
	
	
	// AFP Off
	sfOFF[side] = A2[side][bin];
	sfOFF_err[side] = A2err[side][bin];
	
	// AFP ON
	sfON[side] = A5[side][bin];
	sfON_err[side] = A5err[side][bin];

      }
      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
	0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      superSum[0][0][bin] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      superSumError[0][0][bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
    }
  }

  if (isGoodPair[0][1]) {
    for (unsigned int bin=0; bin<asymmetry[0][1].size(); bin++) {
      for (unsigned int side=0; side<2; side++) {
	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;
        

	// AFP Off
	sfOFF[side] = A10[side][bin];
	sfOFF_err[side] = A10err[side][bin];
	
	// AFP ON
	sfON[side] = A7[side][bin];
	sfON_err[side] = A7err[side][bin];
	
      }
      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
	0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      superSum[0][1][bin] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      superSumError[0][1][bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
      
    }
  }

  if (isGoodPair[1][0]) {
    for (unsigned int bin=0; bin<asymmetry[1][0].size(); bin++) {
      for (unsigned int side=0; side<2; side++) {
	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;
        
	// AFP Off
	sfOFF[side] = B5[side][bin];
	sfOFF_err[side] = B5err[side][bin];
	
	// AFP ON
	sfON[side] = B2[side][bin];
	sfON_err[side] = B2err[side][bin];
	
      }
      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
	0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      superSum[1][0][bin] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      superSumError[1][0][bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
    }
  }

  if (isGoodPair[1][1]) {
    for (unsigned int bin=0; bin<asymmetry[1][1].size(); bin++) {
      
      for (unsigned int side=0; side<2; side++) {
	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;
	
	// AFP Off
	sfOFF[side] = B7[side][bin];
	sfOFF_err[side] = B7err[side][bin];
	
	// AFP ON
	sfON[side] = B10[side][bin];
	sfON_err[side] = B10err[side][bin];
	
      }
      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
	0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      superSum[1][1][bin] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      superSumError[1][1][bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
    }
  }  


  boolSuperSum = true;
};

void PairAsymmetry::writeAsymToFile() {
  if (!isAsymmetry()) {
    std::cout << "No Asymmetry has been calculated. Not writing to file!\n\n";
    return;
  }

  //Setting paths to output files
  std::string outpathA0 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) + 
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Pair_A0.dat";
  std::string outpathA1 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) + 
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Pair_A1.dat";
  std::string outpathB0 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) +
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Pair_B0.dat";
  std::string outpathB1 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) +
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Pair_B1.dat";

  //Open and fill file for quartet A
  std::ofstream outfileA0(outpathA0.c_str()); 
  
  if (isGoodPair[0][0]) {
    for (unsigned int i=0; i<asymmetry[0][0].size(); i++) {
      outfileA0 << binLowerEdge[i] << " " << asymmetry[0][0][i] << " " << asymmetryError[0][0][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
    }
  }
  else outfileA0 << "BAD PAIR";
  outfileA0.close();

  std::ofstream outfileA1(outpathA1.c_str()); 
  
  if (isGoodPair[0][1]) {
    for (unsigned int i=0; i<asymmetry[0][1].size(); i++) {
      outfileA1 << binLowerEdge[i] << " " << asymmetry[0][1][i] << " " << asymmetryError[0][1][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
    }
  }
  else outfileA1 << "BAD PAIR";
  outfileA1.close();

  std::ofstream outfileB0(outpathB0.c_str()); 

  if (isGoodPair[1][0]) {
    for (unsigned int i=0; i<asymmetry[1][0].size(); i++) {
      outfileB0 << binLowerEdge[i] << " " << asymmetry[1][0][i] << " " << asymmetryError[1][0][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
    }
  }
  else outfileB0 << "BAD PAIR";
  outfileB0.close();

  std::ofstream outfileB1(outpathB1.c_str()); 
  
  if (isGoodPair[1][1]) {
    for (unsigned int i=0; i<asymmetry[1][1].size(); i++) {
      outfileB1 << binLowerEdge[i] << " " << asymmetry[1][1][i] << " " << asymmetryError[1][1][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[1][1][i] << " " << asymmetryError[1][1][i] << std::endl;
    }
  }
  else outfileB1 << "BAD PAIR";
  outfileB1.close();

  
  std::cout << "Wrote Asymmetry to file for anaChoice " << analysisChoice << "\n";
};

void PairAsymmetry::writeSuperSumToFile() {
  if (!isSuperSum()) {
    std::cout << "No Super-Sum has been calculated. Not writing to file!\n\n";
    return;
  }

  //Setting paths to output files
  std::string outpathA0 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) + 
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Pair_A0.dat";
  std::string outpathA1 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) + 
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Pair_A1.dat";
  std::string outpathB0 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) +
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Pair_B0.dat";
  std::string outpathB1 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) +
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+analysisChoice+"_Pair_B1.dat";
  
  //Open and fill file for quartet A
  std::ofstream outfileA0(outpathA0.c_str()); 
  
  if (isGoodPair[0][0]) {
    for (unsigned int i=0; i<superSum[0][0].size(); i++) {
      outfileA0 << binLowerEdge[i] << " " << superSum[0][0][i] << " " << superSumError[0][0][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << superSum[0][0][i] << " " << superSumError[0][0][i] << std::endl;
    }
  }
  else outfileA0 << "BAD PAIR";
  outfileA0.close();
  
  std::ofstream outfileA1(outpathA1.c_str()); 
  
  if (isGoodPair[0][1]) {
    for (unsigned int i=0; i<superSum[0][1].size(); i++) {
      outfileA1 << binLowerEdge[i] << " " << superSum[0][1][i] << " " << superSumError[0][1][i] << std::endl;
    }
  }
  else outfileA1 << "BAD PAIR";
  outfileA1.close();
  
  std::ofstream outfileB0(outpathB0.c_str()); 
  
  if (isGoodPair[1][0]) {
    for (unsigned int i=0; i<superSum[1][0].size(); i++) {
      outfileB0 << binLowerEdge[i] << " " << superSum[1][0][i] << " " << superSumError[1][0][i] << std::endl;
    }
  }
  else outfileB0 << "BAD PAIR";
  outfileB0.close();
  
  std::ofstream outfileB1(outpathB1.c_str()); 
  
  if (isGoodPair[1][1]) {
    for (unsigned int i=0; i<superSum[1][1].size(); i++) {
      outfileB1 << binLowerEdge[i] << " " << superSum[1][1][i] << " " << superSumError[1][1][i] << std::endl;
    }
  }
  else outfileB1 << "BAD PAIR";
  outfileB1.close();
  
  
  std::cout << "Wrote Super-Sum to file for anaChoice " << analysisChoice << "\n";
  //std::cout << "Wrote to " << outpathB1 << "\n";
};

////////////////////////////////////////////////////////////////////////////////

void calcBinByBinSuperRatio(std::vector< std::vector<double> > &A2, 
			    std::vector< std::vector<double> > &A5,
			    std::vector< std::vector<double> > &A7, 
			    std::vector< std::vector<double> > &A10,
			    std::vector< std::vector<double> > &B2, 
			    std::vector< std::vector<double> > &B5,
			    std::vector< std::vector<double> > &B7,
			    std::vector< std::vector<double> > &B10,
			    std::vector< std::vector<double> > &A2err, 
			    std::vector< std::vector<double> > &A5err,
			    std::vector< std::vector<double> > &A7err, 
			    std::vector< std::vector<double> > &A10err,
			    std::vector< std::vector<double> > &B2err, 
			    std::vector< std::vector<double> > &B5err,
			    std::vector< std::vector<double> > &B7err, 
			    std::vector< std::vector<double> > &B10err,
			    std::vector<double> &asymmetry, std::vector<double> &asymmetryError) {

//arrays to hold the weighted averages of the rates for each spin state and side
  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  for (unsigned int bin=0; bin<asymmetry.size(); bin++) {

    double R = 0.;
    double deltaR = 0.;

    for (unsigned int side=0; side<2; side++) {

      double weightsum=0.; // Holds the sum of the weights for the denominator of the weighted av
      
      //Start with flipper off... Only use rate if it has an error to avoid dividing by 0
      sfOFF[side] = ( ( A2err[side][bin] > 0. ? A2[side][bin]/power(A2err[side][bin],2)   : 0. ) +
		      ( A10err[side][bin]> 0. ? A10[side][bin]/power(A10err[side][bin],2) : 0. ) +
		      (	B5err[side][bin] > 0. ? B5[side][bin]/power(B5err[side][bin],2)   : 0. ) +
		      ( B7err[side][bin] > 0. ? B7[side][bin]/power(B7err[side][bin],2)   : 0. ) );

      weightsum = ( ( A2err[side][bin] > 0. ? 1./power(A2err[side][bin],2)  : 0. ) + 
		    ( A10err[side][bin]> 0. ? 1./power(A10err[side][bin],2) : 0. ) +
		    ( B5err[side][bin] > 0. ? 1./power(B5err[side][bin],2)  : 0. ) +
		    ( B7err[side][bin] > 0. ? 1./power(B7err[side][bin],2)  : 0. ) );

      sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
      sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;


      weightsum=0.;

      // Now for flipper on
      sfON[side] = ( ( B2err[side][bin] > 0. ? B2[side][bin]/power(B2err[side][bin],2)   : 0. ) +
		      ( B10err[side][bin]> 0. ? B10[side][bin]/power(B10err[side][bin],2) : 0. ) +
		      (	A5err[side][bin] > 0. ? A5[side][bin]/power(A5err[side][bin],2)   : 0. ) +
		      ( A7err[side][bin] > 0. ? A7[side][bin]/power(A7err[side][bin],2)   : 0. ) );

      weightsum = ( ( B2err[side][bin] > 0. ? 1./power(B2err[side][bin],2)  : 0. ) + 
		    ( B10err[side][bin]> 0. ? 1./power(B10err[side][bin],2) : 0. ) +
		    ( A5err[side][bin] > 0. ? 1./power(A5err[side][bin],2)  : 0. ) +
		    ( A7err[side][bin] > 0. ? 1./power(A7err[side][bin],2)  : 0. ) );
      
      sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
      sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      

      
    }
    
    //Calculate super-ratio, R. Note that any of the rates that go into this could be negative 
    R = (sfOFF[1]*sfON[0]) / (sfON[1]*sfOFF[0]) ;
    
    //Only use rates which are real and not zero to avoid nan (couldn't get TMath::IsNaN() to work properly :( ...)
    if ( R!=0. && sfON[1]!=0. && sfOFF[0]!=0.) { 
     
      deltaR = sqrt( power( sfOFF_err[1]*sfON[0]/(sfOFF[0]*sfON[1]), 2 ) + 
		     power( sfON_err[0]*sfOFF[1]/(sfOFF[0]*sfON[1]), 2 ) +
		     power( sfOFF_err[0]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfOFF[0]*sfON[1]), 2 ) +
		     power( sfON_err[1]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfON[1]*sfON[1]), 2 ) );
      //deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      asymmetry[bin] = (1.-sqrt(std::abs(R)))/(1+sqrt(std::abs(R)));
      asymmetryError[bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      asymmetry[bin] = 0.;
      asymmetryError[bin] = 1.;
    }

  } 
};
