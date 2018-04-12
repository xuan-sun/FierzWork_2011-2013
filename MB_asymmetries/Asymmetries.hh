#ifndef ASYMMETRIES_HH
#define ASYMMETRIES_HH

#include <map>
#include <string>
#include <vector>
#include <iostream>



// base class for loading octet information
class AsymmetryBase {
public:
  AsymmetryBase(int oct, std::string anaCh, double enBinWidth, double fidCut, bool ukdata, bool simulation, bool unblind=false);
  ~AsymmetryBase() {}

  void readOctetFile(); //populates the map runType
  int getBGrunNumber(std::string runType); //returns the BG runNumber of current Beta run
  int getBGrunNumber(int run); //returns the BG runNumber of current Beta run
  std::string getBGrunLabel(std::string runType); //returns the BG run label of current Beta run
  int getRunsInOctet() {return runsInOctet;}
  bool isFullOctet(); //Checks whether the octet has all the necessary beta runs
  bool isPair(int quartNum, int pairNum); // There are 4 pairs in every Octet (A0,A1,B0,B1)->((0,0),(0,1),),(1,0),(1,1))
  bool isFullQuartet(int quartNum); //There are 2 quartets in every octet (0,1)
  void loadRates(); //Loads the BG subtracted rates for each run and side
  std::vector <int> makeRunVec(std::string runType); // returns a vector of the runs associated with a run t
  bool checkIfBetaRun(std::string type); //Checks if run of type is a beta run or not
  void writeRatesToFile(); //For now this only uses type0 evts
  std::string getCurrentAnaChoice() {return analysisChoice;}
  std::vector < std::vector<double> > returnBGsubtractedRate(std::string runType); // returns A2, A5, etc below
  std::vector < std::vector<double> > returnBGsubtractedRateError(std::string runType);
  
  std::vector < std::vector<double> > returnBGsubtractedRateQuad(std::string runType, int quad); // returns A2, A5, etc below
  std::vector < std::vector<double> > returnBGsubtractedRateErrorQuad(std::string runType, int quad);
  std::vector < std::vector<double> > returnBGsubtractedRateRad(std::string runType, int rad); // returns A2, A5, etc below
  std::vector < std::vector<double> > returnBGsubtractedRateErrorRad(std::string runType, int rad);
  //std::vector < double > returnRunLength(std::string runType); //Returns both beta and BG run length

  

protected:
  
  bool UKdata; //Boolean which is true for UK data and false if mpm
  bool Simulation; //Boolean to use simulated data

  //*******************************
  //       UNBLINDING BOOLEAN
  bool UNBLIND;

  int octet; // Holds the octet being analyzed
  
  double energyBinWidth;
  unsigned int numEnBins;
  double fiducialCut;
  int runsInOctet;
  std::string analysisChoice;
  std::map < int, std::string > runType;
  //std::map <std::string,int> runType; // the key is the run type, mapped val is run number

  // The following vectors are bin by bin rates for each evt type on each side
  // for example, A2[side][bin] (A2[0-1][0-nBins])

  //The following vectors store the length of the runs where A2len[0][0]=A2East length and A2[1][0] = A2East BG run length(A1)
  std::vector < std::vector < double > > A2len, A5len, A7len, A10len, B2len, B5len, B7len, B10len;

  // The following vectors sum over the appropriate event types for the analysisChoice [side][bin]
  std::vector < std::vector <double> > A2, A5, A7, A10, B2, B5, B7, B10; //BG subtr rates for each beta run 
  std::vector < std::vector <double> > A2err, A5err, A7err, A10err, B2err, B5err, B7err, B10err; //BG subtr rates for each beta run

  // The following vectors sum over the appropriate event types for the analysisChoice [side][bin]
  std::vector < std::vector < std::vector <double> > > A2Quad, A5Quad, A7Quad, A10Quad, B2Quad, B5Quad, B7Quad, B10Quad; //BG subtr rates for each beta run 
  std::vector < std::vector < std::vector <double> > > A2errQuad, A5errQuad, A7errQuad, A10errQuad, B2errQuad, B5errQuad, B7errQuad, B10errQuad; //BG subtr rates for each beta run

  // The following vectors sum over the appropriate event types for the analysisChoice [side][bin]
  std::vector < std::vector < std::vector <double> > > A2Rad, A5Rad, A7Rad, A10Rad, B2Rad, B5Rad, B7Rad, B10Rad; //BG subtr rates for each beta run 
  std::vector < std::vector < std::vector <double> > > A2errRad, A5errRad, A7errRad, A10errRad, B2errRad, B5errRad, B7errRad, B10errRad; //BG subtr rates for each beta run

  std::vector<double> binLowerEdge;
  std::vector<double> binUpperEdge; //Hold the Energy of the upper and lower bin edges
};

/////////////////////////////////////////////////////////////////////////////////

class OctetAsymmetry : public AsymmetryBase {
public:
  OctetAsymmetry(int oct, std::string anaCh, double enBinWidth=10., double fidCut=45., bool ukdata=true, bool simulation=false, bool unblind=false);
  ~OctetAsymmetry() {std::cout << "\n\n\n";}

  void makePlots();
  void calcAsymmetryBinByBin(); //Calculates the raw asymmetry bin by bin to be written to file and plotted
  void calcNCSUSumAsymmetryBinByBin(); //Checking the NCSU analyzer sum method
  void calcTotalAsymmetry(double enWinLow, double enWinHigh); //Returns total raw asymmetry over energy window
  void calcSuperSum(); //Calculates the super sum over the entire octet for spectral comparisons
  void calcSuperSumNCSUstyle(); //Calculates the super sum using the NCSU rate summing method
  void writeAsymToFile();
  void writeSuperSumToFile();
  double returnTotalAsymmetry() {return totalAsymmetry;}
  double returnTotalAsymmetryError() {return totalAsymmetryError;}
  bool isSuperSum() {return boolSuperSum;}
  bool isAsymmetry() {return boolAsymmetry;}

  std::vector <double> returnSuperSum() {return superSum;}
  std::vector <double> returnSuperSumError() {return superSumError;}

private:

  std::vector <double> asymmetry; //Raw asymmetry in bins
  std::vector <double> asymmetryError; //Raw Asymmetry error in bins

  std::vector < std::vector<double> > asymmetryQuad; //Raw asymmetry in bins
  std::vector < std::vector<double> > asymmetryErrorQuad; //Raw Asymmetry error in bins
  
  std::vector < std::vector<double> > asymmetryRad; //Raw asymmetry in bins
  std::vector < std::vector<double> > asymmetryErrorRad; //Raw Asymmetry error in bins

  std::vector <double> superSum; //Raw asymmetry in bins
  std::vector <double> superSumError; //Raw Asymmetry error in bins
  double totalAsymmetry; //Bin summed asymmetry
  double totalAsymmetryError;
  bool boolSuperSum;
  bool boolAsymmetry;
};

/////////////////////////////////////////////////////////////////////////////////

class QuartetAsymmetry : public AsymmetryBase {
public:
  QuartetAsymmetry(int oct, std::string anaCh, double enBinWidth=10., double fidCut=45., bool ukdata=true, bool simulation=false, bool unblind=false);
  ~QuartetAsymmetry() {std::cout << "\n\n\n";}

  void makePlots();
  void calcAsymmetryBinByBin(); //Calculates the raw asymmetry bin by bin to be written to file and plotted
  void calcTotalAsymmetry(double enWinLow, double enWinHigh); //Returns total raw asymmetry over energy window
  void calcSuperSum(); //Calculates the super sum over the entire octet for spectral comparisons
  void writeAsymToFile();
  void writeSuperSumToFile();
  double returnTotalAsymmetry_QuartetA() {return totalAsymmetryA;}
  double returnTotalAsymmetryError_QuartetA() {return totalAsymmetryErrorA;}
  double returnTotalAsymmetry_QuartetB() {return totalAsymmetryB;}
  double returnTotalAsymmetryError_QuartetB() {return totalAsymmetryErrorB;}
  bool boolGoodQuartet(int quart) {return isGoodQuartet[quart];}
  bool isSuperSum() {return boolSuperSum;}
  bool isAsymmetry() {return boolAsymmetry;}

  std::vector < std::vector <double> > returnSuperSum() {return superSum;}
  std::vector < std::vector <double> > returnSuperSumError() {return superSumError;}

private:
  std::vector < std::vector < double > > asymmetry; //Raw asymmetry in bins for quartet 0,1 
                                                   // where 0->Atype and 1->Btype
  std::vector < std::vector < double > > asymmetryError; //Raw Asymmetry error in bins
  std::vector < std::vector <double> > superSum; //Raw asymmetry in bins
  std::vector < std::vector <double> > superSumError; //Raw Asymmetry error in bins
  double totalAsymmetryA; //Bin summed asymmetry for Atype
  double totalAsymmetryErrorA;
  double totalAsymmetryB; //Bin summed asymmetry for Btype
  double totalAsymmetryErrorB;
  std::vector <bool> isGoodQuartet;
  bool boolSuperSum;
  bool boolAsymmetry;
};


/////////////////////////////////////////////////////////////////////////////////
class PairAsymmetry : public AsymmetryBase {
public:
  PairAsymmetry(int oct, std::string anaCh, double enBinWidth=10., double fidCut=45., bool ukdata=true, bool simulation=false, bool unblind=false);
  ~PairAsymmetry() {std::cout << "\n\n\n";}

  void makePlots();
  void calcAsymmetryBinByBin(); //Calculates the raw asymmetry bin by bin to be written to file and plotted
  void calcTotalAsymmetry(double enWinLow, double enWinHigh); //Returns total raw asymmetry over energy window
  void calcSuperSum(); //Calculates the super sum over the entire octet for spectral comparisons
  void writeAsymToFile();
  void writeSuperSumToFile();
  double returnTotalAsymmetry_PairA0() {return totalAsymmetryA0;}
  double returnTotalAsymmetryError_PairA0() {return totalAsymmetryErrorA0;}
  double returnTotalAsymmetry_PairB0() {return totalAsymmetryB0;}
  double returnTotalAsymmetryError_PairB0() {return totalAsymmetryErrorB0;}
  double returnTotalAsymmetry_PairA1() {return totalAsymmetryA1;}
  double returnTotalAsymmetryError_PairA1() {return totalAsymmetryErrorA1;}
  double returnTotalAsymmetry_PairB1() {return totalAsymmetryB1;}
  double returnTotalAsymmetryError_PairB1() {return totalAsymmetryErrorB1;}
  bool boolGoodPair(int quart, int pair) {return isGoodPair[quart][pair];}
  bool isSuperSum() {return boolSuperSum;}
  bool isAsymmetry() {return boolAsymmetry;}
  
  std::vector < std::vector < std::vector <double> > > returnSuperSum() {return superSum;}
  std::vector < std::vector < std::vector <double> > > returnSuperSumError() {return superSumError;}

private:
  std::vector < std::vector < std::vector < double > > > asymmetry; //Raw asymmetry in bins for quartet 0,1 
                                                   // where [0][0]->Atype pair 0 , [1][1]->Btype pair 1 etc.
  std::vector < std::vector < std::vector < double > > > asymmetryError; //Raw Asymmetry error in bins
  std::vector < std::vector < std::vector <double> > > superSum; //Raw asymmetry in bins
  std::vector < std::vector < std::vector <double> > > superSumError; //Raw Asymmetry error in bins
  double totalAsymmetryA0; //Bin summed asymmetry for Atype pair 0
  double totalAsymmetryErrorA0;
  double totalAsymmetryB0; //Bin summed asymmetry for Btype pair 0
  double totalAsymmetryErrorB0;
  double totalAsymmetryA1; //Bin summed asymmetry for Atype pair 1
  double totalAsymmetryErrorA1;
  double totalAsymmetryB1; //Bin summed asymmetry for Btype pair 1
  double totalAsymmetryErrorB1;
  std::vector < std::vector < bool > > isGoodPair;
  bool boolSuperSum;
  bool boolAsymmetry;
};

#endif
