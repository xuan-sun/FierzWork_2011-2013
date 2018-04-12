/*
Class which allows for generic reading of data and construction 
of histograms for rates of each event type. This can be any type of data coming
from UCNA root files. What is stored should be histograms of event rates binned
by energy along with the run time for that particular run.

Also included is a class to handle reading in reverse calibrated 
simulation data
*/

#include "EvtRateHandler.hh"
#include "MBUtils.hh"
#include "calibrationTools.hh"
#include <cstdio>
#include <fstream>

#include <TString.h>

const bool writeRatesToFile = true;

const std::vector<UInt_t> neutronPoisonedBGruns {17274,18209,21433,21890,22311};


EvtRateHandler::EvtRateHandler(std::vector<int> rn, bool fg, std::string anaCh, double enBinWidth, double fidCut, bool ukdata, bool unblind) : runs(rn),FG(fg),analysisChoice(anaCh),fiducialCut(fidCut),UKdata(ukdata),unblinded(unblind),pol(0) {
  
  bool useOldTimes = false;

  numEnergyBins = (int)(1200./enBinWidth);

  pol = polarization(runs[0]); //THIS NEEDS TO BE REWRITTEN TO READ IN MASTER FILE
  
  rateEvec.resize(numEnergyBins,0.); 
  rateWvec.resize(numEnergyBins,0.);
  rateEerr.resize(numEnergyBins,0.); 
  rateWerr.resize(numEnergyBins,0.);

  rateEvecQuad.resize(4, std::vector<double> (numEnergyBins,0.) ); 
  rateWvecQuad.resize(4, std::vector<double> (numEnergyBins,0.) );
  rateEerrQuad.resize(4, std::vector<double> (numEnergyBins,0.) ); 
  rateWerrQuad.resize(4, std::vector<double> (numEnergyBins,0.) );

  rateEvecRad.resize(6, std::vector<double> (numEnergyBins,0.) ); 
  rateWvecRad.resize(6, std::vector<double> (numEnergyBins,0.) );
  rateEerrRad.resize(6, std::vector<double> (numEnergyBins,0.) ); 
  rateWerrRad.resize(6, std::vector<double> (numEnergyBins,0.) );

  refSpectraE.resize(numEnergyBins,0.);
  refSpectraW.resize(numEnergyBins,0.);
  totalRefTime = 0.;
  totalRefCountsE = 0.;
  totalRefCountsW = 0.;

  loadReferenceSpectra(); // Loading the proper reference spectra

  totalCountsE=0.;
  totalCountsW=0.;
  totalCountsEQuad.resize(4, 0.);
  totalCountsWQuad.resize(4, 0.);
  totalCountsERad.resize(6, 0.);
  totalCountsWRad.resize(6, 0.);
  
  runLength.resize(runs.size(), std::vector<double> (2,0.));
  totalRunLengthE = totalRunLengthW = 0.;
  UCNMonIntegral.resize(runs.size(), 0.);
  
  //Loading information from log file

  for (unsigned int i = 0; i<runs.size(); i++) {
    std::string logFilePath = std::string(getenv("RUN_INFO_FILES"))+"runInfo_"+itos(runs[i])+".dat";
    std::ifstream logFile(logFilePath.c_str());
    std::vector <std::string> title(6);
    std::vector <std::string> value(6); // Fix this after all of the Run Info files have been fixed. It should be 7
    for (int i=0; i<6; i++) {      // but the monitor integral is attached to the next time.
      logFile >> title[i] >> value[i];
    }
    logFile.close();
    
    if (unblinded) {
      runLength[i][0] = atof(value[2].c_str());
      runLength[i][1] = atof(value[2].c_str());
    }
    else {
      runLength[i][0] = atof(value[0].c_str());
      runLength[i][1] = atof(value[1].c_str());
    }
    UCNMonIntegral[i] = 100.;//value[3];

    if ( std::find(neutronPoisonedBGruns.begin(),neutronPoisonedBGruns.end(),runs[i])
	   !=neutronPoisonedBGruns.end() ) 
	{
	  //Skipping first 100 seconds of bg runs which are 
	  // contaminated with neutrons
	  runLength[i][0] -= 100.;
	  runLength[i][1] -= 100.;
	   
	}

    /* if ( useOldTimes ) {
      if (unblinded) {
	runLength[i][0] = value[6];
	runLength[i][1] = value[6];
      }
      else {
	runLength[i][0] = value[4];
	runLength[i][1] = value[5];
      }
      }*/

    totalRunLengthE += runLength[i][0];
    totalRunLengthW += runLength[i][1];
    
    std::cout << title[0] << "\t" << runLength[i][0] << std::endl;
    std::cout << title[1] << "\t" << runLength[i][1] << std::endl;
    std::cout << title[3] << "\t" << UCNMonIntegral[i] << std::endl;
    //std::cout << holdTitle << "\t" << runLength[0] << std::endl;
  }

  // Initialize all histograms
  hisCounts[0] = new TH1D("EastCounts","East Event Hist",numEnergyBins,0.,1200.);
  hisCounts[1] = new TH1D("WestCounts","West Event Hist",numEnergyBins,0.,1200.);

  for ( int i=0; i<4; ++i) {
    hisCountsQuad[i][0] = new TH1D(TString::Format("EastQuad%i",i),TString::Format("EastQuad%i",i),
				   numEnergyBins,0.,1200.);
    hisCountsQuad[i][1] = new TH1D(TString::Format("WestQuad%i",i),TString::Format("WestQuad%i",i),
				   numEnergyBins,0.,1200.);
  }

  for ( int i=0; i<6; ++i) {
    int radLow = i*10;
    int radHigh = (1+1)*10;
    hisCountsRad[i][0] = new TH1D(TString::Format("EastRad%i-%i",radLow,radHigh),
				  TString::Format("EastRad%i-%i",radLow,radHigh),
				  numEnergyBins,0.,1200.);
    hisCountsRad[i][1] = new TH1D(TString::Format("WestRad%i-%i",radLow,radHigh),
				  TString::Format("WestRad%i-%i",radLow,radHigh),
				  numEnergyBins,0.,1200.);
  }
};

EvtRateHandler::~EvtRateHandler() {
  
    if (hisCounts[0]) delete hisCounts[0];
    if (hisCounts[1]) delete hisCounts[1]; 

    for (int i=0; i<4; ++i) {
      if (hisCountsQuad[i][0]) delete hisCountsQuad[i][0];
      if (hisCountsQuad[i][1]) delete hisCountsQuad[i][1]; 
    }
    for (int i=0; i<6; ++i) {
      if (hisCountsRad[i][0]) delete hisCountsRad[i][0];
      if (hisCountsRad[i][1]) delete hisCountsRad[i][1]; 
    }
};

int EvtRateHandler::polarization(int run) {

  std::ifstream infile(TString::Format("%s/masterBetaRunList.txt",getenv("OCTET_LIST")).Data());

  int rn = 0;
  std::string type = "";
  bool runFound = false;

  while ( infile >> rn >> type ) {
    if ( rn == run ) { runFound = true; infile.close(); continue; }
  }

  if ( !runFound ) throw "You didn't pick a run in an octet...";
  //Flipper OFF is p = -1
  if ( type=="A1" || type=="A2" || type=="A10" || type=="A12" || type=="B4" || type=="B5" || type=="B7" || type=="B9" ) return -1;
  else if ( type=="B1" || type=="B2" || type=="B10" || type=="B12" || type=="A4" || type=="A5" || type=="A7" || type=="A9" ) return 1;
  
  else throw "You chose a Depol Run.. no polarization";
  
};

void EvtRateHandler::loadReferenceSpectra() {

  // NEED TO REWRITE THE REFERENCE THINGIES TO ALSO SPIT OUT TOTAL EVENTS RATHER THAN RATES
  
  TString dir = ( FG ? TString::Format("%s/reference_rates/foregroundRatesByAnaChoice/",getenv("ANALYSIS_CODE")) :
		  TString::Format("%s/reference_rates/backgroundRatesByAnaChoice/",getenv("ANALYSIS_CODE")) );
    
  TString filename = ( dir + TString("ReferenceSpectra_Octets-") +
		       ( runs[0]<20000 ? TString("0-59") : TString("60-121") ) +
		       ( pol==-1 ? TString("_sfOFF-") : TString("_sfON-") ) +
		       TString::Format("AnaCh-%s.txt",analysisChoice.c_str()) );

  std::ifstream infile(filename.Data());

  std::string label;
  infile >> label >> totalRefTime;
  infile >> label >> totalRefTimeE;
  infile >> label >> totalRefTimeW;

  if ( unblinded ) totalRefTimeE = totalRefTimeW = totalRefTime;

  double binMid=0., eastRef=0., eastErr=0., westRef=0., westErr=0.;
  int inc = 0;

  while ( infile >> binMid >> eastRef >> eastErr >> westRef >> westErr ) {

    refSpectraE[inc] = eastRef;
    refSpectraW[inc] = westRef;

    totalRefCountsE += eastRef;
    totalRefCountsW += westRef;
    
    ++inc;

    if (inc > 9 && inc < 15) std::cout << binMid << "\t" << eastRef << "\t" << westRef << std::endl; 
    
  }

  infile.close();
  std::cout << "Finished Loading Reference Rate data...\n\n";
  
};

double EvtRateHandler::referenceError(int side, int bin, double totalCounts) {

  //total Counts method
  if ( side == 0 ) return totalRefCountsE>0. ? sqrt( refSpectraE[bin] * ( totalCounts / totalRefCountsE ) ) : 0.;
  else if ( side == 1 ) return totalRefCountsW>0. ? sqrt( refSpectraW[bin] * ( totalCounts / totalRefCountsW ) ) : 0.;

  // total time method
  //if ( side == 0 ) return totalRefTimeE>0. ? sqrt( refSpectraE[bin] ) * ( totalRunLengthE / totalRefTimeE )  : 0.;
  //else if ( side == 1 ) return totalRefTimeW>0. ? sqrt( refSpectraW[bin] ) * ( totalRunLengthW / totalRefTimeW )  : 0.;

  else throw "BAD SIDE GIVEN TO referenceError(int side, int bin)";
};

void EvtRateHandler::CalcRates() {
  
  this->dataReader();

  Double_t binContentE = 0., binContentW = 0.;

  //NEED THE TOTAL NUMBER OF EVENTS IN EACH SIDES HISTOGRAM FOR SCALING
  for (unsigned int i=0; i<numEnergyBins; i++) {

    binContentE = hisCounts[0]->GetBinContent(i+1);
    binContentW = hisCounts[1]->GetBinContent(i+1);
    
    rateEvec[i] = binContentE / totalRunLengthE;
    rateWvec[i] = binContentW / totalRunLengthW;

    rateEerr[i] = ( binContentE > 25. || FG ) ? sqrt(binContentE) : referenceError(0,i,totalCountsE);
    rateWerr[i] = ( binContentW > 25. || FG ) ? sqrt(binContentW) : referenceError(1,i,totalCountsW);

    rateEerr[i] /= totalRunLengthE;
    rateWerr[i] /= totalRunLengthW;

    
    for (int j=0; j<4; ++j) {

      binContentE = hisCountsQuad[j][0]->GetBinContent(i+1);
      binContentW = hisCountsQuad[j][1]->GetBinContent(i+1);

      rateEvecQuad[j][i] = binContentE / totalRunLengthE;
      rateWvecQuad[j][i] = binContentW / totalRunLengthW;

      rateEerrQuad[j][i] = ( binContentE > 25. || FG ) ? sqrt(binContentE) : referenceError(0,i,totalCountsEQuad[j]);
      rateWerrQuad[j][i] = ( binContentW > 25. || FG ) ? sqrt(binContentW) : referenceError(1,i,totalCountsWQuad[j]);

      rateEerrQuad[j][i] /= totalRunLengthE;
      rateWerrQuad[j][i] /= totalRunLengthW;
    }

    for (int j=0; j<6; ++j) {

      binContentE = hisCountsRad[j][0]->GetBinContent(i+1);
      binContentW = hisCountsRad[j][1]->GetBinContent(i+1);

      rateEvecRad[j][i] = binContentE / totalRunLengthE;
      rateWvecRad[j][i] = binContentW / totalRunLengthW;

      rateEerrRad[j][i] = ( binContentE > 25. || FG ) ? sqrt(binContentE) : referenceError(0,i,totalCountsERad[j]);
      rateWerrRad[j][i] = ( binContentW > 25. || FG ) ? sqrt(binContentW) : referenceError(1,i,totalCountsWRad[j]);

      rateEerrRad[j][i] /= totalRunLengthE;
      rateWerrRad[j][i] /= totalRunLengthW;
    }

    /*if ( binContentE < 25. || binContentW < 25. ) {

      std::cout << i << "\t" << binContentE << "\t" << rateEerr[i]*totalRunLengthE << "\t" 
		<< binContentW << "\t" << rateWerr[i]*totalRunLengthW << "\n";
		}*/
	
  }
      
};


void EvtRateHandler::dataReader() {

  bool sep23 = false;

  // 2/3 separation
  if ( analysisChoice==std::string("C") ||
       analysisChoice==std::string("E") ||
       analysisChoice==std::string("H") ||
       analysisChoice==std::string("J") ||
       analysisChoice==std::string("K") )
    sep23 = true;

  //Initializing the separator
  BackscatterSeparator sep;
  sep.LoadCutCurve(runs[0]);


  //Event types
  bool Type0=false, Type1=false, Type2=false, Type3=false;
  if ( analysisChoice=="A" || analysisChoice=="B" || analysisChoice=="C" || analysisChoice=="D" || analysisChoice=="E" ) Type0 = true;
  if ( analysisChoice=="A" || analysisChoice=="B" || analysisChoice=="C" || analysisChoice=="E" || analysisChoice=="F" ) Type1 = true;
  if ( analysisChoice=="A" || analysisChoice=="C" || analysisChoice=="E" || analysisChoice=="G" || analysisChoice=="H" || analysisChoice=="J" ) Type2 = true;
  if ( analysisChoice=="A" || analysisChoice=="C" || analysisChoice=="E" || analysisChoice=="G" || analysisChoice=="H" || analysisChoice=="K") Type3 = true;


  
  char temp[100];
  
  std::string infile;
  TFile *input = new TFile();
  TTree *Tin;

  double EmwpcX=0., EmwpcY=0., WmwpcX=0., WmwpcY=0., TimeE=0., TimeW=0.,Time=0., Erecon=0., MWPCEnergyE=0., MWPCEnergyW=0.; //Branch Variables being read in
  float EmwpcX_f=0., EmwpcY_f=0., WmwpcX_f=0., WmwpcY_f=0., TimeE_f=0., TimeW_f=0., Time_f=0., Erecon_f=0., MWPCEnergyE_f=0., MWPCEnergyW_f=0.; // For reading in data from MPM replays
  int xE_nClipped=0, yE_nClipped=0, xW_nClipped=0, yW_nClipped=0;
  int xE_mult=0, yE_mult=0, xW_mult=0, yW_mult=0;
  int badTimeFlag = 0; 

  int xeRC=0, yeRC=0, xwRC=0, ywRC=0; //  Wirechamber response class variables 

  int PID, Side, Type;

  Bool_t EvnbGood, BkhfGood;

  Bool_t skipFirst100s = false;

  for ( UInt_t i=0; i<runs.size(); i++ ) {

    if ( std::find(neutronPoisonedBGruns.begin(),neutronPoisonedBGruns.end(),runs[i])
	   !=neutronPoisonedBGruns.end() ) 
      {
	skipFirst100s = true;
      }
  
    //Set branch addresses
    if (UKdata) {
      sprintf(temp,"%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),runs[i]);
      std::string infile = std::string(temp);
      input = TFile::Open(infile.c_str(), "READ");
      Tin = (TTree*)input->Get("pass3");
      
      Tin->SetBranchAddress("PID", &PID);
      Tin->SetBranchAddress("Type", &Type);
      Tin->SetBranchAddress("Side", &Side); 
      Tin->SetBranchAddress("Erecon",&Erecon);
      Tin->SetBranchAddress("TimeE",&TimeE);
      Tin->SetBranchAddress("TimeW",&TimeW);
      Tin->SetBranchAddress("Time",&Time);
      Tin->SetBranchAddress("badTimeFlag",&badTimeFlag);
      Tin->GetBranch("xE")->GetLeaf("center")->SetAddress(&EmwpcX);
      Tin->GetBranch("yE")->GetLeaf("center")->SetAddress(&EmwpcY);
      Tin->GetBranch("xW")->GetLeaf("center")->SetAddress(&WmwpcX);
      Tin->GetBranch("yW")->GetLeaf("center")->SetAddress(&WmwpcY);
      //Tin->GetBranch("gaus_xE")->GetLeaf("center")->SetAddress(&EmwpcX);
      //Tin->GetBranch("gaus_yE")->GetLeaf("center")->SetAddress(&EmwpcY);
      //Tin->GetBranch("gaus_xW")->GetLeaf("center")->SetAddress(&WmwpcX);
      //Tin->GetBranch("gaus_yW")->GetLeaf("center")->SetAddress(&WmwpcY);
      //Tin->GetBranch("old_xE")->GetLeaf("center")->SetAddress(&EmwpcX);
      //Tin->GetBranch("old_yE")->GetLeaf("center")->SetAddress(&EmwpcY);
      //Tin->GetBranch("old_xW")->GetLeaf("center")->SetAddress(&WmwpcX);
      //Tin->GetBranch("old_yW")->GetLeaf("center")->SetAddress(&WmwpcY);
      Tin->GetBranch("xE")->GetLeaf("nClipped")->SetAddress(&xE_nClipped);
      Tin->GetBranch("yE")->GetLeaf("nClipped")->SetAddress(&yE_nClipped);
      Tin->GetBranch("xW")->GetLeaf("nClipped")->SetAddress(&xW_nClipped);
      Tin->GetBranch("yW")->GetLeaf("nClipped")->SetAddress(&yW_nClipped);
      Tin->GetBranch("xE")->GetLeaf("mult")->SetAddress(&xE_mult);
      Tin->GetBranch("yE")->GetLeaf("mult")->SetAddress(&yE_mult);
      Tin->GetBranch("xW")->GetLeaf("mult")->SetAddress(&xW_mult);
      Tin->GetBranch("yW")->GetLeaf("mult")->SetAddress(&yW_mult);

      Tin->SetBranchAddress("xeRC", &xeRC);
      Tin->SetBranchAddress("yeRC", &yeRC);
      Tin->SetBranchAddress("xwRC", &xwRC);
      Tin->SetBranchAddress("ywRC", &ywRC);

      Tin->SetBranchAddress("EvnbGood", &EvnbGood);
      Tin->SetBranchAddress("BkhfGood", &BkhfGood);

      Tin->SetBranchAddress("EMWPC_E",&MWPCEnergyE);
      Tin->SetBranchAddress("EMWPC_W",&MWPCEnergyW);
      //NEED TO ADD IN WIRECHAMBER ENERGY FOR 2/3 SEPARATION
    }
    else {
      sprintf(temp,"%s/spec_%i.root",getenv("UK_SPEC_REPLAY"),runs[i]);
      std::string infile = std::string(temp);
      input = TFile::Open(infile.c_str(), "READ");
      Tin = (TTree*)input->Get("phys");
      
      Tin->SetBranchAddress("PID", &PID);
      Tin->SetBranchAddress("Type", &Type);
      Tin->SetBranchAddress("Side", &Side); 
      Tin->SetBranchAddress("Erecon",&Erecon_f);
      Tin->SetBranchAddress("TimeE",&TimeE_f);
      Tin->SetBranchAddress("TimeW",&TimeW_f);
      Tin->SetBranchAddress("Time",&Time_f);
      Tin->SetBranchAddress("badTimeFlag",&badTimeFlag);
      Tin->GetBranch("xEmpm")->GetLeaf("center")->SetAddress(&EmwpcX_f);
      Tin->GetBranch("yEmpm")->GetLeaf("center")->SetAddress(&EmwpcY_f);
      Tin->GetBranch("xWmpm")->GetLeaf("center")->SetAddress(&WmwpcX_f);
      Tin->GetBranch("yWmpm")->GetLeaf("center")->SetAddress(&WmwpcY_f);
      Tin->GetBranch("xEmpm")->GetLeaf("nClipped")->SetAddress(&xE_nClipped);
      Tin->GetBranch("yEmpm")->GetLeaf("nClipped")->SetAddress(&yE_nClipped);
      Tin->GetBranch("xWmpm")->GetLeaf("nClipped")->SetAddress(&xW_nClipped);
      Tin->GetBranch("yWmpm")->GetLeaf("nClipped")->SetAddress(&yW_nClipped);
      Tin->GetBranch("xEmpm")->GetLeaf("mult")->SetAddress(&xE_mult);
      Tin->GetBranch("yEmpm")->GetLeaf("mult")->SetAddress(&yE_mult);
      Tin->GetBranch("xWmpm")->GetLeaf("mult")->SetAddress(&xW_mult);
      Tin->GetBranch("yWmpm")->GetLeaf("mult")->SetAddress(&yW_mult);

      Tin->SetBranchAddress("xeRC", &xeRC);
      Tin->SetBranchAddress("yeRC", &yeRC);
      Tin->SetBranchAddress("xwRC", &xwRC);
      Tin->SetBranchAddress("ywRC", &ywRC);

      Tin->SetBranchAddress("EvnbGood", &EvnbGood);
      Tin->SetBranchAddress("BkhfGood", &BkhfGood);

      Tin->SetBranchAddress("EMWPC_E",&MWPCEnergyE_f);
      Tin->SetBranchAddress("EMWPC_W",&MWPCEnergyW_f);
    }
    unsigned int nevents = Tin->GetEntriesFast();
    std::cout << "Number of Events: " << nevents << std::endl;
    
    
    //Run over all events in file, fill output histogram with rate
    
    double r2E = 0.; //position of event squared
    double r2W = 0.; //position of event squared
    
    for (unsigned int i=0; i<nevents; i++) {
      
      Tin->GetEvent(i);

      if ( !(EvnbGood && BkhfGood) ) continue; //Checking quality of header/footer (I think)
      
      //Cast float variables to double if they are from MPM replay
      if (!UKdata) {
	EmwpcX = (double) EmwpcX_f;
	EmwpcY = (double) EmwpcY_f;
	WmwpcX = (double) WmwpcX_f;
	WmwpcY = (double) WmwpcY_f;
	Erecon = (double) Erecon_f;
	TimeE = (double) TimeE_f;
	TimeW = (double) TimeW_f;
	Time = (double) Time_f;
      }
	
      // Skipping event in first 100s if bg run is contaminated by neutrons
      if ( skipFirst100s ) {
	if ( Time<100. ) { 
	  continue;
	}
      }

      
      if ( PID==1 && badTimeFlag==0 ) {       //  Cut on electrons not cut out by beam drops or bursts
	
	//Cut out clipped events and bad Wirechamber signals
	if ( Type!=0 ) {
	  //if ( xE_nClipped>1 || yE_nClipped>1 || xW_nClipped>1 || yW_nClipped>1 ) continue;
	  if ( xeRC>6 || yeRC>6 || xwRC>6 || ywRC>6 ) continue; //Must look at both sides
	  else if ( xE_mult<1 || yE_mult<1 || xW_mult<1 || yW_mult<1 ) continue;
	}
	else {
	  if ( Side==0 ) {
	    //if ( xE_nClipped>1 || yE_nClipped>1 ) continue;
	    if ( xeRC>6 || yeRC>6 ) continue; //only look at MWPC signal on East
	    if ( xE_mult<1 || yE_mult<1 ) continue;
	  }
	  else if ( Side==1 ) {
	    //if ( xW_nClipped>1 || yW_nClipped>1 ) continue;
	    if ( xwRC>6 || ywRC>6 ) continue; //only look at MWPC signal on West
	    if ( xW_mult<1 || yW_mult<1 ) continue;
	  }
	}
	
	
	// Determine radial event position squared for cutting on fiducial volume
	r2E = EmwpcX*EmwpcX + EmwpcY*EmwpcY;
	r2W = WmwpcX*WmwpcX + WmwpcY*WmwpcY;
	

	if ( r2E>3600. || r2W>3600. ) continue;
	  
	// If the analysis choice calls for separation of 2/3, do it here
	if ( sep23 ) {
	  if (Erecon>0. && Type==2) {
	    
	    if (Side==0) {
	      Type = sep.separate23(MWPCEnergyE);
	      Side = Type==2 ? 1 : 0;
	      //std::cout << "Side 0: " << MWPCEnergyE << "\t" << Type << "\t" << Side << std::endl;
	    }
	    else if (Side==1) {
	      Type = sep.separate23(MWPCEnergyW);
	      Side = Type==2 ? 0 : 1;
	    }
	  }
	}
	
	//Determine Quadrant 
	int quad = 0;
	
	if (Side==0) {
	  if (EmwpcX>0.) quad = ( EmwpcY>0. ? 0 : 3 );
	  else quad = ( EmwpcY>0. ? 1 : 2 );
	}
	else if (Side==1) {
	  if (WmwpcX>0.) quad = ( WmwpcY>0. ? 0 : 3 );
	  else quad = ( WmwpcY>0. ? 1 : 2 );
	}

	
	//Determine Radial bin 
	int rad = (int) ( Side==0 ? (sqrt(r2E)/10.) : (sqrt(r2W)/10.) );
	if ( rad>5 ) rad = 5; // Catch the bad guys
	
	//Type0
	if ( Type0 && Type==0 ) { 
	  if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut ) ) {
	    hisCounts[Side]->Fill(Erecon);
	    hisCountsQuad[quad][Side]->Fill(Erecon);
	  }
	  hisCountsRad[rad][Side]->Fill(Erecon);
	}
	//Type1
	else if ( Type1 && Type==1 ) { 
	  if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut ) ) { 
	    hisCounts[Side]->Fill(Erecon);
	    hisCountsQuad[quad][Side]->Fill(Erecon);
	  }
	  hisCountsRad[rad][Side]->Fill(Erecon);
	}
	//Type2
	else if ( Type2 && Type==2 ) { 
	  if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut ) ) {
	    hisCounts[Side]->Fill(Erecon);
	    hisCountsQuad[quad][Side]->Fill(Erecon);
	  }
	  hisCountsRad[rad][Side]->Fill(Erecon);
	}
	//Type3
	else if ( Type3 && Type==3 ) { 
	  if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut ) ) { 
	    hisCounts[Side]->Fill(Erecon);
	    hisCountsQuad[quad][Side]->Fill(Erecon);
	  }
	  hisCountsRad[rad][Side]->Fill(Erecon);
	} 

      }
    }
    input->Close();
  }

  totalCountsE = (double) hisCounts[0]->Integral(1,hisCounts[0]->GetNbinsX());
  totalCountsW = (double) hisCounts[1]->Integral(1,hisCounts[1]->GetNbinsX());

  for (int i=0; i<4; ++i) {
    totalCountsEQuad[i] = (double) hisCountsQuad[i][0]->Integral(1,hisCountsQuad[i][0]->GetNbinsX());
    totalCountsWQuad[i] = (double) hisCountsQuad[i][1]->Integral(1,hisCountsQuad[i][1]->GetNbinsX());
    std::cout << totalCountsEQuad[i] << " " << totalCountsWQuad[i] << std::endl;
  }

  for (int i=0; i<6; ++i) {
    totalCountsERad[i] = (double) hisCountsRad[i][0]->Integral(1,hisCountsRad[i][0]->GetNbinsX());
    totalCountsWRad[i] = (double) hisCountsRad[i][1]->Integral(1,hisCountsRad[i][1]->GetNbinsX());
  }
    
  std::cout << "Beta Events: " << totalCountsE+totalCountsW << std::endl;

  // Output the integrated crap to file
  /*std::ofstream integralsT("integrals_timeNorm.txt",std::ofstream::app);
  integralsT << runs[0] << "\t\tE = " << hisCounts[0]->Integral(18,77)/totalRunLengthE 
	    << "\t\tW = " << hisCounts[1]->Integral(18,77)/totalRunLengthW << std::endl;

  integralsT.close();

  std::ofstream integralsC("integrals_counts.txt",std::ofstream::app);
  integralsC << runs[0] << "\t\tE = " << hisCounts[0]->Integral(18,77) 
	    << "\t\tW = " << hisCounts[1]->Integral(18,77) << std::endl;

  integralsC.close();
  */
  
  if (input) delete input;
  
};



void SimEvtRateHandler::dataReader() {

  bool sep23 = false;

  // 2/3 separation
  if ( analysisChoice==std::string("C") ||
       analysisChoice==std::string("E") ||
       analysisChoice==std::string("H") ||
       analysisChoice==std::string("J") ||
       analysisChoice==std::string("K") )
    sep23 = true;

  //Event types
  bool Type0=false, Type1=false, Type2=false, Type3=false;
  if ( analysisChoice=="A" || analysisChoice=="B" || analysisChoice=="C" || analysisChoice=="D" || analysisChoice=="E") Type0 = true;
  if ( analysisChoice=="A" || analysisChoice=="B" || analysisChoice=="C" || analysisChoice=="E" || analysisChoice=="F" ) Type1 = true;
  if ( analysisChoice=="A" || analysisChoice=="C" || analysisChoice=="E" || analysisChoice=="G" || analysisChoice=="H" || analysisChoice=="J" ) Type2 = true;
  if ( analysisChoice=="A" || analysisChoice=="C" || analysisChoice=="E" || analysisChoice=="G" || analysisChoice=="H" || analysisChoice=="K") Type3 = true;
  
  char temp[100];
  std::string infile;
  TFile *input = new TFile();
  TTree *Tin;

  double TimeE=0., TimeW=0., Erecon=0., MWPCEnergyE=0., MWPCEnergyW=0.; //Branch Variables being read in
  int PID, Side, Type;

  double AsymWeight=1.;
  double scintPosE[3]={0.};
  double scintPosW[3]={0.};
  double cathRespPosE[3]={0.}; //Position reconstructed the same way it is in data
  double cathRespPosW[3]={0.};
  double mwpcPosE[3]={0.};
  double mwpcPosW[3]={0.}; //holds the position of the event in the MWPC for simulated data
  double EmwpcX=0., EmwpcY=0., WmwpcX=0., WmwpcY=0.;
  int xE_nClipped, yE_nClipped, xW_nClipped, yW_nClipped;
  int xE_mult=0, yE_mult=0, xW_mult=0, yW_mult=0;


  xE_nClipped = yE_nClipped = xW_nClipped = yW_nClipped = 0;
  
  for ( unsigned int i=0; i<runs.size(); i++ ) {
    
    //sprintf(temp,"%s/beta/revCalSim_%i_Beta.root",getenv("REVCALSIM"),runs[i]);
    sprintf(temp,"%s/beta_highStatistics/revCalSim_%i_Beta.root",getenv("REVCALSIM"),runs[i]);
    infile = std::string(temp);
    input = TFile::Open(infile.c_str(), "READ");
    Tin = (TTree*)input->Get("revCalSim");
    //std::cout << "made it here in dataReader\n";
    //Set branch addresses
    
    Tin->SetBranchAddress("PID", &PID);
    Tin->SetBranchAddress("type", &Type);
    Tin->SetBranchAddress("side", &Side); 
    Tin->SetBranchAddress("Erecon",&Erecon);
    Tin->GetBranch("time")->GetLeaf("timeE")->SetAddress(&TimeE);
    Tin->GetBranch("time")->GetLeaf("timeW")->SetAddress(&TimeW);
    Tin->GetBranch("ScintPosAdjusted")->GetLeaf("ScintPosAdjE")->SetAddress(scintPosE);
    Tin->GetBranch("ScintPosAdjusted")->GetLeaf("ScintPosAdjW")->SetAddress(scintPosW);
    Tin->GetBranch("cathRespPos")->GetLeaf("cathRespPosE")->SetAddress(cathRespPosE);
    Tin->GetBranch("cathRespPos")->GetLeaf("cathRespPosW")->SetAddress(cathRespPosW);
    Tin->GetBranch("MWPCPos")->GetLeaf("MWPCPosE")->SetAddress(mwpcPosE);
    Tin->GetBranch("MWPCPos")->GetLeaf("MWPCPosW")->SetAddress(mwpcPosW);
    Tin->GetBranch("xE")->GetLeaf("center")->SetAddress(&EmwpcX);
    Tin->GetBranch("yE")->GetLeaf("center")->SetAddress(&EmwpcY);
    Tin->GetBranch("xW")->GetLeaf("center")->SetAddress(&WmwpcX);
    Tin->GetBranch("yW")->GetLeaf("center")->SetAddress(&WmwpcY);
    //Tin->SetBranchAddress("AsymWeight",&AsymWeight);
    //Tin->SetBranchAddress("nClipped_EX",&nClipped_EX);
    //Tin->SetBranchAddress("nClipped_EY",&nClipped_EY);
    //Tin->SetBranchAddress("nClipped_WX",&nClipped_WX);
    //Tin->SetBranchAddress("nClipped_WY",&nClipped_WY);
    //Tin->GetBranch("gaus_xE")->GetLeaf("center")->SetAddress(&EmwpcX);
    //Tin->GetBranch("gaus_yE")->GetLeaf("center")->SetAddress(&EmwpcY);
    //Tin->GetBranch("gaus_xW")->GetLeaf("center")->SetAddress(&WmwpcX);
    //Tin->GetBranch("gaus_yW")->GetLeaf("center")->SetAddress(&WmwpcY);
    //Tin->GetBranch("old_xE")->GetLeaf("center")->SetAddress(&EmwpcX);
    //Tin->GetBranch("old_yE")->GetLeaf("center")->SetAddress(&EmwpcY);
    //Tin->GetBranch("old_xW")->GetLeaf("center")->SetAddress(&WmwpcX);
    //Tin->GetBranch("old_yW")->GetLeaf("center")->SetAddress(&WmwpcY);
    Tin->GetBranch("xE")->GetLeaf("nClipped")->SetAddress(&xE_nClipped);
    Tin->GetBranch("yE")->GetLeaf("nClipped")->SetAddress(&yE_nClipped);
    Tin->GetBranch("xW")->GetLeaf("nClipped")->SetAddress(&xW_nClipped);
    Tin->GetBranch("yW")->GetLeaf("nClipped")->SetAddress(&yW_nClipped);
    Tin->GetBranch("xE")->GetLeaf("mult")->SetAddress(&xE_mult);
    Tin->GetBranch("yE")->GetLeaf("mult")->SetAddress(&yE_mult);
    Tin->GetBranch("xW")->GetLeaf("mult")->SetAddress(&xW_mult);
    Tin->GetBranch("yW")->GetLeaf("mult")->SetAddress(&yW_mult);
    Tin->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyE")->SetAddress(&MWPCEnergyE);
    Tin->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyW")->SetAddress(&MWPCEnergyW);

    // ADD IN WIRECHAMBER ENERGY DEPOSITION
    

    unsigned int nevents = Tin->GetEntriesFast();
    //std::cout << nevents << std::endl;
    
    //Run over all events in file, fill output histogram with rate
    
    double r2E = 0.; //position of event squared
    double r2W = 0.;

    BackscatterSeparator sep;
    sep.LoadCutCurve(runs[0]);
    
    for (unsigned int i=0; i<nevents; ++i) {
      
      Tin->GetEvent(i);

      if ( Type>3  || Side>1 || Erecon<=0. ) continue;
      
      if (PID==1) {

	//Cut out clipped events
	if ( Type!=0 ) {
	  //if (xE_nClipped>0 || yE_nClipped>0 || xW_nClipped>0 || yW_nClipped>0) continue;
	  if (xE_mult<1 || yE_mult<1 || xW_mult<1 || yW_mult<1) continue;
	}
	else {
	  if ( Side==0 ) {
	    //if ( xE_nClipped>0 || yE_nClipped>0 )  continue;
	    if ( xE_mult<1 || yE_mult<1 )  continue;
	  }
	  else if ( Side==1 ) {
	    //if ( xW_nClipped>0 || yW_nClipped>0 )  continue;
	    if ( xW_mult<1 || yW_mult<1 )  continue;
	  }
	}

	//r2E = ( mwpcPosE[0]*mwpcPosE[0] + mwpcPosE[1]*mwpcPosE[1] ) * 0.6 * 10.;
	//r2W = ( mwpcPosW[0]*mwpcPosW[0] + mwpcPosW[1]*mwpcPosW[1] ) * 0.6 * 10.;
	//r2E = ( cathRespPosE[0]*cathRespPosE[0] + cathRespPosE[1]*cathRespPosE[1] ) ; //Transforming to decay trap coords
	//r2W = ( cathRespPosW[0]*cathRespPosW[0] + cathRespPosW[1]*cathRespPosW[1] ) ;
	r2E = ( EmwpcX*EmwpcX + EmwpcY*EmwpcY ) ; //Transforming to decay trap coords
	r2W = ( WmwpcX*WmwpcX + WmwpcY*WmwpcY );
	  
	  
	if ( r2E>3600. || r2W>3600. ) continue;

	if ( sep23 ) {
	  if (Type==2) {
	    
	    if (Side==0) {
	      Type = sep.separate23(MWPCEnergyE);
	      Side = Type==2 ? 1 : 0;
	    }
	    else if (Side==1) {
	      Type = sep.separate23(MWPCEnergyW);
	      Side = Type==2 ? 0 : 1;
	    }
	  }
	}
	  
	//Determine Quadrant 
	int quad = 0;
	
	if (Side==0) {
	  if (EmwpcX>0.) quad = ( EmwpcY>0. ? 0 : 3 );
	  else quad = ( EmwpcY>0. ? 1 : 2 );
	}
	else if (Side==1) {
	  if (WmwpcX>0.) quad = ( WmwpcY>0. ? 0 : 3 );
	  else quad = ( WmwpcY>0. ? 1 : 2 );
	}

	//Determine Radial bin 
	int rad = (int) ( Side==0 ? (sqrt(r2E)/10.) : (sqrt(r2W)/10.) );
	if ( rad>5 ) rad = 5; // catch any bad guys
	

	// Fill histograms according to event types in this analysis choice
	
	//Type0
	if ( Type0 && Type==0 ) { 
	  if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut ) ) {
	    hisCounts[Side]->Fill(Erecon);
	    hisCountsQuad[quad][Side]->Fill(Erecon);
	  }
	  hisCountsRad[rad][Side]->Fill(Erecon);
	}
	//Type1
	else if ( Type1 && Type==1 ) { 
	  if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut ) ) { 
	    hisCounts[Side]->Fill(Erecon);
	    hisCountsQuad[quad][Side]->Fill(Erecon);
	  }
	  hisCountsRad[rad][Side]->Fill(Erecon);
	}
	//Type2
	else if ( Type2 && Type==2 ) { 
	  if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut ) ) {
	    hisCounts[Side]->Fill(Erecon);
	    hisCountsQuad[quad][Side]->Fill(Erecon);
	  }
	  hisCountsRad[rad][Side]->Fill(Erecon);
	}
	//Type3
	else if ( Type3 && Type==3 ) { 
	  if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut ) ) { 
	    hisCounts[Side]->Fill(Erecon);
	    hisCountsQuad[quad][Side]->Fill(Erecon);
	  }
	  hisCountsRad[rad][Side]->Fill(Erecon);
	} 
	
      }
    }
  
    input->Close();
  }
  totalCountsE = (double) hisCounts[0]->Integral();
  totalCountsW = (double) hisCounts[1]->Integral();

  for (int i=0; i<4; ++i) {
    totalCountsEQuad[i] = (double) hisCountsQuad[i][0]->Integral(1,hisCountsQuad[i][0]->GetNbinsX());
    totalCountsWQuad[i] = (double) hisCountsQuad[i][1]->Integral(1,hisCountsQuad[i][1]->GetNbinsX());
  }

  for (int i=0; i<6; ++i) {
    totalCountsERad[i] = (double) hisCountsRad[i][0]->Integral(1,hisCountsRad[i][0]->GetNbinsX());
    totalCountsWRad[i] = (double) hisCountsRad[i][1]->Integral(1,hisCountsRad[i][1]->GetNbinsX());
  }
    
  std::cout << "Beta Events: " << totalCountsE+totalCountsW << std::endl;

  if (input) delete input;
  
};


BGSubtractedRate::BGSubtractedRate(std::vector <int> fgruns, std::vector <int> bgruns, std::string anaCh,
				   double enBin, double fidCut, bool ukdata, bool sim, bool unblind) : FGruns(fgruns), BGruns(bgruns), analysisChoice(anaCh), EnergyBinWidth(enBin), fiducialCut(fidCut), UKdata(ukdata), Simulation(sim), UNBLIND(unblind) {
  
  int numBins = int(1200./double(EnergyBinWidth));

  BetaRateE.resize(numBins,0.);
  BGRateE.resize(numBins,0.);
  FinalRateE.resize(numBins,0.);
  BetaRateW.resize(numBins,0.);
  BGRateW.resize(numBins,0.);
  FinalRateW.resize(numBins,0.);
  BetaRateErrorE.resize(numBins,0.);
  BGRateErrorE.resize(numBins,0.);
  FinalRateErrorE.resize(numBins,0.);
  BetaRateErrorW.resize(numBins,0.);
  BGRateErrorW.resize(numBins,0.);
  FinalRateErrorW.resize(numBins,0.);

  BetaRateEQuad.resize(4,std::vector<double>(numBins,0.));
  BGRateEQuad.resize(4,std::vector<double>(numBins,0.));
  FinalRateEQuad.resize(4,std::vector<double>(numBins,0.));
  BetaRateWQuad.resize(4,std::vector<double>(numBins,0.));
  BGRateWQuad.resize(4,std::vector<double>(numBins,0.));
  FinalRateWQuad.resize(4,std::vector<double>(numBins,0.));
  BetaRateErrorEQuad.resize(4,std::vector<double>(numBins,0.));
  BGRateErrorEQuad.resize(4,std::vector<double>(numBins,0.));
  FinalRateErrorEQuad.resize(4,std::vector<double>(numBins,0.));
  BetaRateErrorWQuad.resize(4,std::vector<double>(numBins,0.));
  BGRateErrorWQuad.resize(4,std::vector<double>(numBins,0.));
  FinalRateErrorWQuad.resize(4,std::vector<double>(numBins,0.));

  BetaRateERad.resize(6,std::vector<double>(numBins,0.));
  BGRateERad.resize(6,std::vector<double>(numBins,0.));
  FinalRateERad.resize(6,std::vector<double>(numBins,0.));
  BetaRateWRad.resize(6,std::vector<double>(numBins,0.));
  BGRateWRad.resize(6,std::vector<double>(numBins,0.));
  FinalRateWRad.resize(6,std::vector<double>(numBins,0.));
  BetaRateErrorERad.resize(6,std::vector<double>(numBins,0.));
  BGRateErrorERad.resize(6,std::vector<double>(numBins,0.));
  FinalRateErrorERad.resize(6,std::vector<double>(numBins,0.));
  BetaRateErrorWRad.resize(6,std::vector<double>(numBins,0.));
  BGRateErrorWRad.resize(6,std::vector<double>(numBins,0.));
  FinalRateErrorWRad.resize(6,std::vector<double>(numBins,0.));


  runLengthBeta.resize(2,0.);
  runLengthBG.resize(2,0.);
};


void BGSubtractedRate::calcBGSubtRates() {

  /*if (Simulation) 
    {
      SimEvtRateHandler *evt = new SimEvtRateHandler(FGruns, true, analysisChoice, EnergyBinWidth, fiducialCut, UNBLIND);
      evt->CalcRates();
      FinalRateE = evt->getRateVectors(0);
      FinalRateW = evt->getRateVectors(1);
      FinalRateErrorE = evt->getRateErrors(0);
      FinalRateErrorW = evt->getRateErrors(1);
      delete evt;
      //for (unsigned int t=0;t<FinalRateE.size();t++) {
      //for (unsigned int i=0;i<FinalRateE[0].size();i++) {
      //  FinalRateErrorE[t][i] = sqrt(FinalRateE[t][i]);
      //  FinalRateErrorW[t][i] = sqrt(FinalRateW[t][i]);
      //}
      //}
      }*/
 
  LoadRatesByBin();
  std::cout << "Loaded Rates for run " << FGruns[0]  << std::endl;
  CalcFinalRate();
  
};

std::vector<double > BGSubtractedRate::returnBGSubtRate(int side) {  

  if (side==0) return FinalRateE;
  else if (side==1) return FinalRateW;
  else throw "Invalid side when returning BG subtracted rate";
  
};

std::vector<double > BGSubtractedRate::returnBGSubtRateError(int side) {
  
  if (side==0) return FinalRateErrorE;
  else if (side==1) return FinalRateErrorW;
  else throw "Invalid side when returning BG subtracted rate error";
  
};

std::vector< double > BGSubtractedRate::returnBGSubtRateQuad(int side, int quad) {  

  if (side==0) return FinalRateEQuad[quad];
  else if (side==1) return FinalRateWQuad[quad];
  else throw "Invalid side when returning BG subtracted rate";
  
};

std::vector< double > BGSubtractedRate::returnBGSubtRateErrorQuad(int side, int quad) {
  
  if (side==0) return FinalRateErrorEQuad[quad];
  else if (side==1) return FinalRateErrorWQuad[quad];
  else throw "Invalid side when returning BG subtracted rate error";
  
};

std::vector< double > BGSubtractedRate::returnBGSubtRateRad(int side, int rad) {  

  if (side==0) return FinalRateERad[rad];
  else if (side==1) return FinalRateWRad[rad];
  else throw "Invalid side when returning BG subtracted rate";
  
};

std::vector< double > BGSubtractedRate::returnBGSubtRateErrorRad(int side, int rad) {
  
  if (side==0) return FinalRateErrorERad[rad];
  else if (side==1) return FinalRateErrorWRad[rad];
  else throw "Invalid side when returning BG subtracted rate error";
  
};


void BGSubtractedRate::LoadRatesByBin() {
  
  if (!Simulation) {
    
    EvtRateHandler *evtBG = new EvtRateHandler(BGruns, false, analysisChoice, EnergyBinWidth, fiducialCut, UKdata, UNBLIND);
    evtBG->CalcRates();
    BGRateE =          evtBG->getRateVectors(0);
    BGRateErrorE =     evtBG->getRateErrors(0);
    BGRateW =          evtBG->getRateVectors(1);
    BGRateErrorW =     evtBG->getRateErrors(1);
    BGRateEQuad =      evtBG->getRateVectorsQuad(0);
    BGRateErrorEQuad = evtBG->getRateErrorsQuad(0);
    BGRateWQuad =      evtBG->getRateVectorsQuad(1);
    BGRateErrorWQuad = evtBG->getRateErrorsQuad(1);
    BGRateERad =       evtBG->getRateVectorsRad(0);
    BGRateErrorERad =  evtBG->getRateErrorsRad(0);
    BGRateWRad =       evtBG->getRateVectorsRad(1);
    BGRateErrorWRad =  evtBG->getRateErrorsRad(1);
    runLengthBG[0] =   evtBG->returnRunLength(0);
    runLengthBG[1] =   evtBG->returnRunLength(1);
    delete evtBG;
    
    EvtRateHandler *evt = new EvtRateHandler(FGruns, true, analysisChoice, EnergyBinWidth, fiducialCut, UKdata, UNBLIND);
    evt->CalcRates();
    BetaRateE =          evt->getRateVectors(0);
    BetaRateErrorE =     evt->getRateErrors(0);
    BetaRateW =          evt->getRateVectors(1);
    BetaRateErrorW =     evt->getRateErrors(1);
    BetaRateEQuad =      evt->getRateVectorsQuad(0);
    BetaRateErrorEQuad = evt->getRateErrorsQuad(0);
    BetaRateWQuad =      evt->getRateVectorsQuad(1);
    BetaRateErrorWQuad = evt->getRateErrorsQuad(1);
    BetaRateERad =       evt->getRateVectorsRad(0);
    BetaRateErrorERad =  evt->getRateErrorsRad(0);
    BetaRateWRad =       evt->getRateVectorsRad(1);
    BetaRateErrorWRad =  evt->getRateErrorsRad(1);
    runLengthBeta[0] =   evt->returnRunLength(0);
    runLengthBeta[1] =   evt->returnRunLength(1);
    delete evt;
  }

  else { 

    //Note that for the BG runs for simulation, we only need the errors. The rates remain initialized to zero
    EvtRateHandler *evtBG = new EvtRateHandler(BGruns, false, analysisChoice, EnergyBinWidth, fiducialCut, true, UNBLIND);
    evtBG->CalcRates();
    BGRateErrorE =     evtBG->getRateErrors(0);
    BGRateErrorW =     evtBG->getRateErrors(1);
    BGRateErrorEQuad = evtBG->getRateErrorsQuad(0);
    BGRateErrorWQuad = evtBG->getRateErrorsQuad(1);
    BGRateErrorERad = evtBG->getRateErrorsRad(0);
    BGRateErrorWRad = evtBG->getRateErrorsRad(1);
    runLengthBG[0] =   evtBG->returnRunLength(0);
    runLengthBG[1] =   evtBG->returnRunLength(1);
    delete evtBG;
    
    // Simulation data
    SimEvtRateHandler *evt = new SimEvtRateHandler(FGruns, true, analysisChoice, EnergyBinWidth, fiducialCut, UNBLIND);
    evt->CalcRates();
    BetaRateE =          evt->getRateVectors(0);
    BetaRateErrorE =     evt->getRateErrors(0);
    BetaRateW =          evt->getRateVectors(1);
    BetaRateErrorW =     evt->getRateErrors(1);
    BetaRateEQuad =      evt->getRateVectorsQuad(0);
    BetaRateErrorEQuad = evt->getRateErrorsQuad(0);
    BetaRateWQuad =      evt->getRateVectorsQuad(1);
    BetaRateErrorWQuad = evt->getRateErrorsQuad(1);
    BetaRateERad =       evt->getRateVectorsRad(0);
    BetaRateErrorERad =  evt->getRateErrorsRad(0);
    BetaRateWRad =       evt->getRateVectorsRad(1);
    BetaRateErrorWRad =  evt->getRateErrorsRad(1);
    runLengthBeta[0] =   evt->returnRunLength(0);
    runLengthBeta[1] =   evt->returnRunLength(1);
    delete evt;

  }

  if ( Simulation || writeRatesToFile ) {
    std::ofstream ofileE, ofileW, ofileEquad, ofileWquad, ofileErad, ofileWrad;
    if ( !Simulation ) {
      
      ofileE.open(TString::Format("BinByBinComparison/%sUK_Erun%i_anaCh%s.dat",(UNBLIND?"UNBLINDED_":""),
				  FGruns[0],analysisChoice.c_str()).Data());
      ofileW.open(TString::Format("BinByBinComparison/%sUK_Wrun%i_anaCh%s.dat",(UNBLIND?"UNBLINDED_":""),
				  FGruns[0],analysisChoice.c_str()).Data());

      ofileEquad.open(TString::Format("BinByBinComparison/%sUK_Erun%i_anaCh%s_Quadrants.dat",(UNBLIND?"UNBLINDED_":""),
				      FGruns[0],analysisChoice.c_str()).Data());
      ofileWquad.open(TString::Format("BinByBinComparison/%sUK_Wrun%i_anaCh%s_Quadrants.dat",(UNBLIND?"UNBLINDED_":""),
				      FGruns[0],analysisChoice.c_str()).Data());
      
      ofileErad.open(TString::Format("BinByBinComparison/%sUK_Erun%i_anaCh%s_Radial.dat",(UNBLIND?"UNBLINDED_":""),
				     FGruns[0],analysisChoice.c_str()).Data());
      ofileWrad.open(TString::Format("BinByBinComparison/%sUK_Wrun%i_anaCh%s_Radial.dat",(UNBLIND?"UNBLINDED_":""),
				     FGruns[0],analysisChoice.c_str()).Data());
    }
    else {
      ofileE.open(TString::Format("BinByBinComparison/%sSIM_Erun%i_anaCh%s.dat",(UNBLIND?"UNBLINDED_":""),
				  FGruns[0],analysisChoice.c_str()).Data());
      ofileW.open(TString::Format("BinByBinComparison/%sSIM_Wrun%i_anaCh%s.dat",(UNBLIND?"UNBLINDED_":""),
				  FGruns[0],analysisChoice.c_str()).Data());
      
      ofileEquad.open(TString::Format("BinByBinComparison/%sSIM_Erun%i_anaCh%s_Quadrants.dat",(UNBLIND?"UNBLINDED_":""),
				      FGruns[0],analysisChoice.c_str()).Data());
      ofileWquad.open(TString::Format("BinByBinComparison/%sSIM_Wrun%i_anaCh%s_Quadrants.dat",(UNBLIND?"UNBLINDED_":""),
				      FGruns[0],analysisChoice.c_str()).Data());
      
      ofileErad.open(TString::Format("BinByBinComparison/%sSIM_Erun%i_anaCh%s_Radial.dat",(UNBLIND?"UNBLINDED_":""),
				     FGruns[0],analysisChoice.c_str()).Data());
      ofileWrad.open(TString::Format("BinByBinComparison/%sSIM_Wrun%i_anaCh%s_Radial.dat",(UNBLIND?"UNBLINDED_":""),
				     FGruns[0],analysisChoice.c_str()).Data());
    }
    
    ofileE << "FG time = " << runLengthBeta[0] << "\tBG time = " << runLengthBG[0] << std::endl;
    ofileW << "FG time = " << runLengthBeta[1] << "\tBG time = " << runLengthBG[1] << std::endl;
    ofileEquad << "FG time = " << runLengthBeta[0] << "\tBG time = " << runLengthBG[0] << std::endl;
    ofileWquad << "FG time = " << runLengthBeta[1] << "\tBG time = " << runLengthBG[1] << std::endl;
    ofileErad << "FG time = " << runLengthBeta[0] << "\tBG time = " << runLengthBG[0] << std::endl;
    ofileWrad << "FG time = " << runLengthBeta[1] << "\tBG time = " << runLengthBG[1] << std::endl;

    for ( unsigned int i = 0; i < BetaRateE.size(); ++i ) {
      ofileE << i*10.+5. << "\t" << BetaRateE[i] << "\t" << BGRateE[i] << std::endl;
      ofileW << i*10.+5. << "\t" << BetaRateW[i] << "\t" << BGRateW[i] << std::endl;

      ofileEquad << i*10.+5. ;
      ofileWquad << i*10.+5. ;

      for ( unsigned int j = 0; j<4; ++j ) {
	ofileEquad << "\t" << BetaRateEQuad[j][i] << "\t" << BGRateEQuad[j][i] ;
	ofileWquad << "\t" << BetaRateWQuad[j][i] << "\t" << BGRateWQuad[j][i] ;
      }

      ofileEquad << "\n";
      ofileWquad << "\n";

      ofileErad << i*10.+5. ;
      ofileWrad << i*10.+5. ;

      for ( unsigned int j = 0; j<6; ++j ) {
	ofileErad << "\t" << BetaRateERad[j][i] << "\t" << BGRateERad[j][i] ;
	ofileWrad << "\t" << BetaRateWRad[j][i] << "\t" << BGRateWRad[j][i] ;
      }

      ofileErad << "\n";
      ofileWrad << "\n";

    }
    ofileE.close(), ofileW.close(), ofileEquad.close(), ofileWquad.close(), ofileErad.close(), ofileWrad.close();
  }
};

std::vector<double> BGSubtractedRate::returnRunLengths(bool beta) {
  if (beta) return runLengthBeta;
  else return runLengthBG;
}

void BGSubtractedRate::CalcFinalRate()  {
  std::cout << "Calculating Final Rates for run " << FGruns[0] << std::endl;
  
  if (BetaRateE.size()==BGRateE.size()) {
    FinalRateE.resize(BetaRateE.size(),0.);
    for (unsigned int i=0; i<BetaRateE.size(); i++)
      {
	FinalRateE[i] = BetaRateE[i]-BGRateE[i];
	FinalRateErrorE[i] = sqrt(power(BetaRateErrorE[i],2)+power(BGRateErrorE[i],2));
	//std::cout << BetaRateE[i] << " " << BGRateE[i] << " " << FinalRateE[i] << std::endl;
      }
  }
  else throw "Number of energy bins do not agree between Beta and BG runs. Can't calculate final rate!";
  
  if (BetaRateW.size()==BGRateW.size()) {
    FinalRateW.resize(BetaRateW.size(),0.);
    for (unsigned int i=0; i<BetaRateW.size(); i++)
      {
	FinalRateW[i] = BetaRateW[i]-BGRateW[i];
	FinalRateErrorW[i] = sqrt(power(BetaRateErrorW[i],2)+power(BGRateErrorW[i],2));
	//std::cout << BetaRate[i] << " " << BGRate[i] << " " << FinalRate[i] << std::endl;
      }
  }
  else throw "Number of energy bins do not agree between Beta and BG runs. Can't calculate final rate!";

  for (int j=0; j<4; ++j) {

    if (BetaRateEQuad[j].size()==BGRateEQuad[j].size()) {
      FinalRateEQuad[j].resize(BetaRateEQuad[j].size(),0.);
      for (unsigned int i=0; i<BetaRateEQuad[j].size(); i++)
	{
	  FinalRateEQuad[j][i] = BetaRateEQuad[j][i]-BGRateEQuad[j][i];
	  FinalRateErrorEQuad[j][i] = sqrt(power(BetaRateErrorEQuad[j][i],2)+power(BGRateErrorEQuad[j][i],2));
	  //std::cout << BetaRateE[i] << " " << BGRateE[i] << " " << FinalRateE[i] << std::endl;
	}
    }
    else throw "Number of energy bins do not agree between Beta and BG runs. Can't calculate final rate!";
    
    if (BetaRateWQuad[j].size()==BGRateWQuad[j].size()) {
      FinalRateWQuad[j].resize(BetaRateWQuad[j].size(),0.);
      for (unsigned int i=0; i<BetaRateWQuad[j].size(); i++)
	{
	  FinalRateWQuad[j][i] = BetaRateWQuad[j][i]-BGRateWQuad[j][i];
	  FinalRateErrorWQuad[j][i] = sqrt(power(BetaRateErrorWQuad[j][i],2)+power(BGRateErrorWQuad[j][i],2));
	  //std::cout << BetaRateE[i] << " " << BGRateE[i] << " " << FinalRateE[i] << std::endl;
	}
    }
    else throw "Number of energy bins do not agree between Beta and BG runs. Can't calculate final rate!";
  }

  for (int j=0; j<6; ++j) {

    if (BetaRateERad[j].size()==BGRateERad[j].size()) {
      FinalRateERad[j].resize(BetaRateERad[j].size(),0.);
      for (unsigned int i=0; i<BetaRateERad[j].size(); i++)
	{
	  FinalRateERad[j][i] = BetaRateERad[j][i]-BGRateERad[j][i];
	  FinalRateErrorERad[j][i] = sqrt(power(BetaRateErrorERad[j][i],2)+power(BGRateErrorERad[j][i],2));
	  //std::cout << BetaRateE[i] << " " << BGRateE[i] << " " << FinalRateE[i] << std::endl;
	}
    }
    else throw "Number of energy bins do not agree between Beta and BG runs. Can't calculate final rate!";
    
    if (BetaRateWRad[j].size()==BGRateWRad[j].size()) {
      FinalRateWRad[j].resize(BetaRateWRad[j].size(),0.);
      for (unsigned int i=0; i<BetaRateWRad[j].size(); i++)
	{
	  FinalRateWRad[j][i] = BetaRateWRad[j][i]-BGRateWRad[j][i];
	  FinalRateErrorWRad[j][i] = sqrt(power(BetaRateErrorWRad[j][i],2)+power(BGRateErrorWRad[j][i],2));
	  //std::cout << BetaRateE[i] << " " << BGRateE[i] << " " << FinalRateE[i] << std::endl;
	}
    }
    else throw "Number of energy bins do not agree between Beta and BG runs. Can't calculate final rate!";
  }

};


