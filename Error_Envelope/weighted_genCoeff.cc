// Standard C and C++ libraries
#include         <vector>
#include         <iostream>
#include         <algorithm>
#include         <functional>
#include         <fstream>
#include         <sstream>
#include         <stdlib.h>
#include         <math.h>
#include         <string.h>
#include         <time.h>
#include         <assert.h>

// Pretty much all the ROOT libraries I have ever used.
#include         <TROOT.h>
#include         <TSystem.h>
#include         <TMath.h>
#include         <TF1.h>
#include         <TGaxis.h>
#include         <TGraph.h>
#include         <TGraphErrors.h>
#include         <TCanvas.h>
#include         <TApplication.h>
#include         <TH1.h>
#include         <TProfile.h>
#include         <TObjArray.h>
#include         <TStyle.h>
#include         <TMarker.h>
#include         <TPaveStats.h>
#include         <TPaveText.h>
#include         <TFile.h>
#include         <TLegend.h>
#include         <TLegendEntry.h>
#include         <TH2F.h>
#include         <TRandom.h>
#include         <TTree.h>
#include         <TChain.h>
#include         <TObjArray.h>
#include         <TFractionFitter.h>
#include         <TLatex.h>
#include         <TMatrixD.h>
#include	 <TRandom3.h>

using            namespace std;

// Fundamental constants that get used
const double m_e = 511.00;                                              ///< electron mass, keV/c^2

// Input and output names and paths used in the code.
// My pseudo version of environment variables.
#define		PARAM_FILE_NAME		"test_weighted_genCoeff.txt"
#define		INPUT_EQ2ETRUE_PARAMS_2010	"/home/xuansun/Documents/MBrown_Work/ParallelAnalyzer/simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat"
#define		INPUT_EQ2ETRUE_PARAMS_2011	"/home/xuansun/Documents/MBrown_Work/ParallelAnalyzer/simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat"
#define		INPUT_EQ2ETRUE_PARAMS_2012	"/home/xuansun/Documents/MBrown_Work/ParallelAnalyzer/simulation_comparison/EQ2EtrueConversion/2012-2013_EQ2EtrueFitParams.dat"

// Plotting functions.
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString command);
void PlotFunc(TCanvas *C, int styleIndex, int canvasIndex, TF1 *fPlot, TString command);

// Michael Brown's SimulationProcessor.cpp uses this to convert to Erecon
vector < vector < vector <double> > > GetEQ2EtrueParams(string geometry);

// Does the Erecon calculation using calibration information.
double CalculateErecon(double totalEvis, vector < vector < vector <double> > > tempEQ2Etrue, int type, int side);

// Michael Mendenhall's 2010 error envelope.
TF1* ErrorEnvelope_2010(double factor);

// implementing the 2011-2012 error envelope
void ErrorEnvelope_2011();
void ErrorEnvelope_2012();
double converter1top2011(double *x, double *par);
double converter2top2011(double *x, double *par);
double converter1bot2011(double *x, double *par);
double converter2bot2011(double *x, double *par);
double converter1top2012(double *x, double *par);
double converter2top2012(double *x, double *par);
double converter1bot2012(double *x, double *par);
double converter2bot2012(double *x, double *par);
void LoadEnvelopeHistogram_2011();
void LoadEnvelopeHistogram_2012();

// variables related to running the 2011-2012 error envelopes.
TH1D *hEnvelope2011 = new TH1D("2011-2012", "2011-2012", 110, 0, 1100);
TF1* errEnv2011_top_1sigma;
TF1* errEnv2011_top_2sigma;
TF1* errEnv2011_bot_1sigma;
TF1* errEnv2011_bot_2sigma;

// same but for 2012-2013 error envelopes
TH1D *hEnvelope2012 = new TH1D("2012-2013", "2012-2013", 110, 0, 1100);
TF1* errEnv2012_top_1sigma;
TF1* errEnv2012_top_2sigma;
TF1* errEnv2012_bot_1sigma;
TF1* errEnv2012_bot_2sigma;

// Perform a single twiddle so we can loop over it in main(), check against a save condition.
// Return whether or not the thrown polynomial passed the save condition.
bool PerformVariation(double a, double b, double c, double d, int numPassed,
                      vector < vector < vector <double> > > EQ2Etrue, TRandom3 *factor,
		      int passNumber, int sideIndex);

// probability weighted on twiddle based on error envelope at source energies
double ProbTwiddleValidity(vector <double> convertedTwiddle, vector <double> energyAxis);

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

// Testing histogram for plotting stuff of interest.
vector <TH1D*> histErecon;

// Counters to allow us to properly sample beyond 1 sigma distribution.
int num1sigma = 0;
int num2sigma = 0;

// coefficients we are saving on first pass. Speeds up program
vector < vector <double> > goodTwiddles;

// Takes an Evis function and converts it to an Erecon function.
struct EreconFunction
{
  EreconFunction(vector < vector < vector <double> > > calibrationCoeffs, int type, int side, TF1 *f_var):
  COEFFS(calibrationCoeffs), TYPE(type), SIDE(side), XAXIS(f_var) {}

  double operator() (double *x, double *p) const
  {
    return CalculateErecon(XAXIS->EvalPar(x,p), COEFFS, TYPE, SIDE);
  }

  vector < vector < vector <double> > > COEFFS;
  int TYPE;
  int SIDE;
  TF1* XAXIS;
};

// Creates a twiddle function in Erecon space. Mostly this just does a subtraction.
struct TwiddleFunctionErecon
{
  TwiddleFunctionErecon(TF1* base, TF1* varied): BASE(base), VARIED(varied) {}

  double operator() (double *x, double *p) const
  {
    return VARIED -> EvalPar(x, p) - BASE -> EvalPar(x, p);
  }

  TF1* BASE;
  TF1* VARIED;
};

// ----------------------------------------------------------------- //
// -------------------- Start of program code. --------------------- //
// ----------------------------------------------------------------- //

int main(int argc, char *argv[])
{
  // Ensures the seed is different for randomizing in ROOT.
  TRandom3* engine = new TRandom3(0);
  gRandom->SetSeed(0);

  // Start the plotting stuff so we can loop and use "SAME" as much as possible.
  TCanvas *C = new TCanvas("canvas", "canvas");
  C -> Divide(3, 2);
  C -> cd(1);
  gROOT->SetStyle("Plain");

  LoadEnvelopeHistogram_2011();
  LoadEnvelopeHistogram_2012();
  ErrorEnvelope_2011();
  ErrorEnvelope_2012();

  errEnv2011_top_2sigma -> GetYaxis() -> SetRangeUser(-15, 15);
  errEnv2011_top_2sigma -> GetYaxis() -> SetTitle("E_{recon} Error (keV)");
  errEnv2011_top_2sigma -> GetXaxis() -> SetTitle("E_{recon} (keV)");
  errEnv2011_top_2sigma -> SetTitle("Non-linearity Polynomial Variations, 2011-2012");
  errEnv2011_top_2sigma -> SetLineStyle(2);
  errEnv2011_top_2sigma -> Draw();

  // Create histograms at fixed Erecon values to look at distribution of polynomials.
  histErecon.push_back(new TH1D("Ce2011-2012", "Erecon = 130", 90, -15, 15));
  histErecon.push_back(new TH1D("Sn2011-2012", "Erecon = 368", 90, -15, 15));
  histErecon.push_back(new TH1D("BiLow2011-2012", "Erecon = 498", 90, -15, 15));
  histErecon.push_back(new TH1D("BiHigh2011-2012", "Erecon = 994", 90, -15, 15));

  // Load the converter to get Erecon from a single EQ value.
  cout << "Using following calibration for 2011-2012 geometry to convert Evis to Erecon..." << endl;
  vector < vector < vector <double> > > converter = GetEQ2EtrueParams("2011-2012");

  int counter, numberSaved;
  counter = 0;
  numberSaved = 0;
  // outer loop, j, is the side index.
  for(int j = 0; j <= 1; j++)
  {
//    for(double a = 2; a <= 2; a = a + 0.1)
    for(double a = -5.0; a <= 5.0; a = a + 0.5)
    {
//      for(double b = 0; b <= 0; b = b + 0.0001)
      for(double b = -0.1; b <= 0.1; b = b + 1e-3)
      {
//        for(double c = 0; c <= 0; c = c + 1e-6)
        for(double c = -1e-4; c <= 1e-4; c = c + 2e-5)
        {
//          for(double d = -1e-7; d <= 1e-7; d = d + 5e-8)
	  for(double d = 0; d <= 0; d++)
          {
            bool save = PerformVariation(a, b, c, d, numberSaved, converter, engine, 1, j);

            // A couple of counters and print-out statements to follow along
	    if(save == true)
	    {
	      numberSaved++;
	    }
            if(counter % 10000 == 0)
            {
     	      cout << "First pass on coefficients. Checking thrown polynomial number... " << counter << endl;
	    }
            counter++;
          }
        }
      }
    }
  }

  ofstream outfile;
  outfile.open(PARAM_FILE_NAME, ios::app);
  numberSaved = 0;
  cout << "Number of good twiddles is " << goodTwiddles.size() << endl;
  for(int j = 0; j <= 1; j++)
  {
    for(unsigned int i = 0; i < goodTwiddles.size(); i++)
    {
      bool save2 = PerformVariation(goodTwiddles[i][0], goodTwiddles[i][1], goodTwiddles[i][2], goodTwiddles[i][3],
				    numberSaved, converter, engine, 2, j);

      if(abs(goodTwiddles[i][0]) < 1e-10)
      {
        goodTwiddles[i][0] = 0;
      }
      if(abs(goodTwiddles[i][1]) < 1e-10)
      {
        goodTwiddles[i][1] = 0;
      }
      if(abs(goodTwiddles[i][2]) < 1e-10)
      {
        goodTwiddles[i][2] = 0;
      }
      if(abs(goodTwiddles[i][3]) < 1e-10)
      {
        goodTwiddles[i][3] = 0;
      }

      if(save2 == true)
      {
        outfile	<< numberSaved << "\t"
		<< goodTwiddles[i][0] << "\t"
		<< goodTwiddles[i][1] << "\t"
		<< goodTwiddles[i][2] << "\t"
		<< goodTwiddles[i][3] << "\t"
		<< goodTwiddles[i][0] << "\t"
		<< goodTwiddles[i][1] << "\t"
		<< goodTwiddles[i][2] << "\t"
		<< goodTwiddles[i][3] << "\n";

        numberSaved++;
      }
    }
  }

  outfile.close();

  cout << "\nNumber of twiddle coefficients thrown: " << counter << endl;
  cout << "Number of twiddle coefficients saved: "<< numberSaved << "\n" << endl;

  cout << "Number of polynomials in first band: " << num1sigma << endl;
  cout << "Number of polynomials in second band: " << num2sigma << endl;

  // Placed here so 1 sigma error envelope goes on top.
  errEnv2011_top_1sigma -> SetLineStyle(2);
  errEnv2011_top_1sigma -> Draw("SAME");
  errEnv2011_bot_1sigma -> SetLineStyle(2);
  errEnv2011_bot_1sigma -> Draw("SAME");
  errEnv2011_bot_2sigma -> SetLineStyle(2);
  errEnv2011_bot_2sigma -> Draw("SAME");
  TLine *line = new TLine(0, 0, 1000, 0);
  line->Draw("SAME");

  // Plot all the additional Erecon slice histograms
  for(unsigned int i = 0; i < histErecon.size(); i++)
  {
    C->cd(i+2);
    histErecon[i]->Draw();
  }

  // Save our plot and print it out as a pdf.
  C -> Print("output_weighted_genCoeff.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

bool PerformVariation(double a, double b, double c, double d, int numPassed,
		      vector < vector < vector <double> > > EQ2Etrue, TRandom3* factor,
		      int passNumber, int sideIndex)
{
  // booleans needed for first pass to count successes
  bool firstBand = true;
  bool secondBand = true;

  bool saveCondition = true;
  bool throwCondition = false;

  double xMin = 0.1;	// For all polynomial ranges, in Evis units.
  double xMax = 1000;

  // Generate twiddle polynomials in EQ space (equivalently, Evis space).
  TF1* pure_Evis = new TF1("pureE", "x", xMin, xMax);
  TF1* twiddle_Evis = new TF1("polyE", "[0] + (1+[1])*x + [2]*x*x + [3]*x*x*x", xMin, xMax);
  twiddle_Evis -> SetParameter(0, a);
  twiddle_Evis -> SetParameter(1, b);
  twiddle_Evis -> SetParameter(2, c);
  twiddle_Evis -> SetParameter(3, d);

  // Calculate Erecon
  int typeIndex = 0;

  // Create our twiddled and untwiddled functions in Erecon space.
  TF1 *Erecon0_East = new TF1("Erecon0", EreconFunction(EQ2Etrue, typeIndex, sideIndex, pure_Evis), xMin, xMax, 0);
  TF1* Erecon_Twiddle_East = new TF1("Erecon_twiddle", EreconFunction(EQ2Etrue, typeIndex, sideIndex, twiddle_Evis), xMin, xMax, 0);
  TF1* delta_Erecon_East = new TF1("Delta_Erecon", TwiddleFunctionErecon(Erecon0_East, Erecon_Twiddle_East), xMin, xMax, 0);

  // Create arrays of Erecon0 and EreconError so we can scatter plot them (hence have error as function of Erecon0).
  vector <double> Evis_axis;
  vector <double> Erecon0_values;
  vector <double> delta_Erecon_values;
  double Evis_min = 1;
  double Evis_max = xMax;
  double Evis_step = 1;
  int nbPoints = 0;
  for(int i = Evis_min; i <= Evis_max; i = i + Evis_step)
  {
    Evis_axis.push_back(i);     // note: i is in whatever units Evis is in.
    Erecon0_values.push_back(Erecon0_East -> Eval(i));
    delta_Erecon_values.push_back(delta_Erecon_East -> Eval(i));
    nbPoints++;
  }

  // Create our scatter plot as a TGraph.
  TGraph* graph = new TGraph(nbPoints, &(Erecon0_values[0]), &(delta_Erecon_values[0]));

  // Get our error envelopes so we can check polynomial values against them.
//  TF1* errEnv1 = ErrorEnvelope_2010(1);
//  TF1* errEnv2 = ErrorEnvelope_2010(2);
  TF1* errEnv1 = errEnv2011_top_1sigma;
  TF1* errEnv2 = errEnv2011_top_2sigma;

  // Check our polynomial (the scatter plot) against a save condition.
  double x, y;
  // These are to save values for histogram plots later.
  double v1 = -10;
  double v2 = -10;
  double v3 = -10;
  double v4 = -10;
  for(int i = 1; i <= graph->GetN(); i++)
  {
    graph->GetPoint(i, x, y);

    // This is to plot a Erecon slice histogram
    if(x > 130 && x < 131)
    {
      v1 = y;
    }
    else if(x > 368 && x < 369)
    {
      v2 = y;
    }
    else if(x > 467.5 && x < 468.5)
    {
      v3 = y;
    }
    else if(x > 993 && x < 994)
    {
      v4 = y;
    }

    // if, at any point, we are over 2 sigma away, exit and don't save and don't throw a number.
    if(abs(y) > errEnv2->Eval(x))
    {
      if(passNumber == 1)
      {
        firstBand = false;
        secondBand = false;
        break;
      }
      else if(passNumber == 2)
      {
        saveCondition = false;
        throwCondition = false;
        break;
      }
    }
    // if we are ever between 1 and 2 sigma and never outside 2 sigma, set flag to throw number to true
    else if(abs(y) > errEnv1->Eval(x) && abs(y) < errEnv2->Eval(x))
    {
      if(passNumber == 1)
      {
        firstBand = false;
        secondBand = true;
      }
      if(passNumber == 2)
      {
        throwCondition = true;
      }
    }
  }

  if(passNumber == 1)
  {
    if(firstBand == true)	// i.e. if we complete the search and never trip either flag, we're in first band
    {
      secondBand = false;
    }
    if(firstBand == true)	// now do the counting.
    {
      num1sigma++;
    }
    else if(secondBand == true)
    {
      num2sigma++;
    }
    if(firstBand == true || secondBand == true)
    {
      // this implements the Gaussian (hopefully) weighting
      if(factor->Rndm() < ProbTwiddleValidity(delta_Erecon_values, Evis_axis))
      {
        vector <double> temp;
        temp.push_back(a);
        temp.push_back(b);
        temp.push_back(c);
        temp.push_back(d);
        goodTwiddles.push_back(temp);
      }
    }
  }
  else if(passNumber == 2)
  {
    // if our condition is tripped, means our curve lies (not exclusively from below) between 1 and 2 sigma
    if(throwCondition == true)
    {
      double percentageToSave = 1;
//      double percentageToSave = (double)num1sigma/num2sigma;
//      percentageToSave = percentageToSave*(0.2718/0.6827);      // this factor accounts for the % difference between 1sigma and 2sigma
      if(factor->Rndm() < (1 - percentageToSave))
      {
        saveCondition = false;
      }
    }
  }

  if(passNumber == 2)
  {
    if(saveCondition == true)
    {
      histErecon[0] -> Fill(v1);
      histErecon[1] -> Fill(v2);
      histErecon[2] -> Fill(v3);
      histErecon[3] -> Fill(v4);
      // Plotting stuff
      graph->SetLineColor(numPassed % 50);
      graph->Draw("SAME");
    }
    // memory management. Delete the left-over pointers. Absolutely necessary or program doesn't run.
    else if(saveCondition == false)
    {
      delete graph;
    }
  }
  else if(passNumber == 1)
  {
    delete graph;
  }
  delete pure_Evis;
  delete twiddle_Evis;
  delete Erecon0_East;
  delete Erecon_Twiddle_East;
  delete delta_Erecon_East;
//  delete errEnv2;
//  delete errEnv1;

  return saveCondition;
}

void PlotFunc(TCanvas *C, int styleIndex, int canvasIndex, TF1 *fPlot, TString command)
{
  C -> cd(canvasIndex);

  fPlot->SetLineColor(styleIndex % 50);	// only 50 colors in set line color.
  fPlot->GetYaxis()->SetRangeUser(-40, 40);
  fPlot->GetYaxis()->SetTitle("Erecon error");
  fPlot->GetXaxis()->SetTitle("Evis");

  fPlot->Draw(command);
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString command)
{
  C -> cd(canvasIndex);
  gPlot->SetLineColor(styleIndex);
  gPlot->GetYaxis()->SetRangeUser(-40, 40);
  gPlot->Draw(command);
}

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Energy (keV)");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("Counts (N)");
  hPlot -> GetYaxis() -> CenterTitle();

  if(styleIndex == 1)
  {
    hPlot -> SetFillColor(46);
    hPlot -> SetFillStyle(3004);
  }
  if(styleIndex == 2)
  {
    hPlot -> SetFillColor(38);
    hPlot -> SetFillStyle(3005);
  }
  if(styleIndex == 3)
  {
    hPlot -> SetFillColor(29);
    hPlot -> SetFillStyle(3001);
  }

  hPlot -> Draw(command);
  C -> Update();
}

//Get the conversion from EQ2Etrue
vector < vector < vector <double> > > GetEQ2EtrueParams(string geometry)
{
  ifstream infile;
  if (geometry=="2010") infile.open(INPUT_EQ2ETRUE_PARAMS_2010);
  else if (geometry=="2011-2012") infile.open(INPUT_EQ2ETRUE_PARAMS_2011);
  else if (geometry=="2012-2013") infile.open(INPUT_EQ2ETRUE_PARAMS_2012);
  else {
    cout << "Bad geometry passed to getEQ2EtrueParams\n";
    exit(0);
  }
  //a check to make sure the file is open
  if(!infile.is_open())
    cout << "Problem opening " << INPUT_EQ2ETRUE_PARAMS_2012 << endl;

  vector < vector < vector < double > > > params;
  params.resize(2,vector < vector < double > > (3, vector < double > (6,0.)));

  char holdType[10];
  int side=0, type=0;
  while (infile >> holdType >> params[side][type][0] >> params[side][type][1] >> params[side][type][2] >> params[side][type][3] >> params[side][type][4] >> params[side][type][5])
  {
    std::cout << holdType << " " << params[side][type][0] << " " << params[side][type][1]
	      << " " << params[side][type][2] << " " << params[side][type][3]
	      << " " << params[side][type][4] << " " << params[side][type][5] << std::endl;
    type+=1;
    if (type==3) {type=0; side=1;}
  }
  return params;
}

double CalculateErecon(double totalEvis, vector < vector < vector <double> > > tempEQ2Etrue, int type, int side)
{
  return tempEQ2Etrue[side][type][0]
	+tempEQ2Etrue[side][type][1]*totalEvis
	+tempEQ2Etrue[side][type][2]/(totalEvis+tempEQ2Etrue[side][type][3])
	+tempEQ2Etrue[side][type][4]/((totalEvis+tempEQ2Etrue[side][type][5])*(totalEvis+tempEQ2Etrue[side][type][5]));;
}

double ProbTwiddleValidity(vector <double> convertedTwiddle, vector <double> energyAxis)
{

  if(convertedTwiddle.size() != energyAxis.size())
  {
    cout << "ERROR. Mismatched vector sizes for input graph." << endl;
    return 0;
  }

  // since both vectors same size, we need to first find index that corresponds to source energies
  int Ce_index = 0;
  int Sn_index = 0;
  int Bi1_index = 0;
  int Bi2_index = 0;
  for(unsigned int i = 0; i < energyAxis.size(); i++)
  {
    if(energyAxis[i] > 129.5 && energyAxis[i] < 130.5)
    {
      Ce_index = i;
    }
    if(energyAxis[i] > 367.5 && energyAxis[i] < 368.5)
    {
      Sn_index = i;
    }
    if(energyAxis[i] > 497.5 && energyAxis[i] < 498.5)
    {
      Bi1_index = i;
    }
    if(energyAxis[i] > 993.5 && energyAxis[i] < 994.5)
    {
      Bi2_index = i;
    }
  }
/*
  cout << "Ce source error: " << convertedTwiddle[Ce_index] << endl;
  cout << "Sn source error: " << convertedTwiddle[Sn_index] << endl;
  cout << "Bi1 source error: " << convertedTwiddle[Bi1_index] << endl;
  cout << "Bi2 source error: " << convertedTwiddle[Bi2_index] << endl;

  cout << "Error envelope at Ce source: " << errEnv2011_top_1sigma->Eval(130.1) << endl;
  cout << "Error envelope at Sn source: " << errEnv2011_top_1sigma->Eval(368.5) << endl;
  cout << "Error envelope at Bi1 source: " << errEnv2011_top_1sigma->Eval(498) << endl;
  cout << "Error envelope at Bi2 source: " << errEnv2011_top_1sigma->Eval(993.8) << endl;
*/
  // need absolute values since we're just counting error bars
  double totalErrorBars = ( abs(convertedTwiddle[Ce_index]) / errEnv2011_top_1sigma->Eval(130.3)
			+ abs(convertedTwiddle[Sn_index]) / errEnv2011_top_1sigma->Eval(368.5)
                        + abs(convertedTwiddle[Bi1_index]) / errEnv2011_top_1sigma->Eval(498)
                        + abs(convertedTwiddle[Bi2_index]) / errEnv2011_top_1sigma->Eval(993.8)
			) / 4.0;

  TF1* gaussian = new TF1("gaus", "TMath::Gaus(x, 0, 1, 1)", -10, 10);

  // because we are doing absolute values, we only care about the half-Gaussian
  // so we want to take the sigma to the endpoint (times 2) as probability of acceptance
  return 1;
//  return 2.0*gaussian->Integral(totalErrorBars, 10);
}

TF1* ErrorEnvelope_2010(double factor)
{
  TF1* fEnv = new TF1("2010_error_envelope", Form("%f*((x <= 200)*2.5 + (x > 200 && x <= 500)*(2.5 + 0.0125*(x-200)) + (x>500)*6.25)", factor), 0, 1050);

  return fEnv;
}

void ErrorEnvelope_2011()
{
  errEnv2011_top_1sigma = new TF1("2011-2012_error_envelope_top1sigma", converter1top2011, 0, 1050, 0);

  errEnv2011_top_2sigma = new TF1("2011-2012_error_envelope_top2sigma", converter2top2011, 0, 1050, 0);

  errEnv2011_bot_1sigma = new TF1("2011-2012_error_envelope_bot1sigma", converter1bot2011, 0, 1050, 0);

  errEnv2011_bot_2sigma = new TF1("2011-2012_error_envelope_bot2sigma", converter2bot2011, 0, 1050, 0);
}

void ErrorEnvelope_2012()
{
  errEnv2012_top_1sigma = new TF1("2012-2013_error_envelope_top1sigma", converter1top2012, 0, 1050, 0);

  errEnv2012_top_2sigma = new TF1("2012-2013_error_envelope_top2sigma", converter2top2012, 0, 1050, 0);

  errEnv2012_bot_1sigma = new TF1("2012-2013_error_envelope_bot1sigma", converter1bot2012, 0, 1050, 0);

  errEnv2012_bot_2sigma = new TF1("2012-2013_error_envelope_bot2sigma", converter2bot2012, 0, 1050, 0);
}

double converter1top2011(double *x, double *par)
{
  double energy = x[0];
  return 1*hEnvelope2011->GetBinContent(hEnvelope2011->FindBin(energy));
}
double converter2top2011(double *x, double *par)
{
  double energy = x[0];
  return 2*hEnvelope2011->GetBinContent(hEnvelope2011->FindBin(energy));
}
double converter1bot2011(double *x, double *par)
{
  double energy = x[0];
  return (-1)*hEnvelope2011->GetBinContent(hEnvelope2011->FindBin(energy));
}
double converter2bot2011(double *x, double *par)
{
  double energy = x[0];
  return (-2)*hEnvelope2011->GetBinContent(hEnvelope2011->FindBin(energy));
}

double converter1top2012(double *x, double *par)
{
  double energy = x[0];
  return 1*hEnvelope2012->GetBinContent(hEnvelope2012->FindBin(energy));
}
double converter2top2012(double *x, double *par)
{
  double energy = x[0];
  return 2*hEnvelope2012->GetBinContent(hEnvelope2012->FindBin(energy));
}
double converter1bot2012(double *x, double *par)
{
  double energy = x[0];
  return (-1)*hEnvelope2012->GetBinContent(hEnvelope2012->FindBin(energy));
}
double converter2bot2012(double *x, double *par)
{
  double energy = x[0];
  return (-2)*hEnvelope2012->GetBinContent(hEnvelope2012->FindBin(energy));
}

void LoadEnvelopeHistogram_2011()
{
  TString fileName = "/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/Error_Envelope/MB_errorEnvelopes_publication/MB_ErrorEnvelope_2011-2012.txt";

  int binNum = 0;
  double energy = 0;
  double error = 0;
  double errorHigh = 0;
  double errorLow = 0;

  // opens the file named above
  string buf1;
  ifstream infile1;
  cout << "The file being opened is: " << fileName << endl;
  infile1.open(fileName);

  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << fileName << endl;

  while(true)
  {
    getline(infile1, buf1);
    istringstream bufstream1(buf1);

    if(!infile1.eof())
    {
      bufstream1 >> binNum
                 >> energy
                 >> error
                 >> errorHigh
                 >> errorLow;

      hEnvelope2011->SetBinContent(hEnvelope2011->FindBin(energy), errorHigh);
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Done filling in data from " << fileName.Data() << endl;

}

void LoadEnvelopeHistogram_2012()
{
  TString fileName = "/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/Error_Envelope/MB_errorEnvelopes_publication/MB_ErrorEnvelope_2012-2013.txt";

  int binNum = 0;
  double energy = 0;
  double error = 0;
  double errorHigh = 0;
  double errorLow = 0;

  // opens the file named above
  string buf1;
  ifstream infile1;
  cout << "The file being opened is: " << fileName << endl;
  infile1.open(fileName);

  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << fileName << endl;

  while(true)
  {
    getline(infile1, buf1);
    istringstream bufstream1(buf1);

    if(!infile1.eof())
    {
      bufstream1 >> binNum
                 >> energy
                 >> error
                 >> errorHigh
                 >> errorLow;

      hEnvelope2012->SetBinContent(hEnvelope2012->FindBin(energy), errorHigh);
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Done filling in data from " << fileName.Data() << endl;

}

