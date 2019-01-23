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
#define		PARAM_FILE_NAME		"test_errEnv2.txt"
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
		      int sideIndex);

// probability weighted on twiddle based on error envelope at source energies
void ProbTwiddleValidity(vector <double> convertedTwiddle, vector <double> energyAxis, TRandom3 *randNum);
bool TwiddleGoodOrBad;
double wCe;
double wSn;
double wEnd;
double wBi1;
double wBi2;
int printIndex = 0;

// Prints the twiddles to file so we don't need to store massive vectors
bool PrintTwiddlesToFile(double a, double b, double c, double d);

// fits our Erecon cross-section histograms
void FitHistogram(TH1D* h);

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

// Testing histogram for plotting stuff of interest.
vector <TH1D*> histErecon;

int GlobalTwiddleCounter = 0;

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

  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (w1) (w2) (w3) (w4)" << endl;
    return 0;
  }
/*
  wCe = atof(argv[1]);
  wSn = atof(argv[2]);
  wBi1 = atof(argv[3]);
  wBi2 = atof(argv[4]);
*/
  wCe = 6.0;
  wSn = 4.0;
  wBi1 = 2.5;
  wEnd = 2.0;
  wBi2 = 1.0;

  // Ensures the seed is different for randomizing in ROOT.
  TRandom3* engine = new TRandom3(0);
  gRandom->SetSeed(0);

  // Start the plotting stuff so we can loop and use "SAME" as much as possible.
  TCanvas *C = new TCanvas("canvas", "canvas");
  C->Divide(3, 2);
  C->cd(1);
  gROOT->SetStyle("Plain");

  LoadEnvelopeHistogram_2011();
  ErrorEnvelope_2011();

  errEnv2011_top_2sigma -> GetYaxis() -> SetRangeUser(-20, 20);
  errEnv2011_top_2sigma -> GetYaxis() -> SetTitle("E_{recon} Error (keV)");
  errEnv2011_top_2sigma -> GetXaxis() -> SetTitle("E_{recon} (keV)");
  errEnv2011_top_2sigma -> SetTitle("Non-linearity Polynomial Variations, 2011-2012");
  errEnv2011_top_2sigma -> SetLineStyle(2);
  errEnv2011_top_2sigma -> Draw();

  // Create histograms at fixed Erecon values to look at distribution of polynomials.
  histErecon.push_back(new TH1D("Ce2011-2012", "Erecon = 130, #sigma_{env} = 3.1", 120, -30, 30));
  histErecon.push_back(new TH1D("Sn2011-2012", "Erecon = 368, #sigma_{env} = 2.1", 120, -30, 30));
  histErecon.push_back(new TH1D("BiLow2011-2012", "Erecon = 498, #sigma_{env} = 2.76", 120, -30, 30));
  histErecon.push_back(new TH1D("endpoint", "Erecon = 782, #sigma_{env} = 3.8", 120, -30, 30));
  histErecon.push_back(new TH1D("BiHigh2011-2012", "Erecon = 994, #sigma_{env} = 7.6", 120, -30, 30));

  // Load the converter to get Erecon from a single EQ value.
  cout << "Using following calibration for 2011-2012 geometry to convert Evis to Erecon..." << endl;
  vector < vector < vector <double> > > converter = GetEQ2EtrueParams("2011-2012");

  int counter, numberSaved;
  counter = 0;
  numberSaved = 0;

  // outer loop, j, is the side index.
  for(int j = 0; j <= 1; j++)
  {
    for(double a = -20.0; a <= 20.0; a = a + 0.5)
    {
      for(double b = -0.1; b <= 0.1; b = b + 1e-3)
      {
        for(double c = -1e-4; c <= 1e-4; c = c + 2e-5)
        {
	  for(double d = 0; d <= 0; d++)
          {
            bool save = PerformVariation(a, b, c, d, numberSaved, converter, engine, j);

            // A couple of counters and print-out statements to follow along
	    if(save == true)
	    {
	      numberSaved++;
	    }
            if(counter % 1000 == 0)
            {
     	      cout << "First pass on coefficients. Checking thrown polynomial number... " << counter << endl;
	    }
            counter++;
          }
        }
      }
    }
  }

  cout << "Number of twiddles saved should be = " << GlobalTwiddleCounter << endl;

  // Placed here so 1 sigma error envelope goes on top.
  errEnv2011_top_1sigma -> SetLineStyle(2);
  errEnv2011_top_1sigma -> Draw("SAME");
  errEnv2011_bot_1sigma -> SetLineStyle(2);
  errEnv2011_bot_1sigma -> Draw("SAME");
  errEnv2011_bot_2sigma -> SetLineStyle(2);
  errEnv2011_bot_2sigma -> Draw("SAME");

  errEnv2011_top_2sigma -> Draw("SAME");

  TLine *line = new TLine(0, 0, 1000, 0);
  line->Draw("SAME");

  // Plot all the additional Erecon slice histograms

  for(unsigned int i = 0; i < histErecon.size(); i++)
  {
    C->cd(i+2);
    printIndex = i + 1;
    histErecon[i]->Draw();
    FitHistogram(histErecon[i]);
  }
  printIndex = 0;

  // Save our plot and print it out as a pdf.
  C -> Print("output_errEnv2.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

bool PerformVariation(double a, double b, double c, double d, int numPassed,
		      vector < vector < vector <double> > > EQ2Etrue, TRandom3* factor,
		      int sideIndex)
{
  bool saveCondition = true;

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

  // Get our error envelope so we can check polynomial values against (multiples of) them.
  TF1* errEnv1 = errEnv2011_top_1sigma;
  TF1* errEnv2 = errEnv2011_top_2sigma;

  // Check our polynomial (the scatter plot) against a save condition.
  double x, y;
  // These are to save values for histogram plots later.
  double v1 = -10;
  double v2 = -10;
  double v3 = -10;
  double vEnd = -10;
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
    else if(x > 487.5 && x < 488.5)
    {
      v3 = y;
    }
    else if(x > 781.5 && x < 782.5)
    {
      vEnd = y;
    }
    else if(x > 993.5 && x < 994.5)
    {
      v4 = y;
    }

    // if, at any point, we are outside error envelope, exit and don't save and don't throw a number.
    if(abs(y) > 1.0*errEnv1->Eval(x))
    {
      saveCondition = false;
      break;
    }

  }

  if(saveCondition == true)
  {
    // The reason this is done as a [void + global bool] is because of a crazy memory leak.
    // Tons of debugging, can't find the cause. But this is a work-around.
    ProbTwiddleValidity(delta_Erecon_values, Evis_axis, factor);

    if(TwiddleGoodOrBad == true)
    {
      GlobalTwiddleCounter++;
      PrintTwiddlesToFile(a, b, c, d);

      histErecon[0] -> Fill(v1);
      histErecon[1] -> Fill(v2);
      histErecon[2] -> Fill(v3);
      histErecon[3] -> Fill(vEnd);
      histErecon[4] -> Fill(v4);

      // Plotting stuff
      graph->SetLineColor(numPassed % 50);
      graph->Draw("SAME");
//      delete graph;
    }
  }
  else if(saveCondition == false)
  {
    delete graph;
  }
  delete pure_Evis;
  delete twiddle_Evis;
  delete Erecon0_East;
  delete Erecon_Twiddle_East;
  delete delta_Erecon_East;

  Evis_axis.clear();
  Erecon0_values.clear();
  delta_Erecon_values.clear();

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

void ProbTwiddleValidity(vector <double> convertedTwiddle, vector <double> energyAxis, TRandom3 *randNum)
{

  if(convertedTwiddle.size() != energyAxis.size())
  {
    cout << "ERROR. Mismatched vector sizes for input graph." << endl;
  }

  // since both vectors same size, we need to first find index that corresponds to source energies
  int Ce_index = 0;
  int Sn_index = 0;
  int End_index = 0;
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
    if(energyAxis[i] > 781.5 && energyAxis[i] < 782.5)
    {
      End_index = i;
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

  double CeErrorBar = abs(convertedTwiddle[Ce_index]) / errEnv2011_top_1sigma->Eval(130.3);
  double SnErrorBar = abs(convertedTwiddle[Sn_index]) / errEnv2011_top_1sigma->Eval(368.5);
  double Bi1ErrorBar = abs(convertedTwiddle[Bi1_index]) / errEnv2011_top_1sigma->Eval(498.0);
  double EndErrorBar = abs(convertedTwiddle[End_index]) / errEnv2011_top_1sigma->Eval(782);
  double Bi2ErrorBar = abs(convertedTwiddle[Bi2_index]) / errEnv2011_top_1sigma->Eval(993.8);

  // weights for 2011-2012 error bars that work best
//  double totalErrorBars = (2.8*CeErrorBar + 1.2*SnErrorBar + 0.8*Bi1ErrorBar + 1.6*Bi2ErrorBar) / 4.0;
//  double totalErrorBars = (wCe*CeErrorBar + wSn*SnErrorBar + wBi1*Bi1ErrorBar + wEnd*EndErrorBar + wBi2*Bi2ErrorBar) / 5.0;


  TF1* sampleGaussian = new TF1("sampleGaussian", "TMath::Gaus(x, 0, 1, 1)", -10, 10);

  // because we are doing absolute values, we only care about the half-Gaussian
  // so we want to take the sigma to the endpoint (times 2) as probability of acceptance

  double sampleGaussianValue = 2.0*sampleGaussian->Integral(wCe*CeErrorBar, 10)
			     * 2.0*sampleGaussian->Integral(wSn*SnErrorBar, 10)
			     * 2.0*sampleGaussian->Integral(wBi1*Bi1ErrorBar, 10)
			     * 2.0*sampleGaussian->Integral(wEnd*EndErrorBar, 10)
			     * 2.0*sampleGaussian->Integral(wBi2*Bi2ErrorBar, 10);

//  double sampleGaussianValue = 2.0*sampleGaussian->Integral(totalErrorBars, 10);

  delete sampleGaussian;

  // We straight up need to do all this because there is a weird memory leak.
  // This is more complicated than it needs to be but DON'T TOUCH IT!
  double randNb = randNum->Rndm();
  if(randNb <= sampleGaussianValue)
  {
    TwiddleGoodOrBad = true;
  }
  else if(randNb > sampleGaussianValue)
  {
    TwiddleGoodOrBad = false;
  }
}

bool PrintTwiddlesToFile(double a, double b, double c, double d)
{
  ofstream outfile;
  outfile.open(PARAM_FILE_NAME, ios::app);

  if(abs(a) < 1e-10)
  {
    a = 0;
  }
  if(abs(b) < 1e-10)
  {
    b = 0;
  }
  if(abs(c) < 1e-10)
  {
    c = 0;
  }
  if(abs(d) < 1e-10)
  {
    d = 0;
  }

  outfile << GlobalTwiddleCounter << "\t"
	  << a << "\t"
	  << b << "\t"
	  << c << "\t"
	  << d << "\t"
          << a << "\t"
	  << b << "\t"
          << c << "\t"
	  << d << "\n";

  outfile.close();

  return true;
}

void FitHistogram(TH1D* h)
{
  h->Fit("gaus");
  TF1* fFitResults = h->GetFunction("gaus");

  TLatex t2;
  t2.SetTextSize(0.03);
  t2.SetTextAlign(13);
  t2.DrawLatex(0.5*(h->GetNbinsX()/2.0)*(h->GetBinWidth(5)), 0.7*(h->GetMaximum()), Form("#mu = %f", fFitResults->GetParameter(1)));
  TLatex t3;
  t3.SetTextSize(0.03);
  t3.SetTextAlign(13);
  t3.DrawLatex(0.5*(h->GetNbinsX()/2.0)*(h->GetBinWidth(5)), 0.6*(h->GetMaximum()), Form("#mu_{err} = %f", fFitResults->GetParError(1)));

  TLatex t4;
  t4.SetTextSize(0.03);
  t4.SetTextAlign(13);
  t4.DrawLatex(0.5*(h->GetNbinsX()/2.0)*(h->GetBinWidth(5)), 0.5*(h->GetMaximum()), Form("#sigma = %f", fFitResults->GetParameter(2)));
  TLatex t5;
  t5.SetTextSize(0.03);
  t5.SetTextAlign(13);
  t5.DrawLatex(0.5*(h->GetNbinsX()/2.0)*(h->GetBinWidth(5)), 0.4*(h->GetMaximum()), Form("#sigma_{err} = %f", fFitResults->GetParError(2)));

/*
  ofstream outfile;
  outfile.open("ErrorBarWeightingCoeff_Scan_2012-2013_pass2.txt", ios::app);

  if(printIndex == 1)
  {
  outfile << wCe << "\t"
          << wSn << "\t"
          << wBi1 << "\t"
          << wBi2 << "\t"
          << fFitResults->GetParameter(2) << "\t";
  }
  else if(printIndex > 1 && printIndex < 5)
  {
    outfile << fFitResults->GetParameter(2) << "\t";
  }
  else if(printIndex == 5)
  {
    outfile << fFitResults->GetParameter(2) << "\n";
  }
  outfile.close();
*/
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

      // creating a "flat" error envelope of 20keV. Further work with weighting.
      hEnvelope2011->SetBinContent(hEnvelope2011->FindBin(energy), 20 /*errorHigh*/);
      cout << "At energy " << energy << ", we have symmetrized 2011-2012 error envelope: " << errorHigh << endl;

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
/*
      if(energy >= 900 && energy < 925)
      {
        hEnvelope2012->SetBinContent(hEnvelope2012->FindBin(energy), 0.95*errorHigh);
      }
      else if(energy >= 925 && energy < 950)
      {
        hEnvelope2012->SetBinContent(hEnvelope2012->FindBin(energy), 0.9*errorHigh);
      }
      else if(energy >= 950 && energy < 975)
      {
        hEnvelope2012->SetBinContent(hEnvelope2012->FindBin(energy), 0.85*errorHigh);
      }
      else if(energy >= 975 && energy < 1000)
      {
        hEnvelope2012->SetBinContent(hEnvelope2012->FindBin(energy), 0.8*errorHigh);
      }
      else if(energy >= 1000 && energy < 1025)
      {
        hEnvelope2012->SetBinContent(hEnvelope2012->FindBin(energy), 0.75*errorHigh);
      }
      else if(energy >= 1025 && energy < 1050)
      {
        hEnvelope2012->SetBinContent(hEnvelope2012->FindBin(energy), 0.7*errorHigh);
      }
      else if(energy >= 1050 && energy < 1075)
      {
        hEnvelope2012->SetBinContent(hEnvelope2012->FindBin(energy), 0.65*errorHigh);
      }
      else if(energy >= 1075 && energy < 1100)
      {
        hEnvelope2012->SetBinContent(hEnvelope2012->FindBin(energy), 0.6*errorHigh);
      }
*/
//      else
//      {
        hEnvelope2012->SetBinContent(hEnvelope2012->FindBin(energy), errorHigh);
//      }


//      hEnvelope2012->SetBinContent(hEnvelope2012->FindBin(energy), errorHigh);

      cout << "At energy " << energy << ", we have symmetrized 2012-2013 error envelope: " << errorHigh << endl;

    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Done filling in data from " << fileName.Data() << endl;

}


