#include	 "BetaSpectrum.hh"

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
#include 	 <TTree.h>
#include	 <TChain.h>
#include	 <TVector.h>
#include	 <vector>
#include	 <utility>
#include	 <TLeaf.h>
#include	 <TLine.h>
#include	 <TLatex.h>
using		 namespace std;

struct Event
{
  double primKE;
  double Erecon;
  int side;
  int type;
  int eventNum;
  double time[2];
  double tE;
  double tW;
  int pid;
  double xEastPos;
  double yEastPos;
  double xWestPos;
  double yWestPos;
};

// forward declarations for useful functions
double CalculatebFromPercentageMixing(TString fileName);
double CalculateAveragemOverE(TH1D* gammaSM, int binMin, int binMax);
TH1D* LoadMBAsymmetry(TString fileName);
TH1D* BlindAsymmetry(TH1D *unblindA, double b_forBlinding, double avg_mE_forBlinding);

// plotting functions
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xTitle, TString yTitle, TString command);

// these are actual beta run indices
const int index_A2 = 0;
const int index_A5 = 1;
const int index_A7 = 2;
const int index_A10 = 3;
const int index_B2 = 4;
const int index_B5 = 5;
const int index_B7 = 6;
const int index_B10 = 7;
double avg_mE = 0;

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

//-------------------------------------------------//
//------------ Start of Program -------------------//
//-------------------------------------------------//

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  // creating canvas for plotting
  TCanvas *C = new TCanvas("canvas", "canvas", 800, 400);

  int octNb = 43;
  double xMin = 170;
  double xMax = 660;

  TFile fMC0(TString::Format("ExtractedHistograms/MC_A_0_b_0/MC_A_0_b_0_Octet_%i_ssHist_%s.root", octNb, "type0"));
  TH1D* mcTheoryHistBeta = (TH1D*)fMC0.Get("Super sum");
  avg_mE = CalculateAveragemOverE(mcTheoryHistBeta, mcTheoryHistBeta->FindBin(xMin), mcTheoryHistBeta->FindBin(xMax));

  cout << "Average m/E for octet " << octNb << " is equal to " << avg_mE
 	<< " over a fit range of " << xMin << " to " << xMax << endl;

  double bMixing = CalculatebFromPercentageMixing("ExtractedHistograms/randomMixingSeeds.txt");

  TString asymmFile = Form("MB_asymmetries/AsymmFilesFromMB/AllCorr_OctetAsymmetries_AnaChD_Octets_0-59_BinByBin.txt");

  TH1D *asymm = LoadMBAsymmetry(asymmFile);

  TH1D *blindedAsymm = BlindAsymmetry(asymm, bMixing, avg_mE);

  TF1 *fit = new TF1("beta fit", Form("( [0]*(1.0 + [1]*(%f)) ) / (1.0 + [1]*(%f)/(x + %f))", avg_mE, m_e, m_e), xMin, xMax);

  fit->SetParName(0, "asymm");
  fit->SetParName(1, "b");
  blindedAsymm->Fit("beta fit", "R");
  TF1 *fitResults = blindedAsymm->GetFunction("beta fit");
  cout << "Chi squared value is " << fitResults->GetChisquare()
	<< " with ndf of " << fitResults->GetNDF()
	<< ". For a final chisquared/ndf = " << fitResults->GetChisquare() / fitResults->GetNDF() << endl;


  PlotHist(C, 1, 1, blindedAsymm, asymmFile, "Reconstructed Energy (keV)", "Blinded Fully Corrected A(E)", "");


  TLine *yLow = new TLine(xMin, gPad->GetUymin(), xMin, gPad->GetUymax());
  TLine *yHigh = new TLine(xMax, gPad->GetUymin(), xMax, gPad->GetUymax());
  yLow->SetLineWidth(2);
  yHigh->SetLineWidth(2);
  yLow->SetLineColor(1);
  yHigh->SetLineColor(1);
  yLow->Draw("SAME");
  yHigh->Draw("SAME");

  TLatex t2;
  t2.SetTextSize(0.03);
  t2.SetTextAlign(13);
  t2.DrawLatex(1000, 0.14, Form("A_{fit} = %f", fitResults->GetParameter(0)));
  TLatex t3;
  t3.SetTextSize(0.03);
  t3.SetTextAlign(13);
  t3.DrawLatex(1000, 0.13, Form("AErr_{fit} = %f", fitResults->GetParError(0)));
  TLatex t4;
  t4.SetTextSize(0.03);
  t4.SetTextAlign(13);
  t4.DrawLatex(1000, 0.12, Form("b_{fit} = %f", fitResults->GetParameter(1)));
  TLatex t5;
  t5.SetTextSize(0.03);
  t5.SetTextAlign(13);
  t5.DrawLatex(1000, 0.11, Form("bErr_{fit} = %f", fitResults->GetParError(1)));
  TLatex t6;
  t6.SetTextSize(0.03);
  t6.SetTextAlign(13);
  t6.DrawLatex(900, 0.10, Form("#frac{#Chi^{2}}{ndf} = #frac{%f}{%i} = %f",
				fitResults->GetChisquare(), fitResults->GetNDF(), fitResults->GetChisquare() / fitResults->GetNDF()));




  // Save our plot and print it out as a pdf.
  C -> Print("ReBLINDed_b_fit_fromAsymmData.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

TH1D* LoadMBAsymmetry(TString fileName)
{
  TH1D* AofE = new TH1D("AofE", "AofE", 120, 0, 1200);

  double energy, asymm, asymmErr;

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
      bufstream1 >> energy >> asymm >> asymmErr;

      AofE->SetBinContent(AofE->FindBin(energy), (-1)*asymm);
      AofE->SetBinError(AofE->FindBin(energy), asymmErr);
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  return AofE;
}


void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xTitle, TString yTitle, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle(xTitle);
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle(yTitle);
  hPlot -> GetYaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetRangeUser(0.08, 0.16);


  if(styleIndex == 1)
  {
    hPlot->SetLineColor(46);
    hPlot -> SetFillColor(46);
    hPlot -> SetFillStyle(3004);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 2)
  {
    hPlot->SetLineColor(38);
    hPlot -> SetFillColor(38);
    hPlot -> SetFillStyle(3005);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 3)
  {
    hPlot->SetLineColor(29);
    hPlot -> SetFillColor(29);
//    hPlot -> SetFillStyle(3005);
    hPlot -> SetFillStyle(3001);
  }

  hPlot -> Draw(command);

  C -> Update();
}

double CalculatebFromPercentageMixing(TString fileName)
{
  double b = 0;

  // read in our random mixing seed so I stay pretty blind.
  double s0, s1, s2, s3, s4, s5, s6, s7, s8, s9;

  string buf1;
  ifstream infile1;
  cout << "The file being opened is: " << fileName.Data() << endl;
  infile1.open(fileName);

  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << fileName.Data() << endl;

  while(true)
  {
    getline(infile1, buf1);
    istringstream bufstream1(buf1);

    if(!infile1.eof())
    {
      bufstream1 >> s0 >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8 >> s9;
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }


  b = 0.05 / ( (1 - 0.05) * (avg_mE) );

  return b;
}

double CalculateAveragemOverE(TH1D* gammaSM, int binMin, int binMax)
{
  double num = 0;
  double denom = 0;

  for(int i = binMin; i < binMax; i++)
  {
    num = num + (m_e*gammaSM->GetBinContent(i)) / (gammaSM->GetBinCenter(i) + m_e);
    denom = denom + gammaSM->GetBinContent(i);
  }

  return num/denom;
}


TH1D* BlindAsymmetry(TH1D *unblindA, double b_forBlinding, double avg_mE_forBlinding)
{
  TH1D *blindA = new TH1D("blinded AofE", "blinded AofE", 120, 0, 1200);

  double blindingFactor = 0;

  for(int i = 0; i <= unblindA->GetNbinsX(); i++)
  {
    blindingFactor = (1.0 + b_forBlinding*avg_mE_forBlinding) / (1.0 + b_forBlinding*m_e/(unblindA->GetBinCenter(i) + m_e));
    blindA->SetBinContent(i, (unblindA->GetBinContent(i))*blindingFactor);
    blindA->SetBinError(i, (unblindA->GetBinError(i))*blindingFactor);
  }

  return blindA;
}
