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
#include         <TRandom3.h>
#include	 <TLegend.h>

#define		TYPE	"type0"
#define		GEOM	"2012-2013"
#define		FITMINBIN	17
#define		FITMAXBIN	65

using            namespace std;

// Fundamental constants that get used
const double m_e = 511.00;                      ///< electron mass, keV/c^2
double NDF = 1;				// taken from the .txt file. All twiddles fit the same range so same NDF

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

struct entry
{
  int indexNb;
  double chi2_ndf;
  double Ea;
  double Eb;
  double Ec;
  double Ed;
  double Wa;
  double Wb;
  double Wc;
  double Wd;
};

// has to go below the struct because I take, as input, the struct
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command, double maxBinContents);
void FillArrays(TString fileName, TH1D *h, int hFillOption);
void PrintChosenTwiddles(entry cevt);
void RepopulateHist(TString fileName, TH1D* hRef, TH1D* hNew);

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

//  int option = atoi(argv[1]);

  TCanvas *C = new TCanvas("canvas", "canvas");
  C->cd();
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  // load the full ch2 histogram first
  TString fPath = Form("chi2_errEnv2_asymm_2012-2013.txt");
  TH1D* h = new TH1D("twiddles", "twiddles w/ 2012-2013 calibration sources", 100, 0.1, 10);
  FillArrays(fPath, h, 1);

  // use previous histogram to resample according to theory chi2 distribution
  TH1D* h2 = new TH1D("twiddlesRe", "twiddles w/ 2012-2013 calibration sources, chi2 sampling", 100, 0.1, 10);
  RepopulateHist(fPath, h, h2);
  int max = h2->GetMaximum();
  PlotHist(C, 2, 1, h2, Form("twiddles, %s", GEOM), "chisquared/ndf", "N", "", max);

  // overlay a (normalized) plot of the theory chi2 distribution.
  TF1 *theoryChiP = new TF1("theoryPlot", Form("-1*(TMath::Prob(x*%f, %f) - TMath::Prob((x-0.1)*%f, %f))", NDF, NDF, NDF, NDF), 0.1, 10);
  TH1D *theoryChiHistP = (TH1D*)(theoryChiP->GetHistogram());
  double hTot = 0;
  double theoryHTot = 0;
  // must do minimum value of 0.1 or else the chisquared function diverges down, messing up normalization.
  for(int i = h2->FindBin(0.1); i <= h2->GetNbinsX(); i++)
  {
    hTot = hTot + h2->GetBinContent(i);
    theoryHTot = theoryHTot + theoryChiHistP->GetBinContent(i);

  }
  theoryChiHistP->Scale(hTot / theoryHTot);

  PlotHist(C, 1, 1, theoryChiHistP, "", "", "", "SAME", max);

  TLegend* leg1 = new TLegend(0.7,0.6,0.9,0.8);
  leg1->AddEntry(h,"twiddle chi2/ndf","f");
  leg1->AddEntry(theoryChiHistP,"theory #frac{#Chi^{2}}{NDF = 1}","f");
  leg1->Draw();


  //prints the canvas with a dynamic TString name of the name of the file
//  C->Print("viewNewXuanFitter_SymmetricTwiddles_finerGrid.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

void FillArrays(TString fileName, TH1D* h, int hFillOption)
{
  entry evt;

  //opens the file that I name in DATA_FILE_IN
  string buf;
  ifstream infile1;
  cout << "The file being opened is: " << fileName << endl;
  infile1.open(fileName);

  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << fileName << endl;

  while(true)
  {
    getline(infile1, buf);
    istringstream bufstream(buf);

    if(!infile1.eof())
    {
      bufstream >> evt.indexNb
		>> evt.chi2_ndf
		>> evt.Ea
		>> evt.Eb
		>> evt.Ec
		>> evt.Ed
		>> evt.Wa
		>> evt.Wb
		>> evt.Wc
		>> evt.Wd;

       h->Fill(evt.chi2_ndf);
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}

void RepopulateHist(TString fileName, TH1D* hRef, TH1D* hNew)
{
  TRandom3 *randNum = new TRandom3(0);
  TF1 *theoryChi = new TF1("theory", Form("-1*(TMath::Prob(x*%f, %f) - TMath::Prob((x-0.1)*%f, %f))", NDF, NDF, NDF, NDF), 0.1, 10);
  TH1D *theoryChiHist = (TH1D*)(theoryChi->GetHistogram());

  entry evt;

  //opens the file that I name in DATA_FILE_IN
  string buf;
  ifstream infile1;
  cout << "The file being opened is: " << fileName << endl;
  infile1.open(fileName);

  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << fileName << endl;

  while(true)
  {
    getline(infile1, buf);
    istringstream bufstream(buf);

    if(!infile1.eof())
    {
      bufstream >> evt.indexNb
                >> evt.chi2_ndf
                >> evt.Ea
                >> evt.Eb
                >> evt.Ec
                >> evt.Ed
                >> evt.Wa
                >> evt.Wb
                >> evt.Wc
                >> evt.Wd;

      double hMax = theoryChiHist->GetMaximum();
      int hMaxBin = theoryChiHist->GetMaximumBin();
      double randSample = randNum->Rndm() * hMax;
      double acceptProb = theoryChiHist->GetBinContent(theoryChiHist->FindBin(evt.chi2_ndf))
		        * ( hRef->GetBinContent(hMaxBin)
		        / hRef->GetBinContent(hRef->FindBin(evt.chi2_ndf)));

      if(randSample <= theoryChiHist->GetBinContent(theoryChiHist->FindBin(evt.chi2_ndf)))
      {
        hNew->Fill(evt.chi2_ndf);
        PrintChosenTwiddles(evt);
      }

    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been read again successfully." << endl;




}


void PrintChosenTwiddles(entry cevt)
{
  ofstream outfile;
  outfile.open("chosenTwiddles_chi2Checked_1NDF_2012-2013_index19.txt", ios::app);
  outfile << cevt.indexNb << "\t"
          << cevt.Ea << "\t"
          << cevt.Eb << "\t"
          << cevt.Ec << "\t"
          << cevt.Ed << "\t"
          << cevt.Wa << "\t"
          << cevt.Wb << "\t"
          << cevt.Wc << "\t"
          << cevt.Wd << "\n";
  outfile.close();
}


void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command, double maxBinContents)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot->GetXaxis()->SetTitle(xAxis);
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle(yAxis);
  hPlot -> GetYaxis() -> CenterTitle();
  hPlot->GetYaxis()->SetRangeUser(0, 1.2*maxBinContents);

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
    hPlot->SetFillColor(30);
    hPlot->SetFillStyle(3003);
  }
  if(styleIndex == 4)
  {
    hPlot->SetFillColor(41);
    hPlot->SetFillStyle(3016);
  }

  hPlot->Draw(command);

  C->Update();
}

