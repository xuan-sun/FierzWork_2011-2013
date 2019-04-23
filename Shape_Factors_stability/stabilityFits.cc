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
#include         <TLegend.h>
#include 	 <TFitResult.h>
#include	 <TMatrixD.h>

#define		 NAMEFLAG	"endpointCorrection"

using            namespace std;

//required later for plot_program
//TApplication plot_program("FADC_readin",0,0,0,0);

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString xAxis, TString yAxis, TString command);
void FillArrays(TString fileName, int flag);



int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  TCanvas *C = new TCanvas("canvas", "canvas");
//  C -> Divide(2,1);
  gROOT -> SetStyle("Plain");   //on my computer this sets background to white, finally!

  int octetLow = 80;
  int octetHigh = 122;

  // for normalizing our spectra in the fit region
  int fitBinMin = 17;
  int fitBinMax = 65;
  double N_data = 0;
  double N_beta = 0;
  double yAxisMin = -0.1;
  double yAxisMax = 0.1;


  // load the sum of all the data octets, with proper error propagation (if ROOT isn't broken).
  TH1D *hTotalData = new TH1D("totalData", "totalData", 120, 0, 1200);
  hTotalData->Sumw2();
  TH1D *hTotalBeta = new TH1D("totalBeta", "totalBeta", 120, 0, 1200);
  hTotalBeta->Sumw2();

  // add lines to our plot in case we want to save the plots.
  TLine *y0 = new TLine(0, 0, 1200, 0);
  TLine *xMin = new TLine(hTotalData->GetBinCenter(fitBinMin), yAxisMin, hTotalData->GetBinCenter(fitBinMin), yAxisMax);
  TLine *xMax = new TLine(hTotalData->GetBinCenter(fitBinMax), yAxisMin, hTotalData->GetBinCenter(fitBinMax), yAxisMax);

  ofstream outfile;
  outfile.open(Form("zero_crossing_shapeFactors_%s_2012-2013.txt", NAMEFLAG), ios::app);

  for(int i = octetLow; i < octetHigh; i++)
  {
    if(i == 9 || i == 59 || i == 91 || i == 93 || i == 101 || i == 107 || i == 121)
    {
      continue;
    }
//    TFile fData(Form("../PositionCuts/radialCut_0-49/Octet_%i_ssDataHist_type0_radialCut_0-49mm.root", i));
    TFile fData(Form("../PositionCuts/radialCut_0-49/Octet_%i_ssDataHist_type0_radialCut_0-49mm_endpointCorrected.root", i));
    TH1D *hTempData = (TH1D*)fData.Get("Super sum");
    hTotalData->Add(hTempData);
    fData.Close();

    // no need to load an endpoint corrected MC simply because the MC is far more stable.
    TFile fBeta(Form("../PositionCuts/radialCut_0-49/MC_A_0_b_0_Octet_%i_ssHist_type0_posCut_0-0.049000m.root", i));
    TH1D *hTempBeta = (TH1D*)fBeta.Get("Super sum");
    hTotalBeta->Add(hTempBeta);
    fBeta.Close();

    // scale MC histogram to data histogram number of events
    for(int j = fitBinMin; j <= fitBinMax; j++)
    {
      N_data = N_data + hTotalData->GetBinContent(j);
      N_beta = N_beta + hTotalBeta->GetBinContent(j);
    }

    hTotalBeta->Scale(N_data / N_beta);

    TH1D *hResidual = new TH1D("residual", "residual", 120, 0, 1200);
    hResidual->Sumw2();
    hResidual->Add(hTotalData, hTotalBeta, 1.0, -1.0);

    hResidual->Divide(hTotalBeta);

    // straight line fit to determine 0 crossing point
    TF1 *fit1 = new TF1("fit1", "[0] + [1]*x", 200, 400);
    TFitResultPtr results = hResidual->Fit(fit1, "RS");
    TMatrixD cor = results->GetCorrelationMatrix();
    TMatrixD cov = results->GetCovarianceMatrix();
    cout << "Correlation matrix: " << endl;
    cor.Print();
    cout << "Covariance matrix: " << endl;
    cov.Print();

    cout << "Covariance = " << cov[1][0] << endl;


    double a = fit1->GetParameter(0);
    double ae = fit1->GetParError(0);
    double b = fit1->GetParameter(1);
    double be = fit1->GetParError(1);
    double sigma_ab = cov[1][0];
    double f = -a/b;

    double sigma_f = abs(f)*sqrt( pow(ae/a, 2.0) + pow(be/b, 2.0) - 2*sigma_ab/(a*b) );

    outfile << i << "\t"
	    << f << "\t"
	    << sigma_f << "\n";


    hResidual->GetYaxis()->SetRangeUser(yAxisMin, yAxisMax);
    PlotHist(C, 2, 1, hResidual, Form("2012-2013, octet %i, %s, (S_{data}-S_{MC}) / S_{MC}", i, NAMEFLAG), "Reconstructed Energy (keV)", "Fractional residual", "");

    y0->Draw();
    xMin->Draw();
    xMax->Draw();

//    C->Print(Form("octet_%i.pdf", i));
    // reset all the values for next octet
    hTotalData->Reset();
    hTotalBeta->Reset();
    N_data = 0;
    N_beta = 0;
    delete hResidual;
  }

  outfile.close();


  //prints the canvas with a dynamic TString name of the name of the file
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString xAxis, TString yAxis, TString command)
{
  C->cd(canvasIndex);
  gPlot->SetTitle(title);
  gPlot->GetXaxis()->SetTitle(xAxis);
  gPlot->GetXaxis()->CenterTitle();
  gPlot->GetYaxis()->SetTitle(yAxis);
  gPlot->GetYaxis()->CenterTitle();
//  C->SetLogy();

  gPlot->SetMarkerStyle(21);
  gPlot->SetMarkerSize(0.75);
  gPlot->SetMarkerColor(styleIndex);
  gPlot->SetLineColor(styleIndex);

  gPlot->Draw(command);

  C->Update();

}

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xAxis, TString yAxis, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot->GetXaxis()->SetTitle(xAxis);
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle(yAxis);
  hPlot -> GetYaxis() -> CenterTitle();

  hPlot->SetFillColor(styleIndex);
  hPlot->SetFillStyle(styleIndex);
  hPlot->SetLineColor(styleIndex);

  hPlot->Draw(command);

  C->Update();
}

