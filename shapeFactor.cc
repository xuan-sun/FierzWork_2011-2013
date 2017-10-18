#include         <iostream>
#include         <fstream>
#include         <TGaxis.h>
#include         <sstream>
#include         <TGraph.h>
#include         <TGraphErrors.h>
#include         <TCanvas.h>
#include         <TApplication.h>
#include         <stdlib.h>
#include         <TF1.h>
#include         <TH1.h>
#include         <TProfile.h>
#include         <TObjArray.h>
#include         <TStyle.h>
#include         <TMarker.h>
#include         <math.h>
#include         <TStyle.h>
#include         <TPaveStats.h>
#include         <TPaveText.h>
#include         <vector>
#include         <string.h>
#include         <fstream>
#include         <TROOT.h>
#include         <TFile.h>
#include         <TLegend.h>
#include         <TLegendEntry.h>
#include         <time.h>
#include         <TH2F.h>
#include         <assert.h>
#include         <string>
#include         <TRandom.h>
#include         <TTree.h>
#include         <TChain.h>
#include         <TVector.h>
#include         <vector>
#include         <utility>
#include         <TLeaf.h>
#include	 <TLatex.h>
using            namespace std;

// Used for visualization, keeps the graph on screen.
//TApplication plot_program("FADC_readin",0,0,0,0);

struct entry
{
  int octetNum;
  double b;
  double avg_mE;
  double chisquared;
  double ndf;
  double bErr;
  double GluckErr;
  int dataEntriesNumber;
  int fitEntriesNumber;
  double numberOfOriginalEvents;
};

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (octet #)" << endl;
    return 0;
  }

  int octNb = atoi(argv[1]);

  TCanvas *C = new TCanvas("canvas", "canvas");

  TFile fData(Form("/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/ExtractedHistograms/Data_Hists/Octet_%i_ssDataHist.root", octNb));
  TFile fMCSM(Form("/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/ExtractedHistograms/MC_A_0_b_0/MC_A_0_b_0_Octet_%i_ssHist.root", octNb));

  TH1D *hData = new TH1D("Data", "Data", 100, 0, 1000);
  TH1D *hMCSM = new TH1D("MC", "MC", 100, 0, 1000);

  // note: ->Sumw2() not needed because I saved histograms with ->Sumw2() in these .root files I'm opening.
  hData = (TH1D*)fData.Get("Super sum");
  hMCSM = (TH1D*)fMCSM.Get("Super sum");

  double hDataNorm = 0;
  double hMCNorm = 0;

  for(int i = 10; i <= 65; i++)
  {
    hDataNorm = hDataNorm + hData->GetBinContent(i);
    hMCNorm = hMCNorm + hMCSM->GetBinContent(i);
  }

  hData->Scale(1.0/hDataNorm);
  hMCSM->Scale(1.0/hMCNorm);

  vector <double> energy;
  vector <double> energyErr;
  vector <double> shape;
  vector <double> shapeErr;
  double estimatedErr = 0;
  double shapeValue = 0;
  int numPoints = 0;

  for(int i = 0; i < 100; i++)
  {
    energy.push_back(hData->GetXaxis()->GetBinCenter(i));
    energyErr.push_back(5);	// error bar of 5 KeV aka half the bin width
    if(hMCSM->GetBinContent(i) == 0 || hData->GetBinContent(i) == 0)
    {
      shape.push_back(0);
      shapeErr.push_back(0);
    }
    else
    {
      shapeValue = (hData->GetBinContent(i) - hMCSM->GetBinContent(i)) / hMCSM->GetBinContent(i);
      shape.push_back(shapeValue);

      estimatedErr = (hData->GetBinContent(i))/(hMCSM->GetBinContent(i))
		     * sqrt( (hData->GetBinError(i)*hData->GetBinError(i))/(hData->GetBinContent(i)*hData->GetBinContent(i))
			   + (hMCSM->GetBinError(i)*hMCSM->GetBinError(i))/(hMCSM->GetBinContent(i)*hMCSM->GetBinContent(i)) );

      shapeErr.push_back(estimatedErr);
    }
    numPoints++;
  }

  TGraphErrors *g = new TGraphErrors(numPoints, &(energy[0]), &(shape[0]), &(energyErr[0]), &(shapeErr[0]));
  g->SetTitle(Form("Shape Factor for Octet %i, Type 0", octNb));
  g->SetMarkerSize(1);
  g->SetMarkerStyle(21);
  g->SetMarkerColor(38);
  g->GetHistogram()->SetMaximum(0.1);
  g->GetHistogram()->SetMinimum(-0.1);
  g->Draw("AP");

  // extract the fitted b value and plot the shape factor with 'b' on the same graph.
  double fierzVal = -1;
  double mOverE = -1;
  entry evt;

  TString fileName = "Results_comparehist_bValues_withEvtNumbers.txt";

  //opens the file that I name in DATA_FILE_IN
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
      bufstream1 >> evt.octetNum
                >> evt.b
                >> evt.avg_mE
                >> evt.chisquared
                >> evt.ndf
                >> evt.bErr
                >> evt.GluckErr
                >> evt.dataEntriesNumber
                >> evt.fitEntriesNumber
                >> evt.numberOfOriginalEvents;
      if(evt.octetNum == octNb)
      {
        fierzVal = evt.b;
	mOverE = evt.avg_mE;
      }
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  if(fierzVal == -1 && mOverE == -1)
  {
    return 0;
  }

  vector <double> Sfactor;
  vector <double> restrictedEnergy;
  int restrictedNumPoints = 0;
  double sFierz = 0;
  for(int i = 10; i <= 65; i++)
  {
    sFierz = ( fierzVal*( (511.0 / (511.0 + energy[i])) - mOverE ) ) / (1 + fierzVal*mOverE);
    Sfactor.push_back(sFierz);
    restrictedEnergy.push_back(energy[i]);
    restrictedNumPoints++;
  }

  TGraph *gb = new TGraph(restrictedNumPoints, &(restrictedEnergy[0]), &(Sfactor[0]));
  gb -> SetLineWidth(3.0);
  gb -> SetLineColor(46);
  gb -> Draw("LSAME");

  TLine *yLow = new TLine(95, -0.1, 95, 0.1);
  yLow -> SetLineStyle(9);
  TLine *yHigh = new TLine(645, -0.1, 645, 0.1);
  yHigh -> SetLineStyle(9);
  yLow->Draw("SAME");
  yHigh->Draw("SAME");

  TLatex t;
  t.SetTextSize(0.03);
  t.SetTextAlign(13);
  t.DrawLatex(900, 0.09, Form("b = %f", fierzVal));

  // calculate the chi-squared by hand from the theory with fit b values to shape factor
  double chisquared = 0;
  double ndf = 0;
  for(unsigned int i = 0; i < Sfactor.size(); i++)
  {
    chisquared = chisquared + ((shape[i+10]-Sfactor[i])*(shape[i+10]-Sfactor[i]))/(shapeErr[i+10]*shapeErr[i+10]);
    ndf = ndf + 1;
  }

  TLatex t2;
  t2.SetTextSize(0.03);
  t2.SetTextAlign(13);
  t2.DrawLatex(900, 0.08, Form("#frac{#chi^{2}}{n} = %f", chisquared/ndf));

  C->Print(Form("ShapeFactor_%i_allTypes.pdf", octNb));

  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

