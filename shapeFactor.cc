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

#define		TYPE	"allTypes"		// allTypes, type0, type1 are acceptable
#define		GEOM	"2011-2012"	// 2011-2012, 2012-2013
#define		FITMINBIN	17
#define		FITMAXBIN	65

// Used for visualization, keeps the graph on screen.
//TApplication plot_program("FADC_readin",0,0,0,0);

struct entry
{
  int octetNum;
  double avg_mE;
  double chisquared_readin;
  double ndf_readin;
  double chisquaredperdf_readin;
  double b;
  double bErr;
  int fitStatus;
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

  TFile fData(Form("/home/xuansun/Documents/Analysis_Code/FierzWork_2011-2013/ExtractedHistograms/Data_Hists/Octet_%i_ssDataHist_%s.root", octNb, TYPE));
  TFile fMCSM(Form("/mnt/Data/xuansun/BLIND_MC_files/%s_geom/BLIND_MC_A_0_b_0_Octet_%i_ssHist_%s.root", GEOM, octNb, TYPE));

  TH1D *hData = new TH1D("Data", "Data", 100, 0, 1000);
  TH1D *hMCSM = new TH1D("MC", "MC", 100, 0, 1000);

  // note: ->Sumw2() not needed because I saved histograms with ->Sumw2() in these .root files I'm opening.
  hData = (TH1D*)fData.Get("Super sum");
  hMCSM = (TH1D*)fMCSM.Get("Super sum");

  double hDataNorm = 0;
  double hMCNorm = 0;

  for(int i = FITMINBIN; i <= FITMAXBIN; i++)
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
  g->SetTitle(Form("Shape Factor for Octet %i, %s", octNb, TYPE));
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
  entry correctOctetEntry;

  TString fileName = Form("BLIND_TMinuitbValues_%s_%s_Bins_%i-%i.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN);

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
      bufstream1 >> evt.octetNum
                >> evt.avg_mE
                >> evt.chisquared_readin
                >> evt.ndf_readin
		>> evt.chisquaredperdf_readin
		>> evt.b
                >> evt.bErr
                >> evt.fitStatus;
      if(evt.octetNum == octNb)
      {
        fierzVal = evt.b;
	mOverE = evt.avg_mE;
	correctOctetEntry.octetNum = evt.octetNum;
	correctOctetEntry.avg_mE = evt.avg_mE;
	correctOctetEntry.chisquared_readin = evt.chisquared_readin;
	correctOctetEntry.ndf_readin = evt.ndf_readin;
	correctOctetEntry.chisquaredperdf_readin = evt.chisquaredperdf_readin;
	correctOctetEntry.b = evt.b;
	correctOctetEntry.bErr = evt.bErr;
	correctOctetEntry.fitStatus = evt.fitStatus;
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
  for(int i = FITMINBIN; i <= FITMAXBIN; i++)
  {
    sFierz = ( fierzVal*( (511.0 / (511.0 + energy[i])) - mOverE ) ) / (1 + fierzVal*mOverE);
    Sfactor.push_back(sFierz);
    restrictedEnergy.push_back(energy[i]);
    restrictedNumPoints++;
  }

  int energyMin = hData->GetBinCenter(FITMINBIN);
  int energyMax = hData->GetBinCenter(FITMAXBIN);

  TGraph *gb = new TGraph(restrictedNumPoints, &(restrictedEnergy[0]), &(Sfactor[0]));
  gb -> SetLineWidth(3.0);
  gb -> SetLineColor(30);
  gb -> Draw("LSAME");

  TLine *yLow = new TLine(energyMin, -0.1, energyMin, 0.1);
  yLow -> SetLineStyle(9);
  TLine *yHigh = new TLine(energyMax, -0.1, energyMax, 0.1);
  yHigh -> SetLineStyle(9);
  yLow->Draw("SAME");
  yHigh->Draw("SAME");

  TLatex t;
  t.SetTextSize(0.03);
  t.SetTextAlign(13);
  t.DrawLatex(900, 0.09, Form("b_{Minuit} = %f", fierzVal));

  // calculate the chi-squared by hand from the theory with TMinuit fit b values to shape factor
  double chisquared_minuitShape = 0;
  double ndf_minuitShape = 0;
  for(unsigned int i = 0; i < Sfactor.size(); i++)
  {
    chisquared_minuitShape = chisquared_minuitShape + ((shape[i+10]-Sfactor[i])*(shape[i+10]-Sfactor[i]))/(shapeErr[i+10]*shapeErr[i+10]);
    ndf_minuitShape = ndf_minuitShape + 1;
  }
  ndf_minuitShape = ndf_minuitShape - 1;	// because 1 degree of freedom used in the 'b' fit.

  TLatex t2;
  t2.SetTextSize(0.03);
  t2.SetTextAlign(13);
  t2.DrawLatex(900, 0.08, Form("#frac{#chi^{2}_{MinuitShape}}{n} = %f", chisquared_minuitShape/ndf_minuitShape));

  // Now we want to create a "shape factor" function where b is a variable
  // and fit our existing shape factor points.
  TF1 *fShapeFit = new TF1("Shape factor b fit", Form("([0]*( (511.0 / (511.0 + x)) - %f )) / (1 + [0]*%f)", mOverE, mOverE)
				, energyMin, energyMax);
  g->Fit(fShapeFit, "R", "", energyMin, energyMax);
  TF1 *fitFunc = g->GetFunction("Shape factor b fit");
  double bFit = fitFunc->GetParameter(0);
  double bFitErr = fitFunc->GetParError(0);

  TLatex t3;
  t3.SetTextSize(0.03);
  t3.SetTextAlign(13);
  t3.DrawLatex(850, 0.06, Form("b_{HF} = %f", bFit));

  TLatex t4;
  t4.SetTextSize(0.03);
  t4.SetTextAlign(13);
  t4.DrawLatex(850, 0.05, Form("bErr_{HF} = %f", bFitErr));

  // get a second measure of chi-squared using the TF1 hist fit.
  // note: we use the restrictedEnergy vector above because that corresponds to actual x values.
  // so there is some "cross-over" between these "sections" of code.
  double chisquared_hf = 0;
  double ndf_hf = 0;
  for(unsigned int i = 0; i < restrictedEnergy.size(); i++)
  {
    chisquared_hf = chisquared_hf + ((shape[i+10]-fitFunc->Eval(restrictedEnergy[i]))*(shape[i+10]-fitFunc->Eval(restrictedEnergy[i])))
		/ (shapeErr[i+10]*shapeErr[i+10]);
    ndf_hf = ndf_hf + 1;
  }
  ndf_hf = ndf_hf - 1;	// again, because 1 dof used in 'b' param fit.

  TLatex t5;
  t5.SetTextSize(0.03);
  t5.SetTextAlign(13);
  t5.DrawLatex(900, 0.04, Form("#frac{#chi^{2}_{HF}}{n} = %f", chisquared_hf/ndf_hf));

  // print both chi-squareds (TFF and HF) to file so we can plot it in other code.
  // and the new fit values using TH1::Fit()
  ofstream outfile;
  TString chisquaredFileName = Form("Minuit_HF_chisquared_%s_%s_Bins_%i-%i_shapeFactor.txt", TYPE, GEOM, FITMINBIN, FITMAXBIN);
  outfile.open(chisquaredFileName.Data(), ios::app);
  outfile << octNb << "\t"
	  << correctOctetEntry.fitStatus << "\t"
          << chisquared_minuitShape << "\t"
          << ndf_minuitShape << "\t"
          << chisquared_minuitShape/ndf_minuitShape << "\t"
	  << chisquared_hf << "\t"
	  << ndf_hf << "\t"
	  << chisquared_hf/ndf_hf << "\t"
	  << correctOctetEntry.chisquared_readin << "\t"
	  << correctOctetEntry.ndf_readin << "\t"
	  << correctOctetEntry.chisquaredperdf_readin << "\t"
	  << correctOctetEntry.b << "\t"
	  << correctOctetEntry.bErr << "\t"
	  << bFit << "\t"
	  << bFitErr << "\n";
  outfile.close();

  C->Print(Form("BLIND_ShapeFactor_%03i_%s_%i-%iKeV.pdf", octNb, TYPE, energyMin, energyMax));

  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

  return 0;
}

