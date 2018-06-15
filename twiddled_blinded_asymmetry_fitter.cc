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
  double Erecon;
  int PID;
  int type;
  int side;
};

// forward declarations for useful functions
double CalculatebFromPercentageMixing(TString fileName);
double CalculateAveragemOverE(TH1D* gammaSM, int binMin, int binMax);
TH1D* LoadMBAsymmetry(TString fileName);
TH1D* BlindAsymmetry(TH1D *unblindA, double b_forBlinding, double avg_mE_forBlinding);
TH1D* LoadTwiddledPairAsymm(TString pathBase, int twiddle);

// plotting functions
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xTitle, TString yTitle, TString command);

double avg_mE = 0;

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

//-------------------------------------------------//
//------------ Start of Program -------------------//
//-------------------------------------------------//

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (twiddle #)" << endl;
    return 0;
  }

  // creating canvas for plotting
  TCanvas *C = new TCanvas("canvas", "canvas", 800, 400);

  int twiddleIndex = atoi(argv[1]);

  int octNb = 43;
  double xMin = 170;
  double xMax = 660;

  TFile fMC0(TString::Format("ExtractedHistograms/MC_A_0_b_0/MC_A_0_b_0_Octet_%i_ssHist_%s.root", octNb, "type0"));
  TH1D* mcTheoryHistBeta = (TH1D*)fMC0.Get("Super sum");
  avg_mE = CalculateAveragemOverE(mcTheoryHistBeta, mcTheoryHistBeta->FindBin(xMin), mcTheoryHistBeta->FindBin(xMax));

  cout << "Average m/E for octet " << octNb << " is equal to " << avg_mE
 	<< " over a fit range of " << xMin << " to " << xMax << endl;

  double bMixing = CalculatebFromPercentageMixing("ExtractedHistograms/randomMixingSeeds.txt");

  TString dataAsymmFile = Form("MB_asymmetries/AsymmFilesFromMB/UnCorr_OctetAsymmetries_AnaChD_Octets_0-59_BinByBin.txt");

  TH1D *dataAsymm = LoadMBAsymmetry(dataAsymmFile);

  TH1D *blindedDataAsymm = BlindAsymmetry(dataAsymm, bMixing, avg_mE);

  TH1D *twiddledAsymm = LoadTwiddledPairAsymm("/mnt/Data/xuansun/analyzed_files", twiddleIndex);

  TH1D* blindedDivision = new TH1D("blindedDivision", "blindedDivision", 120, 0, 1200);

  for(int i = 0; i <= blindedDataAsymm->GetNbinsX(); i++)
  {
    if(twiddledAsymm->GetBinContent(i) == 0 || blindedDataAsymm->GetBinContent(i) == 0)
    {
      blindedDivision->SetBinContent(i, 1);
      blindedDivision->SetBinError(i, 1);
    }
    else
    {
      blindedDivision->SetBinContent(i, twiddledAsymm->GetBinContent(i) / blindedDataAsymm->GetBinContent(i));
      blindedDivision->SetBinError(i, sqrt(pow(twiddledAsymm->GetBinContent(i), 2.0) + pow(blindedDataAsymm->GetBinContent(i), 2.0)));
    }
  }


  TF1 *fit = new TF1("beta fit", Form("( (1.0 + [0]*(%f))/(x + %f) ) / (1.0 + [0]*(%f))", m_e, m_e, avg_mE), xMin, xMax);

  fit->SetParName(0, "b");
  blindedDivision->Fit("beta fit", "R");
  TF1 *fitResults = blindedDivision->GetFunction("beta fit");
  cout << "Chi squared value is " << fitResults->GetChisquare()
	<< " with ndf of " << fitResults->GetNDF()
	<< ". For a final chisquared/ndf = " << fitResults->GetChisquare() / fitResults->GetNDF() << endl;


  PlotHist(C, 1, 1, blindedDivision, "", "Reconstructed Energy (keV)", "#frac{A_{twiddled}}{A_{data}}", "");



  TLine *yLow = new TLine(xMin, gPad->GetUymin(), xMin, gPad->GetUymax());
  TLine *yHigh = new TLine(xMax, gPad->GetUymin(), xMax, gPad->GetUymax());
  yLow->SetLineWidth(2);
  yHigh->SetLineWidth(2);
  yLow->SetLineColor(1);
  yHigh->SetLineColor(1);
  yLow->Draw("SAME");
  yHigh->Draw("SAME");

/*
  TLatex t2;
  t2.SetTextSize(0.03);
  t2.SetTextAlign(13);
  t2.DrawLatex(1000, 0.14, Form("A_{fit} = %f", fitResults->GetParameter(0)));
  TLatex t3;
  t3.SetTextSize(0.03);
  t3.SetTextAlign(13);
  t3.DrawLatex(1000, 0.13, Form("AErr_{fit} = %f", fitResults->GetParError(0)));
*/
  TLatex t4;
  t4.SetTextSize(0.03);
  t4.SetTextAlign(13);
  t4.DrawLatex(1000, 1.5, Form("b_{fit} = %f", fitResults->GetParameter(0)));
  TLatex t5;
  t5.SetTextSize(0.03);
  t5.SetTextAlign(13);
  t5.DrawLatex(1000, 1.0, Form("bErr_{fit} = %f", fitResults->GetParError(0)));
  TLatex t6;
  t6.SetTextSize(0.03);
  t6.SetTextAlign(13);
  t6.DrawLatex(900, 0.5, Form("#frac{#Chi^{2}}{ndf} = #frac{%f}{%i} = %f",
				fitResults->GetChisquare(), fitResults->GetNDF(), fitResults->GetChisquare() / fitResults->GetNDF()));


  ofstream outfile;
  outfile.open(Form("Xuan_asymmetries/Twiddled_Asymmetries/TwiddledbValues_asymmetryFitter_twiddleIndex_%i.txt", twiddleIndex), ios::app);
  outfile << twiddleIndex << "\t"
          << avg_mE << "\t"
          << fitResults->GetChisquare() << "\t"
          << fitResults->GetNDF() << "\t"
          << fitResults->GetChisquare() / fitResults->GetNDF() << "\t"
          << fitResults->GetParameter(0) << "\t"
          << fitResults->GetParError(0) << "\n";
  outfile.close();


  // Save our plot and print it out as a pdf.
  C -> Print(Form("Xuan_asymmetries/Twiddled_Asymmetries/b_fit_toTwiddledAsymm_twiddleIndex_%i.pdf", twiddleIndex));
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

TH1D* LoadTwiddledPairAsymm(TString pathBase, int twiddle)
{

  // load files that have already been twiddled
//  TFile fDown(Form("%s/TwiddledSimFiles_A_-1_b_0/SimAnalyzed_2011-2012_Beta_paramSet_%i_0.root", pathBase.Data(), twiddle));
//  TFile fUp(Form("%s/TwiddledSimFiles_A_1_b_0/SimAnalyzed_2011-2012_Beta_paramSet_%i_0.root", pathBase.Data(), twiddle));

//  TTree *tDown = (TTree*)fDown.Get("SimAnalyzed");
//  TTree *tUp = (TTree*)fUp.Get("SimAnalyzed");

  // You NEED to use the loop method over events.
  // For some reason, TTree->Draw into histogram is causing a seg fault.
  TChain* chainDown = new TChain("SimAnalyzed");
  chainDown->Add(Form("%s/TwiddledSimFiles_A_-1_b_0/SimAnalyzed_2011-2012_Beta_paramSet_%i_0.root", pathBase.Data(), twiddle));
  TChain* chainUp = new TChain("SimAnalyzed");
  chainUp->Add(Form("%s/TwiddledSimFiles_A_1_b_0/SimAnalyzed_2011-2012_Beta_paramSet_%i_0.root", pathBase.Data(), twiddle));

  cout << "Finished making chains of our events... " << endl;

  Event* evtDown = new Event;
  Event* evtUp = new Event;

  chainDown->SetBranchAddress("PID", &evtDown->PID);
  chainDown->SetBranchAddress("Erecon", &evtDown->Erecon);
  chainDown->SetBranchAddress("side", &evtDown->side);
  chainDown->SetBranchAddress("type", &evtDown->type);
  chainUp->SetBranchAddress("PID", &evtUp->PID);
  chainUp->SetBranchAddress("Erecon", &evtUp->Erecon);
  chainUp->SetBranchAddress("side", &evtUp->side);
  chainUp->SetBranchAddress("type", &evtUp->type);

  TH1D* hEastDown = new TH1D("hEastDown", "hEastDown", 120, 0, 1200);
  TH1D* hWestDown = new TH1D("hWestDown", "hWestDown", 120, 0, 1200);
  TH1D* hEastUp = new TH1D("hEastUp", "hEastUp", 120, 0, 1200);
  TH1D* hWestUp = new TH1D("hWestUp", "hWestUp", 120, 0, 1200);

  // put everything into east/west, spin up/down histograms.
  for(int i = 0; i < chainDown->GetEntries(); i++)
  {
    chainDown->GetEntry(i);
    if(evtDown->PID == 1 && evtDown->type == 0 && evtDown->Erecon > 0 && evtDown->side < 2)
    {
      if(evtDown->side == 0)
      {
        hEastDown->Fill(evtDown->Erecon);
      }
      if(evtDown->side == 1)
      {
        hWestDown->Fill(evtDown->Erecon);
      }
    }
  }
  for(int i = 0; i < chainUp->GetEntries(); i++)
  {
    chainUp->GetEntry(i);
    if(evtUp->PID == 1 && evtUp->type == 0 && evtUp->Erecon > 0 && evtUp->side < 2)
    {
      if(evtUp->side == 0)
      {
        hEastUp->Fill(evtUp->Erecon);
      }
      if(evtUp->side == 1)
      {
        hWestUp->Fill(evtUp->Erecon);
      }
    }
  }

  cout << "Constructing a super ratio... " << endl;

  // construct a super ratio
  TH1D* hSR = new TH1D("hSuperRatio", "hSuperRatio", 120, 0, 1200);
  vector <double> mySetErrors;

  double r1up, r1down, r2up, r2down;
  double er1up, er1down, er2up, er2down;
  for(int i = 0; i <= hSR->GetNbinsX(); i++)
  {
    r1up = hEastUp->GetBinContent(i);
    r1down = hEastDown->GetBinContent(i);
    r2up = hWestUp->GetBinContent(i);
    r2down = hWestDown->GetBinContent(i);
    er1up = hEastUp->GetBinError(i);
    er1down = hEastDown->GetBinError(i);
    er2up = hWestUp->GetBinError(i);
    er2down = hWestDown->GetBinError(i);


    if(r1up <= 0 || r1down <= 0 || r2up <= 0 || r2down <= 0 )
    {
      hSR->SetBinContent(i, 0);
      mySetErrors.push_back(0);
    }
    else
    {
      hSR->SetBinContent(i, (r1down*r2up) / (r1up*r2down));

      mySetErrors.push_back( ( (r1down*r2up) / (r1up*r2down) )
                           * sqrt( (er1down/r1down)*(er1down/r1down)
                                 + (er2up/r2up)*(er2up/r2up)
                                 + (er1up/r1up)*(er1up/r1up)
                                 + (er2down/r2down)*(er2down/r2down) ) );
    }
  }

  hSR->SetError(&(mySetErrors[0]));

  cout << "Created a super ratio.." << endl;

  // make the asymmetry histogram.
  cout << "Making asymmetry now... " << endl;
  TH1D* hTwiddledAsymm = new TH1D("hTwiddledAsymm", "hTwiddledAsymm", 120, 0, 1200);

  double asymm = 0;
  double asymmErr = 0;

  for(int i = 0; i <= hSR->GetNbinsX(); i++)
  {
    if(hSR->GetBinContent(i) == 0)
    {
      asymm = 0;
      asymmErr = 0;
    }
    else
    {
      asymm = (1 - sqrt(hSR->GetBinContent(i)) ) / ( 1 + sqrt(hSR->GetBinContent(i)) );
      asymmErr = (hSR->GetBinError(i))
                 / ( sqrt(hSR->GetBinContent(i)) * (sqrt(hSR->GetBinContent(i)) + 1) * (sqrt(hSR->GetBinContent(i)) + 1) );
    }

    hTwiddledAsymm->SetBinContent(i, asymm);
    hTwiddledAsymm->SetBinError(i, asymmErr);
  }

  return hTwiddledAsymm;
}






void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString xTitle, TString yTitle, TString command)
{
  cout << "1" << endl;
  C -> cd(canvasIndex);
  cout << "2" << endl;
  hPlot -> SetTitle(title);
  cout << "3" << endl;
  hPlot -> GetXaxis() -> SetTitle(xTitle);
  cout << "4" << endl;
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle(yTitle);
  hPlot -> GetYaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetRangeUser(0, 2);


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

//  C -> Update();
}
