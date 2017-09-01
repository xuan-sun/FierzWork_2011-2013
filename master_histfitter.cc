#include	"comparehist.hh"

#define		HIST_IMAGE_PRINTOUT_NAME	"Test_master_histfitter"
#define		OUTPUT_ANALYSIS_FILE		"AnalyzedTextFiles/b_0_SimProcessed_noTwiddles_100keV-650keV_firstPass_forStatError.txt"

int main()
{
  TString treeName = Form("SimAnalyzed");
  TChain *MCTheoryChainBeta = MakeTChain("Data/BigSims_b_xsunCode/SimAnalyzed_2010_Beta_paramSet", treeName, 0, 100, 42);
  TChain *MCTheoryChainFierz = MakeTChain("Data/BigSims_b_xsunCode/SimAnalyzed_2010_Beta_fierz_paramSet", treeName, 0, 100, 42);
  TChain *dataChain = MakeTChain("Data/BigSims_b_xsunCode/SimAnalyzed_2010_Beta_paramSet", treeName
				, ReplaceWithIndexLow, ReplaceWithIndexHigh, 42);

  TString variableName = Form("Erecon");
  TString cutsUsed = Form("type != 4 && side != 2");
  TH1D* mcTheoryHistBeta = ExtractHistFromChain(variableName, cutsUsed, MCTheoryChainBeta,
                                      "mcBeta", "Beta", 100, 0, 1000);
  TH1D* mcTheoryHistFierz = ExtractHistFromChain(variableName, cutsUsed, MCTheoryChainFierz,
                                      "mcFierz", "Fierz", 100, 0, 1000);
  TH1D* dataHist = ExtractHistFromChain(variableName, cutsUsed, dataChain, "myHist", "Data", 100, 0, 1000);

  TObjArray *MCTheory = new TObjArray(2);
  MCTheory -> Add(mcTheoryHistBeta);
  MCTheory -> Add(mcTheoryHistFierz);
  TFractionFitter* fit = new TFractionFitter(dataHist, MCTheory, "Q");  // initialise
  TVirtualFitter* vfit = fit->GetFitter();

  int fitMin = 10;
  int fitMax = 65;
  fit -> SetRangeX(fitMin, fitMax);

  int status = fit->Fit();

  int fitPassNumber = 1;
  double value = 1.19;
  int entries = -100;
  while(status != 0)
  {
    cout << "Fit unsuccessful at attempt number " << fitPassNumber << ". Trying again..." << endl;

    vfit->SetParameter(0, "a", value, 0.01, -10, 10);
    vfit->SetParameter(1, "c", 1-value, 0.01, -10, 10);

    status = fit->Fit();

    if(status == 0)
    {
      TH1D* resultHist = (TH1D*)fit->GetPlot();     // extract the plot from the fit.
      entries = 0;
      for(int i = fitMin; i < fitMax; i++)
      {
        entries = entries + resultHist->GetBinContent(i);
      }
      if(entries > 1000)
      {
        break;	// get out of our fitter loop
      }
      else if(entries <= 0)
      {
        status = 1;	// set our status flag back to != 0, continue loop.
      }
    }
    value = value - 0.05;	// try again with a seed value slightly lower
    fitPassNumber++;

    if(value < 0.5)
    {
      cout << "Fit not successful at lowest boundary of fraction values. Exiting." << endl;
      break;
    }
  }


  double avg_mE = CalculateAveragemOverE(mcTheoryHistBeta, fitMin, fitMax);

  // Get valuable numbers for later
  double chisquared = fit->GetChisquare();
  int ndf = fit->GetNDF();
  double frac0Val, frac0Err, frac1Val, frac1Err;
  fit->GetResult(0, frac0Val, frac0Err);
  fit->GetResult(1, frac1Val, frac1Err);

  ofstream outfile;
  outfile.open(OUTPUT_ANALYSIS_FILE, ios::app);
  outfile << frac1Val/(frac0Val*avg_mE) << "\t"
	  << avg_mE << "\t"
	  << Fierz_b_Error(frac0Val, frac0Err, frac1Val, frac1Err, avg_mE,
                              vfit->GetCovarianceMatrixElement(0,0), vfit->GetCovarianceMatrixElement(0,1),
                              vfit->GetCovarianceMatrixElement(1,0), vfit->GetCovarianceMatrixElement(1,1) ) << "\t"
          << fitMin << "\t" << fitMax << "\t"
	  << entries << "\t"
          << "SimAnalyzed_2010_Beta_paramSet_" << 42 << "_" << ReplaceWithIndexLow << ".root" << "\t"
	  << ReplaceWithIndexLow << "\t"
	  << chisquared << "\t"
	  << ndf << "\t"
	  << chisquared/ndf << "\n";
  outfile.close();

  return 0;
}
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Energy (keV)");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("Counts (N)");
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(0, 0.000004);

  if(styleIndex == 1)
  {
    hPlot -> SetFillColor(46);
    hPlot -> SetFillStyle(3004);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 2)
  {
    hPlot -> SetFillColor(38);
    hPlot -> SetFillStyle(3005);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 3)
  {
    hPlot -> SetFillColor(29);
//    hPlot -> SetFillStyle(3005);
    hPlot -> SetFillStyle(3001);
  }

  hPlot -> Draw(command);
  C -> Update();
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors* gPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  gPlot -> SetTitle(title);
  gPlot -> GetXaxis() -> SetTitle("Energy (keV)");
  gPlot -> GetXaxis() -> SetRangeUser(0, 1000);
  gPlot -> GetXaxis() -> CenterTitle();
  gPlot -> GetYaxis() -> SetTitle("Hits (N)");
  gPlot -> GetYaxis() -> CenterTitle();

  if(styleIndex == 1)
  {
    gPlot -> SetMarkerColor(46);
  }
  if(styleIndex == 2)
  {
    gPlot -> SetMarkerColor(38);
  }
  if(styleIndex == 3)
  {
    gPlot -> SetMarkerColor(30);
  }
  gPlot -> SetMarkerStyle(7);
  gPlot -> SetMarkerSize(1);

  // add a flat line at y = 0
  TLine *line = new TLine(0, 0, 1000, 0);

  gPlot -> Draw(command);
  line -> Draw("SAME");

}

TH1D* ExtractHistFromChain(TString varName, TString cutsUsed, TChain* chain,
                           TString name, TString title, int nbBins, double minX, double maxX)
{
  TH1D* hist = new TH1D(name.Data(), title.Data(), nbBins, minX, maxX);

  chain -> Draw(Form("%s >> %s",varName.Data(),name.Data()), cutsUsed.Data());

  cout << "Completed storing histogram of variable " << varName.Data() << " into histogram " << name.Data() << endl;

  return hist;
}

TChain* MakeTChain(TString baseName, TString treeName, int fileNumMin, int fileNumMax, int paramIndex)
{
  TChain* chain = new TChain(treeName.Data());

  for(int i = fileNumMin; i < fileNumMax; i++)
  {
    chain -> AddFile(Form("%s_%i_%i.root", baseName.Data(), paramIndex, i));
  }

  cout << "Loaded trees from files identified by the template: " << baseName.Data() << "_#.root \n"
       << "Lower index = " << fileNumMin << " and upper index = " << fileNumMax << "\n"
       << "Parameter (i.e. twiddle) index = " << paramIndex << endl;

  return chain;
}

double CalculateChiSquared(TH1D* hdat, TH1D* hthBeta, TH1D* hthFierz, double frac0, double frac1, double xBinMin, double xBinMax)
{
  TH1D* hmc = new TH1D("MCHist", "MCHist", 100, 0, 1000);
  hmc -> Add(hthBeta, hthFierz, frac0, frac1);
  double norm = hdat->GetEntries() / hmc -> GetEntries();
  hmc -> Scale(norm);

  double chisquare = 0;
  for(int i = xBinMin; i <= xBinMax; i++)
  {
    chisquare = chisquare + ((hdat->GetBinContent(i) - hmc->GetBinContent(i))*
			     (hdat->GetBinContent(i) - hmc->GetBinContent(i)))/
			    hdat->GetBinContent(i);
  }

  return chisquare;
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

double Fierz_b_Error(double f0v, double f0e, double f1v, double f1e, double avgInverseW
		    , double cov00, double cov01, double cov10, double cov11)
{
  double errb = 0;

  errb = abs(f1v/(f0v*avgInverseW))
	 * sqrt((f0e/f0v)*(f0e/f0v) + (f1e/f1v)*(f1e/f1v) - (2*cov01*avgInverseW*f0v*avgInverseW*f1v));


  return errb;
}
