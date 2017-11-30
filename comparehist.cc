#include	"comparehist.hh"

#define		HIST_IMAGE_PRINTOUT_NAME	"Test_comparehist"
#define		TYPE	"allTypes"

//required later for plot_program
//TApplication plot_program("FADC_readin",0,0,0,0);

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (octet #)" << endl;
    return 0;
  }

  int octNb = atoi(argv[1]);

/*  TString treeName = Form("Evts");
  TChain *MCTheoryChainBeta = MakeTChain("/home/xuansun/Documents/Analysis_Code/ucna_g4_2.1/UCN/UK_EventGen_2016/Evts_Files/b_0_300mill/Evts", treeName, 0, 100, 42);
  TChain *MCTheoryChainFierz = MakeTChain("/home/xuansun/Documents/Analysis_Code/ucna_g4_2.1/UCN/UK_EventGen_2016/Evts_Files/b_inf_100mill/Evts", treeName, 0, 100, 42);
  TChain *dataChain = MakeTChain("/home/xuansun/Documents/Analysis_Code/ucna_g4_2.1/UCN/UK_EventGen_2016/Evts_Files/b_1/Evts", treeName, 30, 31);
  TString variableName = Form("KE");
  TString cutsUsed = Form("");

  // Create a TChain
  TString treeName = Form("SimAnalyzed");
  TChain *MCTheoryChainBeta = MakeTChain("Data/BigSims_b_xsunCode/SimAnalyzed_2010_Beta_paramSet", treeName, 0, 100, 42);
  TChain *MCTheoryChainFierz = MakeTChain("Data/BigSims_b_xsunCode/SimAnalyzed_2010_Beta_fierz_paramSet", treeName, 0, 100, 42);
  TChain *dataChain = MakeTChain("/mnt/Data/xuansun/analyzed_files/SimAnalyzed_2010_Beta_paramSet", treeName, 0, 1, 7);

  // Get the Erecon histogram out with appropriate cuts
  TString variableName = Form("Erecon");
  TString cutsUsed = Form("type != 4 && side != 2");

  TH1D* dataHist = ExtractHistFromChain(variableName, cutsUsed, dataChain,
				      "myHist", "Test of comparehist code", 100, 0, 1000);
  TH1D* mcTheoryHistBeta = ExtractHistFromChain(variableName, cutsUsed, MCTheoryChainBeta,
                                      "mcBeta", "Test of comparehist code", 100, 0, 1000);
  TH1D* mcTheoryHistFierz = ExtractHistFromChain(variableName, cutsUsed, MCTheoryChainFierz,
                                      "mcFierz", "Test of comparehist code", 100, 0, 1000);
*/

  TFile fData(TString::Format("ExtractedHistograms/Data_Hists/Octet_%i_ssDataHist_%s.root", octNb, TYPE));
  TFile fMC0(TString::Format("/mnt/Data/xuansun/BLIND_MC_files/2011-2012_geom/BLIND_MC_A_0_b_0_Octet_%i_ssHist_%s.root", octNb, TYPE));
//  TFile fMC0(TString::Format("ExtractedHistograms/MC_A_0_b_0/MC_A_0_b_0_Octet_%i_ssHist%s.root", octNb, ""));
  TFile fMCinf(TString::Format("ExtractedHistograms/MC_A_0_b_inf/MC_A_0_b_inf_Octet_%i_ssHist%s.root", octNb, ""));

  TH1D* dataHist = (TH1D*)fData.Get("Super sum");
  TH1D* mcTheoryHistBeta = (TH1D*)fMC0.Get("Super sum");
  TH1D* mcTheoryHistFierz = (TH1D*)fMCinf.Get("Super sum");


  // Create a TFractionFitter and do the fit.
  TObjArray *MCTheory = new TObjArray(2);
  MCTheory -> Add(mcTheoryHistBeta);
  MCTheory -> Add(mcTheoryHistFierz);
  TFractionFitter* fit = new TFractionFitter(dataHist, MCTheory, "Q");	// initialise
  TVirtualFitter* vfit = fit->GetFitter();
  int fitMin = 10;
  int fitMax = 65;
  fit -> SetRangeX(fitMin, fitMax);	// Set range in bin numbers

  int status = fit->Fit();

  int fitPassNumber = 1;
  double value = 1.698;
  int entries = -100;
  int entriesData = -100;
  while(status != 0 || isnan(fit->GetChisquare()))
  {
    cout << "Fit unsuccessful at attempt number " << fitPassNumber << ". Trying again..." << endl;

    vfit->SetParameter(0, "a", value, 0.01, -10, 10);
    vfit->SetParameter(1, "c", 1-value, 0.01, -10, 10);

    status = fit->Fit();

    if(status == 0)
    {
      TH1D* resultHist = (TH1D*)fit->GetPlot();     // extract the plot from the fit.
      entries = 0;
      entriesData = 0;
      for(int i = fitMin; i < fitMax; i++)
      {
        entries = entries + resultHist->GetBinContent(i);
        entriesData = entriesData + dataHist->GetBinContent(i);
      }
    }
    value = value - 0.05;     // try again with a seed value slightly lower
    fitPassNumber++;

    if(fitPassNumber >= 30)
    {
      cout << "Fit attempted 30 times and failed. Exiting..." << endl;
      break;
    }
  }

  cout << "Octet number we are currently on: " << octNb << endl;

  cout << "Number events in fitted histogram (blue): " << entries << endl;
  cout << "Number events in data histogram (red): " << entriesData << endl;

  double avg_mE = CalculateAveragemOverE(mcTheoryHistBeta, fitMin, fitMax);
  cout << "Our average value and hence scaling factor is: " << avg_mE << endl;

  // Get valuable numbers for later
  double chisquared = fit->GetChisquare();
  int ndf = fit->GetNDF();
  double frac0Val, frac0Err, frac1Val, frac1Err;
  fit->GetResult(0, frac0Val, frac0Err);
  fit->GetResult(1, frac1Val, frac1Err);

  double bErr = Fierz_b_Error(frac0Val, frac0Err, frac1Val, frac1Err, avg_mE,
			      vfit->GetCovarianceMatrixElement(0,0), vfit->GetCovarianceMatrixElement(0,1),
			      vfit->GetCovarianceMatrixElement(1,0), vfit->GetCovarianceMatrixElement(1,1) );
  cout << "b error calculated by hand: " << bErr << endl;
  cout << "b value: " << frac1Val/(frac0Val*avg_mE) << endl;
  cout << "To compare, the limitation from 100KeV and up is: " << 10.1/sqrt(dataHist->GetEntries()) << endl;
  cout << "Entries used in theoretical limit: " << dataHist->GetEntries() << endl;

  ofstream outfile;
  outfile.open(Form("BLIND_ExtractedbValues_%s_comparehist.txt", TYPE), ios::app);
  outfile << octNb << "\t"
	  << frac1Val/(frac0Val*avg_mE) << "\t"
          << avg_mE << "\t"
	  << chisquared << "\t"
	  << ndf << "\t"
	  << bErr << "\t"
	  << 10.1/sqrt(dataHist->GetEntries()) << "\t"
	  << dataHist->GetEntries() << "\t"
	  << entries << "\n";
  outfile.close();


  // plot everything and visualize
/*  TCanvas *C = new TCanvas("canvas", "canvas");
  gROOT->SetStyle("Plain");	//on my computer this sets background to white, finally!
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat("en");
  gStyle->SetStatH(0.45);
  gStyle->SetStatW(0.45);
  PlotHist(C, 1, 1, dataHist, "Data histogram", "");
//  PlotHist(C, 2, 1, resultHist, "Fit Result histogram", "SAME");


  // update the legend to include valuable variables.
  TPaveStats *ps = (TPaveStats*)C->GetPrimitive("stats");
  ps->SetName("mystats");
  TList *listOfLines = ps->GetListOfLines();
  TLatex *myText1 = new TLatex(0,0,Form("#Chi^{2} = %f", chisquared));
  listOfLines->Add(myText1);
  TLatex *myText2 = new TLatex(0,0,Form("NDF = %d", ndf));
  listOfLines->Add(myText2);
  TLatex *myText3 = new TLatex(0,0,Form("#frac{#Chi^{2}}{NDF} = %f", chisquared/ndf));
  listOfLines->Add(myText3);
  TLatex *myText4 = new TLatex(0,0,Form("frac0 = %f #pm %f", frac0Val, frac0Err));
  listOfLines->Add(myText4);
  TLatex *myText5 = new TLatex(0,0,Form("frac1 = %f #pm %f", frac1Val, frac1Err));
  listOfLines->Add(myText5);
  // the following line is needed to avoid that the automatic redrawing of stats
  dataHist->SetStats(0);
*/
  // prints the canvas with a dynamic TString name of the name of the file
//  C -> Print(Form("%s.pdf", HIST_IMAGE_PRINTOUT_NAME));
  cout << "-------------- End of Program ---------------" << endl;
//  plot_program.Run();

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

  cout << "Loaded trees from files identified by the template: " << baseName.Data() << "_#.root" << endl;

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
