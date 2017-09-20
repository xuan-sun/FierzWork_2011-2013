#include "comparehist.hh"

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString command);
void FillArrays(TString fileName, TH1D *hist);

struct entry
{
  int octNb;
  double avg_mE;
  double chisquared;
  double ndf;
  double bErr;
  double b;
  double GluckErr;
  int entryNum;
  double numberOfOriginalEvents;
};

vector <double> octets;
vector <double> FierzValues;
vector <double> errorBarsEntries;

int main()
{
  TCanvas *C = new TCanvas("canvas", "canvas");
  C -> Divide(2,1);
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  TH1D *h1 = new TH1D("myhist", "myhist", 50, -0.5, 0.5);

  FillArrays("ResultsOfAllTFractionFitters_Fierz2011-2012_noCuts.txt", h1);

  vector <double> xErr;
  vector <double> yErr;
  for(unsigned int i = 0; i < errorBarsEntries.size(); i++)
  {
    xErr.push_back(0.5);	// half an octet number
    yErr.push_back(sqrt(2)*10.1 / sqrt(errorBarsEntries[i]));	// converting using Gluck formula for 100KeV and up fit.
  }

  TGraphErrors *g1 = new TGraphErrors(octets.size(), &(octets[0]), &(FierzValues[0]), &(xErr[0]), &(yErr[0]));

  PlotHist(C, 1, 1, h1, "Extracted b values.", "");
  PlotGraph(C, 1, 2, g1, "b values by octet", "AP");

  //prints the canvas with a dynamic TString name of the name of the file
  C -> Print(Form("%s.png", "plotFierz"));
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Extracted b");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("N(b)");
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(0, 0.000004);

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

  hPlot -> Draw(command);
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors *gPlot, TString title, TString command)
{
  C->cd(canvasIndex);
  gPlot->SetTitle(title);
  gPlot->GetXaxis()->SetTitle("Octet Number");
  gPlot->GetXaxis()->CenterTitle();
  gPlot->GetYaxis()->SetTitle("Fitted b Value");
  gPlot->GetYaxis()->CenterTitle();

  if(styleIndex == 1)
  {
    gPlot->SetMarkerStyle(21);
    gPlot->SetMarkerSize(0.5);
    gPlot->SetMarkerColor(2);
  }

  gPlot->Draw(command);

}

void FillArrays(TString fileName, TH1D* hist)
{

  entry evt;
  int counter = 0;

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
      bufstream1 >> evt.octNb
		>> evt.avg_mE
		>> evt.chisquared
		>> evt.ndf
		>> evt.bErr
		>> evt.b
		>> evt.GluckErr
		>> evt.entryNum
		>> evt.numberOfOriginalEvents;
      {
	counter++;
        hist -> Fill(evt.b);
	octets.push_back(evt.octNb);
	FierzValues.push_back(evt.b);
	errorBarsEntries.push_back(evt.numberOfOriginalEvents);
      }
    }

    if(infile1.eof() == true)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}
