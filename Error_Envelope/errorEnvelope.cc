
//Should make an array of all the calibration periods and the
// runs they apply to, and then put vertical lines at each of these
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>

#include <TString.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <TF1.h>
#include <TList.h>
#include <TPaveStats.h>
#include <TLegend.h>

using namespace	std;

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

const static double peakCe = 130.3;// 131.5956;//80.5;
const static double peakIn = 174.35;
const static double peakSn = 368.4938;//317.8;
const static double peakBiLow = 498.;//501.420;//448.8;
const static double peakBiHigh = 993.789;//926.;

bool color = true;

vector < vector <Double_t > > GetGeomEnvelope(TString geom);
void ErrorEnvelope();

int main(int argc, char* argv[])
{
  if(argc < 1)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable)" << endl;
    return 0;
  }

  ErrorEnvelope();

  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}




//This returns mean, sigma, and rms, where the order is Ce,In,Sn,Bi2,Bi1
std::vector < std::vector < Double_t > > GetGeomEnvelope(TString geom) {

  
  std::vector < std::vector <Double_t> > ret(3,std::vector<Double_t>(0));
  Int_t lowCal, highCal;
  
  if ( geom==TString("2011-2012") ) lowCal=1, highCal=12;
  else if ( geom==TString("2012-2013") ) lowCal=16, highCal=24;
  
  else std::cout << "Bad geometry passed to getGeomEnvelope\n", exit(0);

  TString path = TString::Format("/home/xuansun/Documents/MBrown_Work/ParallelAnalyzer-1/error_envelope/error_envelope_calPeriods_%i-%i.dat",lowCal,highCal);
  std::cout << path << std::endl;

  std::ifstream infile(path.Data());  
  std::string txt = "";

  if (infile.is_open()) {
    for ( UInt_t i = 1; i<16; ++i ) {
      infile >> txt >> txt >> txt;
      if ( i%3 == 0 && i!=6 ) ret[2].push_back(atof(txt.c_str())); // RMS
      if ( (i+1)%3 == 0 && i!=5 ) ret[1].push_back(atof(txt.c_str())); // Err
      if ( (i+2)%3 == 0 && i!=4 ) ret[0].push_back(atof(txt.c_str())); // mean
     
    }
  }
  else {
    std::cout << "Couldn't open error envelope file!!\n";
    exit(0);
  }
  infile.close();
  return ret;
}


void ErrorEnvelope() {

  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleYOffset(0.8);
  gStyle->SetTitleXSize(0.06);
  //gStyle->SetTitleXOffset(0.7);
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetOptFit(1111);
  gStyle->SetStatY(0.85);
  gStyle->SetStatX(0.975);
  gStyle->SetStatW(.09);

  gStyle->SetFillStyle(0000); 
  //gStyle->SetStatStyle(0); 
  //gStyle->SetTitleStyle(0); 
  //gStyle->SetCanvasBorderSize(0); 
  //gStyle->SetFrameBorderSize(0); 
  gStyle->SetLegendBorderSize(0); 
  //gStyle->SetStatBorderSize(0); 
  //gStyle->SetTitleBorderSize(0);


  Double_t xmin = 0.0;
  Double_t xmax = 1050.;
  TLine *zeroLine = new TLine(xmin,0.,xmax,0.);

  Int_t col2011 = color?9:14;
  Int_t col2012 = color?8:17;
  Int_t fill2011 = 1001;
  Int_t fill2012 = 1001;
  Int_t marker2011 = 20;
  Int_t marker2012 = 24;

  Double_t sep = 20.;


  TCanvas *c2 = new TCanvas("c2","c2");
  
  double En2011[4];
  En2011[0] = peakCe;//98.2;
  //  En[1] = peakIn;
  En2011[1] = peakSn;//331.2;
  En2011[2] = peakBiLow; //443.0;
  En2011[3] = peakBiHigh;//928.0;

  double En2012[4];
  En2012[0] = peakCe+sep;//98.2;
  //  En[1] = peakIn;
  En2012[1] = peakSn+sep;//331.2;
  En2012[2] = peakBiLow+sep; //443.0;
  En2012[3] = peakBiHigh+sep;//928.0;
 
  double dEn[4] = {0.};
  
  double rmsY[4] = {0.};

  std::vector < std::vector<Double_t> > env2011 = GetGeomEnvelope("2011-2012");
  std::vector < std::vector<Double_t> > env2012 = GetGeomEnvelope("2012-2013");

 

  /////////////////////////////////////////////////// 2011-2012 /////////////////////////////////////////////////////////


  TGraphErrors *RMS2011 = new TGraphErrors(4,En2011,rmsY,dEn,&env2011[2][0]);
  RMS2011->SetTitle(TString::Format("RMS"));
  RMS2011->SetMarkerColor(8);
  RMS2011->SetLineColor(8);
  RMS2011->SetLineWidth(15);
  RMS2011->SetMarkerStyle(21);
  RMS2011->SetMarkerSize(0);
  RMS2011->SetFillStyle(0);

  TGraphErrors *gr2011 = new TGraphErrors(4,En2011,&env2011[0][0],dEn,&env2011[1][0]);
  gr2011->SetTitle(TString::Format("2011-2012"));
  gr2011->SetMarkerColor(1);
  gr2011->SetLineColor(1);
  gr2011->SetLineWidth(2);
  gr2011->SetMarkerStyle(marker2011);
  gr2011->SetMarkerSize(1.25);
  gr2011->SetFillStyle(fill2011);
  gr2011->SetFillColorAlpha(col2011,1.0);  


  //Read in the error envelope from Kevin Hickerson
  std::ifstream errEnv2011(TString::Format("/home/xuansun/Documents/MBrown_Work/ParallelAnalyzer-1/systematics/EnergyUncertainty/envolopeValues-%i-deg2-cal4-curves1000.tsv",2011));
  std::vector <Double_t> en2011;
  std::vector<Double_t> maxEnv2011;
  std::vector <Double_t> low2011;
  std::vector <Double_t> high2011;
  std::vector <Double_t> yenv2011;
  std::vector <Double_t> xerr2011;
  
  Double_t e,l,h;

  while ( errEnv2011 >> e >> l >> h ) {
    en2011.push_back(e);
    std::cout << e << " " << l << " " << h << std::endl;
    low2011.push_back(TMath::Abs(l));
    high2011.push_back(h);
    maxEnv2011.push_back(h>TMath::Abs(l)? h : TMath::Abs(l));
    yenv2011.push_back(0.);
    xerr2011.push_back(5.);
  }
/*
  TGraphAsymmErrors *finenv2011 = new TGraphAsymmErrors(high2011.size(),&en2011[0],&yenv2011[0],&xerr2011[0],&xerr2011[0],&maxEnv2011[0],&maxEnv2011[0]);
  finenv2011->SetFillColorAlpha(col2011,0.8);
  finenv2011->SetFillStyle(fill2011);
*/
  TGraphAsymmErrors *finenv2011 = new TGraphAsymmErrors(high2011.size(),&en2011[0],&yenv2011[0],&xerr2011[0],&xerr2011[0],&low2011[0],&high2011[0]);
  finenv2011->SetFillColorAlpha(col2011,0.8);
  finenv2011->SetFillStyle(fill2011);

  // This little bit is all my own (Xuan) creation in order to print out the values of the error envelopes.
  TCanvas *c3 = new TCanvas("c3", "c3");
  c3->cd();
  finenv2011->Draw("AP");
  c2->cd();
  double x2011 = -1;
  double y2011 = -1;
  ofstream outfile2011;
  outfile2011.open("MB_ErrorEnvelope_2011-2012.txt", ios::app);
  for(int i = 0; i < finenv2011->GetN(); i++)
  {
    finenv2011->GetPoint(i, x2011, y2011);
    outfile2011 << i << "\t"
	        << x2011 << "\t"
	        << y2011 << "\t"
	        << finenv2011->GetErrorYhigh(i) << "\t"
	        << finenv2011->GetErrorYlow(i) << "\n";
  }
  outfile2011.close();


  /////////////////////////////////////////////////// 2012-2013 /////////////////////////////////////////////////////////

  TGraphErrors *RMS2012 = new TGraphErrors(4,En2012,rmsY,dEn,&env2012[2][0]);
  RMS2012->SetTitle(TString::Format("RMS"));
  RMS2012->SetMarkerColor(8);
  RMS2012->SetLineColor(8);
  RMS2012->SetLineWidth(15);
  RMS2012->SetMarkerStyle(21);
  RMS2012->SetMarkerSize(0);
  RMS2012->SetFillStyle(0);

  TGraphErrors *gr2012 = new TGraphErrors(4,En2012,&env2012[0][0],dEn,&env2012[1][0]);
  gr2012->SetTitle(TString::Format("2012-2013"));
  gr2012->SetMarkerColor(1);
  gr2012->SetLineColor(1);
  gr2012->SetLineWidth(2);
  gr2012->SetMarkerStyle(marker2012);
  gr2012->SetMarkerSize(1.25);
  gr2012->SetFillStyle(fill2012);
  gr2012->SetFillColor(col2012);  


  //Read in the error envelope from Kevin Hickerson
  std::ifstream errEnv2012(TString::Format("/home/xuansun/Documents/MBrown_Work/ParallelAnalyzer-1/systematics/EnergyUncertainty/envolopeValues-%i-deg2-cal4-curves1000.tsv",2012));
  std::vector <Double_t> en2012;
  std::vector<Double_t> maxEnv2012;
  std::vector <Double_t> low2012;
  std::vector <Double_t> high2012;
  std::vector <Double_t> yenv2012;
  std::vector <Double_t> xerr2012;
  
  while ( errEnv2012 >> e >> l >> h ) {
    en2012.push_back(e);
    std::cout << e << " " << l << " " << h << std::endl;
    low2012.push_back(TMath::Abs(l));
    high2012.push_back(h);
    maxEnv2012.push_back(h>TMath::Abs(l)? h : TMath::Abs(l));
    yenv2012.push_back(0.);
    xerr2012.push_back(5.);
  }

  TGraphAsymmErrors *finenv2012 = new TGraphAsymmErrors(high2012.size(),&en2012[0],&yenv2012[0],&xerr2012[0],&xerr2012[0],&maxEnv2012[0],&maxEnv2012[0]);
  finenv2012->SetFillColorAlpha(col2012,0.7);
  finenv2012->SetFillStyle(fill2012);

  // This little bit is all my own (Xuan) creation in order to print out the values of the error envelopes.
  double x2012 = -1;
  double y2012 = -1;
  ofstream outfile2012;
  outfile2012.open("MB_ErrorEnvelope_2012-2013.txt", ios::app);
  for(int i = 0; i < finenv2012->GetN(); i++)
  {
    finenv2012->GetPoint(i, x2012, y2012);
    outfile2012 << i << "\t"
                << x2012 << "\t"
                << y2012 << "\t"
                << finenv2012->GetErrorYhigh(i) << "\t"
                << finenv2012->GetErrorYlow(i) << "\n";
  }
  outfile2012.close();


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  TMultiGraph *mg = new TMultiGraph();
  mg->Add(finenv2011,"3");
  mg->Add(finenv2012,"3");
  //mg->Add(RMS,"PZ");
  //mg->Draw("A");
  mg->Add(gr2011,"P");
  mg->Add(gr2012,"P");
 

  mg->Draw("A");
  mg->GetXaxis()->SetTitle("E_{recon} [keV]");
  mg->GetXaxis()->SetTitleOffset(1.2);
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetTitle("Calibration Residual [keV]");
  // mg->GetYaxis()->SetTitleOffset(1.2);
  mg->GetYaxis()->CenterTitle();
  

  mg->GetXaxis()->SetLimits(xmin,xmax);
  mg->SetMinimum(-20.0);
  mg->SetMaximum( 20.0);

  zeroLine->SetLineWidth(2);
  zeroLine->SetLineStyle(2);
  zeroLine->Draw("SAME");

  TLegend *leg = new TLegend(0.65,0.70,0.875,0.875);
  leg->AddEntry(gr2011,0,"fp");
  leg->AddEntry(gr2012,0,"fp");
  leg->Draw();


 
  /*TGraph *env_upper = new TGraph(high.size(),&en[0],&high[0]);
  env_upper->Draw("Same");
  env_upper->SetLineWidth(3);
  env_upper->SetLineColor(2);
  env_upper->SetLineStyle(8);
  
  TGraph *env_lower = new TGraph(nn,x2,y_lower);
  env_lower->Draw("Same");
  env_lower->SetLineWidth(3);
  env_lower->SetLineColor(2);
  env_lower->SetLineStyle(8); 
  */
  c2->Update();
  TString filename = TString::Format("energyErrorEnvelope%s.pdf",(color?"_color":"")); 
  c2->Print("EnergyErrorEnvelope_2011-2013.pdf");
 
}
