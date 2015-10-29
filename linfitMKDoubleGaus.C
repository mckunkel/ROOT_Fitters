#include "TROOT.h"
#include "TRint.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TVector3.h"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TGStatusBar.h"
#include "TSystem.h"
#include "TXMLEngine.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TLine.h"
#include "TLatex.h"

#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include<iomanip>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

//This is in the root directory and as of May 15 2012, is the final version

void linfitMKDoubleGaus(TH1F *h33 , Double_t low, Double_t high, Double_t p0, Double_t p1, Double_t initialPar, Double_t width, Double_t factor, Int_t draw_opt){

  double nEnt = h33->GetEntries();
  TF1 *fitter = new TF1("fitter","gaus(0) + gaus(3) + [6] + [7]*x",low,high);
  fitter->SetParameters(100, initialPar, width, 100, initialPar, width/2., p0, p1);
  fitter->SetParLimits(0,10,nEnt); fitter->SetParLimits(1,initialPar-width,initialPar+width);
  fitter->SetParLimits(3,10,nEnt); fitter->SetParLimits(4,initialPar-width/2.,initialPar+width/2.);
h33->Fit("fitter","REM","same");
  //h33->Draw("E");
TF1 *backFcn = new TF1("backFcn", "[0] + [1]*x",low,high);
TF1 *signalFcn = new TF1("signalFcn", "gaus(0) + gaus(3)",low,high);
signalFcn->SetLineColor(2);
signalFcn->SetLineWidth(2);
Double_t par[8];
fitter->GetParameters(par);
  
  signalFcn->SetParameters(par);
  signalFcn->SetLineStyle(2);
  signalFcn->SetLineColor(4);
  signalFcn->SetLineWidth(1);
  
backFcn->SetParameters(&par[6]);
backFcn->SetLineStyle(2);
backFcn->SetLineColor(6);
backFcn->SetLineWidth(1);



backFcn->Draw("same");
signalFcn->Draw("same");
Double_t Intg = abs(signalFcn->Integral(par[1]-factor*fabs(par[2]),par[1]+factor*fabs(par[2])));
Double_t Intb = abs(backFcn->Integral(par[1]-factor*fabs(par[2]),par[1]+factor*fabs(par[2])));


Double_t binw = h33->GetBinWidth(1);
Int_t yield = Intg/binw;
Int_t bckgd = Intb/binw;
Double_t ratio = double(yield)/(double(yield+bckgd)); //(double(yield+bckgd)
cout << yield << "\t" <<bckgd<<"\t"<< ratio << endl;
TAxis *x=h33->GetXaxis();
TAxis *y=h33->GetYaxis();

Double_t startx=x->GetXmin()+0.75*(x->GetXmax()-x->GetXmin());
Double_t starty=0.65*h33->GetMaximum();
Double_t starty1=0.55*h33->GetMaximum();
Double_t starty0=0.75*h33->GetMaximum();
Double_t starty2=0.95*h33->GetMaximum();
Double_t starty3=0.85*h33->GetMaximum();

double meanError = fitter->GetParError(4);
double sigmaError = fitter->GetParError(5);
  if (draw_opt ==1) {

TLatex *sum = new TLatex(startx*0.83, starty,Form("Yield: %i",yield));
TLatex *sum12 = new TLatex(startx*0.83, starty*0.84,Form("Background: %i",bckgd));
TLatex *sum0=new TLatex(startx*0.83, starty0,Form("Range: #pm %2.1f #sigma",factor));
TLatex *sum2=new TLatex(startx*0.83, starty2,Form("Mean:%4.4f #pm %.4f GeV",par[4], meanError));
TLatex *sum3=new TLatex(startx*0.83, starty3,Form("#sigma:%5.4f #pm %.4f GeV",par[5], sigmaError));
TLatex *ra = new TLatex(startx*0.83, starty*0.6,Form("#frac{S}{S+B}= %.3f", ratio));
sum->SetTextSize(0.04);
sum->SetTextColor(2);
sum->Draw("same");
sum12->SetTextSize(0.04);
sum12->SetTextColor(6);
sum12->Draw("same");

ra->Draw("same");


sum0->SetTextSize(0.04);
sum0->SetTextColor(2);
sum0->Draw("same");
sum2->SetTextSize(0.04);
sum2->SetTextColor(4);
sum2->Draw("same");
sum3->SetTextSize(0.04);
sum3->SetTextColor(4);
sum3->Draw("same");
  }
  
//  h33->GetXaxis()->SetTitle("M(e^{+}e^{-}#gamma) [GeV]");
//  h33->GetYaxis()->SetTitleOffset(1.2);
//
//  h33->GetYaxis()->SetTitle("Entries / 6 MeV");
  
  
  //can1->Print("/Volumes/Mac_Storage/Physics_Papers/Annual_Review_Paper/ODU_Group_Pres/PLOTS_for_Review/New_Mass_Plots/etaprime.pdf")

}

