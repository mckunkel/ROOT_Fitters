#include "TComplex.h"
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
#include "TLegend.h"
#include "TPaveStats.h"
#include "TList.h"
#include "TMatrixT.h"
#include "TLatex.h"

//This is in the root directory and as of May 15 2012, is the final version

class Chebyshev {
public:
  Chebyshev(int n, double xmin, double xmax) :
  fA(xmin), fB(xmax),
  fT(std::vector<double>(n) )  {}
  
  double operator() (const double * xx, const double *p) {
    double x = (xx[0] - fA -fB)/(fB-fA);
    int order = fT.size();
    if (order == 1) return p[0];
    if (order == 2) return p[0] + x*p[1];
    // build the polynomials
    fT[0] = 1;
    fT[1] = x;
    for (int i = 1; i< order; ++i) {
      fT[i+1] =  2 *x * fT[i] - fT[i-1];
    }
    double sum = p[0]*fT[0];
    for (int i = 1; i<= order; ++i) {
      sum += p[i] * fT[i];
    }
    return sum;
  }
  
private:
  double fA;
  double fB;
  std::vector<double> fT; // polynomial
  std::vector<double> fC; // coefficients
};


double CBshape (double *x, double *p ) {
  //x[0] x
  //p[0] mean
  //p[1] sigma
  //p[2] alpha
  //p[3] n
  //p[4] normalization
  Double_t t = (x[0]-p[0])/p[1];
  if (p[2] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)p[2]);
  
  if (t >= -absAlpha) {
    return p[4]*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(p[3]/absAlpha,p[3])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= p[3]/absAlpha - absAlpha;
    
    return p[4]*a/TMath::Power(b - t, p[3]);
  }
  
}


// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];//+ par[3]*x[0]*x[0]*x[0]; // + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0]
  

  
  
}
// Lorentzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
  return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10,(x[0]-par[2])*(x[0]-par[2])+ .25*par[1]*par[1]);
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  //return background(x,par) + lorentzianPeak(x,&par[3]);
//  Int_t n =4;
//  Chebyshev * cheb = new Chebyshev(n,low,high);
//  TF1 * f1 = new TF1("f1",cheb,low,high,n+1,"Chebyshev");
  return CBshape(x,par) + background(x,&par[5]);
  //return CBshape(x,par) + cheb(x,&par[5]);

}




void CBShape(TH1 *h33 , Double_t low, Double_t high, Double_t initialPar, Double_t width, Double_t alpha, Double_t n, Double_t factor, Int_t draw_opt){
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  
  //void CBShape(){

  //[0] -> m0
  //[1] -> width
  //[2] -> alpha
  //[3] -> n
  //[4] -> Normalization

//  Double_t low = -0.1;
//  Double_t high = 0.1;
//  TF1 * f1 = new TF1("f1",CBshape,low,high,5);
//  f1->SetParameters(0.0,0.01,1.,2,1000);
//  TF1 * f2 = new TF1("f2",CBshape,low,high,5);
//  f2->SetParameters(0.0,0.01,-1.,2,1000);  //these parameters give tail to positive side if mean is 0
//  f2->SetLineColor(kGreen);
//  
//  f1->Draw();
//  f2->Draw("same");
  
  
  
  TF1 *CB = new TF1("CB",CBshape,low, high, 5);
  //CB->SetParameters(initialPar,width,alpha,n);

//  double nEnt = h33->GetEntries();
//  
//  cout<<nEnt<<"###################### ENTRIES ##############"<<endl;
//  double nEnt_val;
//  if (nEnt >100000) {
//    nEnt_val = nEnt*3;
//  }
//  else if(nEnt < 100000 && nEnt > 50000){
//    nEnt_val = nEnt*3;
//  }
//  else if(nEnt < 50000 && nEnt > 25000){
//    nEnt_val = nEnt*2;
//  }
//  else if(nEnt < 25000){
//    nEnt_val = nEnt;
//  }
//  else{nEnt_val = 0.;}
  
  
  double nEnt = h33->GetMaximum();
  double nEnt_val = 10*nEnt;

  double newlow = low;//0.0;
  for (int yval = 1 ; yval<=h33->GetNbinsX()/2; yval++) {
    if (h33->GetBinContent(yval) == 0. && h33->GetBinContent(yval+1) != 0.){ // && hslice->GetBinContent(yval+1) <0.025
      cout<< h33->GetBinCenter(yval-1)<<"  LIMITS "<<endl;
      newlow = h33->GetBinCenter(yval-1);
      
      break;
    }
  }

  
  Double_t c_val = n/fabs(alpha)*exp(-0.5*fabs(alpha)*fabs(alpha));
  Double_t d_val = sqrt(TMath::Pi()/2.)*(1.+TMath::Erf(fabs(alpha)/sqrt(2)));
  //nEnt_val = nEnt_val/(width*(c_val*d_val));
  
  TF1 *fitter = new TF1("fitter",fitFunction,newlow,high,9);//  [4] +[5]*x + [6]*x*x + [6]*x*x + [7]*x*x*x
  //TF1 *fitter = new TF1("fitter",CBshape,newlow,high,5);//  [4] +[5]*x + [6]*x*x + [6]*x*x + [7]*x*x*x
  
  fitter->SetParameters(initialPar,width,alpha,n,nEnt);  fitter->SetParLimits(0,initialPar-2.*width,initialPar+2.*width); fitter->SetParLimits(4,nEnt_val/10.,nEnt_val);
  fitter->SetParLimits(3,-1000,1000);
  fitter->SetParLimits(6,1,100);
  //fitter->SetParLimits(7,1,100);
  //fitter->SetParLimits(8,1,100);
  //fitter->SetParLimits(9,3,100);

  
  h33->Fit("fitter","REM","same");
  h33->Draw("E");
  
  
  TF1 *backFcn = new TF1("backFcn", background,newlow,high,3);
  TF1 *signalFcn = new TF1("signalFcn", CBshape,newlow,high,5);
  signalFcn->SetLineColor(2);
  signalFcn->SetLineWidth(2);
  Double_t par[8];
  fitter->GetParameters(par);
  backFcn->SetParameters(&par[5]);
  backFcn->SetLineStyle(2);
  backFcn->SetLineColor(6);
  backFcn->SetLineWidth(1);
  
  
  signalFcn->SetParameters(par);
  signalFcn->SetLineStyle(2);
  signalFcn->SetLineColor(4);
  signalFcn->SetLineWidth(1);
  
  
  //backFcn->Draw();
  //signalFcn->Draw("same");
  
  
  Double_t Intg = abs(signalFcn->Integral(par[0]-factor*par[1],par[0]+factor*par[1]));
  Double_t Intb = abs(backFcn->Integral(par[0]-factor*par[1],par[0]+factor*par[1]));
  
  
  cout<<"$$$$$$$$$$$$$$$$$$$$"<<par[0]-factor*par[1]<<"  "<<par[0]+factor*par[1]<<"  "<<par[0]<<"  "<<par[1]<<endl;
  //h->Integral(h->FindBin(-1), h->FindBin(0) - 1)
  
  TF1 *sigIntegralFcn = new TF1("sigIntegralFcn", CBshape,newlow,high,5);
  sigIntegralFcn->SetParameters(par);
  //sigIntegralFcn->SetFillStyle(3010);
  //sigIntegralFcn->SetFillColor(kGreen);
  //sigIntegralFcn->SetLineWidth(0.01);


  sigIntegralFcn->Draw("same");
  backFcn->Draw("same");

  
  //Double_t Intg = abs(signalFcn->Integral(newlow,high));
  double num_Intg = 0.0;
  double num_Intb = 0.0;

  double frac = 100000.;
  double divisor = frac*h33->GetNbinsX(); //1000.;
  Double_t binw = h33->GetBinWidth(1);
  
  //for (double interval = low; interval<=high; interval+=(high-low)/divisor) {
  for (double interval = par[0]-factor*par[1]; interval<=par[0]+factor*par[1]; interval+=binw) {

    double s_value = signalFcn->Eval(interval);
    double b_value = backFcn->Eval(interval);

    num_Intg = num_Intg + s_value;
    num_Intb = num_Intb + b_value;

    //cout<<value<<"   "<<interval<<"  "<<num_Intg<<"  "<<divisor<<endl;
    //cout<<"   "<<interval<<"  "<<low<<"  "<<high<<endl;
    
  }
  cout<<ceil(num_Intg)<<"################# NUM INTG"<<endl;
  cout<<abs(ceil(num_Intb))<<"################# NUM INTB"<<endl;
  cout<<binw<<"################# BINW"<<endl;

  //Double_t Intb = abs(backFcn->Integral(newlow,high));
  

  //Int_t yield = Intg/binw;
  //Int_t bckgd = Intb/binw;
  Int_t yield = ceil(num_Intg);
  Int_t bckgd = abs(ceil(num_Intb));
  cout<<yield+bckgd<<"################# INTG +  BINW"<<endl;

  //Double_t ratio = double(yield)/TMath::Sqrt(double(bckgd));
  Double_t ratio = double(yield)/(double(yield+bckgd));
  
  cout << yield << "\t" << ratio << endl;
  TAxis *x=h33->GetXaxis();
  TAxis *y=h33->GetYaxis();
  
  Double_t startx=x->GetXmin()+0.55*(x->GetXmax()-x->GetXmin());
  
  Double_t starty0=0.35*h33->GetMaximum();
  Double_t starty1=0.5*h33->GetMaximum();
  Double_t starty2=0.6*h33->GetMaximum();
  Double_t starty3=0.7*h33->GetMaximum();
  Double_t starty4=0.8*h33->GetMaximum();
  Double_t starty5=0.9*h33->GetMaximum();
  
  
  double meanError = fitter->GetParError(0);
  double sigmaError = fitter->GetParError(1);
  if (draw_opt ==1) {
    
    TLatex *sum = new TLatex(startx*0.73, starty5,Form("Yield: %i",yield));
    TLatex *sum12 = new TLatex(startx*0.73, starty4,Form("Background: %i",bckgd));
    TLatex *sum0=new TLatex(startx*0.73, starty3,Form("Range: #pm %2.1f #sigma",factor));
    TLatex *sum2=new TLatex(startx*0.73, starty2,Form("Mean:%4.4f #pm %.2e GeV^{2}",par[0], meanError));
    TLatex *sum3=new TLatex(startx*0.73, starty1,Form("#sigma:%5.4f #pm %.2e GeV^{2}",par[1], sigmaError));
    TLatex *ra = new TLatex(startx*0.73, starty0,Form("#frac{S}{S+B}= %.3e", ratio));
    
    sum->SetTextSize(0.05);
    sum->SetTextColor(2);
    sum->Draw("same");
    sum12->SetTextSize(0.05);
    sum12->SetTextColor(6);
    sum12->Draw("same");
    
    ra->SetTextSize(0.05);
    ra->Draw("same");
    
    
    sum0->SetTextSize(0.05);
    sum0->SetTextColor(2);
    sum0->Draw("same");
    sum2->SetTextSize(0.05);
    sum2->SetTextColor(4);
    sum2->Draw("same");
    sum3->SetTextSize(0.05);
    sum3->SetTextColor(4);
    sum3->Draw("same");
    
    
  }
  
  
}
