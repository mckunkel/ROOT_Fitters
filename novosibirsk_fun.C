#include "TComplex.h"

//This is in the root directory and as of May 15 2012, is the final version

Double_t Novosibirsk(Double_t *x, Double_t *par)
{
  //From RooNovosibirsk
  double qa=0,qb=0,qc=0,qx=0,qy=0;
  double tail = par[2];
  double width = par[1];
  double peak = par[0];
  
  //  if(TMath::Abs(tail) < 1.e-7)
  //    qc = 0.5*TMath::Power(((x[0]-peak)/width),2);
  //  else {
  qa = tail*sqrt(log(4.));
  qb = sinh(qa)/qa;
  qx = (x[0]-peak)/width*qb;
  qy = 1.+tail*qx;
  
  //    //---- Cutting curve from right side
  //
  //    if( qy > 1.E-7)
  //      qc = 0.5*(TMath::Power((log(qy)/tail),2) + tail*tail);
  //    else
  //      qc = 15.0;
  //  }
  
  //---- Normalize the result
  
  return par[3]*exp(-qc);
  
  //return ret;
}


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
    if (order == 3) return p[0] + x*p[1] + x*x*p[2];
    if (order == 4) return p[0] + x*p[1] + x*x*p[2] + x*x*x*p[3];
    if (order == 5) return p[0] + x*p[1] + x*x*p[2] + x*x*x*p[3] + x*x*x*x*p[4];

    
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


void novosibirsk_fun(TH1D *h33 , Double_t low, Double_t high, Double_t initialPar, Double_t width, Double_t tail, Double_t factor, Int_t draw_opt){

//  double qa=0,qb=0,qc=0,qx=0,qy=0;
//  qa = tail*sqrt(log(4.));
//  qb = sinh(qa)/qa;
//  qx = (x[0]-peak)/width*qb;
//  qy = 1.+tail*qx;
//  qc = 0.5*(TMath::Power((log(qy)/tail),2) + tail*tail);
  
  //[0] -> Normalization
  //[1] -> lambda
  //[2] -> x0  i.e. peak
  //[3] -> sigma

  TF1 *nov_fun = new TF1("nov_fun","[0]*exp(-0.5*pow(log(1.+ [1]*(x-[2])/[3]*(sinh([1]*log(sqrt(4)))/[1]*sqrt(log(4)))),2)/([1]*[1]) +([1]*[1]))");
  Int_t n =4;
  Chebyshev * cheb = new Chebyshev(n,low,high);
  TF1 * f1 = new TF1("f1",cheb,low,high,n+1,"Chebyshev");
  
  double nEnt = h33->GetEntries();
  TF1 *fitter = new TF1("fitter","[0]*exp(-0.5*pow(log(1.+ [1]*(x-[2])/[3]*(sinh([1]*log(sqrt(4)))/[1]*sqrt(log(4)))),2)/([1]*[1]) +([1]*[1])) +[4] +[5]*x + [6]*x*x ",low,high);//  [4] +[5]*x + [6]*x*x + [6]*x*x + [7]*x*x*x
fitter->SetParameters(100, tail, initialPar, width,3); fitter->SetParLimits(0,10,nEnt*5.); fitter->SetParLimits(2,initialPar-width,initialPar+width);
  //fitter->SetParLimits(4,-10,100);
  //fitter->SetParLimits(5,-10,100);
  //fitter->SetParLimits(6,-10,100);
h33->Fit("fitter","REM","same");
h33->Draw("E");
  
  
TF1 *backFcn = new TF1("backFcn", "pol2",low,high);
TF1 *signalFcn = new TF1("signalFcn", "nov_fun",low,high);
signalFcn->SetLineColor(2);
signalFcn->SetLineWidth(2);
Double_t par[8];
fitter->GetParameters(par);
backFcn->SetParameters(&par[4]);
backFcn->SetLineStyle(2);
backFcn->SetLineColor(6);
backFcn->SetLineWidth(1);


signalFcn->SetParameters(par);
signalFcn->SetLineStyle(2);
signalFcn->SetLineColor(4);
signalFcn->SetLineWidth(1);
  
  
  
  cout<<par[0]<<"  "<<par[1]<<"  "<<par[2]<<"  "<<par[3]<<"  "<<endl;
  double newlow = 0.0;
  for (double interval = low; interval<high; interval+=(high-low)/10000.) {
    double value = signalFcn->Eval(interval);
    if (TMath::IsNaN(value)) {
      newlow = interval;

    }
    //cout<<value<<"   "<<interval<<"  "<<low<<"  "<<high<<endl;
    //cout<<"   "<<interval<<"  "<<low<<"  "<<high<<endl;
    
  }
  double newbcklow = 0.0;
  for (double interval = low; interval<high; interval+=(high-low)/10000.) {
    double value = backFcn->Eval(interval);
    if (TMath::IsNaN(value) || value < 0.0) {
      newbcklow = interval;
      
    }
    //cout<<value<<"   "<<interval<<"  "<<low<<"  "<<high<<endl;
    //cout<<"   "<<interval<<"  "<<low<<"  "<<high<<endl;
    
  }
  
  signalFcn->SetRange(newlow,high);
  backFcn->SetRange(newbcklow,high);

backFcn->Draw("same");
signalFcn->Draw("same");

  
  TF1 *sigIntegralFcn = new TF1("sigIntegralFcn", "nov_fun",par[2]-factor*par[3],par[2]+factor*par[3]/(2*par[1]));
  sigIntegralFcn->SetParameters(par);
  sigIntegralFcn->SetFillStyle(1001);
  sigIntegralFcn->SetFillColor(kGreen);
  sigIntegralFcn->SetLineColor(kCyan);

  sigIntegralFcn->Draw("same");
  
  
  TF1 *backIntegralFcn = new TF1("backIntegralFcn", "backFcn",newbcklow,par[2]+factor*par[3]/par[1]);
  backIntegralFcn->SetParameters(&par[4]);
  backIntegralFcn->SetFillStyle(3010);
  backIntegralFcn->SetFillColor(kBlue);
  backIntegralFcn->SetLineColor(kBlue);

  backIntegralFcn->Draw("same");
//  TH1 * h1 = new TH1D(*((TH1D*) signalFcn->GetHistogram()) );
//  for (int i = 1; i <= h1->GetNbinsX(); ++i){
//    //h1->SetBinContent(i, signalFcn->Integral( h1->GetXaxis()->GetBinLowEdge(i), h1->GetXaxis()->GetBinUpEdge(i) )/h1->GetBinWidth(i) );
//    h1->SetBinContent(i,signalFcn->Integral(par[2]-factor*par[3],par[2]+factor*par[3]));
//}
//  h1->Draw();
  //signalFcn->Draw();
  //backFcn->Draw();

Double_t Intg = abs(signalFcn->Integral(par[2]-factor*par[3],par[2]+factor*par[3]/(2*par[1])));
Double_t Intb = abs(backFcn->Integral(newbcklow,par[2]+factor*par[3]/(2*par[1])));


Double_t binw = h33->GetBinWidth(1);
Int_t yield = Intg/binw;
Int_t bckgd = Intb/binw;
//Double_t ratio = double(yield)/TMath::Sqrt(double(bckgd));
  Double_t ratio = double(yield)/(double(yield+bckgd));

cout << yield << "\t" << ratio << endl;
TAxis *x=h33->GetXaxis();
TAxis *y=h33->GetYaxis();

Double_t startx=x->GetXmin()+0.75*(x->GetXmax()-x->GetXmin());
  
  Double_t starty0=0.35*h33->GetMaximum();
  Double_t starty1=0.45*h33->GetMaximum();
  Double_t starty2=0.55*h33->GetMaximum();
  Double_t starty3=0.65*h33->GetMaximum();
  Double_t starty4=0.75*h33->GetMaximum();
  Double_t starty5=0.85*h33->GetMaximum();


  double meanError = fitter->GetParError(2);
  double sigmaError = fitter->GetParError(3);
  if (draw_opt ==1) {

    TLatex *sum = new TLatex(startx*0.93, starty5,Form("Yield: %i",yield));
    TLatex *sum12 = new TLatex(startx*0.93, starty4,Form("Background: %i",bckgd));
    TLatex *sum0=new TLatex(startx*0.93, starty3,Form("Range: #pm %2.1f #sigma",factor));
    TLatex *sum2=new TLatex(startx*0.93, starty2,Form("Mean:%4.4f #pm %.4f GeV",par[2], meanError));
    TLatex *sum3=new TLatex(startx*0.93, starty1,Form("#sigma:%5.4f #pm %.4f GeV",par[3], sigmaError));
    TLatex *ra = new TLatex(startx*0.93, starty0,Form("#frac{S}{S+B}= %.3f", ratio));
    
    sum->SetTextSize(0.04);
    sum->SetTextColor(2);
    sum->Draw("same");
    sum12->SetTextSize(0.04);
    sum12->SetTextColor(6);
    sum12->Draw("same");
    
    ra->SetTextSize(0.04);
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


}
