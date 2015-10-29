
//This is in the root directory and as of May 15 2012, is the final version

void MKGausPol2(TH1D *h33 , Double_t low, Double_t high, Double_t p0, Double_t p1, Double_t p2, Double_t initialPar, Double_t width, Double_t factor, Int_t draw_opt){

  double nEnt = h33->GetEntries();
  TF1 *fitter = new TF1("fitter","gaus + [3] + [4]*x + [5]*x*x",low,high);
fitter->SetParameters(100, initialPar, width, p0, p1, p2); fitter->SetParLimits(0,10,nEnt); fitter->SetParLimits(1,initialPar-width,initialPar+width);
h33->Fit("fitter","REM","same");
  //h33->Draw("E");
TF1 *backFcn = new TF1("backFcn", "pol2",low,high);
TF1 *signalFcn = new TF1("signalFcn", "gaus",low,high);
signalFcn->SetLineColor(2);
signalFcn->SetLineWidth(2);
Double_t par[6];
fitter->GetParameters(par);
backFcn->SetParameters(&par[3]);
backFcn->SetLineStyle(2);
backFcn->SetLineColor(6);
backFcn->SetLineWidth(1);


signalFcn->SetParameters(par);
signalFcn->SetLineStyle(2);
signalFcn->SetLineColor(4);
signalFcn->SetLineWidth(1);
backFcn->Draw("same");
signalFcn->Draw("same");
Double_t Intg = signalFcn->Integral(par[1]-factor*par[2],par[1]+factor*par[2]);
Double_t Intb = backFcn->Integral(par[1]-factor*par[2],par[1]+factor*par[2]);


Double_t binw = h33->GetBinWidth(1);
Int_t yield = Intg/binw;
Int_t bckgd = Intb/binw;
Double_t ratio = double(yield)/TMath::Sqrt(double(bckgd));
cout << yield << "\t" << ratio << endl;
TAxis *x=h33->GetXaxis();
TAxis *y=h33->GetYaxis();

Double_t startx=x->GetXmin()+0.85*(x->GetXmax()-x->GetXmin());
  
  Double_t starty0=0.35*h33->GetMaximum();
  Double_t starty1=0.45*h33->GetMaximum();
  Double_t starty2=0.55*h33->GetMaximum();
  Double_t starty3=0.65*h33->GetMaximum();
  Double_t starty4=0.75*h33->GetMaximum();
  Double_t starty5=0.85*h33->GetMaximum();


  double meanError = fitter->GetParError(1);
  double sigmaError = fitter->GetParError(2);
  if (draw_opt ==1) {

    TLatex *sum = new TLatex(startx*0.93, starty5,Form("Yield: %i",yield));
    TLatex *sum12 = new TLatex(startx*0.93, starty4,Form("Background: %i",bckgd));
    TLatex *sum0=new TLatex(startx*0.93, starty3,Form("Range: #pm %2.1f #sigma",factor));
    TLatex *sum2=new TLatex(startx*0.93, starty2,Form("Mean:%4.4f #pm %.4f GeV",par[1], meanError));
    TLatex *sum3=new TLatex(startx*0.93, starty1,Form("#sigma:%5.4f #pm %.4f GeV",par[2], sigmaError));
    TLatex *ra = new TLatex(startx*0.93, starty0,Form("#frac{S}{#sqrt{B}}= %.1f", ratio));
    
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
  
//  h33->GetXaxis()->SetTitle("M(e^{+}e^{-}#gamma) [GeV]");
//  h33->GetYaxis()->SetTitleOffset(1.2);
//
//  h33->GetYaxis()->SetTitle("Entries / 6 MeV");
  
  
  //can1->Print("/Volumes/Mac_Storage/Physics_Papers/Annual_Review_Paper/ODU_Group_Pres/PLOTS_for_Review/New_Mass_Plots/etaprime.pdf")

}
