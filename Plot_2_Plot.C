void Plot_2_Plot(TH1D *gr1, TH1D *gr2)
{
  
  gROOT->SetStyle("Bold");
  gStyle->SetOptStat(kFALSE);

  
  gr1->SetDirectory(0);
  gr1->SetMarkerStyle(4);
  gr1->SetTitle("");
  
  gr2->SetDirectory(0);
  gr2->SetMarkerStyle(23);
  
  
  TCanvas *c1 = new TCanvas("c1","trans pad",200,10,1050,700);
  
  gr1->Draw("P9 E1");
  c1->Update();
  //scale gr2 to the pad coordinates
  Float_t rightmax = 1.1*gr2->GetMaximum();
  Float_t scale = gPad->GetUymax()/rightmax;
  gr2->SetLineColor(kRed);
  gr2->Scale(scale);
  gr2->Draw("same P9 E1");

  
  //draw an axis on the right side
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->Draw();


}

