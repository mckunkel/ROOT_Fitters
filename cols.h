#include <iostream>
#include <cmath>

using namespace std;

void printcols()
{
  cout << "usage: colpalette(a, l, i, n)" << endl
  << "    a = [-1.0, 1.0]" << endl
  << "      (scaling factor for colors. Default: 0.0)" << endl
  << "    l = [0.0, 1.0]" << endl
  << "      (lightness factor for color scale. Default: 1.0)" << endl
  << "    i = [true|false]" << endl
  << "      (invert color scale. Default: false)" << endl
  << "    n = [1, 256]" << endl
  << "      (number of color divisions. Default: 256)" << endl;
  cout << "available color palettes:"
  " ucla"
  " jet"
  " flame"
  " rainbow"
  " rainbow2"
  " pastel"
  " bone"
  " land"
  " land2"
  " earth"
  " bear"
  << endl;
}

void adjust_stops(double* stops, int n, double a = 0.0)
{
  /// force the limits on a:
  /// a = [-1.0, 1.0]
  if (a < -1.0) { a = -1.0; }
  else if (a >  1.0) { a =  1.0; }
  
  double x;
  for (int i=0; i<n; i++)
  {
    x = stops[i];
    if (a > 0)
    {
      stops[i] = (x - 1.) * exp(-10. *  a * x) + 1.;
    }
    else
    {
      stops[i] = x * exp(-10. * a * (x - 1.));
    }
  }
}

void lighten(double* red, double* green, double* blue, const int& n, const double& l)
{
  if (0.0 < l && l < 1.0)
  {
    int i;
    for(i=0; i<n; i++)
    {
      red[i] += (1.0 - red[i]) * l;
      green[i] += (1.0 - green[i]) * l;
      blue[i] += (1.0 - blue[i]) * l;
    }
  }
}

void invertf(double* stops, double* red, double* green, double* blue, const int& n)
{
  double invstops[n];
  double invred[n];
  double invgreen[n];
  double invblue[n];
  int i;
  for (i=0; i<n; i++)
  {
    invstops[i] = 1.0 - stops[n-i-1];
    invred[i] = red[n-i-1];
    invgreen[i] = green[n-i-1];
    invblue[i] = blue[n-i-1];
  }
  for (i=0; i<n; i++)
  {
    stops[i] = invstops[i];
    red[i] = invred[i];
    green[i] = invgreen[i];
    blue[i] = invblue[i];
  }
}

void ucla(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 4;
  
  double stops[nstops] = { 0.00, 0.50, 0.50, 1.00 };
  double red[nstops]   = { 1.00, 0.00, 1.00, 0.00 };
  double green[nstops] = { 1.00, 0.00, 1.00, 0.00 };
  double blue[nstops]  = { 1.00, 0.75, 0.75, 0.50 };
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void jet(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 5;
  
  double stops[nstops] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  double red[nstops]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  double green[nstops] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  double blue[nstops]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void flame(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 5;
  
  double stops[nstops] = { 0.00, 0.10, 0.40, 0.70, 1.00 };
  double red[nstops]   = { 1.00, 0.00, 0.87, 1.00, 0.51 };
  double green[nstops] = { 1.00, 0.70, 1.00, 0.20, 0.00 };
  double blue[nstops]  = { 1.00, 1.00, 0.12, 0.00, 0.00 };
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void rainbow(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 6;
  
  double stops[nstops] = { 0.00, 0.20, 0.40, 0.60, 0.80, 1.00 };
  double red[nstops]   = { 0.00, 0.00, 0.00, 0.80, 1.00, 0.50 };
  double green[nstops] = { 0.00, 0.80, 0.90, 1.00, 0.20, 0.00 };
  double blue[nstops]  = { 0.50, 1.00, 0.50, 0.10, 0.00, 0.00 };
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void rainbow2(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 5;
  
  double stops[nstops] = { 0.00, 0.20, 0.40, 0.65, 1.00 };
  double red[nstops]   = { 0.00, 0.00, 0.00, 1.00, 1.00 };
  double green[nstops] = { 0.00, 1.00, 1.00, 1.00, 0.00 };
  double blue[nstops]  = { 1.00, 1.00, 0.00, 0.00, 0.00 };
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void pastel(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 5;
  
  double stops[nstops] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  double red[nstops]   = { 0.60, 0.25, 0.90, 1.00, 1.00 };
  double green[nstops] = { 0.40, 0.65, 1.00, 0.60, 0.40 };
  double blue[nstops]  = { 0.70, 0.00, 0.00, 0.50, 0.30 };
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void bone(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 3;
  
  double stops[nstops] = { 0.00, 0.30, 1.00 };
  double red[nstops]   = { 1.00, 0.50, 0.00 };
  double green[nstops] = { 1.00, 0.70, 0.00 };
  double blue[nstops]  = { 1.00, 0.60, 0.00 };
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void land(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 5;
  
  double stops[nstops] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  double red[nstops]   = { 0.00, 0.50, 1.00, 0.75, 0.50 };
  double green[nstops] = { 0.50, 0.75, 1.00, 0.20, 0.20 };
  double blue[nstops]  = { 0.00, 0.10, 0.60, 0.10, 0.00 };
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void land2(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 3;
  
  double stops[nstops] = { 0.00, 0.33, 1.00 };
  double red[nstops]   = { 1.00, 0.00, 0.50 };
  double green[nstops] = { 1.00, 0.50, 0.10 };
  double blue[nstops]  = { 0.90, 0.00, 0.00 };
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void earth(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 4;
  
  double stops[nstops] = { 0.00, 0.33, 0.67, 1.00 };
  double red[nstops]   = { 0.05, 1.00, 0.10, 0.70 };
  double green[nstops] = { 0.10, 1.00, 0.40, 0.20 };
  double blue[nstops]  = { 0.50, 0.80, 0.20, 0.00 };
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void bear(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 3;
  
  double stops[nstops] = { 0.00, 0.50, 1.00 };
  double red[nstops]   = { 1.00, 0.50, 0.00 };
  double green[nstops] = { 1.00, 0.30, 0.00 };
  double blue[nstops]  = { 0.90, 0.00, 0.00 };
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void obiwan(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 6;
  
  Double_t stops[nstops] = { 0.00, 0.1, 0.4, 0.6, 0.70, 1.00 };
  Double_t red[nstops]   = { 0.00, 0.10, 0.2, 0.10, 0.10, 1.0 };
  Double_t green[nstops] = { 0.00, 0.10, 0.4, 0.50, 0.40, 1.0 };
  Double_t blue[nstops]  = { 0.00, 0.20, 0.5, 0.8,  1.00, 1.0 };
  
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}


void luke(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 6;
  
  Double_t stops[nstops] = { 0.00, 0.1, 0.4, 0.6, 0.70, 1.00 };
  Double_t red[nstops]   = { 0.00, 0.10, 0.2, 0.10, 0.10, 1.0 };
  Double_t blue[nstops] = { 0.00, 0.10, 0.2, 0.30, 0.40, 1.0 };
  Double_t green[nstops]  = { 0.00, 0.20, 0.6, 0.8,  1.00, 1.0 };
  
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void vader(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 6;
  
  Double_t stops[nstops] = { 0.00, 0.1, 0.4, 0.6, 0.70, 1.00 };
  Double_t blue[nstops]   = { 0.00, 0.0, 0.1, 0.10, 0.10, 1.0 };
  Double_t green[nstops]   = { 0.00, 0.00, 0.0, 0.10, 0.10, 1.0 };
  Double_t red[nstops]  = { 0.00, 0.20, 0.5, 0.8,  1.00, 1.0 };
  
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}

void mace(double a = 0.0, double l = 0.0, bool invert = false, int n = 256)
{
  const int nstops = 6;
  
  Double_t stops[nstops] = { 0.00, 0.1, 0.4, 0.6, 0.70, 1.00 };
  Double_t blue[nstops]   = { 0.00, 0.20, 0.5, 0.8,  1.00, 1.0 };
  Double_t green[nstops]   = { 0.00, 0.1, 0.2, 0.10, 0.10, 1.0 };
  Double_t red[nstops]  = { 0.00, 0.20, 0.5, 0.8,  1.00, 1.0 };
  
  
  if (invert) { invertf(stops, red, green, blue, nstops); }
  lighten(red, green, blue, nstops, l);
  adjust_stops(stops, nstops, a);
  TColor::CreateGradientColorTable(nstops, stops, red, green, blue, n);
  gStyle->SetNumberContours(n);
}
