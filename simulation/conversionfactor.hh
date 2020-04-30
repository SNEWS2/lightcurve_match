#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>

extern TGraph *gsigmatab;
Double_t funcsigma(Double_t *x, Double_t *par);
Double_t funcnudistrwrapper(Double_t *x, Double_t *par);
Double_t funcnuintdistrwrapper(Double_t *x, Double_t *par);
Double_t funcnudistr(Double_t Enu, Double_t norm, Double_t alpha, Double_t avEnu);
double conversionfactor();
