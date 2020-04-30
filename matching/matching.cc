#include "matching.hh"
#include <string>
#include "TGraph.h"
#include "TRandom3.h"
#include <iostream>

using namespace std;

double getT0(std::string method, TH1D *hdet1, TH1D *hdet2, double &minchi2, double rebinwidth, double chi2step, double windowmin, double windowmax, double shift_min, double shift_max, double shift_step){

  minchi2=1e20; 
  double Tfit;

  TGraph *gchi2dt;
  if (debug >= 1) {
      gchi2dt = new TGraph();
      gchi2dt->SetName("gchi2dt");
  }

  double nev_det1, nev_det2, err2_det1, err2_det2;
  double chi2_value;
  double mean_det1, sq_det1;
  double mean_det2, sq_det2;
  //double nsteps;

    for(double dt=shift_min;dt<=shift_max; dt+=shift_step){
      chi2_value = 0;
      mean_det1  = 0;
      sq_det1    = 0;
      mean_det2  = 0;
      sq_det2    = 0;
      //nsteps = 0;
      for (double b = windowmin; b<=windowmax; b+=chi2step){

        //double nev_det1 = hdet1->Interpolate(b-dt);
        //double nev_det2 = hdet2->Interpolate(b);
        //double nev_det1 = VInterpolate(hdet1,b-dt);
        //double nev_det2 = VInterpolate(hdet2,b);
        nev_det1 = VRebin(hdet1,b-dt,rebinwidth*10, err2_det1);
        nev_det2 = VRebin(hdet2,b,rebinwidth*10, err2_det2);

        //WARNING case with 0 background in both detectors needs to be threated somehow
        if (method=="chi2") {
            if (err2_det1+err2_det2 > 0)
                chi2_value += pow(nev_det1-nev_det2,2)/(err2_det1+err2_det2);
        } else if (method=="corr") {
            chi2_value -= nev_det1*nev_det2; //we do it negatively so the function has minimum
        } else if (method=="corrslide"){
            chi2_value += nev_det1*nev_det2; //we do it negatively so the function has minimum
            mean_det1 += nev_det1;
            sq_det1 += nev_det1*nev_det1;
            mean_det2 += nev_det2;
            sq_det2 += nev_det2*nev_det2;
            //nsteps += 1;
        }
        else {
            cout << "ERROR: the method " << method << " is unknown" << endl;
        }
      }
      if (method=="corrslide"){
          double nsteps = int((windowmax-windowmin)/chi2step)+1;
          sq_det1 = sqrt((sq_det1-mean_det1*mean_det1/nsteps)/double(nsteps-1));
          sq_det2 = sqrt((sq_det2-mean_det2*mean_det2/nsteps)/double(nsteps-1));
          mean_det1 = mean_det1/nsteps;
          mean_det2 = mean_det2/nsteps;
          chi2_value = -(chi2_value-mean_det1*mean_det2*nsteps)/sq_det1/sq_det2;
      }
      
      if (debug >= 1) gchi2dt->SetPoint(gchi2dt->GetN(),dt,chi2_value);
      if (chi2_value<minchi2) {
          minchi2=chi2_value;
          Tfit=dt;
      }
    }

    if (debug >= 1) gchi2dt->Write();
    return Tfit;
}

TH1D *corrhist(TH1D *hdet, int bgopt, int sgopt, double Meff, double Rbg, double tleft, double tright, double rebinwidth){
    string hname = hdet->GetName();
    hname += "_corr";
    hdet = (TH1D *)hdet->Clone(hname.c_str());

    double delta;

//BEGIN background subraction
    double Rbgcorr = 0;
    if (bgopt == 0) {
        Rbgcorr = 0;
    } else if (bgopt == 1) {
        Rbgcorr = Rbg*hdet->GetBinWidth(0)/1000.;
    } else if (bgopt == 2) {
        Rbgcorr = hdet->Integral(1,1000)/1000.; //first 1000 bins have no signal
    } else if (bgopt == 3) {
        Rbgcorr = hdet->Integral(hdet->FindBin(tleft),hdet->FindBin(tright))/(hdet->FindBin(tright)-hdet->FindBin(tleft)+1); //mean of the lightcurve for zero-normalized cross-correlation
    } else cerr << "ERROR: option bgopt = " << bgopt << " is not supported" << endl;
    
    for (int ii = 1; ii <= hdet->GetNbinsX(); ii++) {
        hdet->SetBinError(ii,hdet->GetBinError(ii));
        hdet->SetBinContent(ii,hdet->GetBinContent(ii)-Rbgcorr); //we keep stat error since SetBinError tried to rewrite it if it was never set
    }

//END background subtraction
    bool samecorrection = true; //to use same systematics for signal and background

//BEGIN signal level normalisation
     double scale = 1;
     if (sgopt == 0) {
         scale = 1;
     } else if (sgopt == 1) {
         scale = Meff;
         //scale = Meff*111.36; //the value 111.36 is the expected rate in 1 Mt detector (model dependent)
     } else if (sgopt == 2) {
         scale = hdet->Integral(hdet->FindBin(tleft),hdet->FindBin(tright));
     } else if (sgopt == 3) {
        if (bgopt != 3) cout << "WARNING: you are using wrong bgopt=" << bgopt << " for zero-normalization cross-correlation" << endl;
        scale = 0;
        TH1D *hdetrebin = (TH1D *)hdet->Clone("hdetbin"); //we need this rebin histogram because sigma should be calculated for the proper bin lenght
        hdetrebin->Rebin(rebinwidth*10);
        for (int ibin = hdetrebin->FindBin(tleft); ibin <= hdetrebin->FindBin(tright); ibin++) {
            scale += pow(hdetrebin->GetBinContent(ibin),2);  //WARNING  - it is n_i-<n> only if bgopt = 3 is used
        }
        scale = sqrt(scale/double(hdetrebin->FindBin(tright)-hdetrebin->FindBin(tleft)));
        hdetrebin->Delete();
     } else cerr << "ERROR: option sgopt = " << sgopt << " is not supported" << endl;
     if (debug >= 1) cout << "debug: Rbgcorr " << Rbgcorr << " bgopt " << bgopt << " scale: " << scale << " Meff " << Meff << " sgopt " << sgopt << endl;
     hdet->Scale(1./scale);
//END signal level normalisation

     return hdet;
}
