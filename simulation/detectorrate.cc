#include "TROOT.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TRandom3.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "conversionfactor.hh"

static int debug = 0;
Double_t funcrate(Double_t *x, Double_t *par) {
    if (x[0]>par[6])  
        return par[5]*exp(-pow((par[0]/(x[0]-par[6])),par[1])) / pow((1+pow(((x[0]-par[6])/par[2]),par[3])),par[4]/par[3]);
    else return 0;

}

using namespace std;

/*
 * Function that generates experimental detection rate curve
 * Inputs:
 * Meff - Detector effective mass in kT (77kT is the default value for ORCA)
 * Rbg background rate in detector in Hz    (500Hz*2070 is is the default value for ORCA)
 * tdelay of the signal respect to t=0 in s
 * binwidth is the time duration of the bin in ms
 * drawopt is to draw interactively a canvas with histograms
 * drawopt is to save detector light curve to the txt file
 * fillopt : 1 - FillRandom method, 2 - Poisson sampling in each bin. First method shold be efficient for lower number of total events while second for low number of bins
 * seed - set random seed
 */
TH1D* detectorrate(double Meff, double Rbg, double tdelay, double binwidth, bool drawopt, bool printopt, int fillopt, int seed){
    double tleft = -2; //time of the lightcurve histo start
    double tright = 8; //time of the lightcurve histo end
    int nbins = (tright-tleft)/binwidth*1000; //1000 is since time for histo is in s, but bin width is in ms
    
    if ((tright-tleft)/nbins*1000 != binwidth) cerr << "WARNING: bin width is not multple of 0.0001ms" << endl;

    /*
     * Conversion factor that is multiplied by Luminosity in 1e51 ergs per detector effective mas      s in Kt gives detector rate per s
     * obtained by running conversionfactor.cc
     */

    double cfactor;
#ifdef IFACTOR
     cfactor = conversionfactor();
     cout << "Conversion factor (I): " << cfactor << endl;
#else
     cfactor = 4.28725; //Conversion factor that is multiplied by Luminosity in 1e51 ergs per detector effective mass in Kt gives detector rate per s
#endif
    

    //parametrisation from "Constraints on neutrino masses from a galacticsupernova neutrino signal at present and future detectors"
    //Nuclear Physics B 731 (2005) 140–163
    //http://bibliotecadigital.udea.edu.co/dspace/bitstream/10495/8533/1/ZuluagaJorge_2005_ConstraintsNeutrinoMasses.pdf
    TF1 *fflux = new TF1("fflux","[5]*exp(-pow(([0]/x),[1])) / pow((1+pow((x/[2]),[3])),[4]/[3])",0,3);
    fflux->SetParNames("ta","na","tc","np","nc","norm");
    //fflux->SetParameters(0.035,2,0.15,8,0.80,1); //initial from Fig.1 of [1]
    //fflux->SetParameters(0.035,2,0.15,8,1.5,1); //the decay is forced so integral converges fast
    fflux->SetParameters(0.035,2,0.2,20,1.5,1);
    double totalnorm = 3e53/6.;
    fflux->SetParameter(5,totalnorm/fflux->Integral(0,1e9));
    fflux->SetTitle("luminosity;t [s]; L_{#nu} [erg/s]");

    TCanvas *crates;
    if (drawopt) {
        crates = new TCanvas();
        crates->Divide(2,2);
        crates->cd(1);
        TGaxis::SetMaxDigits(3);
        fflux->SetNpx(1000);
        fflux->Draw();
    }
    
    if (debug > 0) {
        cout << "Integral: (0,10)  " << fflux->Integral(0,10)  << " erg "<< endl;
        cout << "Integral: (0,600) " << fflux->Integral(0,1e2) << " erg "<< endl;
        cout << "Integral: (0,1e3) " << fflux->Integral(0,1e3) << " erg "<< endl;
        cout << "Integral: (0,1e4) " << fflux->Integral(0,1e4) << " erg "<< endl;
        cout << "Integral: (0,1e5) " << fflux->Integral(0,1e5) << " erg "<< endl;
        cout << "Integral: (0,1e6) " << fflux->Integral(0,1e6) << " erg "<< endl;
    }
        
    TF1 *frate = new TF1("frate",funcrate,tdelay,tdelay+0.6,7);
    totalnorm = 3e53/6./1e51*Meff*cfactor;  //total luminosity is reported to 1e51erg units
    frate->SetParNames("ta","na","tc","np","nc","norm","tdelay"); 
    frate->SetParameters(0.035,2,0.2,20,1.5,1,0);
    frate->SetParameter(5,totalnorm/frate->Integral(0,1e9));
    frate->SetParameter(6,tdelay);
    frate->SetTitle("event rate;t [s]; R [events/s]");

    if (debug > 0){
        cout << "Number of events: (peak-0.3s,peak+0.3) " << frate->Integral(frate->GetMaximumX()-0.3,frate->GetMaximumX()+0.3) << endl;
    }
    

    if (drawopt){
        crates->cd(2);
        frate->Draw();
    }

    

    TH1I *hsignal = new TH1I("hsignal","",nbins,tleft,tright);
    //cout << "Bin size: " << hsignal->GetBinWidth(0) << endl;
    //TRandom3 *r3 = new TRandom3(0);
    TRandom3 *r3 = new TRandom3();
    if(gRandom) delete gRandom;
    gRandom = r3;
    r3->SetSeed(seed);

    if (fillopt == 1) {
        int totalnormsmeared = r3->Poisson(totalnorm*frate->Integral(hsignal->GetXaxis()->GetXmin(),hsignal->GetXaxis()->GetXmax())/frate->Integral(tdelay,1e9+tdelay));  //correction to integral is added since number of events will only be placed in the histogram range
        if (debug > 0) {
            cout << "Check norm: " << frate->Integral(tdelay,0.4+tdelay)/frate->Integral(tdelay,1e9+tdelay) << endl;
            cout << "Check norm: " << frate->Integral(tdelay+0.1,0.4+tdelay+0.1)/frate->Integral(tdelay,1e9+tdelay) << endl;
            cout << "Check norm: " << frate->Integral(tdelay+0.05,0.4+tdelay+0.05)/frate->Integral(tdelay,1e9+tdelay) << endl;
            cout << "Check norm: " << frate->Integral(tdelay,1+tdelay)/frate->Integral(tdelay,1e9+tdelay) << endl;
            cout << "Check norm: " << frate->Integral(tdelay,0.2+tdelay)/frate->Integral(tdelay,1e9+tdelay) << endl;
        }
        hsignal->FillRandom("frate",totalnormsmeared);  //This will create Poisson number of events in each bin
    } else if (fillopt == 2) {
        double meanval;
        double lowedge;
        for (int ii = 1; ii <= hsignal->GetNbinsX(); ii++) {
            lowedge =  hsignal->GetBinLowEdge(ii);
            meanval = frate->Integral(lowedge,lowedge+binwidth/1000.);
            if (meanval > 0) hsignal->SetBinContent(ii,r3->Poisson(meanval));
        }
    } else if (fillopt == 3) {
        double meanval;
        double lowedge;
        for (int ii = 1; ii <= hsignal->GetNbinsX(); ii++) {
            lowedge =  hsignal->GetBinLowEdge(ii);
            meanval = frate->Eval(lowedge+binwidth/2/1000.)*(binwidth/1000.);
            if (meanval > 0) hsignal->SetBinContent(ii,r3->Poisson(meanval));
        }
    }
    else cerr << "ERROR: fillopt should be 1, 2 or 3and it is instead: " << fillopt << endl;
    if (drawopt) {
        hsignal->SetTitle("signal events per bin;t [s]; N [events]");
        crates->cd(3);
        hsignal->SetStats(0);
        hsignal->Draw();
    }

    TH1D *hsignalbg = new TH1D("hsignalbg","",nbins,tleft,tright);
    if (fillopt == 1) {
        int totalbg = r3->Poisson(Rbg*10); //histogram duration is 10s hardcoded
        hsignalbg->FillRandom("pol0",totalbg);
    } else if (fillopt == 2 || fillopt == 3) {
        double meanval = Rbg*binwidth/1000.;
        for (int ii = 1; ii <= hsignalbg->GetNbinsX(); ii++) {
            hsignalbg->SetBinContent(ii,r3->Poisson(meanval));
        }
    } else cerr << "ERROR: fillopt should be 1-3 and it is instead: " << fillopt << endl;
    hsignalbg->Add(hsignal);
    hsignalbg->SetTitle("signal and background events per bin;t [s]; N [events]");
    ostringstream ssfname;
    ssfname << "fluxparametrisation_" << Meff << "kT_" << Rbg << "Hz_" << tdelay << "s";
    hsignalbg->SetName(ssfname.str().c_str());

    if (drawopt) {
        crates->cd(4);
        hsignalbg->SetStats(0);
        hsignalbg->Draw();
        ssfname.str("");
        ssfname << "fluxparametrisation_" << Meff << "kT_" << Rbg << "Hz_";
        ios oldState(nullptr);
        oldState.copyfmt(ssfname);
        ssfname << fixed << setprecision(1) << tdelay*1000;
        ssfname.copyfmt(oldState);
        ssfname<< "msT0_" << binwidth << "ms" << "bin.pdf";
        crates->Print(ssfname.str().c_str());
    }

    if (printopt) {
        ssfname.str("");
        ssfname << "fluxparametrisation_" << Meff << "kT_" << Rbg << "Hz_";
        ios oldState(nullptr);
        oldState.copyfmt(ssfname);
        ssfname << fixed << setprecision(1) << tdelay*1000;
        ssfname.copyfmt(oldState);
        ssfname<< "msT0_" << binwidth << "ms" << "bin.txt";
        ofstream fout(ssfname.str().c_str());
        for (int ii = 1; ii <= hsignalbg->GetNbinsX(); ii++) {
            fout << hsignalbg->GetXaxis()->GetBinLowEdge(ii) << "\t" <<  hsignalbg->GetBinContent(ii) << endl;
        }
    }

    if (! drawopt) {
        delete hsignal;
        delete frate;
    }

    return hsignalbg;
}

