#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <sstream>
#include <cstring>
//#include <cstdlib> //for system
#include <fstream>
#include <TTree.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TNtuple.h>
#include "../simulation/detectorrate.hh"
#include "performance.hh"
#include "matching.hh"
#include <ctime> //for clock

#include "TSystem.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Script to compute the uncertainty deltat using a chi2 method //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void getDeltaT(string method, int bgopt1, int sgopt1, int bgopt2, int sgopt2, double rebinwidth, double chi2step, double windowmin, double windowmax, double _tdelay, double Meff1, double Rbg1, double Meff2, double Rbg2, double searchwindow){
    ostringstream oss;
    oss << "deltaT_bg" << bgopt1 << "_sg" << sgopt1 << "_bg" << bgopt2 << "_sg" << sgopt2 << "_rbn" << rebinwidth  << "_cs" << chi2step << "_l" << windowmin << "_r" << windowmax<< "_t" << _tdelay << "_M" << Meff1 << "_R" << Rbg1 << "_M" << Meff2 << "_R"     << Rbg2 << "_sw" << searchwindow << ".root";
    TFile *fout = new TFile(oss.str().c_str(), "RECREATE");

    const int ntoys = 1000;  //Number of Toy MC to simulate
    //These should be the official values for scan in ms, the maximum time of ths scan, driven by the distance between experiments / c_light
    double shift_min;
    double shift_max;
    const double shift_step = 0.1; //in ms
    int nscan = (shift_max-shift_min)/shift_step;
    TRandom3 *rd = new TRandom3(20001); //fix seed to have reproducible light curves

    // Output histogram with deltat uncertainty
    TH1D* hdeltat = new TH1D("hdt", "", 200/shift_step+1, -100-shift_step/2.,100+shift_step/2.);
    hdeltat->SetDirectory(0);
    hdeltat->SetXTitle("T0_{fit}-T0_{true} [ms]");

    TNtuple *nt = new TNtuple("fitresults", "","t0:tfit:q");

    //TH1D* hdeltachi2 = new TH1D("hdchi2", "", 600,500,1900);
    //hdeltachi2->SetDirectory(0);
    //hdeltachi2->SetXTitle("#chi^{2} min");

    TCanvas *c1 = new TCanvas();

    // Loop over nunmber of Toy MC realizations
    for(int itoy=0; itoy<ntoys; itoy++){
        double t0;
        if (fabs(_tdelay) <= 30) t0 = _tdelay;
        else t0=rd->Uniform(-30,30);  

        if (searchwindow >= 100) { //the correct blind choise
            shift_min = -100;
            shift_max = 100;  //these should be the official values for scan
        } else {
            //values to speed up the tests
            shift_min = t0 - searchwindow; // searchwindow should be some >3 RMS of the final histo
            shift_max = t0 + searchwindow; 
        }

        double startloop = std::clock();

        TH1D *hdet1;
        TH1D *hdet2;

        hdet1 = detectorrate(Meff1,Rbg1,0,0.1,false,false,3,itoy+1);
        hdet2 = detectorrate(Meff2,Rbg2,t0/1000,0.1,false,false,3,(itoy+1)+10000);

        hdet1->GetXaxis()->SetLimits(-2000,8000);
        hdet2->GetXaxis()->SetLimits(-2000,8000);

        hdet1->SetDirectory(0);
        hdet2->SetDirectory(0);

        double startchi2 = std::clock();

        //we define chi2 calculation window respect to the center of the first histogram (which also has better sensitivity to have better performance)
        double light_curve_max = hdet1->GetBinLowEdge(hdet1->GetMaximumBin());

        bool samecorrection = false; //to use same systematics for signal and background
        double delta1,delta2;
        double Rbg1corr = Rbg1;
        double Rbg2corr = Rbg2;
        double Meff1corr = Meff1;
        double Meff2corr = Meff2;
        if (bgopt1 > 100 && bgopt1 < 200) {
            delta1 = (bgopt1-100)/100.;
            TRandom3 *rd = new TRandom3(30000+itoy); //fix seed to have reproducible results
            delta1 = rd->Uniform(1-delta1,1+delta1);
            Rbg1corr = Rbg1*hdet1->GetBinWidth(0)/1000.*delta1;

        }
        if (sgopt1 > 100 && sgopt1 < 200) { 
            TRandom3 *rd = new TRandom3(40000+itoy); //fix seed to have reproducible results
            double Meffdelta = (sgopt1-100)/100.;    
            if (!samecorrection) delta1 = rd->Uniform(1-Meffdelta,1+Meffdelta);
            Meff1corr = Meff1*delta1;
        }
         if (bgopt2 > 100 && bgopt2 < 200) {
            delta2 = (bgopt2-100)/100.;
            TRandom3 *rd = new TRandom3(50000+itoy); //fix seed to have reproducible results
            delta2 = rd->Uniform(1-delta2,1+delta2);
            Rbg2corr = Rbg2*hdet2->GetBinWidth(0)/1000.*delta2;

        }
        if (sgopt2 > 100 && sgopt2 < 200) { 
            TRandom3 *rd = new TRandom3(60000+itoy); //fix seed to have reproducible results
            double Meffdelta = (sgopt2-100)/100.; ;
            if (!samecorrection) delta2 = rd->Uniform(1-Meffdelta,1+Meffdelta);
            Meff2corr = Meff2*delta2;
        }
        
        TH1D *hdet1corr = corrhist(hdet1, bgopt1, sgopt1, Meff1corr, Rbg1corr, light_curve_max+windowmin, light_curve_max+windowmax,rebinwidth);
        TH1D *hdet2corr = corrhist(hdet2, bgopt2, sgopt2, Meff2corr, Rbg2corr, light_curve_max+windowmin, light_curve_max+windowmax,rebinwidth);

        if (debug > 0) cout << "debug: light curve max " << light_curve_max << endl;
        double minchi2;
        double Tfit = getT0(method, hdet1corr,hdet2corr,minchi2,rebinwidth,chi2step,light_curve_max+windowmin,light_curve_max+windowmax,shift_min,shift_max,shift_step);
        cout << "toy #" << itoy << " T0:  " << t0 << " T0fit " << Tfit << " chi2 " << minchi2 << endl;
        nt->Fill(t0,Tfit,minchi2);

        hdeltat->Fill(t0-Tfit);

        double stoploop =  std::clock();
        cout << "Full loop, chi2 calc (s): " << ( stoploop - startloop ) / (double) CLOCKS_PER_SEC << " " << ( stoploop - startchi2 ) / (double) CLOCKS_PER_SEC << endl;

        hdet1corr->Delete();
        hdet2corr->Delete();
        hdet1->Delete();
        hdet2->Delete();
    }

    TCanvas* canv1 = new TCanvas("canv1", "", 800, 700);
    canv1->cd();
    hdeltat->Draw();
    if (hdeltat->GetRMS() > 10) hdeltat->Rebin(10);
    hdeltat->Fit("gaus");
    gStyle->SetOptStat("neMR");
    gStyle->SetOptFit(1);
    hdeltat->GetXaxis()->SetRangeUser(hdeltat->GetMean()-5*hdeltat->GetRMS(), hdeltat->GetMean()+5*hdeltat->GetRMS());


    oss.str("");
    oss << "deltaT_bg" << bgopt1 << "_sg" << sgopt1 << "_bg" << bgopt2 << "_sg" << sgopt2 << "_rbn" << rebinwidth  << "_cs" << chi2step << "_l" << windowmin << "_r" << windowmax<< "_t" << _tdelay << "_M" << Meff1 << "_R" << Rbg1 << "_M" << Meff2 << "_R" << Rbg2 << "_sw" << searchwindow << ".pdf";
    canv1->Print(oss.str().c_str());
    oss.str("");
    hdeltat->Write();
    nt->Write();
    fout->Close();

    cout << "Mean " << hdeltat->GetMean() << "+-" << hdeltat->GetMeanError() << " RMS " << hdeltat->GetRMS() << "+-" << hdeltat->GetRMSError() << endl;
}
int main(int argc, char *argv[]) {
    if (argc == 16) {

        string method;
        int bgopt1, sgopt1, bgopt2, sgopt2;
        double rebinwidth, chi2step, windowmin, windowmax, tdelay, Meff1, Rbg1,  Meff2, Rbg2, searchwindow;

        unsigned int iarg = 0;
        
        istringstream issmethod(argv[++iarg]);  
        issmethod >> method;

        istringstream iss1(argv[++iarg]);
        iss1 >> bgopt1; 

        istringstream iss2(argv[++iarg]);
        iss2 >> sgopt1; 

        istringstream issbgopt2(argv[++iarg]);
        issbgopt2 >> bgopt2;

        istringstream isssgopt2(argv[++iarg]);   
        isssgopt2 >> sgopt2;

        istringstream iss3(argv[++iarg]);
        iss3 >> rebinwidth;

        istringstream iss4(argv[++iarg]);
        iss4>> chi2step;

        istringstream isswindowmin(argv[++iarg]);
        isswindowmin >> windowmin;

        istringstream isswindowmax(argv[++iarg]);
        isswindowmax >> windowmax;

        istringstream isstdelay(argv[++iarg]);
        isstdelay >> tdelay;

        istringstream iss5(argv[++iarg]);
        iss5 >> Meff1;

        istringstream iss6(argv[++iarg]);
        iss6 >> Rbg1;

        istringstream iss7(argv[++iarg]);
        iss7 >> Meff2;

        istringstream iss8(argv[++iarg]);
        iss8 >> Rbg2;

        istringstream iss9(argv[++iarg]);
        iss9 >> searchwindow;

        getDeltaT(method, bgopt1, sgopt1, bgopt2, sgopt2, rebinwidth, chi2step, windowmin, windowmax, tdelay, Meff1, Rbg1,  Meff2, Rbg2, searchwindow);

        return 0;
    }
    else {
        cerr << "run with: " << endl;
        cerr << argv[0] << " chi2/corrslide bgopt1[1-2] sgopt1[1-2] bgopt2[1-2] sgopt2[1-2] rebinwidth[ms] chi2step[ms] windowmin[ms] windowmax[ms] T0[-30,30]s|100for_random Meff1[kT] Rbg1[Hz] Meff2[kT] Rbg2[Hz] searchwindow[ms]" << endl;

    }

    return 1;
}
