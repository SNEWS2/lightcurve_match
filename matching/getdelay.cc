#include "matching.hh"
#include <sstream>
#include <iostream>
#include <fstream>
#include "TFile.h"
using namespace std;
int main(int argc, char *argv[]) {
    if (argc == 11) {

        string f1, f2;
        string method;
        double rebinwidth, chi2step, windowmin, windowmax, shiftmin, shiftmax, shift_step;

        unsigned int iarg = 0;

        istringstream issf1(argv[++iarg]);
        issf1 >> f1;

        istringstream issf2(argv[++iarg]);
        issf2 >> f2;

        istringstream issmethod(argv[++iarg]);  
        issmethod >> method;

        istringstream issrebinwidth(argv[++iarg]);
        issrebinwidth >> rebinwidth;

        istringstream isschi2step(argv[++iarg]);
        isschi2step>> chi2step;

        istringstream isswindowmin(argv[++iarg]);
        isswindowmin >> windowmin;

        istringstream isswindowmax(argv[++iarg]);
        isswindowmax >> windowmax;

        istringstream issshiftmin(argv[++iarg]);
        issshiftmin >> shiftmin;

        istringstream issshiftmax(argv[++iarg]);
        issshiftmax >> shiftmax;

        istringstream issshift_step(argv[++iarg]);
        issshift_step >> shift_step;

        TFile *fout = new TFile("getdelay.root","RECREATE");

	TH1D *hdet1 = new TH1D("h1", "", 100000, -2000, 8000);
	TH1D *hdet2 = new TH1D("h2", "", 100000, -2000, 8000);
	ifstream ifs(f1);
	double _t, _val;
	int bin = 0;
	while (ifs >> _t >> _val) {
	    bin++;
	    hdet1->SetBinContent(bin,_val);
	}
	ifstream ifs2(f2);
	bin = 0;
	while (ifs2 >> _t >> _val) {
	    bin++;
	    hdet2->SetBinContent(bin,_val);
	}
	
        int bgopt1, sgopt1, bgopt2, sgopt2;

        if (method == "chi2") {
            bgopt1 = 2;
            sgopt1 = 2;
            bgopt2 = 2;
            sgopt2 = 2;
        } else if (method == "corrslide") {
            bgopt1 = 0;
            sgopt1 = 0;
            bgopt2 = 0;
            sgopt2 = 0;
        } else {
            cerr << "Option " << method << " is not supported" << endl;
            return 1;
        }

        double light_curve_max = hdet1->GetBinLowEdge(hdet1->GetMaximumBin());
	TH1D *hdet1corr = corrhist(hdet1, bgopt1, sgopt1, 0, 0, light_curve_max+windowmin, light_curve_max+windowmax,rebinwidth);
	TH1D *hdet2corr = corrhist(hdet2, bgopt2, sgopt2, 0, 0, light_curve_max+windowmin, light_curve_max+windowmax,rebinwidth);

        double minchi2;
        double Tfit = getT0(method, hdet1corr,hdet2corr,minchi2,rebinwidth,chi2step,light_curve_max+windowmin,light_curve_max+windowmax,shiftmin,shiftmax,shift_step);
        cout << "T0match: " << Tfit << endl;

        fout->Close();


        return 0;
    }
    else {
        cerr << "run with: " << endl;
        cerr << argv[0] << " filename1 filename2 chi2/corrslide rebinwidth[ms] chi2step[ms] windowmin[ms] windowmax[ms] shift_min[ms] shift_max[ms] shift_step[ms]" << endl;

    }

    return 1;
}
