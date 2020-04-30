#ifndef _INC_MATCHING
#define _INC_MATCHING
#include <string>
#include "TH1D.h"
const int debug = 1; //to save debug info
/* 
 * method: chi2 (chi-square) corr (not to use) corrslide (cross-corellation)
 *
 * rebinwidth width of the bins after rebin in ms. Optimal value is a tradeoff between statistics and curve features size
 *
 * chi2step[ms] - step used for calculating chi2/corr sum (it is calculated at some fixed points with a step 
 * correspondent to integer number of original bin widths.
 *
 * It can be set to rebinwidth to speed up, but technically, improvement (more information) can be obtained
 * if smaller chi2step up to 0.1ms is used (which corresponds to bin size in the original histo).
 *
 * windowmin, windowmax  - define xragne of histogram [ms] (respect to the first histogram) where chi2 is calculated
 *
 * shift_min/shift_max range of the matching time search [ms]
 *
 * Meff1, Meff2 are detector effective masses [Mt]
 *
 * Rbg1, Rbg2 are detector backgrounds [Hz]
 *
 * shift_step - step for matching time
 *
 */
double getT0(std::string method, TH1D *hdet1, TH1D *hdet2, double &minchi2, double rebinwidth = 50, double chi2step = 50, double windowmin = -100, double windowmax = 500, double shift_min = -30, double shift_max = 30, double shift_step = 0.1);

/*
 * Correct detector light curve histogram.
 *
 * bgopt - correct background:
 * 1       Lightcurve - known apriori background level (value that you add in the simulation).
 * 2       Lightcurve - background level as a median value of the bins in (tleft, tright)
 *
 * sgopt - correct signal level (after backround subraction):
 * 1       rescale to known apriori detector effective mass
 * 2       rescale lightcurves from to have integral of curves in as 1
 *
 * bgopt=3 an sgopt=3 - "zero-normalisation" correction for cross-correllation
 *
 * Rbg - backgorund rate in the detector in Hz
 * 
 * Meff - detector effective mass (in Mt), 1 if no scale is needed.
 *
 * Default options are chosen in the way that histogram is unchanged.
 *
 * tleft, tright - time window for signal integration (background normalization)
 *
 * rebinwidth is needed for cross-correlation correction
 *
 */
TH1D *corrhist(TH1D *hdet, int bgopt = 1, int sgopt = 1, double Meff = 1, double Rbg = 0, double tleft = -100, double tright=500, double rebinwidth = 1);

//function that sums nbin next bins starting from the bin with x
inline double VRebin(TH1D *hdet, double x,int nbin, double &err2){
    int bin0 = hdet->FindBin(x+1e-5); //shift is added to avoid values on the bin borders
    //cout << "DEBBB: " << bin0 << " " << bin0+nbin << " " << x << endl;
    double result = 0;
    err2 = 0;
    for (int ii = bin0; ii<bin0+nbin; ii++) { //we don't care that the right bins are taken only since this procedure is done both for h1 and h2
        result += hdet->GetBinContent(ii);
        err2 += pow(hdet->GetBinError(ii),2);
    }
    return result;
}
//interpolation function - not used
inline double VInt(double x, double x0, double x1, double y0, double y1) {
    return (x-x0)/(x1-x0)*(y1-y0)+y0;
}

//improved interpolation that should have similar uncertainty for each interpolated value - not used
inline double VInterpolate(TH1D *hdet, double x) {
    int bin0 = hdet->FindBin(x);
    double xbin[] = {hdet->GetBinCenter(bin0-1),hdet->GetBinCenter(bin0),hdet->GetBinCenter(bin0+1)};
    double ybin[] = {hdet->GetBinContent(bin0-1),hdet->GetBinContent(bin0),hdet->GetBinContent(bin0+1)};
    double yfar  = VInt(x,xbin[0],xbin[2],ybin[0],ybin[2]);
    double ynear;
    if (x < xbin[1]) {
        ynear = VInt(x,xbin[0],xbin[1],ybin[0],ybin[1]); 
    } else if (x == xbin[1]) {
        ynear = xbin[1];
    } else {
        ynear = VInt(x,xbin[1],xbin[2],ybin[1],ybin[2]);
    }
    return (yfar+ynear)/2.;
}


#endif
