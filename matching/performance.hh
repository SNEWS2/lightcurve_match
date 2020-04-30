#ifndef _INC_PERFORMANCE
#define _INC_PERFORMANCE 
/* 
 * bgopt - correct background:
 * 1 Lightcurve - known apriori background level (Rbg1&Rbg2).
 * 100-200 Lightcurve - known apriori background level +- 10% random variation.
 * 2 Lightcurve - background level as a median value of the bins.  
 *
 * sgopt - correct signal level (after backround subraction):
 * 1 rescale to known apriori detector effective mass (Meff1&Meff2)
 * 100-200 rescale to known apriori detector effective mass +- 10% random variation
 * 2 rescale lightcurves from to have integral of curves in as 1
 * if samecorrection=true is set than the same correction scale from bg is used (now it is hardcoded to false in performance.cc)
 *
 * rebinwidth width of the bins after rebin in ms. Optimal value is a tradeoff between statistics and curve features size
 *
 * chi2step[ms] - step used for calculating chi2 sum (it is calculated at some fixed points with a step 
 * correspondent to integer number of original bin widths.
 *
 * It can be set to rebinwidth to speed up, but technically, improvement (more information) can be obtained
 * if smaller chi2step up to 0.1ms is used (which corresponds to bin size in the original histo).
 *
 * windowmin, windowmax  - define xragne of histogram (respect to the first histogram) where chi2 is calculated
 *
 * _tdelay if set [-30,30]ms it is used as a fixed delay, otherwise it is random (-30,30)
 *
 * Meff1, Meff2 are detector effective masses in Mt
 *
 * Rbg1, Rbg2 are detector backgrounds in Hz
 *
 * searchwindow: if >= 100, fair bling search of T0 in the range -100, 100 is performed
 * if < 100, then T0 is searched in the range (T0_true - searchwindow, T0_true + searchwindow) to speed up the calculation
 *
 */
#include <string>
void getDeltaT(std::string method, int bgopt1 = 1, int sgopt1 = 1, int bgopt2 = 1, int sgopt2 = 1, double rebinwidth = 50, double chi2step = 50, double windowmin = -100, double windowmax = 500, double _tdelay = 7.5, double Meff1 = 3290, double Rbg1 = 0, double Meff2 = 3290, double Rbg2 = 0, double searchwindow = 100);
#endif
