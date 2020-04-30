# Matching codes and performance estimation.
## Description

1. Chi2 with:

    - ready to use "real" normalization (parameters are taken from the light curves)
    - "ideal" normalization (Meff and Rbg are taken as in the detector simulations)

2. Cross-corellation with zero-level normalization.

For the performance optimisation number of light curve realizations is set to 1000 (hardcoded in `chi2code.cc` as `ntoys`).

## Running performance estimation

To run all the performance simulations simply type:

    make

Use `-j <N>` to parallelize the calculation.

The outputs are stored in log files. The lines in the end contain the mean value and the resolution.
For example, for IceCube/Hyper-K with chi2 and "ideal" correction:

     Mean 0.017+-0.0164442 RMS 0.520011+-0.0116278

Parameters description can be found in `getDeltaT` description in `chi2code.cc`.

## Fire drill

A separate application `getdelay` is developed to get the delay of two input histograms.
Histograms should be providied as text files. 
Each row should have time (units are not important) and number of events.
The step should be fixed to 0.1s.

Run with:

```
getdelay filename1 filename2 chi2/corrslide rebinwidth[ms] chi2step[ms] windowmin[ms] windowmax[ms] shift_min[ms] shift_max[ms] shift_step[ms]
```

See a mini-test with 100 alerts (to compare performance) on Icecube/Hyper-K with:

```
make testgetdelay
```
Or:
```
./testgetdelay.sh t1[ms] t2[ms]
```
