# Supplementary materials for "Combining neutrino experimental light-curves for pointing to the next Galactic Core-Collapse Supernova"

Authors: A. Coleiro, D. Dornic, M. Colomer Molla, M. Lincetto, V. Kulikovskiy

## Content

- `Makefile` runs all analyses codes with `make -j N`... takes about 1--2 hours with 4 cores (N=4).
- `det.csv` data file with simplified detector characteristics which are used in the `simulation` and `matching` codes
- `simulation` detected neutrino light curve simulation codes (`README.md` there for more details) 
- `matching` chi-square and cross-correlation matching codes (`README.md` there for more details)
    - for printing the results after `make` one can use a bash loop:
```
    for i in `ls log*`; do echo $i; cat $i | grep "^Mean"; done
```
- `skymap` HEALPix skymap creation code and error box area calculation (`README.md` there for more details) 
    - for the code ready-to use matching code that takes two txt files with histograms see `getdelay.cc`
