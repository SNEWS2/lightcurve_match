# Skymap script 

## Installation 

```
yum/apt-get install python3
pip install numpy matplotlib random healpy scipy
```
Note that if this error is encounted: `ModuleNotFoundError: No module named 'numpy.testing.decorators` one can try to fix with:
```
pip install numpy==1.16.4
```

## Running

```
make -j <N>
```
About 3 hours are needed per map per core (the HEALPix resolution and number of realisations can be reduced in MAkefile to speed it up).
