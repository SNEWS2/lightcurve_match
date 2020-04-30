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
./makeskymap.py
```

For other combinations: 

1. comment/decomment the chosen possibility on the code
2. put the correct number of degrees of freedom for pvalue calculation:
    -  6 for 4 detector combination
    - 3 for 3 detectors combination
