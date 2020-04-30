# Detector neutrino light curve simulation

- Uses nu-bar cross-sections from A. Strumia, F. Vissani, Phys. Lett. B 564, 42 (2003).
- Uses neutrino light curve parametrisation from E. Nardia, I. Zuluaga, Nuclear Physics B 731, 140-163 (2005)

## Running standalone
     
```
make detectorrate
./detectorrate Meff[kt] Fbck[Hz] tdelay[s]
```

Or for pretabulated detectors:

```
make [IC|ARCA|ORCA|SK|HK|JUNO]
```

## Conversion factor

Conversion factor I is hardcoded to `detectorrate.cc`.
If this factor should be recalculated -DIFACTOR should be added to the code compillation:

```
g++ main.cc conversionfactor.cc detectorrate.cc -o detectorrate -I. `root-config --glibs --cflags`  -DIFACTOR
make IC
```
