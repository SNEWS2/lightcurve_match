#!/bin/bash

if [ ! -f ../simulation/detectorrate ]; then
    cd ../simulation
    make detectorrate
    cd -
fi

DELAY1_MS=${1:--1}
DELAY2_MS=${2:-2}

DELAY1_S=`echo "$DELAY1_MS/1000" | bc -l`
DELAY2_S=`echo "$DELAY2_MS/1000" | bc -l`

DELAY1_MS=`printf %0.1f $DELAY1_MS`
DELAY2_MS=`printf %0.1f $DELAY2_MS`

echo "Running with T1=$DELAY1_MS ms and T2=$DELAY2_MS ms"

for i in {1..100}; do
    ../simulation/detectorrate 3500 1548000 $DELAY1_S 0.1 3 &> /dev/null
    ../simulation/detectorrate 560 0 $DELAY2_S 0.1 3 &> /dev/null
    echo getdelay \#$i
    { time ./getdelay fluxparametrisation_3500kT_1.548e+06Hz_${DELAY1_MS}msT0_0.1msbin.txt fluxparametrisation_560kT_0Hz_${DELAY2_MS}msT0_0.1msbin.txt chi2 50 50 -300 300 -30 30 0.1 ; } &> getdelay.txt
    tail -n5 getdelay.txt | head -n1
    tail -n3 getdelay.txt | head -n1
done

med=0
var=0
counter=0
for i in `cat getdelay.txt | grep T0match | cut -f 2 -d " "`; do
    let counter=counter+1
    med=`echo "$med+$i" | bc -l`
    var=`echo "$var+$i*$i" | bc -l`
done
var=`echo "sqrt(($counter*$var - $med*$med)/$counter/($counter-1))" | bc -l`
med=`echo "$med/$counter" | bc -l`

offset=`echo "$med-($DELAY2_MS-$DELAY1_MS)" | bc -l`

echo "Mean: $med Offset: $offset RMS: $var"
