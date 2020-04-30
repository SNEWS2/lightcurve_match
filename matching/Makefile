MEFF_IC   = $(shell cat ../det.csv | grep IC | cut -f 2)
RBG_IC    = $(shell cat ../det.csv | grep IC | cut -f 3)
MEFF_ARCA = $(shell cat ../det.csv | grep ARCA | cut -f 2)   
RBG_ARCA  = $(shell cat ../det.csv | grep ARCA | cut -f 3)    
MEFF_HK   = $(shell cat ../det.csv | grep HK | cut -f 2)
RBG_HK    = $(shell cat ../det.csv | grep HK | cut -f 3) 
MEFF_SK   = $(shell cat ../det.csv | grep SK | cut -f 2)
RBG_SK    = $(shell cat ../det.csv | grep SK | cut -f 3)
MEFF_JUNO = $(shell cat ../det.csv | grep JUNO | cut -f 2)
RBG_JUNO  = $(shell cat ../det.csv | grep JUNO | cut -f 3)

.PHONY: clean testgetdelay all

all: performance ICARCA ICSK ICHK SKARCA HKARCA JUNOARCA SKJUNO HKJUNO ICJUNO

performance: performance.cc matching.cc performance.hh matching.hh ../simulation/conversionfactor.cc ../simulation/detectorrate.cc
	g++ -O3 performance.cc matching.cc ../simulation/conversionfactor.cc ../simulation/detectorrate.cc -o performance -I. `root-config --glibs --cflags`

getdelay: getdelay.cc matching.cc matching.hh ../simulation/conversionfactor.cc ../simulation/detectorrate.cc
	g++ -O3 getdelay.cc matching.cc ../simulation/conversionfactor.cc ../simulation/detectorrate.cc -o getdelay -I. `root-config --glibs --cflags`

testgetdelay: getdelay
	./testgetdelay.sh 10 2

ICARCA: performance
	{ time ./performance chi2 1 1 1 1 50 50 -300 300 100 $(MEFF_IC) $(RBG_IC) $(MEFF_ARCA) $(RBG_ARCA) 100 ; } &> logICARCA_50_50_ideal.txt
	{ time ./performance chi2 2 2 2 2 50 50 -300 300 100 $(MEFF_IC) $(RBG_IC) $(MEFF_ARCA) $(RBG_ARCA) 100 ; } &> logICARCA_50_50_real.txt
	{ time ./performance corrslide 0 0 0 0 10 10 -300 300 100 $(MEFF_IC) $(RBG_IC) $(MEFF_ARCA) $(RBG_ARCA) 100 ; } &> logICARCA_10_10_corr.txt

ICSK: performance
	{ time ./performance chi2 1 1 1 1 50 50 -300 300 100 $(MEFF_IC) $(RBG_IC) $(MEFF_SK) $(RBG_SK) 100 ; } &> logICSK_50_50_ideal.txt
	{ time ./performance chi2 2 2 2 2 50 50 -300 300 100 $(MEFF_IC) $(RBG_IC) $(MEFF_SK) $(RBG_SK) 100 ; } &> logICSK_50_50_real.txt
	{ time ./performance corrslide 0 0 0 0 10 10 -300 300 100 $(MEFF_IC) $(RBG_IC) $(MEFF_SK) $(RBG_SK) 100 ; } &> logICSK_10_10_corr.txt

ICHK: performance
	{ time ./performance chi2 1 1 1 1 50 50 -300 300 100 $(MEFF_HK) $(RBG_HK) $(MEFF_IC) $(RBG_IC) 100 ; } &> logHKIC_50_50_ideal.txt
	{ time ./performance chi2 2 2 2 2 50 50 -300 300 100 $(MEFF_HK) $(RBG_HK) $(MEFF_IC) $(RBG_IC) 100 ; } &> logHKIC_50_50_real.txt
	{ time ./performance corrslide 0 0 0 0 10 10 -300 300 100 $(MEFF_HK) $(RBG_HK) $(MEFF_IC) $(RBG_IC) 100 ; } &> logHKIC_10_10_corr.txt

ICJUNO: performance
	{ time ./performance chi2 1 1 1 1 50 50 -300 300 100 $(MEFF_IC) $(RBG_IC) $(MEFF_JUNO) $(RBG_JUNO) 100 ; } &> logICJUNO_50_50_ideal.txt
	{ time ./performance chi2 2 2 2 2 50 50 -300 300 100 $(MEFF_IC) $(RBG_IC) $(MEFF_JUNO) $(RBG_JUNO) 100 ; } &> logICJUNO_50_50_real.txt
	{ time ./performance corrslide 0 0 0 0 10 10 -300 300 100 $(MEFF_IC) $(RBG_IC) $(MEFF_JUNO) $(RBG_JUNO) 100 ; } &> logICJUNO_10_10_corr.txt

SKARCA: performance
	{ time ./performance chi2 1 1 1 1 50 50 -300 300 100 $(MEFF_SK) $(RBG_SK) $(MEFF_ARCA) $(RBG_ARCA) 100 ; } &> logSKARCA_50_50_ideal.txt
	{ time ./performance chi2 2 2 2 2 50 50 -300 300 100 $(MEFF_SK) $(RBG_SK) $(MEFF_ARCA) $(RBG_ARCA) 100 ; } &> logSKARCA_50_50_real.txt
	{ time ./performance corrslide 0 0 0 0 10 10 -300 300 100 $(MEFF_SK) $(RBG_SK) $(MEFF_ARCA) $(RBG_ARCA) 100 ; } &> logSKARCA_10_10_corr.txt

HKARCA: performance
	{ time ./performance chi2 1 1 1 1 50 50 -300 300 100 $(MEFF_HK) $(RBG_HK) $(MEFF_ARCA) $(RBG_ARCA) 100 ; } &> logHKARCA_50_50_ideal.txt
	{ time ./performance chi2 2 2 2 2 50 50 -300 300 100 $(MEFF_HK) $(RBG_HK) $(MEFF_ARCA) $(RBG_ARCA) 100 ; } &> logHKARCA_50_50_real.txt
	{ time ./performance corrslide 0 0 0 0 10 10 -300 300 100 $(MEFF_HK) $(RBG_HK) $(MEFF_ARCA) $(RBG_ARCA) 100 ; } &> logHKARCA_10_10_corr.txt

JUNOARCA: performance
	{ time ./performance chi2 1 1 1 1 50 50 -300 300 100 $(MEFF_JUNO) $(RBG_JUNO) $(MEFF_ARCA) $(RBG_ARCA) 100 ; } &> logJUNOARCA_50_50_ideal.txt
	{ time ./performance chi2 2 2 2 2 50 50 -300 300 100 $(MEFF_JUNO) $(RBG_JUNO) $(MEFF_ARCA) $(RBG_ARCA) 100 ; } &> logJUNOARCA_50_50_real.txt
	{ time ./performance corrslide 0 0 0 0 10 10 -300 300 100 $(MEFF_JUNO) $(RBG_JUNO) $(MEFF_ARCA) $(RBG_ARCA) 100 ; } &> logJUNOARCA_10_10_corr.txt

SKJUNO: performance
	{ time ./performance chi2 1 1 1 1 50 50 -300 300 100 $(MEFF_SK) $(RBG_SK) $(MEFF_JUNO) $(RBG_JUNO) 100 ; } &> logSKJUNO_50_50_ideal.txt
	{ time ./performance chi2 2 2 2 2 50 50 -300 300 100 $(MEFF_SK) $(RBG_SK) $(MEFF_JUNO) $(RBG_JUNO) 100 ; } &> logSKJUNO_50_50_real.txt
	{ time ./performance corrslide 0 0 0 0 10 10 -300 300 100 $(MEFF_SK) $(RBG_SK) $(MEFF_JUNO) $(RBG_JUNO) 100 ; } &> logSKJUNO_10_10_corr.txt

HKJUNO: performance
	{ time ./performance chi2 1 1 1 1 50 50 -300 300 100 $(MEFF_HK) $(RBG_HK) $(MEFF_JUNO) $(RBG_JUNO) 100 ; } &> logHKJUNO_50_50_ideal.txt
	{ time ./performance chi2 2 2 2 2 50 50 -300 300 100 $(MEFF_HK) $(RBG_HK) $(MEFF_JUNO) $(RBG_JUNO) 100 ; } &> logHKJUNO_50_50_real.txt
	{ time ./performance corrslide 0 0 0 0 10 10 -300 300 100 $(MEFF_HK) $(RBG_HK) $(MEFF_JUNO) $(RBG_JUNO) 100 ; } &> logHKJUNO_10_10_corr.txt

clean:
	rm deltaT*pdf performance
