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

.PHONY: IC clean

all: detectorrate IC

detectorrate: main.cc detectorrate.cc conversionfactor.cc conversionfactor.hh
	g++ main.cc conversionfactor.cc detectorrate.cc -o detectorrate -I. `root-config --glibs --cflags`

IC: detectorrate
	./detectorrate $(MEFF_IC) $(RBG_IC) 0

ORCA: detectorrate
	./detectorrate $(MEFF_ORCA) $(RBG_ORCA) 0

ARCA: detectorrate
	./detectorrate $(MEFF_ARCA) $(RBG_ARCA) 0

SK: detectorrate
	./detectorrate $(MEFF_SK) $(RBG_SK) 0

HK: detectorrate
	./detectorrate $(MEFF_HK) $(RBG_HK) 0

JUNO: detectorrate
	./detectorrate $(MEFF_JUNO) $(RBG_JUNO) 0
clean:
	rm fluxparametrisation*pdf fluxparametrisation*txt detectorrate
