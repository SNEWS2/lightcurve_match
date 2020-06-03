NSIDE = 256
NTOY = 100000

GClon = -94.4
GClat = -28.92

Blon = 88.8
Blat = 7.4

Clon = -45
Clat = 40

.PHONY: clean all

all: GC_IC_HK_ARCA_JUNO GC_IC_HK_ARCA GC_IC_HK_JUNO GC_IC_ARCA_JUNO GC_HK_ARCA_JUNO Cygnus_IC_HK_ARCA_JUNO Betelgeuse_HK_ARCA_JUNO Cygnus_HK_ARCA_JUNO GC_IC_SK_JUNO GC_IC_SK_ARCA

GC_IC_HK_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py $(GClon) $(GClat) $(NSIDE) $(NTOY) IC HK ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
GC_IC_HK_ARCA: makeskymap.py
	{ time python makeskymap.py $(GClon) $(GClat) $(NSIDE) $(NTOY) IC HK ARCA; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
GC_IC_HK_JUNO: makeskymap.py
	{ time python makeskymap.py $(GClon) $(GClat) $(NSIDE) $(NTOY) IC HK JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
GC_IC_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py $(GClon) $(GClat) $(NSIDE) $(NTOY) IC ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
GC_HK_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py $(GClon) $(GClat) $(NSIDE) $(NTOY) HK ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
Betelgeuse_IC_HK_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py $(Blon) $(Blat) $(NSIDE) $(NTOY) IC HK ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
Cygnus_IC_HK_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py $(Clon) $(Clat) $(NSIDE) $(NTOY) IC HK ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
Betelgeuse_HK_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py $(Blon) $(Blat) $(NSIDE) $(NTOY) HK ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
Cygnus_HK_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py $(Clon) $(Clat) $(NSIDE) $(NTOY) HK ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt
GC_IC_SK_JUNO: makeskymap.py
	{ time python makeskymap.py $(GClon) $(GClat) $(NSIDE) $(NTOY) IC SK JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt
GC_IC_SK_ARCA: makeskymap.py
	{ time python makeskymap.py $(GClon) $(GClat) $(NSIDE) $(NTOY) IC SK ARCA; } &> log_$(NSIDE)_$(NTOY)_$@.txt



clean:
	rm *png *txt
