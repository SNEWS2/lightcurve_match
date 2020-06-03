NSIDE = 256
NTOY = 100000

.PHONY: clean all

all: GC_IC_HK_ARCA_JUNO GC_IC_HK_ARCA GC_IC_HK_JUNO GC_IC_ARCA_JUNO GC_HK_ARCA_JUNO Cygnus_IC_HK_ARCA_JUNO Betelgeuse_HK_ARCA_JUNO Cygnus_HK_ARCA_JUNO GC_IC_SK_JUNO GC_IC_SK_ARCA

GC_IC_HK_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py -94.4 -28.94 $(NSIDE) $(NTOY) IC HK ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
GC_IC_HK_ARCA: makeskymap.py
	{ time python makeskymap.py -94.4 -28.94 $(NSIDE) $(NTOY) IC HK ARCA; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
GC_IC_HK_JUNO: makeskymap.py
	{ time python makeskymap.py -94.4 -28.94 $(NSIDE) $(NTOY) IC HK JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
GC_IC_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py -94.4 -28.94 $(NSIDE) $(NTOY) IC ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
GC_HK_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py -94.4 -28.94 $(NSIDE) $(NTOY) HK ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
Betelgeuse_IC_HK_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py 88.8 7.4 $(NSIDE) $(NTOY) IC HK ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
Cygnus_IC_HK_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py -45 40 $(NSIDE) $(NTOY) IC HK ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
Betelgeuse_HK_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py 88.8 7.4 $(NSIDE) $(NTOY) HK ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt 
Cygnus_HK_ARCA_JUNO: makeskymap.py
	{ time python makeskymap.py -45 40 $(NSIDE) $(NTOY) HK ARCA JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt
GC_IC_SK_JUNO: makeskymap.py
	{ time python makeskymap.py -94.4 -28.94 $(NSIDE) $(NTOY) IC SK JUNO; } &> log_$(NSIDE)_$(NTOY)_$@.txt
GC_IC_SK_ARCA: makeskymap.py
	{ time python makeskymap.py -94.4 -28.94 $(NSIDE) $(NTOY) IC SK ARCA; } &> log_$(NSIDE)_$(NTOY)_$@.txt



clean:
	rm *png *txt
