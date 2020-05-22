NSIDE = 32
NTOY = 10000

.PHONY: clean testgetdelay all

all: IC_HK_ARCA_JUNO



GC_IC_HK_ARCA_JUNO: makeskymapadvanced.py
	{ time python makeskymapadvanced.py -94.4 -28.94 $(NSIDE) $(NTOY) IC HK ARCA JUNO; } &> log$@.txt
GC_IC_HK_ARCA: makeskymapadvanced.py
	{ time python makeskymapadvanced.py -94.4 -28.94 $(NSIDE) $(NTOY) IC HK ARCA; } &> log$@.txt
GC_IC_HK_JUNO: makeskymapadvanced.py
	{ time python makeskymapadvanced.py -94.4 -28.94 $(NSIDE) $(NTOY) IC HK JUNO; } &> log$@.txt
GC_IC_ARCA_JUNO: makeskymapadvanced.py
	{ time python makeskymapadvanced.py -94.4 -28.94 $(NSIDE) $(NTOY) IC ARCA JUNO; } &> log$@.txt
GC_HK_ARCA_JUNO: makeskymapadvanced.py
	{ time python makeskymapadvanced.py -94.4 -28.94 $(NSIDE) $(NTOY) HK ARCA JUNO; } &> log$@.txt
Betelgeuse_HK_ARCA_JUNO: makeskymapadvanced.py
	{ time python makeskymapadvanced.py -88.8 7.4 $(NSIDE) $(NTOY) HK ARCA JUNO; } &> log$@.txt
Cygnus_HK_ARCA_JUNO: makeskymapadvanced.py
	{ time python makeskymapadvanced.py -45 40 $(NSIDE) $(NTOY) HK ARCA JUNO; } &> log$@.txt
clean:
	rm deltaT*png *txt
