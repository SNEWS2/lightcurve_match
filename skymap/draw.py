#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib import ticker, cm
import healpy as hp
from matplotlib.colors import LogNorm
import scipy.stats
import itertools
import sys

if len(sys.argv)>1:
    name = sys.argv[1]
    ra_deg = float(sys.argv[2])
    dec_deg = float(sys.argv[3])
#    NSIDE = int(sys.argv[4]) 
else:
    name = "skymap_fitpos_-45.0_40.0_256_100_HK_ARCA_JUNO.npy"
    ra_deg = -45.0
    dec_deg = 40.0
#    NSIDE = 256

#NPIX = hp.nside2npix(NSIDE)

ra_true = np.radians(ra_deg)
dec_true = np.radians(dec_deg)

print("making plot...")
#fitted_distr = np.zeros(NPIX)
with open(name, 'rb') as f:
    fitted_distr = np.load(f)

# plot CL controus and true position
hp.mollview(fitted_distr, norm=None, min=0, max=np.amax(fitted_distr), unit='Fitted values distribution, %', cmap='coolwarm', title='', flip='geo')
hp.graticule()
hp.projscatter(np.pi/2.-dec_true,ra_true, color='black') #colatitude and longitude in radian

#add axis labels
plt.text(2.0,0., r"$0^\circ$", ha="left", va="center")
plt.text(1.9,0.45, r"$30^\circ$", ha="left", va="center")
plt.text(1.4,0.8, r"$60^\circ$", ha="left", va="center")
plt.text(1.9,-0.45, r"$-30^\circ$", ha="left", va="center")
plt.text(1.4,-0.8, r"$-60^\circ$", ha="left", va="center")
plt.text(1.333, -0.15, r"$120^\circ$", ha="center", va="center")
plt.text(.666, -0.15, r"$60^\circ$", ha="center", va="center")
plt.text(0.0, -0.15, r"$0^\circ$", ha="center", va="center")
plt.text(-.666, -0.15, r"$-60^\circ$", ha="center", va="center")
plt.text(-1.333, -0.15, r"$-120^\circ$", ha="center", va="center")
plt.text(-2.0, -0.15, r"$-180^\circ$", ha="center", va="center")

plt.show()

name = name+".png"
plt.savefig(name)
