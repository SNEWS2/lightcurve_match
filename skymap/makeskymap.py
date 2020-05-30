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
    ra_deg = float(sys.argv[1])
    dec_deg = float(sys.argv[2]) 
    NSIDE = int(sys.argv[3])
    ntoy = int(sys.argv[4])
    activedetectors = sys.argv[5:]
    print(activedetectors)
else:
    ra_deg = -94.4
    dec_deg = -28.94
    NSIDE = 64
    ntoy = 1000
    activedetectors = ["IC","ARCA","HK","JUNO"]

ra_true = np.radians(ra_deg)
dec_true = np.radians(dec_deg) 
xsource = np.cos(ra_true)*np.cos(dec_true)
ysource = np.sin(ra_true)*np.cos(dec_true)
zsource = np.sin(dec_true)
#unit vector that indicates the emission direction
sourcedir = -np.array([xsource, ysource, zsource])

c = 3e8 #speed of light in [m/s]
Rearth = 6.4e6 #earth radius in [m]

cls = [0.68,0.9] #confidence levels in growing order

#Galactic scancoordinates of the detectors
lonKM3 = 16*(np.pi/180)
latKM3 = 0.632973
lonIC = -63.453056*(np.pi/180)
latIC = -89.99*(np.pi/180)
lonSK = 129*(np.pi/180)
latSK = 36*(np.pi/180)
lonJUNO = 112.51867*(np.pi/180)
latJUNO = 22.11827*(np.pi/180)

detscancoord = {
        "ARCA": [lonKM3,latKM3],
        "IC": [lonIC,latIC],
        "SK": [lonSK,latSK],
        "HK": [lonSK,latSK],
        "JUNO": [lonJUNO,latJUNO]
        }

#estimated uncertainties in [s]
sigmatICKM3 = 0.00665
sigmatICHK = 0.00055
sigmatHKKM3 = 0.0067
sigmatKM3SK = 0.0074
sigmatICSK = 0.00195
sigmatJUNOIC = 0.00195
sigmatJUNOHK = 0.00199
sigmatJUNOKM3 = 0.0074
sigmatSKJUNO = 0.00275

detsigma = {
        frozenset(["IC","ARCA"]): sigmatICKM3,
        frozenset(["HK","IC"]): sigmatICHK,
        frozenset(["HK","ARCA"]): sigmatHKKM3,
        frozenset(["JUNO","IC"]): sigmatJUNOIC,
        frozenset(["JUNO","HK"]): sigmatJUNOHK,
        frozenset(["JUNO","ARCA"]): sigmatJUNOKM3,
        frozenset(["SK","ARCA"]): sigmatKM3SK,
        frozenset(["SK","IC"]): sigmatICSK,
        frozenset(["SK","JUNO"]): sigmatSKJUNO
        }

detvec = {}

for detname in activedetectors:
    detvec[detname] = [Rearth*np.cos(detscancoord[detname][0])*np.cos(detscancoord[detname][1]),
            Rearth*np.sin(detscancoord[detname][0])*np.cos(detscancoord[detname][1]),
            Rearth*np.sin(detscancoord[detname][1])]

posdiff = {} #vector connecting detector positions in [(m,m,m)]
truedelay = {} #expected delay for the source in [s]
for detpair in itertools.combinations(activedetectors, 2):
    posdiff[detpair] = np.array(detvec[detpair[0]])-np.array(detvec[detpair[1]])
    truedelay[detpair] = np.dot( posdiff[detpair], sourcedir )/c

#define the resolution of the pixel of the healpix map
NPIX = hp.nside2npix(NSIDE)
APIX = hp.nside2pixarea(NSIDE, degrees=True)
print("pixels:",NPIX,"size:",APIX,"deg^2")

areaarr = np.zeros([len(cls),ntoy])
ninside = np.zeros(len(cls)) #counter of number of toys for which the pixel containing the source is included to the confidence area (to check coverage)

fitted_distr = np.zeros(NPIX)

scancoord = np.array(hp.pix2ang(nside=NSIDE, ipix=range(0,NPIX))).transpose()
scancoord[:,0]=np.pi/2.-scancoord[:,0] #declination = 90-colatitude
scandirX = -np.cos(scancoord[:,1])*np.cos(scancoord[:,0])
scandirY = -np.sin(scancoord[:,1])*np.cos(scancoord[:,0])
scandirZ = -np.sin(scancoord[:,0])


true_pixel = hp.ang2pix(nside=NSIDE, theta=np.pi/2-dec_true, phi=ra_true+2*np.pi)
deltachi2 = scipy.stats.distributions.chi2.ppf(cls,2)

ploteachtoy = False

for itoy in range(0,ntoy):
    print("toy:",itoy,'/',ntoy,end='\r')
    measdelay = {} #measured delay in [s]
    scandelay = {} #vector of delays for each pixel
    chi2pair = {}
    chi2 = np.zeros(NPIX)
    for detpair in itertools.combinations(activedetectors, 2):
        measdelay[detpair] = random.gauss(truedelay[detpair],detsigma[frozenset(detpair)])
        scandelay[detpair] = (posdiff[detpair][0]*scandirX + 
                posdiff[detpair][1]*scandirY +
                posdiff[detpair][2]*scandirZ)/c
        chi2pair[detpair] = ( (scandelay[detpair] - measdelay[detpair])/detsigma[frozenset(detpair)] )**2
        chi2 = chi2 + chi2pair[detpair]

    chi2min = np.amin(chi2)
    fitted_distr[np.where(chi2 == chi2min)[0]] += 1
    chi2 = chi2 - chi2min

    #compute the 90% and 68% CL area                                                                              
    for icl in range(0,len(cls)):
        area = np.sum(chi2 <= deltachi2[icl]) * APIX
        areaarr[icl][itoy] = area
        if chi2[true_pixel] <=deltachi2[icl]: ninside[icl] += 1

    if ploteachtoy is True:
        cdf = scipy.stats.distributions.chi2.cdf(chi2,2)*100
        hp.mollview(cdf, norm=None, min=0, max=100, unit='Confidence level, %', cmap='tab10', title='', flip='geo')
        hp.graticule()
        hp.projscatter(np.pi/2.-dec_true,ra_true, color='black') #colatitude and longitude in radian
        name = "skymap_chi2_"+str(ra_deg)+"_"+str(dec_deg)+"_toy"+str(itoy)
        for detname in activedetectors:
            name+="_"+detname
        name+=".png"
        plt.savefig(name)



print("--------------------Results---------------")
for icl in range(0,len(cls)):
    print(cls[icl]*100,"% conf average area:",areaarr.mean(axis=1)[icl],"+-",areaarr.std(axis=1)[icl],"real coverage: ", ninside[icl]/float(ntoy)*100,"+-",np.sqrt(ninside[icl])/float(ntoy)*100, "%")

fitted_distr = fitted_distr/float(ntoy)  
data_sorted = np.sort(fitted_distr)[::-1]
integral = 0
area = 0
icl = 0
for i in data_sorted:
    #print("calc:",i)
    integral += i
    area += hp.nside2pixarea(NSIDE, degrees=True)
    if integral > cls[icl]: 
        print(cls[icl]*100,"% conf area from fitted positions distribution:",area)
        icl = icl + 1
        if icl == len(cls): break

#---------------------------------threating measured delays as true delays----------------------------------
measdelay = {} #measured delay in [s]
scandelay = {} #vector of delays for each pixel
chi2pair = {}
chi2 = np.zeros(NPIX)

for detpair in itertools.combinations(activedetectors, 2):
    measdelay[detpair] = truedelay[detpair]
    scandelay[detpair] = (posdiff[detpair][0]*scandirX + 
            posdiff[detpair][1]*scandirY +
            posdiff[detpair][2]*scandirZ)/c
    chi2pair[detpair] = ( (scandelay[detpair] - measdelay[detpair])/detsigma[frozenset(detpair)] )**2
    chi2 = chi2 + chi2pair[detpair]

#compute the 90% and 68% CL area                                                                              
for icl in range(0,len(cls)):
    area = np.sum(chi2 <= deltachi2[icl]) * APIX
    print(cls[icl]*100,'% CL area with true delays: ', area, 'deg2')


#----------------------------------------plotting---------------------------------------------------------
print("making plot...")
# plot CL controus and true position
fitted_distr=fitted_distr*100
fitted_distr = fitted_distr+1./float(ntoy)*10. #conversion to percents / 10
scale = 1./float(ntoy)
while np.amax(fitted_distr) > scale:
    scale *= 10
hp.mollview(fitted_distr, norm='log', min=np.amin(fitted_distr), max=scale, unit='Fitted values distribution, %', cmap='coolwarm', title='', flip='geo')
hp.graticule()
hp.projscatter(np.pi/2.-dec_true,ra_true, edgecolors='black', facecolors='none') #colatitude and longitude in radian
for detname in activedetectors:
    hp.projscatter(np.pi/2.-detscancoord[detname][1], detscancoord[detname][0], color="black", marker="s")

name = "skymap_fitpos_"+str(ra_deg)+"_"+str(dec_deg)+"_"+str(NSIDE)+"_"+str(ntoy)
for detname in activedetectors:
    name+="_"+detname
name+=".npy"

with open(name, 'wb') as f:
    np.save(f, fitted_distr)

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

#plt.show()

name = "skymap_fitpos_"+str(ra_deg)+"_"+str(dec_deg)+"_"+str(NSIDE)+"_"+str(ntoy)
for detname in activedetectors:
    name+="_"+detname
name+=".png"
plt.savefig(name)

cdf = scipy.stats.distributions.chi2.cdf(chi2,2)*100
hp.mollview(cdf, norm=None, min=0, max=100, unit='Confidence level, %', cmap='tab10', title='', flip='geo')
hp.graticule()
hp.projscatter(np.pi/2.-dec_true,ra_true, edgecolors='black', facecolors='none') #colatitude and longitude in radian
for detname in activedetectors:
    hp.projscatter(np.pi/2.-detscancoord[detname][1], detscancoord[detname][0], color="black", marker="s")

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

#plt.show()
name = "skymap_chi2_"+str(ra_deg)+"_"+str(dec_deg)+"_"+str(NSIDE)+"_"+str(ntoy)
for detname in activedetectors:
    name+="_"+detname
name+=".png"
plt.savefig(name)
exit()
