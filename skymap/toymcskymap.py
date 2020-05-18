#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib import ticker, cm
import healpy as hp
#from mpl_toolkits.basemap import Basemap,shiftgrid
from matplotlib.colors import LogNorm
import scipy.stats

# Defnition of the True position of GC
ra_true = (-94.4)*(np.pi/180.)
dec_true = (-28.94)*(np.pi/180.)

c = 3e8 #speed of light
Rearth = 6.4e6 #earth radius in m

#Galactic coordinates of the detectors
lonKM3 = 16*(np.pi/180)
latKM3 = 0.632973
lonIC = -63.453056*(np.pi/180)
latIC = -89.99*(np.pi/180)
latSK = 36*(np.pi/180)
lonSK = 129*(np.pi/180)
latJUNO = 22.11827*(np.pi/180)
lonJUNO = 112.51867*(np.pi/180)

#positions of the detectors and true source
xKM3 = Rearth*np.cos(lonKM3)*np.cos(latKM3)
yKM3 = Rearth*np.sin(lonKM3)*np.cos(latKM3)
zKM3 = Rearth*np.sin(latKM3)

xIC = Rearth*np.cos(lonIC)*np.cos(latIC)
yIC = Rearth*np.sin(lonIC)*np.cos(latIC)
zIC = Rearth*np.sin(latIC)

xSK = Rearth*np.cos(lonSK)*np.cos(latSK)
ySK = Rearth*np.sin(lonSK)*np.cos(latSK)
zSK = Rearth*np.sin(latSK)

xJUNO = Rearth*np.cos(lonJUNO)*np.cos(latJUNO)
yJUNO = Rearth*np.sin(lonJUNO)*np.cos(latJUNO)
zJUNO = Rearth*np.sin(latJUNO)

xsource = np.cos(ra_true)*np.cos(dec_true)
ysource = np.sin(ra_true)*np.cos(dec_true)
zsource = np.sin(dec_true)

posKM3 = [xKM3, yKM3, zKM3]
posIC = [xIC, yIC, zIC]
posSK = [xSK, ySK, zSK]
posJUNO = [xJUNO,yJUNO,zJUNO]
sourcepos = np.array([xsource, ysource, zsource])

#distance between different detectors
posdiffICKM3 = np.array([xIC-xKM3, yIC-yKM3, zIC-zKM3])
posdiffHKKM3 = np.array([xSK-xKM3, ySK-yKM3, zSK-zKM3])
posdiffSKM3 = posdiffHKKM3
posdiffICHK = np.array([xIC-xSK, yIC-ySK, zIC-zSK])
posdiffICSK = posdiffICHK
posdiffJUNOKM3 = np.array([xJUNO-xKM3, yJUNO-yIC, zJUNO-zKM3])
posdiffJUNOIC = np.array([xJUNO-xIC, yJUNO-yIC, zJUNO-zIC])
posdiffJUNOHK = np.array([xJUNO-xSK, yJUNO-ySK, zJUNO-zSK])
posdiffJUNOSK = posdiffJUNOHK

#theoretical time delay for each detector pair from the true source
tdelay_trueICKM3= np.dot( posdiffICKM3, sourcepos )/c
tdelay_trueHKKM3 = np.dot( posdiffHKKM3, sourcepos )/c
tdelay_trueICHK = np.dot( posdiffICHK, sourcepos )/c
tdelay_trueJUNOIC = np.dot( posdiffJUNOIC, sourcepos )/c
tdelay_trueJUNOHK = np.dot( posdiffJUNOHK, sourcepos )/c
tdelay_trueJUNOKM3 = np.dot( posdiffJUNOKM3, sourcepos )/c

#estimated uncertainties
sigmatICKM3 = 0.00665
sigmatICHK = 0.00055
sigmatHKKM3 = 0.0067
sigmatKM3SK = 0.0074
sigmatICSK = 0.00195
sigmatJUNOIC = 0.00195
sigmatJUNOHK = 0.00199
sigmatJUNOKM3 = 0.0074
sigmatSKJUNO = 0.00275

#define the resolution of the pixel of the healpix map
#NSIDE = 256
NSIDE = 128
NPIX = hp.nside2npix(NSIDE)
print("pixels:",NPIX)

chi2_stat = np.zeros(NPIX)

ntoy=10000

coord = np.array([hp.pix2ang(nside=NSIDE, ipix=pixel) for pixel in range(0, NPIX)])  #dec, ra
coord[:,0]=coord[:,0]-np.pi/2.

for itoy in range(0,ntoy):
    print("toy:",itoy)
    tdelay_ICKM3 = random.gauss(tdelay_trueICKM3,sigmatICKM3)
    tdelay_HKKM3 = random.gauss(tdelay_trueHKKM3,sigmatHKKM3)
    tdelay_ICHK = random.gauss(tdelay_trueICHK,sigmatICHK)
    tdelay_JUNOIC = random.gauss(tdelay_trueJUNOIC,sigmatJUNOIC)
    tdelay_JUNOHK = random.gauss(tdelay_trueJUNOHK,sigmatJUNOHK)
    tdelay_JUNOKM3 = random.gauss(tdelay_trueJUNOKM3,sigmatJUNOKM3)

    #compute chi2 for each pixel, scanning all the sky
    chi2_list_min = 1e20
    bestpixel=-1

    findposX = np.cos(coord[:,1])*np.cos(coord[:,0])
    findposY = np.sin(coord[:,1])*np.cos(coord[:,0])
    findposZ = np.sin(coord[:,0])

    tdelay_obsICKM3 = (posdiffICKM3[0]*findposX+posdiffICKM3[1]*findposY+posdiffICKM3[2]*findposZ)/c
    tdelay_obsHKKM3 = (posdiffHKKM3[0]*findposX+posdiffHKKM3[1]*findposY+posdiffHKKM3[2]*findposZ)/c
    tdelay_obsICHK = (posdiffICHK[0]*findposX+posdiffICHK[1]*findposY+posdiffICHK[2]*findposZ)/c
    tdelay_obsJUNOIC = (posdiffJUNOIC[0]*findposX+posdiffJUNOIC[1]*findposY+posdiffJUNOIC[2]*findposZ)/c
    tdelay_obsJUNOHK = (posdiffJUNOHK[0]*findposX+posdiffJUNOHK[1]*findposY+posdiffJUNOHK[2]*findposZ)/c
    tdelay_obsJUNOKM3 = (posdiffJUNOKM3[0]*findposX+posdiffJUNOKM3[1]*findposY+posdiffJUNOKM3[2]*findposZ)/c

    chi2_ICKM3 = ( (tdelay_obsICKM3 - tdelay_ICKM3)/sigmatICKM3 )**2
    chi2_HKKM3 = ( (tdelay_obsHKKM3 - tdelay_HKKM3)/sigmatHKKM3 )**2
    chi2_ICHK = ( (tdelay_obsICHK - tdelay_ICHK)/sigmatICHK )**2
    chi2_JUNOIC = ( (tdelay_obsJUNOIC - tdelay_JUNOIC)/sigmatJUNOIC )**2
    chi2_JUNOHK = ( (tdelay_obsJUNOHK - tdelay_JUNOHK)/sigmatJUNOHK )**2
    chi2_JUNOKM3 = ( (tdelay_obsJUNOKM3 - tdelay_JUNOKM3)/sigmatJUNOKM3 )**2

    chi2 = chi2_ICKM3 + chi2_HKKM3 + chi2_ICHK + chi2_JUNOIC + chi2_JUNOHK + chi2_JUNOKM3
    chi2_stat[np.where(chi2 == np.amin(chi2))[0]] += 1


#Conversion to have good RA,dec convention
chi2_stat = [float(value)/float(ntoy)*100 for value in reversed(chi2_stat)]
data=np.array(chi2_stat)

data_sorted = np.sort(data)[::-1]
integral = 0
area = 0
for i in data_sorted:
    print("calc:",i)
    integral += i
    area += hp.nside2pixarea(NSIDE, degrees=True)
    if integral > 68: break
print("68% conf area:",area)

#area = np.sum(data <= 90) * hp.nside2pixarea(NSIDE, degrees=True)
#print('90% CL area: ', area, 'deg2')

print("making plot...")
# plot CL controus and true position
hp.mollview(data, norm=None, min=0, max=20, unit='Fitted values distribution, %', cmap='coolwarm', title='')
hp.graticule()
hp.projscatter(dec_true-np.pi/2.,-ra_true+np.pi, color='black')

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

plt.savefig("skymap.png")
#exit()
