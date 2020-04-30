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
sourcepos = [xsource, ysource, zsource]

#distance between different detectors
posdiff1 = [xIC-xKM3, yIC-yKM3, zIC-zKM3]
posdiff2 = [xSK-xKM3, ySK-yKM3, zSK-zKM3]
posdiff3 = [xIC-xSK, yIC-ySK, zIC-zSK]

posdiff4 = [xJUNO-xKM3, yJUNO-yIC, zJUNO-zKM3]
posdiff5 = [xJUNO-xIC, yJUNO-yIC, zJUNO-zIC]
posdiff6 = [xJUNO-xSK, yJUNO-ySK, zJUNO-zSK]

#theoretical time delay for each detector pair from the true source
tdelay_true1 = np.dot( np.array(posdiff1), np.array(sourcepos) )/c
tdelay_true2 = np.dot( np.array(posdiff2), np.array(sourcepos) )/c
tdelay_true3 = np.dot( np.array(posdiff3), np.array(sourcepos) )/c
tdelay_true4 = np.dot( np.array(posdiff4), np.array(sourcepos) )/c
tdelay_true5 = np.dot( np.array(posdiff5), np.array(sourcepos) )/c
tdelay_true6 = np.dot( np.array(posdiff6), np.array(sourcepos) )/c

#estimated uncertainties
sigmatICKM3 = 0.00665
sigmatICHK = 0.00055
sigmatKM3HK = 0.0067
sigmatKM3SK = 0.0074
sigmatICSK = 0.00195
sigmatICJUNO = 0.00195
sigmatHKJUNO = 0.00199
sigmatKM3JUNO = 0.0074
sigmatSKJUNO = 0.00275

#define the resolution of the pixel of the healpix map
NSIDE = 256
NPIX = hp.nside2npix(NSIDE)
pixels = np.arange(NPIX)

#compute chi2 for each pixel, scanning all the sky
print("computing chi2...")
chi2_list_total = []
for pixel in range(0,NPIX):

        dec, ra = hp.pix2ang(nside=NSIDE, ipix=pixel)
        ra_test = ra*(180./np.pi)
        dec_test = (dec-np.pi/2.)*(180./np.pi)

        xfind = np.cos(ra_test*(np.pi/180))*np.cos(dec_test*(np.pi/180))
        yfind = np.sin(ra_test*(np.pi/180))*np.cos(dec_test*(np.pi/180))
        zfind = np.sin(dec_test*(np.pi/180))
        findpos = [xfind,yfind,zfind]

        tdelay_obs1 = np.dot( np.array(posdiff1), np.array(findpos) )/c
        tdelay_obs2 = np.dot( np.array(posdiff2), np.array(findpos) )/c
        tdelay_obs3 = np.dot( np.array(posdiff3), np.array(findpos) )/c
        tdelay_obs4 = np.dot( np.array(posdiff4), np.array(findpos) )/c
        tdelay_obs5 = np.dot( np.array(posdiff5), np.array(findpos) )/c
        tdelay_obs6 = np.dot( np.array(posdiff6), np.array(findpos) )/c

        chi2_ICKM3 = ( (tdelay_obs1 - tdelay_true1)/sigmatICKM3 )**2
        chi2_ICSK = ( (tdelay_obs3 - tdelay_true3)/sigmatICHK )**2
        chi2_KM3SK = ( (tdelay_obs2 - tdelay_true2)/sigmatKM3HK )**2
        chi2_JUNO_KM3 = ( (tdelay_obs4 - tdelay_true4)/sigmatKM3JUNO )**2
        chi2_JUNO_IC = ( (tdelay_obs5 - tdelay_true5)/sigmatICJUNO )**2
        chi2_JUNO_SK = ( (tdelay_obs6 - tdelay_true6)/sigmatHKJUNO )**2

        chi2 = chi2_KM3SK + chi2_ICSK + chi2_ICKM3 + chi2_JUNO_KM3 + chi2_JUNO_SK + chi2_JUNO_IC
        #chi2 = chi2_ICSK + chi2_JUNO_IC + chi2_JUNO_SK
        #chi2 = chi2_ICKM3 + chi2_KM3SK + chi2_ICSK
        #chi2 = chi2_ICKM3 + chi2_JUNO_IC + chi2_JUNO_KM3
        #chi2 = chi2_KM3SK + chi2_JUNO_SK + chi2_JUNO_KM3
        chi2_list_total.append(chi2)

# convert chi2 to p-value
ndof = 6 # ndof = number of detector pairs in the tot chi2 sum
         #(6 for comibnation of 4 detectors, 3 when combining 3 detectors)
chi2_stat = scipy.stats.distributions.chi2.cdf(chi2_list_total,ndof)

# Trey to add Galactic Plane map from Planck
#fermi = hp.read_map('../HFI_Mask_GalPlane_2048_R1.10.fits')

#Conversion needed from p-value to confidence level and to have good RA,dec convention
chi2_stat = [value*100 for value in reversed(chi2_stat)]
data=np.array(chi2_stat)

#compute the 90% and 68% CL area
contour = np.ma.masked_values(data, data <= 90)
area = np.sum(data <= 90) * hp.nside2pixarea(NSIDE, degrees=True)
contour2 = np.ma.masked_values(data, data <= 68)
area2 = np.sum(data <= 68) * hp.nside2pixarea(NSIDE, degrees=True)
print('90% CL area: ', area, 'deg2')

print("making plot...")
# plot CL controus and true position
hp.mollview(data, norm=None, min=0, max=100, unit='Confidence Level %', cmap='tab10', title='')
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
exit()
