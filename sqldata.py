import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os
from pylab import *

import commands
import pyfits
import sys
directory=commands.getoutput("pwd")
sys.path.append(directory)
import sqlcl

objid=[]
teff=[]
feh=[]
extinction=[]
names=[]
fiber=[]
mjd=[]
plate=[]
camcol=[]
run=[]
ids=[]

alldata= sqlcl.query("SELECT top 10 p.objID, \
p.extinction_g, s.elodieTEff, s.elodieFeH, s.elodieObject, p.camcol, p.run, p.field, \
p.obj, s.plate, s.fiberID, s.mjd \
FROM PhotoObj AS p \
JOIN SpecObj as s ON s.specobjID=p.specobjID \
WHERE psfMag_r<19  \
and psfMag_r>15 and type=6 \
and  dbo.fPhotoStatus('PRIMARY')>0 and dbo.fPhotoFlags('STATIONARY')>0 \
and calibStatus_r=1 \
and s.elodieTEff!=0 and s.elodieFeH!=0 and s.elodieLogG!=0 \
and ((flags&dbo.fPhotoFlags('BLENDED')) \
+(flags&dbo.fPhotoFlags('DEBLEND_TOO_MANY_PEAKS')) + \
(flags&dbo.fPhotoFlags('SATURATED')) \
+(flags&dbo.fPhotoFlags('BADSKY'))+ \
(flags&dbo.fPhotoFlags('COSMIC_RAY')) \
+(flags&dbo.fPhotoFlags('PEAKS_TOO_CLOSE'))+ \
(flags&dbo.fPhotoFlags('NOTCHECKED_CENTER')) \
+(flags&dbo.fPhotoFlags('SATUR_CENTER'))+ \
(flags&dbo.fPhotoFlags('INTERP_CENTER')) \
+(flags&dbo.fPhotoFlags('INTERP'))+ \
(flags&dbo.fPhotoFlags('PSF_FLUX_INTERP')))=0 \
and sqrt((power(psfMag_u-psfmag_g-0.82,2)+ \
power(psfMag_g-psfMag_r-0.3,2)+ \
power(psfMag_r-psfMag_i-0.09,2)+\
power(psfMag_i-psfMag_z-0.02,2)))<0.08").read()
interim=alldata.replace("\n",",")
compiled=interim.split(",")
for i in range(13,len(compiled)-1,12):
    extinction.append(compiled[i])
for i in range(14,len(compiled)-1,12):
    teff.append(compiled[i])
for i in range(15,len(compiled)-1,12):
    feh.append(compiled[i])
for i in range(16,len(compiled)-1,12):
    names.append(compiled[i])
for i in range(17,len(compiled)-1,12):
    objid.append(compiled[i])
    
for i in range(18,len(compiled)-1,12):
    camcol.append(compiled[i])
for i in range(19,len(compiled)-1,12):
    run.append(compiled[i])
for i in range(20,len(compiled)-1,12):
    ids.append(compiled[i])
for i in range(21,len(compiled)-1,12):
    plate.append(compiled[i])
for i in range(22,len(compiled)-1,12):
    fiber.append(compiled[i])
for i in range(23,len(compiled)-1,12):
    mjd.append(compiled[i])

tabs=[] #this will contain each fits file in one super-array
fluxes=[]
sn2s=[]

for i in range(len(camcol)):
    plateid=plate[i]
    mjdid=mjd[i]
    fiberid=fiber[i]
    commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate='+plateid+'&fiber='+fiberid+'&mjd='+mjdid+'"')
    
    tab = pyfits.open(commands.getoutput("pwd")+'/spec-'+plateid.zfill(4)+'-'+mjdid+'-'+fiberid.zfill(4)+'.fits')
    tabs.append(tab)
    flux=tabs[i][1].data.field(0) #flux
    fluxes.append(flux)
    sn2=tabs[i][2].data.field(6)[0]+tabs[i][2].data.field(7)[0] #one of these entries is 0 always
    sn2s.append(sn2)
    errormag=1/sn2
