import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os
import commands
import sys
directory=commands.getoutput("pwd")
sys.path.append(directory)
import sqlcl

objid=[]
teff=[]
feh=[]
extinction=[]
names=[]
alldata= sqlcl.query("SELECT top 1000 p.objID, \
p.extinction_g, s.elodieTEff, s.elodieFeH, s.elodieObject \
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
for i in range(6,len(compiled)-1,5):
    extinction.append(compiled[i])
for i in range(7,len(compiled)-1,5):
    teff.append(compiled[i])
for i in range(8,len(compiled)-1,5):
    feh.append(compiled[i])
for i in range(9,len(compiled)-1,5):
    names.append(compiled[i])
for i in range(5,len(compiled)-1,5):
    objid.append(compiled[i])

plt.show()
