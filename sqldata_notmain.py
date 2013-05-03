import numpy as np
import os
from pylab import *

import commands
import pyfits
import sys
os.chdir("/Users/admin/Desktop/Maser files")
directory=commands.getoutput("pwd")
sys.path.append(directory)
import sqlcl
def sqldata(n):
    objids=[] #this one is left empty, ignore

    objs=[] #8
    teff=[] #2
    feh=[] #3
    extinction=[] #1
    names=[] #4
    fiber=[] #10
    mjd=[] #11
    plate=[] #9
    camcol=[] #5
    run=[] #6
    ids=[] #7
    ras=[] #12
    decs=[] #13
    magu=[] #14
    magg=[]#15
    magr=[] #16
    magi=[] #17
    magz=[] #18

    array1=[objids,extinction,teff,feh,names,camcol,run,ids,objs,plate,fiber,mjd,\
            ras,decs,magu,magg,magr,magi,magz]

    alldata= sqlcl.query("SELECT top 1000 p.objID, \
    p.extinction_g, s.elodieTEff, s.elodieFeH, s.elodieObject, p.camcol, p.run, p.field, \
    p.obj, s.plate, s.fiberID, s.mjd, p.ra, p.dec, \
    p.psfMag_u, p.psfMag_g, psfMag_r, psfMag_i, psfMag_z \
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
    for i in range(20,len(compiled)-1,19):
        extinction.append(float(compiled[i]))
    for i in range(21,len(compiled)-1,19):
        teff.append(float(compiled[i]))
    for i in range(22,len(compiled)-1,19):
        feh.append(float(compiled[i]))
    for i in range(23,len(compiled)-1,19):
        names.append((compiled[i]))
    for i in range(24,len(compiled)-1,19):
        objs.append(float(compiled[i]))
    for i in range(25,len(compiled)-1,19):
        camcol.append(float(compiled[i]))
    for i in range(26,len(compiled)-1,19):
        run.append(int(compiled[i]))
    for i in range(27,len(compiled)-1,19):
        ids.append(int(compiled[i]))
    for i in range(28,len(compiled)-1,19):
        plate.append(int(compiled[i]))
    for i in range(29,len(compiled)-1,19):
        fiber.append(int(compiled[i]))
    for i in range(30,len(compiled)-1,19):
        mjd.append(int(compiled[i]))
    for i in range(31,len(compiled)-1,19):
        ras.append(float(compiled[i]))
    for i in range(32,len(compiled)-1,19):
        decs.append(float(compiled[i]))
    for i in range(33,len(compiled)-1,19):
        magu.append(float(compiled[i]))
    for i in range(34,len(compiled)-1,19):
        magg.append(float(compiled[i]))
    for i in range(35,len(compiled)-1,19):
        magr.append(float(compiled[i]))
    for i in range(36,len(compiled)-1,19):
        magi.append(float(compiled[i]))
    for i in range(37,len(compiled)-1,19):
        magz.append(float(compiled[i]))
        '''
    tabs=[] #this will contain each fits file in one super-array
    fluxes=[]
    sn2s=[]
    errormags=[]

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
        errormags.append(errormag)'''
    return array1

