#the below code just takes two spectra for now but can be easily modified.

import os
import commands
import pyfits
from pylab import *
from scipy import *
import matplotlib.pyplot as plt

#below section is just because my computer refuses to find the sqlcl module 
#import sys
#directory=commands.getoutput("pwd")
#sys.path.append(directory)
import sqlcl

#def __main__():
for i in range(1):

    #below 2 lines just downloads the fits files
    commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=2236&fiber=16&mjd=53729"')
    commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=2231&fiber=10&mjd=53816"')

    tab1a = pyfits.open(commands.getoutput("pwd")+'/spec-2231-53816-0010.fits')
    tab2a = pyfits.open(commands.getoutput("pwd")+'/spec-2236-53729-0016.fits')


    tab=[tab1a,tab2a]
    fluxes=[]   #this and the below lists will contain all the 
    redshifts=[] #needed values for graphing spectra and for
    temps=[]    # the linear discriminant analysis
    sn2s=[]
    wavelengths=[]
    correctedwave=[]
    objids=[]
    logg=[]
    teff=[]
    feh=[]
    extvalues=[]
    

    #below section sorts the data into the lists from above
    for i in range(0,2):
        if type(tab[i][2].data.field(63)[0])==float32: #distinguish SDSS,BOSS
                zm= tab[i][2].data.field(63)[0]
                print  "SDSS - OK"            
        elif type(tab[i][2].data.field(37)[0])==float32:
                zm=	tab[i][2].data.field(37)[0]
                print  "BOSS - error" #BOSS has different order of elements in table
        redshifts.append(zm)

        temp=tab[i][2].data.field(106)[0] #effective temperature
        temps.append(temp)

        flux=tab[i][1].data.field(0) #flux
        fluxes.append(flux)

        sn2=tab[i][2].data.field(6)[0]+tab1a[2].data.field(7)[0] #one of these entries is 0 always
        sn2s.append(sn2)
        errormag=1/sn2

        loglam=tab[i][1].data.field(1) #log(wavelength)
        lam=10**loglam
        wavelengths.append(lam)
        
        objid=tab[i][2].data.field(21)[0] #objID for linking with CAS
        objids.append(objid)
        
        lamcor = zeros(len(flux)) #array for corrected wavelengths
        for j in range(len(flux)): #redshift correction
            lamcor[j]=float(lam[j])/float((1+zm))
        correctedwave.append(lamcor)

            
        alldata= sqlcl.query("SELECT top 100 p.objID, \
        p.extinction_r, s.elodieTEff, s.elodieFeH, s.ElodieLogG \
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
        for i in range(7,len(compiled)-1,5):
            teff.append(compiled[i])
        for i in range(8,len(compiled)-1,5):
            feh.append(compiled[i])
        for i in range(9,len(compiled)-1,5):
            logg.append(compiled[i])
            
        extinctions= sqlcl.query("SELECT top 1 objID, extinction_u, extinction_r, \
        extinction_g, extinction_z, extinction_i FROM PhotoObj WHERE objID="+objid).read()
        extinctions_split=extinctions.split("\n")
        extinctions_values=extinctions_split[1].split(",")
        if extinctions_values[0]!=objid: #check that splitting was successful
            print "Error in splitting"
            sys.exit()
        else:
            ext=extinctions_values[1:6]
            extvalues.append(ext) #five values in each array; UGRIZ

#__main__()

