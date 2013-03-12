## splines

import os
import commands
import pyfits
from pylab import *
from scipy import *
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

print "Let's try a star."

def __main__():
    
    plateid=input("Plate ID (e.g. 1022): ")
    mjd=input("MJD (e.g. 52524): ")
    fiber=input("Fiber (e.g. 5): ")

    commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate='+str(plateid)+'&fiber='+str(fiber)+'&mjd='+str(mjd)+'"')

    tab1 = pyfits.open(commands.getoutput("pwd")+'/spec-'+str(plateid).zfill(4)+'-'+str(mjd)+'-'+str(fiber).zfill(4)+'.fits')
    tab=tab1[1].data
    tab3=tab1[3].data
    z1=tab3.field(5) #redshift
    z1=[x for x in z1 if x!=0] #option to correct out zero values, should be redundant
    z=mean(z1) #redflux value
    flux=tab.field(0)
    loglam=tab.field(1)
    ivar=tab.field(2)
    lam=10**loglam
    sd=sqrt(1/ivar)
    
    
    axvline(x=6563, linewidth=2, color='r') #alpha Balmer line
    axvline(x=4861, linewidth=2, color='r')
    axvline(x=4341, linewidth=2, color='r')
    axvline(x=4102, linewidth=2, color='r')
    axvline(x=3970, linewidth=2, color='r')
    axvline(x=3889, linewidth=2, color='r')
    axvline(x=3835, linewidth=2, color='r')
    axvline(x=3646, linewidth=2, color='r')
    lamcor = zeros(len(flux)) #array for corrected wavelengths
    for j in range(len(flux)):
        lamcor[j]=float(lam[j])/float((1+z))
    s = UnivariateSpline(lamcor, flux, s=3)
    fluxspline=s(flux)
    plt.step(lamcor,fluxspline, color='red')
    plt.step(lamcor,flux+(2*sd), color='blue')
    plt.step(lamcor,flux-(2*sd), color='yellow')
    plt.scatter(lamcor,flux, color='green')
   # plt.vlines(#insert x, ymin, ymax, color='k', linestyles='solid'
    plt.xlabel('Corrected Wavelength [Angstroms]')
    plt.ylabel('Flux [10^-17 erg/cm^2/s/Angstrom]')
    savefig(str(plateid).zfill(4)+'-'+str(mjd)+'-'+str(fiber).zfill(4)+'_corrected.png')
    show()
    
    
__main__()
    
