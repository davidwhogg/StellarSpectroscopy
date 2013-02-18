###this code will pull out the fits file from api, then
### use the data in it to generate flux values then graph it

import os
import commands
import pyfits
from pylab import *
from scipy import *
import matplotlib.pyplot as plt
print "Let's try a star."

def __main__():
    
    plateid=input("Plate ID (e.g. 1022): ")
    mjd=input("MJD (e.g. 52524): ")
    fiber=input("Fiber (e.g. 5): ")

    commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate='+str(plateid)+'&fiber='+str(fiber)+'&mjd='+str(mjd)+'"')

    tab = pyfits.getdata(commands.getoutput("pwd")+'/spec-'+str(plateid).zfill(4)+'-'+str(mjd)+'-'+str(fiber).zfill(4)+'.fits')

    flux=tab.field('flux')
    loglam=tab.field('loglam')
    ivar=tab.field('ivar')
    lam=10**loglam
    sd=sqrt(1/ivar)
    
    
    axvline(x=6563, linewidth=2, color='r')
    axvline(x=4861, linewidth=2, color='r')
    axvline(x=4341, linewidth=2, color='r')
    axvline(x=4102, linewidth=2, color='r')
    axvline(x=3970, linewidth=2, color='r')
    axvline(x=3889, linewidth=2, color='r')
    axvline(x=3835, linewidth=2, color='r')
    axvline(x=3646, linewidth=2, color='r')
    
   
    
    plt.step(lam,flux+(2*sd), color='blue')
    plt.step(lam,flux-(2*sd), color='yellow')
    plt.step(lam,flux, color='green')
   # plt.vlines(#insert x, ymin, ymax, color='k', linestyles='solid'
    plt.xlabel('Wavelength [Angstroms]')
    plt.ylabel('Flux [10^-17 erg/cm^2/s/Angstrom]')
    savefig(str(plateid).zfill(4)+'-'+str(mjd)+'-'+str(fiber).zfill(4)+'_lines.png')
    show()
    
__main__()
    
