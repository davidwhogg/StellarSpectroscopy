###this code will pull out the fits file from api, then
### use the data in it to generate flux values then graph it

import os
import commands
import pyfits
from pylab import *
from scipy import *
import matplotlib.pyplot as plt

plateid=input("Plate ID (e.g. 276): ")
mjd=input("mjd (e.g. 51909): ")
fiber=input("fiber (e.g. 627): ")

commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate='+str(plateid)+'&fiber='+str(fiber)+'&mjd='+str(mjd)+'"')

tab = pyfits.getdata(commands.getoutput("pwd")+'/spec-'+str(plateid).zfill(4)+'-'+str(mjd)+'-'+str(fiber).zfill(4)+'.fits')

flux=tab.field('flux')
loglam=tab.field('loglam')

lam=10**loglam


plt.plot(lam,flux,'.')
plt.xlabel('Wavelength [Angstroms]')
plt.ylabel('Flux [10^-17 erg/cm^2/s/Angstrom]')
savefig(str(plateid).zfill(4)+'-'+str(mjd)+'-'+str(fiber).zfill(4)+'.png')
show()
