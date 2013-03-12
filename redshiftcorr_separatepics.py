## splines
## with 8 spectra

import os
import commands
import pyfits
from pylab import *
from scipy import *
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
'''
#If you want to download all the spectra use below code
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=2231&fiber=10&mjd=53816"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=2236&fiber=16&mjd=53729"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=2236&fiber=5&mjd=53729"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=1955&fiber=2&mjd=53442"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=4532&fiber=20&mjd=55559"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=3253&fiber=20&mjd=54941"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=2151&fiber=9&mjd=54523"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=2251&fiber=3&mjd=53557"')
'''
def __main__():
    tab1a = pyfits.open(commands.getoutput("pwd")+'/spec-2231-53816-0010.fits')
    tab2a = pyfits.open(commands.getoutput("pwd")+'/spec-2236-53729-0016.fits')
    tab3a = pyfits.open(commands.getoutput("pwd")+'/spec-2236-53729-0005.fits')
    tab4a = pyfits.open(commands.getoutput("pwd")+'/spec-1955-53442-0002.fits')
    tab5a = pyfits.open(commands.getoutput("pwd")+'/spec-4532-55559-0020.fits')
    tab6a = pyfits.open(commands.getoutput("pwd")+'/spec-3253-54941-0020.fits')
    tab7a = pyfits.open(commands.getoutput("pwd")+'/spec-2151-54523-0009.fits')
    tab8a = pyfits.open(commands.getoutput("pwd")+'/spec-2251-53557-0003.fits')

    taba=[tab1a,tab2a,tab3a,tab4a,tab5a,tab6a,tab7a,tab8a] #tabarray for simpler coding

    for i in range(0,8):
        z=taba[i][3].data.field(5)
        z=[x for x in z if x!=0] #remove nonzero options
        zm=mean(z)

        flux=taba[i][1].data.field(0)
        loglam=taba[i][1].data.field(1)
        lam=10**loglam
       
        #ivar=taba[i][1].data.field(2)

        axvline(x=6563, linewidth=2, color='r') #alpha Balmer line
        axvline(x=4861, linewidth=2, color='r')
        #axvline(x=4341, linewidth=2, color='r')
        #axvline(x=4102, linewidth=2, color='r')
        #axvline(x=3970, linewidth=2, color='r')
        #axvline(x=3889, linewidth=2, color='r')
        #axvline(x=3835, linewidth=2, color='r')
        #axvline(x=3646, linewidth=2, color='r')

        lamcor = zeros(len(flux)) #array for corrected wavelengths


        for j in range(len(flux)): #redshift correction
            lamcor[j]=float(lam[j])/float((1+zm))
    
        
        s = UnivariateSpline(lamcor, flux, k=3) #fit spline, k=3 means cubic
        fluxspline=s(lamcor) #flux modeled by spline

        lamcora = zeros(len(flux)) #array for +redshift in other direction

        for j in range(len(flux)):
            lamcora[j]=float(lam[j])/float((1-zm))
    
        print zm,i
        plt.plot(lamcora,flux, color='g', label='corrected in wrong direction') #plot mis-corrected lambda
   
        plt.plot(lamcor,flux, color='b', label='corrected') #plot corrected lambda
    
        plt.plot(lam,flux, color='k',label='uncorrected') #plot uncorrected lambda
        
        #plt.step(lamcor,flux+(2*sd), color='blue')
        #plt.step(lamcor,flux-(2*sd), color='yellow')
        # plt.vlines(#insert x, ymin, ymax, color='k', linestyles='solid'
        #savefig(str(plateid).zfill(4)+'-'+str(mjd)+'-'+str(fiber).zfill(4)+'_corrected.png')
        yax=[] #to make the computer pick a reasonable axes
        for i in range(len(lamcora)):
            if lamcora[i]>4845 and lamcora[i]<4880:
                yax.append(flux[i])
                
        plt.xlim(4845,4880)
        plt.ylim(0.85*min(yax),1.15*max(yax))
        plt.xlabel('Wavelength [Angstroms]')
        plt.ylabel('Flux [10^-17 erg/cm^2/s/Angstrom]')
        
        plt.legend(bbox_to_anchor=(1,1), loc=1, borderaxespad=0.)
        #show()

__main__()
