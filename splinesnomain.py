## splines
## with 8 spectra

import os
import commands
import pyfits
from pylab import *
from scipy import *
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

#If you want to download all the spectra use below code
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=2231&fiber=10&mjd=53816"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=2236&fiber=16&mjd=53729"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=2236&fiber=5&mjd=53729"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=1955&fiber=2&mjd=53442"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=4532&fiber=20&mjd=55559"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=3253&fiber=20&mjd=54941"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=2151&fiber=9&mjd=54523"')
commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate=2251&fiber=3&mjd=53557"')

def __main__():    
    tab1a = pyfits.open(commands.getoutput("pwd")+'/spec-2231-53816-0010.fits')
    tab2a = pyfits.open(commands.getoutput("pwd")+'/spec-2236-53729-0016.fits')
    tab3a = pyfits.open(commands.getoutput("pwd")+'/spec-2236-53729-0005.fits')
    tab4a = pyfits.open(commands.getoutput("pwd")+'/spec-1955-53442-0002.fits')
    tab5a = pyfits.open(commands.getoutput("pwd")+'/spec-4532-55559-0020.fits')
    tab6a = pyfits.open(commands.getoutput("pwd")+'/spec-3253-54941-0020.fits')
    tab7a = pyfits.open(commands.getoutput("pwd")+'/spec-2151-54523-0009.fits')
    tab8a = pyfits.open(commands.getoutput("pwd")+'/spec-2251-53557-0003.fits')

    tab1b=tab1a[1].data
    tab1c=tab1a[3].data
    tab2b=tab2a[1].data
    tab2c=tab2a[3].data
    tab3b=tab3a[1].data
    tab3c=tab3a[3].data
    tab4b=tab4a[1].data
    tab4c=tab4a[3].data
    tab5b=tab5a[1].data
    tab5c=tab5a[3].data
    tab6b=tab6a[1].data
    tab6c=tab6a[3].data
    tab7b=tab7a[1].data
    tab7c=tab7a[3].data
    tab8b=tab8a[1].data
    tab8c=tab8a[3].data

    z1=tab1c.field(5) #redshift
    z1=[x for x in z1 if x!=0] #option to correct out zero values, should be redundant
    z1m=mean(z1) #redflux value
    z2=tab2c.field(5)
    z2=[x for x in z2 if x!=0]
    z2m=mean(z2)
    z3=tab3c.field(5)
    z3=[x for x in z3 if x!=0]
    z3m=mean(z3)
    z4=tab4c.field(5)
    z4=[x for x in z4 if x!=0]
    z4m=mean(z4)
    z5=tab5c.field(5)
    z5=[x for x in z5 if x!=0]
    z5m=mean(z5)
    z6=tab6c.field(5)
    z6=[x for x in z6 if x!=0]
    z6m=mean(z6)
    z7=tab7c.field(5)
    z7=[x for x in z7 if x!=0]
    z7m=mean(z7)
    z8=tab8c.field(5)
    z8=[x for x in z8 if x!=0]
    z8m=mean(z8)


    flux1=tab1b.field(0)
    loglam1=tab1b.field(1)
    ivar1=tab1b.field(2)
    lam1=10**loglam1
    sd1=sqrt(1/ivar1)

    flux2=tab2b.field(0)
    loglam2=tab2b.field(1)
    ivar2=tab2b.field(2)
    lam2=10**loglam2
    sd2=sqrt(1/ivar2)

    flux3=tab3b.field(0)
    loglam3=tab3b.field(1)
    ivar3=tab3b.field(2)
    lam3=10**loglam3
    sd3=sqrt(1/ivar3)

    flux4=tab4b.field(0)
    loglam4=tab4b.field(1)
    ivar4=tab4b.field(2)
    lam4=10**loglam4
    sd4=sqrt(1/ivar4)

    flux5=tab5b.field(0)
    loglam5=tab5b.field(1)
    ivar5=tab5b.field(2)
    lam5=10**loglam5
    sd5=sqrt(1/ivar5)

    flux6=tab6b.field(0)
    loglam6=tab6b.field(1)
    ivar6=tab6b.field(2)
    lam6=10**loglam6
    #sd6=sqrt(1/ivar6)

    flux7=tab7b.field(0)
    loglam7=tab7b.field(1)
    ivar7=tab7b.field(2)
    lam7=10**loglam7
    sd7=sqrt(1/ivar7)

    flux8=tab8b.field(0)
    loglam8=tab8b.field(1)
    ivar8=tab8b.field(2)
    lam8=10**loglam8
    sd8=sqrt(1/ivar8)



    axvline(x=6563, linewidth=2, color='r') #alpha Balmer line
    axvline(x=4861, linewidth=2, color='r')
    axvline(x=4341, linewidth=2, color='r')
    axvline(x=4102, linewidth=2, color='r')
    axvline(x=3970, linewidth=2, color='r')
    axvline(x=3889, linewidth=2, color='r')
    axvline(x=3835, linewidth=2, color='r')
    axvline(x=3646, linewidth=2, color='r')

    plt.step(lam1,flux1, marker='', color='blue')
    plt.step(lam2,flux2, marker='', color='blue')
    plt.step(lam3,flux3, marker='', color='blue')
    plt.step(lam4,flux4, marker='', color='blue')
    plt.step(lam5,flux5, marker='', color='blue')
    plt.step(lam6,flux6, marker='', color='blue')
    plt.step(lam7,flux7, marker='', color='blue')
    plt.step(lam8,flux8, marker='', color='blue')


    lamcor1 = zeros(len(flux1)) #array for corrected wavelengths
    lamcor2 = zeros(len(flux2))
    lamcor3 = zeros(len(flux3))
    lamcor4 = zeros(len(flux4))
    lamcor5 = zeros(len(flux5))
    lamcor6 = zeros(len(flux6))
    lamcor7 = zeros(len(flux7))
    lamcor8 = zeros(len(flux8))

    for j in range(len(flux1)):
        lamcor1[j]=float(lam1[j])/float((1+z1m))
    for j in range(len(flux2)):
        lamcor2[j]=float(lam2[j])/float((1+z2m))
    for j in range(len(flux3)):
        lamcor3[j]=float(lam3[j])/float((1+z3m))
    for j in range(len(flux4)):
        lamcor4[j]=float(lam4[j])/float((1+z4m))
    for j in range(len(flux5)):
        lamcor5[j]=float(lam5[j])/float((1+z5m))
    for j in range(len(flux6)):
        lamcor6[j]=float(lam6[j])/float((1+z6m))
    for j in range(len(flux7)):
        lamcor7[j]=float(lam7[j])/float((1+z7m))
    for j in range(len(flux8)):
        lamcor8[j]=float(lam8[j])/float((1+z8m))
        
    s1 = UnivariateSpline(lamcor1, flux1, k=3)
    s2 = UnivariateSpline(lamcor2, flux2, k=3)
    s3 = UnivariateSpline(lamcor3, flux3, k=3)
    s4 = UnivariateSpline(lamcor4, flux4, k=3)
    s5 = UnivariateSpline(lamcor5, flux5, k=3)
    s6 = UnivariateSpline(lamcor6, flux6, k=3)
    s7 = UnivariateSpline(lamcor7, flux7, k=3)
    s8 = UnivariateSpline(lamcor8, flux8, k=3)

    fluxspline1=s1(lamcor1)
    fluxspline2=s2(lamcor2)
    fluxspline3=s3(lamcor3)
    fluxspline4=s4(lamcor4)
    fluxspline5=s5(lamcor5)
    fluxspline6=s6(lamcor6)
    fluxspline7=s7(lamcor7)
    fluxspline8=s8(lamcor8)

    plt.step(lamcor1,flux1, color='green')
    plt.step(lamcor2,flux2, color='green')
    plt.step(lamcor3,flux3, color='green')
    plt.step(lamcor4,flux4, color='green')
    plt.step(lamcor5,flux5, color='green')
    plt.step(lamcor6,flux6, color='green')
    plt.step(lamcor7,flux7, color='green')
    plt.step(lamcor8,flux8, color='green')


    plt.step(lamcor1,fluxspline1, color='black')
    plt.step(lamcor2,fluxspline2, color='black')
    plt.step(lamcor3,fluxspline3, color='black')
    plt.step(lamcor4,fluxspline4, color='black')
    plt.step(lamcor5,fluxspline5, color='black')
    plt.step(lamcor6,fluxspline6, color='black')
    plt.step(lamcor7,fluxspline7, color='black')
    plt.step(lamcor8,fluxspline8, color='black')


    #plt.step(lamcor,flux+(2*sd), color='blue')
    #plt.step(lamcor,flux-(2*sd), color='yellow')
    # plt.vlines(#insert x, ymin, ymax, color='k', linestyles='solid'
    #savefig(str(plateid).zfill(4)+'-'+str(mjd)+'-'+str(fiber).zfill(4)+'_corrected.png')


    plt.xlim(4845,4880)
    plt.ylim(0,300)
    plt.xlabel('Wavelength [Angstroms]')
    plt.ylabel('Flux [10^-17 erg/cm^2/s/Angstrom]')

    show()

__main__()
