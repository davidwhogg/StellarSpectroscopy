import scipy.optimize
from scipy import optimize
from scipy.interpolate import UnivariateSpline
from pylab import *
from scipy import stats
import numpy as np
import sys
import os
import commands
import pickle
import pyfits
os.chdir("/Users/admin/Desktop/Maser files")
directory=commands.getoutput("pwd")
sys.path.append(directory)

f=open("sorted","rb")
array1=pickle.load(f)
f.close()
f=open("datanewdr8bb","rb")
data=pickle.load(f)
f.close()

def deredshift(wls, fluxes, zm, badpoints):
        lamcorb=wls[0]/(1.0-zm)
        s = UnivariateSpline(wls[0], fluxes[0], k=3, s=0)
        xs=wls[0]
        ys=s(lamcorb)
        bady=s(badpoints)
        badx=badpoints
        return xs, ys, badx, bady
def getdata(plate,mjd,fiber):
        
    tabs=[0] #we only get one entry each time, we overwrite each time
    fluxes=[0]
    sn2s=[0]
    errormags=[0]
    wls=[0]
    ivars=[0]
    tab = pyfits.open(commands.getoutput("pwd")+'/FITS_files/spec-'+str(plate).zfill(4)+'-'+str(mjd)+'-'+str(fiber).zfill(4)+'.fits')
    #print "tab success", i
    tabs[0]=tab
    j=0 #-2700 ###########
    if type(tabs[0][2].data.field(63)[0])==float32: #distinguish SDSS,BOSS
        zm= tabs[0][2].data.field(63)[0] #redshift
        #print i, "63"            
    elif type(tabs[0][2].data.field(37)[0])==float32:
        zm= tabs[0][2].data.field(37)[0]
        print i, "37"
    else:
        print "error"
        
    flux=tabs[0][1].data.field(0) #flux
    fluxes[0]=flux
    loglam=tabs[0][1].data.field(1)
    ivar=tabs[0][1].data.field(2)
    ivars[0]=ivar
    
    loglam=np.array(loglam)
    lam=10**loglam
    lamcor=zeros(len(lam))
        
    wls[0]=lam
    sn2=tabs[j][2].data.field(6)[0]+tabs[j][2].data.field(7)[0] #one of these entries is 0 always
    sn2s[0]=sn2
    errormag=1/sn2
    errormags[0]=errormag
    badpoints=[]
    for v in range(len(ivars[0])):
        if ivars[0][v]!=0:
            ivars[0][v]=1/sqrt(ivars[0][v])
        elif ivars[0][v]==0:
            ivars[0][v]=np.inf
            badpoints.append(wls[0][v])
    del tab[0].data #free up memoery
    tab.close()
    return wls, fluxes, sn2s, ivars, badpoints, zm

def LinearFit(p,xsIn,ysIn,zsIn):#fit both
    return p[0]*xsIn+p[1]*ysIn + p[2]*zsIn
fitfunc=lambda p,f,fit_bright,fit_hdew,fit_ext,fit_flux: fabs(f(p,fit_bright,fit_hdew,fit_ext)-fit_flux)

    
ext=np.array(array1[1])
magg=np.array(array1[15])

plates=np.array(array1[9])
mjds=np.array(array1[11])
fibers=np.array(array1[10])


hdew=[]
ext2=[]
for i in range(len(magg)):
    ext2.append(data[0][i][5])
    hdew.append(data[0][i][4])
    
if ext2[0]!=ext[0]:
    print "Error! Data arrays do not align!"
    sys.exit()
       
bright= 10**(-0.4*(magg-22.5))

fit_flux_3900=[]
fit_flux_8000=[]
fit_hdew=[]
fit_ext=[]
fit_bright=[]
for i in range(len(plates)):
    plate=plates[i]
    mjd=mjds[i]
    fiber=fibers[i]
    wls, fluxes, sn2s, ivars, badpoints, zm = getdata(plate,mjd,fiber)
    xs, ys, badx, bady = deredshift(wls, fluxes, zm, badpoints)
    b=(wls[0]>3900)
    c=(wls[0]<8000)
    wls_used = np.trim_zeros(c*b*wls[0])
    flux_used = np.trim_zeros(c*b*ys)
    if wls_used[0]< 3900.4 and wls_used[-1]>7998.0 :
        fit_flux_3900.append(flux_used[0])
        fit_flux_8000.append(flux_used[-1])
        fit_hdew.append(hdew[i])
        fit_ext.append(ext[i])
        fit_bright.append(bright[i])
        print "yay"
    else:
        print "missing flux values", wls_used[0], wls_used[-1]
        pass


x0i=[0.5,0.5,0.5] #initial guess
x0=array(x0i) #optimize requires array
fit_flux_3900 = np.array(fit_flux_3900)
fit_flux_8000=np.array(fit_flux_8000)
fit_hdew=np.array(fit_hdew)
fit_ext=np.array(fit_ext)
fit_bright=np.array(fit_bright)
a1=[]
b1=[]
c1=[]
a2=[]
b2=[]
c2=[]
for j in range(len(fit_bright)):
    
    x=scipy.optimize.leastsq(fitfunc, x0, args=(LinearFit,fit_bright[j],fit_hdew[j],fit_ext[j], fit_flux_3900[j]), Dfun=None, full_output=1, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
    x2=scipy.optimize.leastsq(fitfunc, x0, args=(LinearFit,fit_bright[j],fit_hdew[j],fit_ext[j], fit_flux_8000[j]), Dfun=None, full_output=1, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
    a1.append(x[0][0])
    b1.append(x[0][1])
    c1.append(x[0][2])
    a2.append(x2[0][0])
    b2.append(x2[0][1])
    c2.append(x2[0][2])
