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
import matplotlib.pyplot as plt
os.chdir("/Users/admin/Desktop/Maser files")
directory=commands.getoutput("pwd")
sys.path.append(directory)

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


if __name__=="__main__":
                
        f=open("sorted","rb")
        array1=pickle.load(f)
        f.close()
        f=open("datanewdr8bb","rb")
        data=pickle.load(f)
        f.close()
        
        #Get values from pickled document
        ext=np.array(array1[1])
        magg=np.array(array1[15])
        bright= 10**(-0.4*(magg-22.5))

        plates=np.array(array1[9])
        mjds=np.array(array1[11])
        fibers=np.array(array1[10])


        hdew=[]

        #verify that data are ordered correctly
        ext2=[]
        for i in range(len(magg)):
            ext2.append(data[0][i][5])
            hdew.append(data[0][i][4])
            
        if ext2[0]!=ext[0]:
            print "Error! Data arrays do not align!"
            sys.exit()


        #Standardize flux array length
        fit_flux= [[] for x in xrange(3120)]
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
            if wls_used[0]< 3900.4 and wls_used[-1]>7998.0:
                for j in range(len(flux_used)):
                        fit_flux[j].append(flux_used[j])
                
                fit_hdew.append(hdew[i])
                fit_ext.append(ext[i])
                fit_bright.append(bright[i])
            else:
                print "missing flux values", wls_used[0], wls_used[-1]
                pass


        x0i=[0.1,0.1,0.1] #initial guess
        x0=array(x0i) #optimize requires array
        fit_flux_array = np.array(fit_flux)
        fit_hdew=np.array(fit_hdew)
        fit_ext=np.array(fit_ext)
        fit_bright=np.array(fit_bright)
        aa=fit_ext
        store_values=[]

        #Fitting
        for i in range(len(wls_used)): #=3120
                x=scipy.optimize.leastsq(fitfunc, x0, args=(LinearFit,fit_bright,fit_hdew,fit_ext, fit_flux[i]), Dfun=None, full_output=1, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
                store_values.append(x[0])


        #Save coefficients to file
        h=open("storedvalues","wb")
        pickle.dump(store_values,h)
        h.close()

        #Group and save fitted coefficient to file
        '''
        bcoeff=[]
        hcoeff=[]
        ecoeff=[]

        for i in range(3120):
                bcoeff.append(store_values[i][0])
                hcoeff.append(store_values[i][1])
                ecoeff.append(store_values[i][2])
        coeffs=np.array([bcoeff,hcoeff,ecoeff])
        l=open("coeffs","wb")
        pickle.dump(coeffs,l)
        l.close()'''

        #Plot
        for i in range(3):
                plt.figure()
                if i ==0:
                        plt.title("Coefficient of brightness")
                elif i==1:
                        plt.title("Coefficient of H-D EW")
                else:
                        plt.title("Coefficient of Extinction")
                plt.xlabel("Wavelengths, A")
                plt.ylabel("Coefficient")
                plt.plot(wls_used,coeffs[i])
                plt.savefig("coeffs"+str(i))
                plt.clf()
                
                
        #x2=scipy.optimize.leastsq(fitfunc, x0, args=(LinearFit,fit_bright,fit_hdew,fit_ext, fit_flux_8000), Dfun=None, full_output=1, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
        #plt.scatter(fit_ext*x2[0][2]+fit_bright*x2[0][0]+fit_hdew*x2[0][1],fit_flux_8000,s=10,linewidths=0 )
        #plt.plot([0,max(fit_flux_8000)],[0,max(fit_flux_8000)],'k')
        #plt.plot([min(aa),max(aa)],[min(aa)*x[0][2],max(aa)*x[0][2]],color='k')
        #plt.plot([0,0],[max(fit_hdew),max(fit_hdew)*x[0][0]],color='r')

        #plt.xlabel(str(round(x2[0][0],3))+"*Brightness+"+str(round(x2[0][1],3))+"*H-D EW+"+str(round(x2[0][2],2))+"*Extinction")
        #plt.ylabel("Flux at 8000A")
        #plt.show()
