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
        #if type(tabs[0][2].data.field(63)[0])==float32: #distinguish SDSS,BOSS
        zm= tabs[0][2].data.field("Z") #redshift
        #zm= tabs[0][2].data.field("Z_NOQSO")

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
                    #ivars[0][v]=1/sqrt(ivars[0][v]) convert to Sigma
                        pass
                elif ivars[0][v]==0:
                    #ivars[0][v]=np.inf
                    badpoints.append(wls[0][v])
        del tab[0].data #free up memoery
        tab.close()
        return wls, fluxes, sn2s, ivars, badpoints, zm

def LinearFit(p,xsIn,ysIn,zsIn):#fit both
    return p[0]*xsIn+p[1]*ysIn + p[2]*zsIn
fitfunc=lambda p,f,fit_bright,fit_hdew,fit_ext,fit_flux: np.sqrt(fit_ivar)*fabs(f(p,fit_bright,fit_hdew,fit_ext)-fit_flux)
def ExpFit(p, xsIn, ysIn, zsIn):
        return (p[0]*xsIn+p[1]*ysIn) * np.exp(p[2]*zsIn)
#fitfunc2=lambda p,f,fit_bright,fit_hdew,fit_ext,fit_flux,fit_ivar: fit_ivar*(f(p,fit_bright,fit_hdew,fit_ext)-fit_flux)**2

if __name__=="__main__":
                
        f=open("sorted2","rb")
        array1=pickle.load(f)
        f.close()
        f=open("datanewdr9","rb")
        data=pickle.load(f)
        f.close()
        
        #Get values from pickled document
        ext=np.array(array1[1])
        magg=np.array(array1[15])
        bright= 10**(-0.4*(magg-22.5))

        plates=np.array(array1[9])
        mjds=np.array(array1[11])
        fibers=np.array(array1[10])
        print len(plates)

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
        fit_ivar=[[] for x in xrange(3120)]
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
            ivars_tran = np.trim_zeros(c*b*(ivars[0]+1))
            ivars_used=ivars_tran-1
            if wls_used[0]< 3900.4 and wls_used[-1]>7998.0:
                #for j in [0, -1]:
                for j in range(len(flux_used)):
                        fit_flux[j].append(flux_used[j])
                        fit_ivar[j].append(ivars_used[j])
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
        fit_ivar=np.array(fit_ivar)
        store_values=[]
        store_values2=[]
        store_values3=[]
        sigs=[]
        #Fitting
        #for i in [0,-1]: 
        for i in range(len(wls_used)): #=3120
                #x=scipy.optimize.leastsq(fitfunc, x0, args=(LinearFit,fit_bright,fit_hdew,fit_ext, np.array(fit_flux[i])), Dfun=None, full_output=1, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
                #store_values.append(x[0])
                x2=scipy.optimize.leastsq(fitfunc, x0, args=(ExpFit,fit_bright,fit_hdew,fit_ext, np.array(fit_flux[i])), Dfun=None, full_output=1, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
                store_values2.append(x2[0])
                fit_sig=np.zeros(len(fit_ivar[i]))

                #use Np.LinAlg.LstSq
                '''A=np.vstack([fit_bright, fit_hdew, fit_ext]).T
                b=fit_flux_array[i]
                a=np.linalg.lstsq(A,b)[0]
                store_values3.append(a)'''
                
                for k in range(len(fit_ivar[i])):
                        if fit_ivar[i][k]==0:
                                fit_sig[k]=np.inf
                        else:
                                fit_sig[k]=1.0/fit_ivar[i][k]
                sigs.append(fit_sig)

        sigs=np.array(sigs) 
        
        #Save coefficients to file
        h=open("storedvalues_exp_dr9","wb")
        pickle.dump(store_values2,h)
        h.close()
        #dd=open("storedvaluescheck","wb")
        #pickle.dump(store_values,dd)

        #Group and save fitted coefficient to file
        bcoeff=[]
        hcoeff=[]
        ecoeff=[]

        for i in range(len(store_values2)):
                bcoeff.append(store_values2[i][0])
                hcoeff.append(store_values2[i][1])
                ecoeff.append(store_values2[i][2])
        coeffs=np.array([bcoeff,hcoeff,ecoeff])
        
        #l=open("coeffsAnotherFit","rb")
        #coeffs2=pickle.load(l)
        #
        #pickle.dump(coeffs,l)
        #l.close()

        
        '''
        #Plot residuals
        titles=[3900,8000]
        inds=[0,-1]
        for r in range(2):
                a=inds[r]
                plt.figure()
                #plt.scatter(fit_ext*ecoeff[a]+fit_bright*bcoeff[a]+fit_hdew*hcoeff[a],fit_flux[a],s=10,linewidths=0 )
                plt.errorbar(fit_ext*coeffs2[2][a]+fit_bright*coeffs2[0][a]+fit_hdew*coeffs2[1][a],fit_flux[a], sigs[a], xerr=None, ls='none')
                plt.plot([0,max(fit_flux[a])],[0,max(fit_flux[a])],'k')
                plt.xlabel(str(round(coeffs2[0][a],3))+"*Brightness+"+str(round(coeffs2[1][a],3))+"*H-D EW+"+str(round(coeffs2[2][a],2))+"*Extinction")
                plt.ylabel("Measured Flux")
                plt.title("Measured vs Calculated flux at "+str(titles[r])+"A")
                plt.show()
                plt.savefig(str(titles[r])+"residuals")
                plt.clf()'''
        #Plot coefficients
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
                plt.savefig("coeffsdr9exp"+str(i))
                plt.clf()
                
                '''
        #x2=scipy.optimize.leastsq(fitfunc, x0, args=(LinearFit,fit_bright,fit_hdew,fit_ext, fit_flux_8000), Dfun=None, full_output=1, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
        #plt.scatter(fit_ext*x2[0][2]+fit_bright*x2[0][0]+fit_hdew*x2[0][1],fit_flux_8000,s=10,linewidths=0 )
        #plt.plot([0,max(fit_flux_8000)],[0,max(fit_flux_8000)],'k')
        #plt.plot([min(aa),max(aa)],[min(aa)*x[0][2],max(aa)*x[0][2]],color='k')
        #plt.plot([0,0],[max(fit_hdew),max(fit_hdew)*x[0][0]],color='r')

        #plt.xlabel(str(round(x2[0][0],3))+"*Brightness+"+str(round(x2[0][1],3))+"*H-D EW+"+str(round(x2[0][2],2))+"*Extinction")
        #plt.ylabel("Flux at 8000A")
        #plt.show()'''
