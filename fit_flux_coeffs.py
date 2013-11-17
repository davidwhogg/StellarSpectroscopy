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


def deredshiftone(wls, fluxes, zm, badpoints):
        lamcorb=wls[0]/(1.0-zm)
        s = UnivariateSpline(wls[0], fluxes[0], k=3, s=0)
        xs=wls[0]
        ys=s(lamcorb)
        bady=s(badpoints)
        badx=badpoints
        return xs, ys, badx, bady


def getdataone(plate,mjd,fiber):
        
        tabs=[0] #we only get one entry each time, we overwrite each time
        fluxes=[0]
        sn2s=[0]
        errormags=[0]
        wls=[0]
        ivars=[0]
        tab = pyfits.open(commands.getoutput("pwd")+'/FITS_files/spec-'+str(int(plate)).zfill(4)+'-'+str(int(mjd))+'-'+str(int(fiber)).zfill(4)+'.fits')
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
                if ivars[0][v]==0:
                    #ivars[0][v]=np.inf
                    badpoints.append(wls[0][v])
        del tab[0].data #free up memoery
        tab.close()
        return wls, fluxes, sn2s, ivars, badpoints, zm

# LinearModel no longer used in favor of ExpModel
#def LinearModel(p,xsIn,ysIn,zsIn):#fit both
#    return p[0]*xsIn+p[1]*ysIn + p[2]*zsIn


objfunc=lambda p,f,fit_bright,fit_hdew,fit_ext,fit_flux, fit_ivar: np.sqrt(fit_ivar)*fabs(f(p,fit_bright,fit_hdew,fit_ext)-fit_flux)

def ExpModel(p, xsIn, ysIn, zsIn):
        return (p[0]*xsIn+p[1]*ysIn) * np.exp(-p[2]*zsIn)

if __name__=="__main__":

        #sorted2 is the output of the SDSS query, sorted into lists
        f=open("sorted3","rb")
        array1=np.array(pickle.load(f))
        f.close()

        #datanewdr9 is the output file of eqw_newer with calculated EWs 
        f=open("datanewdr9","rb")
        data=np.array(pickle.load(f))
        f.close()
        
        #Get values from pickle files
        extg=np.array(array1[1])
        magg=np.array(array1[15])
        bright= 10**(-0.4*(magg-22.5-extg))

        plates=np.array(array1[9])
        mjds=np.array(array1[11])
        fibers=np.array(array1[10])
        print len(plates)

        dataline=data[0].ravel()
        datasort=np.reshape(dataline,(10,-1),order='F')
        
        ext2=datasort[5]
        hdew=datasort[4]
        # datasort[6], [7], [8] are plate/mjd/fiber, if needed for assert
        
        #verify that data are ordered correctly
        assert ext2.all()==extg.all()


        #Standardize flux array length
        fit_flux= [[] for x in xrange(3716)]
        fit_hdew=[]
        fit_ext=[]
        fit_bright=[]
        fit_ivar=[[] for x in xrange(3716)]
        for i in range(len(plates)):
            plate=plates[i]
            mjd=mjds[i]
            fiber=fibers[i]
            wls, fluxes, sn2s, ivars, badpoints, zm = getdataone(plate,mjd,fiber)
            xs, ys, badx, bady = deredshiftone(wls, fluxes, zm, badpoints)
            b=(wls[0]>3900)
            c=(wls[0]<9175)
            wls_used = np.trim_zeros(c*b*wls[0])
            flux_used = np.trim_zeros(c*b*ys)
            ivars_tran = np.trim_zeros(c*b*(ivars[0]+1))
            ivars_used=ivars_tran-1
            if wls_used[0]< 3900.4 and wls_used[-1]>9173.0:
                #for j in [0, -1]:
                for j in range(len(flux_used)):
                        fit_flux[j].append(flux_used[j])
                        fit_ivar[j].append(ivars_used[j])
                fit_hdew.append(hdew[i])
                fit_ext.append(extg[i])
                fit_bright.append(bright[i])
            else:
                print "missing flux values", wls_used[0], wls_used[-1]
                pass

        sys.exit()
        x0i=[0.1,0.1,0.1] #initial guess
        x0=array(x0i) #optimize requires array
        fit_flux_array = np.array(fit_flux)
        fit_hdew=np.array(fit_hdew)
        fit_ext=np.array(fit_ext)
        ##############Fitbright updated
        fit_bright=np.array(fit_bright) #########recentered
        fit_ivar=np.array(fit_ivar)
        store_values=[]
        store_values2=[]
        store_values3=[]
        sigs=[]
        #Fitting
        #for i in [0,-1]: 
        for i in range(len(wls_used)): #=3716
                #x=scipy.optimize.leastsq(objfunc, x0, args=(LinearModel,fit_bright,fit_hdew,fit_ext, np.array(fit_flux[i])), Dfun=None, full_output=1, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
                #store_values.append(x[0])
                x2=scipy.optimize.leastsq(objfunc, x0, args=(ExpModel,fit_bright,fit_hdew,fit_ext, np.array(fit_flux[i]), np.array(fit_ivar[i])), Dfun=None, full_output=1, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
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
        h=open("storedvalues_exp_dr9_errs_rebright_expand","wb")
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

        
        
        #Plot residuals
        '''
        inds=[110,622,1493,2218]
        for r in range(4):
                a=inds[r]
                plt.figure()
                plt.errorbar(fit_flux[a],(fit_ext*store_values2[a][2]+fit_bright*store_values2[a][0]+fit_hdew*store_values2[a][1]),sigs[a], xerr=None, ls='none',alpha=0.3)
                plt.scatter(fit_flux[a],fit_ext*store_values2[a][2]+fit_bright*store_values2[a][0]+fit_hdew*store_values2[a][1],c='k',s=10,linewidths=0 )
                #plt.axhline(y=1)
                plt.plot([0,max(fit_flux[a])],[0,max(fit_flux[a])],'k')
                plt.ylabel("Predicted flux")
                #plt.ylabel(str(round(store_values2[a][0],3))+"*Brightness+"+str(round(store_values2[a][1],3))+"*H-D EW+"+str(round(store_values2[a][2],2))+"*Extinction")
                plt.xlabel("Measured Flux")
                plt.title("Measured vs Calculated flux at "+str(wls_used[a])+"A")
                
                plt.savefig(str(wls_used[a])+"residuals_exp.png")
                plt.clf()'''
        #Plot coefficients
        for i in range(3):
                plt.figure()
                if i ==0:
                        plt.title("Coefficient of brightness recentered")
                elif i==1:
                        plt.title("Coefficient of H-D EW")
                else:
                        plt.title("Coefficient of Extinction")
                plt.xlabel("Wavelengths, A")
                plt.ylabel("Coefficient")
                plt.plot(wls_used,coeffs[i])
                dibs=np.array([4430,5449,6284,5780,5778,4727,5382,5535,6177,6005,6590,6613,7224])
                for k in dibs:
                        plt.axvline(x=k, c='b', lw=0.3)
                
                
                plt.axvline(x=6563, c='gray', lw=0.3)
                plt.axvline(x=4861, c='gray', lw=0.3)
                plt.axvline(x=4341, c='gray', lw=0.3)
                plt.axvline(x=4102, c='gray', lw=0.3)
                
                plt.savefig("coeffsdr9experr_newbright_lines_expand"+str(i))
                plt.clf()
                
        '''#x2=scipy.optimize.leastsq(objfunc, x0, args=(LinearFit,fit_bright,fit_hdew,fit_ext, fit_flux_8000), Dfun=None, full_output=1, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
        #plt.scatter(fit_ext*x2[0][2]+fit_bright*x2[0][0]+fit_hdew*x2[0][1],fit_flux_8000,s=10,linewidths=0 )
        #plt.plot([0,max(fit_flux_8000)],[0,max(fit_flux_8000)],'k')
        #plt.plot([min(aa),max(aa)],[min(aa)*x[0][2],max(aa)*x[0][2]],color='k')
        #plt.plot([0,0],[max(fit_hdew),max(fit_hdew)*x[0][0]],color='r')

        #plt.xlabel(str(round(x2[0][0],3))+"*Brightness+"+str(round(x2[0][1],3))+"*H-D EW+"+str(round(x2[0][2],2))+"*Extinction")
        #plt.ylabel("Flux at 8000A")
        #plt.show()'''
