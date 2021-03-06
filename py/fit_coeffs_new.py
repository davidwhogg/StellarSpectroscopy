import scipy.optimize
from scipy import optimize
from scipy.optimize import curve_fit
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
import matplotlib.gridspec as gridspec
os.chdir("/Users/admin/Desktop/Maser files")
directory=commands.getoutput("pwd")
sys.path.append(directory)


def setplotsize(a,b):
        rcParams['figure.figsize'] = a,b
        
def deredshiftone(wls, fluxes, zm, badpoints):
        lamcorb=wls/(1.0-zm)
        s = UnivariateSpline(wls, fluxes, k=3, s=0)
        xs=wls
        ys=s(lamcorb)
        bady=s(badpoints)
        badx=badpoints
        return xs, ys, badx, bady

def replace_badpix(ivars,tiny):#relies on getdataone() output
        #(can change)
        newivars=1.0*ivars
        for i in range(2,len(ivars)-2):
                if ivars[i]<tiny:
                        newivars[i-2:i+3]=0.0
                
                else:
                        pass
        return newivars


def getdataone(plate,mjd,fiber):
        
        tabs=[0] #we only get one entry each time, we overwrite each time
        tab = pyfits.open(commands.getoutput("pwd")+'/FITS_files/spec-'+str(int(plate)).zfill(4)+'-'+str(int(mjd))+'-'+str(int(fiber)).zfill(4)+'.fits')
        #print "tab success", i
        tabs[0]=tab
        j=0 #-2700 ###########
        #if type(tabs[0][2].data.field(63)[0])==float32: #distinguish SDSS,BOSS
        zm= tabs[0][2].data.field("Z") #redshift
        #zm= tabs[0][2].data.field("Z_NOQSO")

        fluxes=tabs[0][1].data.field(0) #flux
        
        loglam=tabs[0][1].data.field(1)
        ivars=tabs[0][1].data.field(2)

        loglam=np.array(loglam)
        wls=10**loglam
        sn2s=tabs[j][2].data.field(6)[0]+tabs[j][2].data.field(7)[0] #one of these entries is 0 always

        badpoints=[]
        for v in range(len(ivars)):
                if ivars[v]==0:
                    #ivars[0][v]=np.inf
                    badpoints.append(wls[v])
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

    #wls_ideal is now my ideal rest wl grid
    wls_ideal, fluxes, sn2s, ivars, badpoints, zm = getdataone(plates[0],mjds[0],fibers[0])
    xs, ys, badx, bady = deredshiftone(wls_ideal, fluxes, zm, badpoints)
    b=(wls_ideal>3900)
    c=(wls_ideal<9175)
    wls_ideal = np.trim_zeros(c*b*wls_ideal)

    
    #Standardize flux array length
    flux_array = np.zeros((len(plates),len(wls_ideal)))
    ivar_array = np.zeros((len(plates),len(wls_ideal)))
    fit_hdew=np.zeros(len(plates))
    fit_ext=np.zeros(len(plates))
    fit_bright=np.zeros(len(plates))
    
    
    for i in range(len(plates)):
        plate=plates[i]
        mjd=mjds[i]
        fiber=fibers[i]
        wls, fluxes, sn2s, ivars, badpoints, zm = getdataone(plate,mjd,fiber)
        ivars=replace_badpix(ivars,0.000001)
        xs, ys, badx, bady = deredshiftone(wls, fluxes, zm, badpoints)
        b=(wls>3900)
        c=(wls<9175)
        wls_used = np.trim_zeros(c*b*wls)
        flux_used = np.trim_zeros(c*b*ys)
        ivars_tran = np.trim_zeros(c*b*(ivars+1))
        ivars_used=ivars_tran-1
        differences = abs(wls_ideal - wls_used[0]) #wls_used[0] is also min(wls_used)
        differences2 = abs(wls_ideal - wls_used[-1])
        minarg = argmin(differences) #gives index of minimum relative to ideal wl. will be 0 if ideal
        maxarg = argmin(differences2)
        
        
        ###Pad zeros to end in case wls incomplete
        flux_used=np.append(flux_used,[0]*(len(wls_ideal)-maxarg-1))
        flux_used=np.insert(flux_used,[0]*minarg,0)
        ivars_used=np.append(ivars_used,[0]*(len(wls_ideal)-maxarg-1))
        ivars_used=np.insert(ivars_used,[0]*minarg,0)
        assert len(flux_used)==len(wls_ideal)
        
        #if wls_used[0]< 3900.4 and wls_used[-1]>9173.0: #equivalently, if maxarg-minarg==3715
        flux_array[i]=flux_used
        ivar_array[i]=ivars_used
        fit_hdew[i]=hdew[i]
        fit_ext[i]=extg[i]
        fit_bright[i]=bright[i]
        
    x0i=[0.1,0.1,0.1] #initial guess
    x0=array(x0i) #optimize requires array

    ##restack
    flux_unravel=flux_array.ravel()
    flux_sort = np.reshape(flux_unravel,(len(wls_ideal),-1),'F')
    ivar_unravel=ivar_array.ravel()
    ivar_sort= np.reshape(ivar_unravel,(len(wls_ideal),-1),'F')
    
    store_values=np.zeros((len(wls_ideal),3))
    b_errors=np.zeros(len(wls_ideal))
    hdew_errors=np.zeros(len(wls_ideal))
    ext_errors=np.zeros(len(wls_ideal))

    store_values2=np.zeros((len(wls_ideal),3))
    b_errors2=np.zeros(len(wls_ideal))
    hdew_errors2=np.zeros(len(wls_ideal))
    ext_errors2=np.zeros(len(wls_ideal))
    
    store_popt=np.zeros((len(wls_ideal),3))
    #Fitting
    for i in range(len(wls_ideal)): #=3716
        x2=scipy.optimize.leastsq(objfunc, x0, args=(ExpModel,fit_bright,fit_hdew,fit_ext, np.array(flux_sort[i]), np.array(ivar_sort[i])), Dfun=None, full_output=1, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
        popt, pcov =scipy.optimize.curve_fit(ExpModel,(fit_bright,fit_hdew,fit_ext),flux_sort[i],sigma=1/sqrt(ivar_sort[i]))

        store_popt[i]=popt
        b_errors2[i]=sqrt(pcov[0][0])
        hdew_errors2[i]=sqrt(pcov[1][1])
        ext_errors2[i]=sqrt(pcov[2][2])

        store_values[i]=x2[0]
        b_errors[i]=sqrt(x2[1][0][0])
        hdew_errors[i]=sqrt(x2[1][1][1])
        ext_errors[i]=sqrt(x2[1][2][2])

    errors=np.array([b_errors,hdew_errors,ext_errors])
    errors2=np.array([b_errors2,hdew_errors2,ext_errors2])
    #Save errors to file
    h=open("fitted_errors_leastsq_ivar","wb")
    pickle.dump(errors,h)
    h.close()
    h=open("fitted_errors_curvefit_ivar","wb")
    pickle.dump(errors2,h)
    h.close()
    
    #group coefficients by type
    coeffs_unravel=store_values.ravel()
    coeffs=np.reshape(coeffs_unravel,(3,-1),'F')

    coeffs_unravel2=store_popt.ravel()
    coeffs2=np.reshape(coeffs_unravel2,(3,-1),'F')
    
    h=open("fitted_coefficients_leastsq_ivar","wb")
    pickle.dump(coeffs,h)
    h.close()
    h=open("fitted_coefficients_curvefit_ivar","wb")
    pickle.dump(coeffs2,h)
    h.close()
    
    
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
    '''for i in range(3):
        plt.figure()
        if i ==0:
            plt.title("Coefficient of brightness recentered")
        elif i==1:
            plt.title("Coefficient of H-D EW")
        else:
            plt.title("Coefficient of Extinction")
        plt.xlabel("Wavelengths, A")
        plt.ylabel("Coefficient")
        plt.plot(wls_ideal,coeffs[i])
        dibs=np.array([4430,5449,6284,5780,5778,4727,5382,5535,6177,6005,6590,6613,7224])
        for k in dibs:
            plt.axvline(x=k, c='b', lw=0.3)
        
        
        plt.axvline(x=6563, c='gray', lw=0.3)
        plt.axvline(x=4861, c='gray', lw=0.3)
        plt.axvline(x=4341, c='gray', lw=0.3)
        plt.axvline(x=4102, c='gray', lw=0.3)
        plt.plot(wls_ideal,errors[i],'k')
        setplotsize(10,5)
        plt.savefig("test"+str(i))
        #plt.savefig("coeffsdr9experr_newbright_lines_expand_restack"+str(i))
        plt.clf()'''



    
    
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
    #Plot coefficients with error in separate subplot
    '''for i in range(2,3):
        fig=plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        if i ==0:
            ax1.set_title("Coefficient of brightness recentered")
        elif i==1:
            ax1.set_title("Coefficient of H-D EW")
        else:
            ax1.set_title("Coefficient of Extinction")
        ax1.plot(wls_ideal,coeffs[i])
        ax1.set_ylabel("Coefficient")
        ax2.set_ylabel("Error")
        ax2.plot(wls_ideal,errors[i])
        ax2.set_xlabel("Wavelengths, A")
        #plt.plot(wls_ideal,coeffs[i])
        dibs=np.array([4430,5449,6284,5780,5778,4727,5382,5535,6177,6005,6590,6613,7224])
        for k in dibs:
            ax1.axvline(x=k, c='b', lw=0.3)
        ax1.grid(b=True)
        ax2.grid(b=True)
        ax1.axvline(x=6563, c='gray', lw=0.3)
        ax1.axvline(x=4861, c='gray', lw=0.3)
        ax1.axvline(x=4341, c='gray', lw=0.3)
        ax1.axvline(x=4102, c='gray', lw=0.3)
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.tight_layout()
        setplotsize(10,2)
        fig.savefig("test2"+str(i))
        
        #plt.savefig("coeffsdr9experr_newbright_lines_expand_restack"+str(i))
        fig.clf()
        #residuals below

        #x2=scipy.optimize.leastsq(objfunc, x0, args=(LinearFit,fit_bright,fit_hdew,fit_ext, fit_flux_8000), Dfun=None, full_output=1, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
        #plt.scatter(fit_ext*x2[0][2]+fit_bright*x2[0][0]+fit_hdew*x2[0][1],fit_flux_8000,s=10,linewidths=0 )
        #plt.plot([0,max(fit_flux_8000)],[0,max(fit_flux_8000)],'k')
        #plt.plot([min(aa),max(aa)],[min(aa)*x[0][2],max(aa)*x[0][2]],color='k')
        #plt.plot([0,0],[max(fit_hdew),max(fit_hdew)*x[0][0]],color='r')

        #plt.xlabel(str(round(x2[0][0],3))+"*Brightness+"+str(round(x2[0][1],3))+"*H-D EW+"+str(round(x2[0][2],2))+"*Extinction")
        #plt.ylabel("Flux at 8000A")
        #plt.show()'''
