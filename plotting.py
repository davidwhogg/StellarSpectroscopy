# -*- coding: utf-8 -*-
from pylab import *
from scipy import stats
import numpy as np
import sys
import os
import commands
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pickle
import scipy.interpolate

os.chdir("/Users/admin/Desktop/Maser files")
directory=commands.getoutput("pwd")
sys.path.append(directory)
def importdata():
        f=open("wavelengths",'rb')
        wls_ideal=pickle.load(f)
        f.close()

        f=open("fitted_coefficients_leastsq_ivar",'rb')
        coeffs=pickle.load(f)
        f.close()

        f=open("fitted_coefficients_curvefit_ivar","rb")
        coeffs2=pickle.load(f)
        f.close()

        f=open("fitted_errors_curvefit_ivar",'rb')
        errors2=pickle.load(f)
        f.close()

        f=open("fitted_errors_leastsq_ivar",'rb')
        errors=pickle.load(f)
        f.close()
        return wls_ideal, coeffs, coeffs2, errors, errors2

def make_spline(i, wls_ideal, coefficients):
        s=scipy.interpolate.UnivariateSpline(wls_ideal, coeffs2[i], w=None, bbox=[None, None], k=3, s=None)
        return s

def setplotsize(a,b):
        rcParams['figure.figsize'] = a,b

# these numbers are typed in from Whittet D.C.B., 1992, \
# Dust in the Galactic Environment, Inst. of Physics Publishing, Bristol, p. 67
knownwls=np.array(  [1000, 1050, 1110, 1180, 1250, 1390, 1490, 1600, 1700, 1800, 1900, 2000, 2100, 2190, 2300, 2400, 2500, 2740, 3440, 4000, 4400, 5500, 7000, 9000, 12500, 16500, 22200, 35000, 48000])
knownal_av=np.array([4.70, 4.21, 3.77, 3.44, 3.15, 2.77, 2.66, 2.59, 2.56, 2.52, 2.61, 2.81, 3.04, 3.15, 2.89, 2.61, 2.37, 2.02, 1.59, 1.43, 1.33, 1.00, 0.74, 0.48, 0.26,  0.15,  0.09,  0.04,  0.02]) #a_L/a_V
ag_av=1.161 #from Schlegel et. al, 1998, ApJ, \

# these numbers are typed in from Cardelli et al., 1989,
# "The relationship between infrared, optical, and ultraviolet extinction." ApJ 345, p. 245-256.
knownwls2=np.array(  [1250,  1500,   1800,  2175, 2500,  3125,  3600, 4400,   5500,  7000,  9000,  12500, 16000, 22000, 34000])
#using R_v=3.1:
knownal_av2=np.array([3.339, 2.667, 2.521, 3.191, 2.317, 1.758, 1.569, 1.337, 1.000, 0.751, 0.479, 0.282, 0.190, 0.114, 0.056])
#Maps of Dust IR Emission for Use in Estimation of Reddening and CMBR Foregrounds
## a_L/a_V / a_g/a_v = a_L/a_g
theoretical_ag = knownal_av/ag_av
theoretical_ag2=knownal_av2/ag_av

wls_ideal, coeffs, coeffs2, errors, errors2 = importdata()


#Plot coefficients. i=0 is brightness, i=1 is H-D EW, i=2 is ext
def plot_coefficients(i, dibs=True, balmer=True):
        s=make_spline(i, wls_ideal, coeffs2)
        setplotsize(10,5)
        median=abs(np.median(coeffs2[i]))
        plt.figure()
        if i ==0:
                plt.title("Coefficient of brightness recentered")
        elif i==1:
                plt.title("Coefficient of H-D EW")
        else:
                plt.title("Coefficient of Extinction")
        plt.xlabel("Wavelengths, A")
        plt.ylabel("Coefficient")
        plt.plot(wls_ideal,coeffs2[i],'k')
        #dibs=np.array([4430,5449,6284,5780,5778,4727,5382,5535,6177,6005,6590,6613,7224])
        if dibs==True:
                #label dibs
                plus2=np.array([0.2, 0.3,  0.2,  0.25, 0.2,  0.2,  0.2])*median
                dibs=np.array([4430, 5449, 6284, 5780, 4727, 6613, 7224])
                for k in range(len(dibs)):
                        #plt.axvline(x=dibs[k], c='b', lw=0.3)
                        plt.annotate(str(dibs[k]), xy=(dibs[k],s(dibs[k])), xycoords='data',
                          xytext=(dibs[k],s(dibs[k])+plus2[k]),
                          va="bottom", ha="center",
                          bbox=dict(boxstyle="round", fc="w"),
                          arrowprops=dict(arrowstyle="->"))
        else:
                pass
        if balmer==True:
                #label Balmer lines
                balmer=np.array([6563,4861,4341,4102])
                minus=np.array([0.05,0.05,0.05,0.05])*median #space between arrow and curve
                minus2=np.array([0.2,0.2,0.3,0.2])*median #space between text and arrow

                for m in range(len(balmer)):
                        sort=abs(wls_ideal-balmer[m])
                        index=nonzero((sort==min(sort)))[0][0]
                        
                        #plt.axvline(x=balmer[m], c='gray', lw=0.3)
                        plt.annotate(str(balmer[m]), xy=(balmer[m],coeffs2[i][index]-minus[m]), xycoords='data',
                                  xytext=(balmer[m],coeffs2[i][index]-minus2[m]),
                                  va="top", ha="center",
                                  bbox=dict(boxstyle="round", fc="w"),
                                  arrowprops=dict(arrowstyle="->"))
                          
        plt.plot(wls_ideal,errors2[i],'k')
        plt.plot(knownwls,theoretical_ag,'k',linestyle='dashed',marker='o'\
             ,label="Theory from Schlegel/Whittet")
        plt.plot(knownwls2,theoretical_ag2,'k',linestyle='dotted',marker='s'\
                 ,label="Theory from Cardelli et al.")
        plt.legend(loc=0)
        
        plt.ylim(min(-0.2*median,min(coeffs[i])),max(2.2*median,max(coeffs[i])))
        plt.xlim(3700,9500)
        plt.axhline(y=0,c='k',lw=0.1)
        plt.savefig("test_curvefit_coefficients_ivar"+str(i))
        #plt.savefig("coeffsdr9experr_newbright_lines_expand_restack"+str(i))
        plt.clf()
        return

#Plot coefficients with error as subplot
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
        fig.savefig("test"+str(i))

        #plt.savefig("coeffsdr9experr_newbright_lines_expand_restack"+str(i))
        fig.clf()'''

