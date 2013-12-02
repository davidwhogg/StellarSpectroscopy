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

f=open("wavelengths",'rb')
wls_ideal=pickle.load(f)
f.close()

f=open("fitted_coefficients_leastsq",'rb')
coeffs=pickle.load(f)
f.close()

f=open("fitted_coefficients_curvefit","rb")
coeffs2=pickle.load(f)
f.close()

f=open("fitted_errors_curvefit",'rb')
errors2=pickle.load(f)
f.close()

f=open("fitted_errors_leastsq",'rb')
errors=pickle.load(f)
f.close()

s0=scipy.interpolate.UnivariateSpline(wls_ideal, coeffs2[0], w=None, bbox=[None, None], k=3, s=None)
s1=scipy.interpolate.UnivariateSpline(wls_ideal, coeffs2[1], w=None, bbox=[None, None], k=3, s=None)
s2=scipy.interpolate.UnivariateSpline(wls_ideal, coeffs2[2], w=None, bbox=[None, None], k=3, s=None)
splines=np.array([s0,s1,s2])

def setplotsize(a,b):
        rcParams['figure.figsize'] = a,b

knownwls=np.array([3440,4000,4400,5500,7000,9000])
knownal_av=np.array([1.59,1.43,1.33,1.00,0.74,0.48]) #a_L/a_V
ag_av=1.161
# a_L/a_V / a_g/a_v = a_L/a_g
theoretical_ag = knownal_av/ag_av

#Plot coefficients
for i in range(3):

    setplotsize(10,5)
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

    #label dibs
    plus2=np.array([0.1,0.3,0.1,0.25,0.1,0.1,0.1])
    dibs=np.array([4430,5449,6284,5780,4727,6613,7224])
    for k in range(len(dibs)):
        #plt.axvline(x=dibs[k], c='b', lw=0.3)
        plt.annotate(str(dibs[k]), xy=(dibs[k],splines[i](dibs[k])), xycoords='data',
                  xytext=(dibs[k],splines[i](dibs[k])+plus2[k]),
                  va="bottom", ha="center",
                  bbox=dict(boxstyle="round", fc="w"),
                  arrowprops=dict(arrowstyle="->"))
    

    #label Balmer lines
    balmer=np.array([6563,4861,4341,4102])
    minus=np.array([0.04,0.04,0.04,0.04]) #space between arrow and curve
    minus2=np.array([0.15,0.15,0.2,0.15]) #space between text and arrow
    
    for m in range(len(balmer)):
        sort=abs(wls_ideal-balmer[m])
        index=nonzero((sort==min(sort)))[0][0]
        print coeffs2[i][index]
        #plt.axvline(x=balmer[m], c='gray', lw=0.3)
        plt.annotate(str(balmer[m]), xy=(balmer[m],coeffs2[i][index]-minus[m]), xycoords='data',
                  xytext=(balmer[m],coeffs2[i][index]-minus2[m]),
                  va="top", ha="center",
                  bbox=dict(boxstyle="round", fc="w"),
                  arrowprops=dict(arrowstyle="->"))
                  
    plt.plot(wls_ideal,errors2[i],'k')
    plt.plot(knownwls,theoretical_ag,'k',linestyle='dashed',marker='o'\
             ,label="Theory from Schlegel/Whittet")
    plt.legend(loc=0)
    plt.ylim(-0.1,1.2)
    plt.xlim(3700,9500)
    plt.axhline(y=0,c='k',lw=0.1)
    plt.savefig("test_curvefit_coefficients"+str(i))
    #plt.savefig("coeffsdr9experr_newbright_lines_expand_restack"+str(i))
    plt.clf()

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

