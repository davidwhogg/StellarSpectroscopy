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

f=open("coefficients",'rb')
coeffs=pickle.load(f)
f.close()

f=open("coefficient_errors",'rb')
errors=pickle.load(f)
f.close()

s0=scipy.interpolate.UnivariateSpline(wls_ideal, coeffs[0], w=None, bbox=[None, None], k=3, s=None)
s1=scipy.interpolate.UnivariateSpline(wls_ideal, coeffs[1], w=None, bbox=[None, None], k=3, s=None)
s2=scipy.interpolate.UnivariateSpline(wls_ideal, coeffs[2], w=None, bbox=[None, None], k=3, s=None)
splines=np.array([s0,s1,s2])

def setplotsize(a,b):
        rcParams['figure.figsize'] = a,b

#Plot coefficients
for i in range(2,3):
    plt.figure()
    if i ==0:
        plt.title("Coefficient of brightness recentered")
    elif i==1:
        plt.title("Coefficient of H-D EW")
    else:
        plt.title("Coefficient of Extinction")
    plt.xlabel("Wavelengths, A")
    plt.ylabel("Coefficient")
    plt.plot(wls_ideal,coeffs[i],'k')
    #dibs=np.array([4430,5449,6284,5780,5778,4727,5382,5535,6177,6005,6590,6613,7224])
    dibs=np.array([4430,5449,6284,5780,4727,6613,7224])
    for k in range(len(dibs)):
        plt.axvline(x=dibs[k], c='b', lw=0.3)
        plt.annotate(str(dibs[k]), xy=(dibs[k],splines[i](dibs[k])), xycoords='data',
                  xytext=(dibs[k],splines[i](dibs[k])+0.1),
                  va="bottom", ha="center",
                  bbox=dict(boxstyle="round", fc="w"),
                  arrowprops=dict(arrowstyle="->"))
    
    balmer=np.array([6563,4861,4341,4102])
    minus=np.array([0.1,0.1,0.2,0.14])
    for m in range(len(balmer)):
        sort=abs(wls_ideal-balmer[m])
        index=nonzero((sort==min(sort)))[0][0]
        print coeffs[i][index]
        plt.axvline(x=balmer[m], c='gray', lw=0.3)
        plt.annotate(str(balmer[m]), xy=(balmer[m],coeffs[i][index]-minus[m]), xycoords='data',
                  va="top", ha="center",
                  bbox=dict(boxstyle="round", fc="w"))

    plt.plot(wls_ideal,errors[i],'k')

    plt.ylim(-0.1,1.2)
    plt.axhline(y=0,c='g',lw=0.1)
    setplotsize(10,5)
    plt.show()

    #plt.savefig("test22"+str(i))
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

