from pylab import *
from scipy import stats
import numpy as np
import sys
import os
import commands
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pickle
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


def setplotsize(a,b):
        rcParams['figure.figsize'] = a,b

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
    plt.plot(wls_ideal,coeffs[i],'k')
    dibs=np.array([4430,5449,6284,5780,5778,4727,5382,5535,6177,6005,6590,6613,7224])
    for k in dibs:
        plt.axvline(x=k, c='b', lw=0.3)
    
    
    plt.axvline(x=6563, c='gray', lw=0.3)
    plt.axvline(x=4861, c='gray', lw=0.3)
    plt.axvline(x=4341, c='gray', lw=0.3)
    plt.axvline(x=4102, c='gray', lw=0.3)
    plt.plot(wls_ideal,errors[i],'k')
    plt.ylim(-0.1,1.2)
    plt.axhline(y=0,c='g',lw=0.1)
    setplotsize(10,5)
    plt.savefig("test2"+str(i))
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

