
#import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import UnivariateSpline
import triangle
import math
#rom pylab import *
#import commands
#import sys
#os.chdir("/Users/admin/Desktop/Maser files")
#directory=commands.getoutput("pwd")
#sys.path.append(directory)
#import sqlcl
import sqldata_notmain

#ns=np.linspace(0.08,0.20,25)
ns=[0.2,0.5]
numbers=[]
for i in range(len(ns)):
    array10=sqldata_notmain.sqldata(ns[i]) #the number is the radius of the 4-Sphere

    #wls=array10[1][1][0] #array of wavelengths for one particular star
    #fluxes=array10[1][0][0] #array of fluxes for one particular star

    #below are lists for graphing from sqlcode
    objs=array10[0][8]
    teff=array10[0][2]
    feh=array10[0][3]
    extinction=array10[0][1]
    ras=array10[0][12]
    decs=array10[0][13]
    magu=array10[0][14]
    magg=array10[0][15]
    magr=array10[0][16]
    magi=array10[0][17]
    magz=array10[0][18]
    ## make one point red
    c=['b']*len(objs)
    c[12]='r'
    ca=np.array(c)

    # make that point larger as well
    s=[5]*len(objs)
    s[12]=100
    sa=np.array(s)

    maggr=np.array(magg)-np.array(magr)
    magug=np.array(magu)-np.array(magg)
    magri=np.array(magr)-np.array(magi)

    data=[0]*(len(maggr))
    for j in range(len(maggr)):
        data[j]=np.array([maggr[j],magri[j],magug[j]])
    print len(maggr)
   
    numbers.append(len(maggr))
    figure=triangle.corner(data, labels=["g-r","r-i","u-g"],plot_ellipse=True)
    name=str(ns[i])+".png"

    ###extinction percentiles
    '''extinction.sort()
    size=len(extinction)
    p1=extinction[int(math.ceil((size*1)/100)-1)]
    p5=extinction[int(math.ceil((size*5)/100)-1)]
    p25=extinctionn
    p50=extinction[int(math.ceil((size*50)/100)-1)]
    p75=extinction[int(math.ceil((size*75)/100)-1)]
    p95=extinction[int(math.ceil((size*95)/100)-1)]
    p99=extinction[int(math.ceil((size*99)/100)-1)]
    '''
   
    figure.savefig(name)
    
    
    ##### SPLINE
    '''
    s = UnivariateSpline(wls, fluxes, k=3, s=0)
    xs=linspace(min(wls),max(wls),len(wls)*10)
    ys=s(xs)
    '''
    #below section graphs the spectrum
    '''
    plt.subplot(212)
    plt.step(xs,ys)
    plt.xlabel("wavelengths (A)")
    plt.ylabel("flux (E-17 ergs/s/cm^2/A)")
    fig = plt.gcf()
    fig.set_size_inches(6,3)
    '''
plt.figure(2)
plt.scatter(ns,np.array(numbers))
plt.xlabel("width of 4-box of i-z")
plt.title("Varying i-z and holding others constant at 0.08")
plt.ylabel("number of stars")
plt.savefig("n--r_for_i-z.png")
