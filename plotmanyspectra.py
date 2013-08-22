import numpy as np
import os
from pylab import *
from scipy.interpolate import UnivariateSpline
from scipy import stats
import commands
import pyfits
import sys
os.chdir("/Users/admin/Desktop/Maser files")
directory=commands.getoutput("pwd")
sys.path.append(directory)
import sqlcl
import pickle
import random



def deredshift(wls, fluxes, zm, badpoints):
    lamcorb=wls[0]/(1.0-zm)
    s = UnivariateSpline(wls[0], fluxes[0], k=3, s=0)
    xs=wls[0]    
    ys=s(lamcorb)
    bady=s(badpoints)
    badx=badpoints
    return xs, ys, badx, bady


########take flux, wl data from fits files##########
def getdata():
        
    tabs=[0] #we only get one entry each time, we overwrite each time
    fluxes=[0]
    sn2s=[0]
    errormags=[0]
    wls=[0]
    ivars=[0]
    ### extinction verification, for checking anomalies. Not default
    '''for i in range(len(plate)):########## cf line 132, 175, 193
        if extinction[i]==0.077685:
            plateid=plate[i]
            mjdid=mjd[i]
            fiberid=fiber[i]'''
    #if extinction[i]!=data[0][i][5]:
        #print "mismatched extinction values!"
        #sys.exit()
        
    ### use this by default
    #for q in range(1): #cf 132, 175, 198
    #plateid=plate[i]
    #mjdid=mjd[i]
    #fiberid=fiber[i]

    #plateid=data[0][i][6]
    #mjdid=data[0][i][7]
    #fiberid=data[0][i][8]
    
    #tab = pyfits.open(commands.getoutput("pwd")+'/FITS_files/spec-'+str(plateid).zfill(4)+'-'+str(mjdid)+'-'+str(fiberid).zfill(4)+'.fits')
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
            
    tab.close()
    return wls, fluxes, sn2s, ivars, badpoints, zm



########### Plot spectra ######## (need grouped_data) - can only plot 
def plot(n):
    #if data[0][n][7]!=mjd[n]:
        #print "mjd mismatch!"
        #sys.exit()
    xs, ys, badx, bady = deredshift(wls, fluxes, zm, badpoints)
    # peaks are [hd, hg, hb,  TiO,  Na,  He-I,   K,    H ] [6306 O-I]
    peaks=[4102, 4340, 4861, 5582, 5896, 3890, 3934, 3970] 
    for w in range(0,1): #for each peak location..

        cont1=data[w][ind[n]][0]
        cont_err=data[w][ind[n]][1]
        flux=data[w][ind[n]][2]
        flux_err=data[w][ind[n]][3]
        eqw=data[w][ind[n]][4]
        peak_loc=peaks[w]
        ext=data[w][ind[n]][5]

        flux_domain = [x for x in xs if x<peak_loc+10 and x>peak_loc-10] #neighborhood of 20A

        ys = s(flux_domain)
        
        #smarty=[s(x) for x in xs if x>3800 and x<7000]
        #plt.fill_between(flux_domain, ys, cont1, color='gray', alpha=0.5)
        #plt.axvline(x=peak_loc+10, color='b')
        #plt.axvline(x=peak_loc-10, color='b')
        #plt.axvline(x=peak_loc, color='r')
        #plt.axvline(x=peak_loc+100, color='k')
        #plt.axvline(x=peak_loc-100, color='k')
        #plt.axvline(x=peak_loc-20, color='k')
        #plt.axvline(x=peak_loc+20, color='k')


        plt.plot(np.array([peak_loc-100,peak_loc-20]),np.array([contv]*2), color='k')
        plt.plot(np.array([peak_loc+20,peak_loc+100]),np.array([contv]*2), color='k')

        #plt.plot(np.array([3850,3880]),np.array([cont1]*2), color='k')
        #plt.plot(np.array([3900,3920]),np.array([cont1]*2), color='k')

        #countdown = [x for x in badpoints if x<peak_loc and x>peaks[w]-100]
        #countup = [x for x in badpoints if x>peak_loc and x<peaks[w]+100]
        #print peaks[w], "number of non-zero ivar pixels in LHS: ", len(countdown)
        #print peaks[w], "number of non-zero ivar pixels in RHS: ", len(countup)
    k=0 #[hd, hg, hb,  TiO,  Na,  He-I,   K,    H ]
    s1="Cont: "+str(round(contv,1)) #h-delta
    s2="Err: "+str(round(conterrv,3))
    s3="Flux: "+str(round(fluxv,1))
    s4="Err: "+str(round(fluxerrv,1))
    s5="EW: "+str(round(hdewv,2))

    plateid=data[0][ind[n]][6] #ignore, use plate below
    mjdid=data[0][ind[n]][7]
    fiberid=data[0][ind[n]][8]
    
    s6=str(plate)+"-"+str(mjd)+"-"+str(fiber)+".fits"
    s7="Ext: "+str(round(extinctionv,4))
    #s8="g-i: "+str(round(gi[n],5))
    plt.text(0.10,0.05, s1, fontsize=11, transform = ax.transAxes)
    plt.text(0.25,0.05, s2, fontsize=11, transform = ax.transAxes)
    plt.text(0.4,0.05, s3, fontsize=11, ha='right', transform = ax.transAxes)
    plt.text(0.5,0.05, s4, fontsize=11, ha='right', transform = ax.transAxes)
    plt.text(0.6,0.05, s5, fontsize=11, ha='right', transform = ax.transAxes)
    plt.text(0.7,0.05,s7,fontsize=11,ha='right', transform=ax.transAxes)
    #plt.text(0.8,0.05,s8,fontsize=11,ha='right', transform=ax.transAxes)
    plt.text(0.5,0.8, s6, fontsize=14, ha='center', transform=ax.transAxes)

    plt.step(xs,ys+sigmas[0], 'g', linewidth=0.4, alpha=1)
    plt.step(xs,ys-sigmas[0],  'g', linewidth=0.4, alpha=1)
    plt.step(xs,ys,'b', linewidth=0.5, alpha=1)
    #plot ivar=0 points
    plt.scatter(np.array(badpoints), np.array(s(badpoints)), c='r', marker='o')    
        
    #plt.xlabel("wavelengths (A)")
    #plt.ylabel("flux (E-17 ergs/s/cm^2/A)")

    plt.tight_layout()

    #plt.title(str(plateid)+"-"+str(mjdid)+"-"+str(fiberid)+".fits")
    plt.grid(True)
    plt.xlim(3800,9000)
    plt.ylim(0, 1.98*np.median(ys))
    return
'''
exts=[]
for i in range(5204):
    exts.append(data[0][i][5])

low=[]
high=[]
   
for i in range(5200):
    if data[0][i][5]<stats.scoreatpercentile(exts,20):
        low.append(i)
    elif data[0][i][5]>stats.scoreatpercentile(exts,75):
        high.append(i)


'''
##################################################
#######################################

def sort(hdew, gi, mjd, plate, fiber, extinction, conts, conterrs, fluxs, fluxerrs)
    ind=argsort(hdew)
    sorthdew=zeros(len(gi))
    sortgi=zeros(len(gi))
    sortplate=zeros(len(gi),int)
    sortmjd=zeros(len(gi),int)
    sortfiber=zeros(len(gi),int)
    sortext=zeros(len(gi),float)
    sortconts=zeros(len(gi),float)
    sortconterrs=zeros(len(gi),float)
    sortfluxs=zeros(len(gi),float)
    sortfluxerrs=zeros(len(gi),float)
    for i in range(len(gi)):
        sorthdew[i]=hdew[ind[i]]
        sortgi[i]=gi[ind[i]]
        sortmjd[i]=int(mjd[ind[i]])
        sortplate[i]=int(plate[ind[i]])
        sortfiber[i]=int(fiber[ind[i]])
        sortext[i]=extinction[ind[i]]
        sortconts[i]=conts[ind[i]]
        sortconterrs[i]=conterrs[ind[i]]
        sortfluxs[i]=fluxs[ind[i]]
        sortfluxerrs[i]=fluxerrs[ind[i]]
    plate=sortplate
    mjd=sortmjd
    fiber=sortfiber
    gi=sortgi
    hdew=sorthdew
    ext=sortext
    cont=sortconts
    conterr=sortconterrs
    flux=sortfluxs
    fluxerr=sortfluxerrs
    array=[plate,mjd,fiber,hdew,ext,cont,conterr,flux,fluxerr]
    alldata=[[],[],[],[],[],[],[],[],[]]
    for i in range(9): #9
        indexstart=377+i*754
        zone=[]
        for j in range(len(array)): #len=5
            zone.append(array[j][indexstart:indexstart+754])

        ind2=argsort(zone[4])
        sortplate2=zeros(754,int)
        sortmjd2=zeros(754,int)
        sortfiber2=zeros(754,int)
        sortext2=zeros(754,float)
        sorthdew2=zeros(754,float)
        sortconts2=zeros(754,float)
        sortconterrs2=zeros(754,float)
        sortfluxs2=zeros(754,float)
        sortfluxerrs2=zeros(754,float)
        for k in range(754):
            sortmjd2[k]=int(zone[1][ind2[k]])
            sortplate2[k]=int(zone[0][ind2[k]])
            sortfiber2[k]=int(zone[2][ind2[k]])
            sortext2[k]=zone[4][ind2[k]]
            sorthdew2[k]=zone[3][ind2[k]]
            sortconts2[k]=zone[5][ind2[k]]
            sortconterrs2[k]=zone[6][ind2[k]]
            sortfluxs2[k]=zone[7][ind2[k]]
            sortfluxerrs2[k]=zone[8][ind2[k]]
            
        array2=[sortplate2,sortmjd2,sortfiber2,sortext2,sorthdew2,sortconts2,sortconterrs2,sortfluxs2,sortfluxerrs2]

        fig, axs = plt.subplots(nrows=9, ncols=1, sharex=True)
        u=np.array([0])
        avgsall=[]
        extsall=[]

        for l in range(9):
            
            indstart=37+75*l
            zone2=[]
            for m in range(len(array2)): #len=5
                zone2.append(array2[m][indstart:indstart+75])

            #wlsall=[]
            xsall=[]
            ysall=[]
            for o in range(len(zone2[0])): #75
                plate=zone2[0][o]
                mjd=zone2[1][o]
                fiber=zone2[2][o]
                wls, fluxes, sn2s, sigmas, badpoints, zm = getdata()
                xs, ys, badx, bady = deredshift(wls, fluxes, zm, badpoints)
                
                #wlsall.append(wls)
                a=min(np.where(xs>3900)[0])
                b=max(np.where(xs<7500)[0])
                if b-a == 2839:
                    xsall.append(xs[a:b])
                    ysall.append(ys[a:b])
                else:
                    pass
            extsall.append([zone2[3][0],zone2[3][-1]])
            avgs=[]
            for p in range(len(ysall[0])):
                avgassist=[]
                for q in range(len(ysall)):
                    avgassist.append(ysall[q][p])
                avg=np.mean(avgassist)
                avgs.append(avg)
            avgsall.append(np.array(avgs))
        avgavg=(avgsall[0]+avgsall[1]+avgsall[2]+avgsall[3]+avgsall[4]+avgsall[5]+avgsall[6]+avgsall[7]+avgsall[8])/9

        for m in range(9):
            j=911+m
            ax=plt.subplot(j)        

            #plt.step(xs[a:b],avgsall[m], 'b')
            #plt.step(xs[a:b],avgavg, 'k')
            plt.step(xs[a:b],avgsall[m]/avgavg, 'b')
            plt.xlim(3890,7500)
            plt.ylim(0,2*np.median(avgsall[m]/avgavg))
            plt.tight_layout()
            plt.grid(True)

            plt.title("Ext range: "+str(extsall[m][0])+" to "+str(extsall[m][1]))
            if m in range(0,8):
                plt.setp(ax.get_xticklabels(), visible=False)
            else:
                plt.setp(ax.get_xticklabels(), visible=True)
            if m==8:
                plt.xlabel("Wavelengths, Ang")
            else:
                pass
        plt.suptitle("Ratio of spectra to average spectra of EWs between "+str(i*10+5)+"th to "+str(i*10+15)+"th percentile in H-D")
        #plt.suptitle("Average Spectra for "+str(i*10+5)+str("th")+" to "+str(i*10+15)+str("th")+ " percentile in HD-EW")
        fig.subplots_adjust(left=0.05, right=0.95, wspace = 0.15, hspace=0)
        fig.set_size_inches(18.0,12.0)
        plt.savefig("spectra_avg_ratio_"+str(i)+".png")
        print i
        return
    
if __name__=="__main__":
    f=open("datanewdr8b","rb")
    data=pickle.load(f)
    f.close()

    f=open("sorted","rb")
    array1=pickle.load(f)
    f.close()
    extinction=np.array(array1[1])
    magu=np.array(array1[14])
    magg=np.array(array1[15])
    magr=np.array(array1[16])
    magi=np.array(array1[17])
    magz=np.array(array1[18])
    plate=np.array(array1[9])
    mjd=np.array(array1[11])
    fiber=np.array(array1[10])

    #platesort=[x for (y,x) in sorted(zip(magg-magi,plate))]
    hdew=[]
    conts=[]
    conterrs=[]
    fluxs=[]
    fluxerrs=[]
    for i in range(len(magu)):
        hdew.append(data[0][i][4])
        conts.append(data[0][i][0])
        conterrs.append(data[0][i][1])
        fluxs.append(data[0][i][2])
        fluxerrs.append(data[0][i][3])
    hdew=np.array(hdew)
    conts=np.array(conts)
    conterrs=np.array(conterrs)
    fluxs=np.array(fluxs)
    fluxerrs=np.array(fluxerrs)


    gi=magg-magi

sys.exit()
'''## work with fluxavg, globalwls
for i in range(1):
    for j in range(1):
        
        

        #plt.step(wls[0],fluxavg, lw=0.3)
        #plt.xlim(3800,9000)
        if l in range(0,8):
            plt.setp(ax.get_xticklabels(), visible=False)
        else:
            plt.setp(ax.get_xticklabels(), visible=True)
        if l==8:
            plt.xlabel("Wavelengths, Ang")
        else:
            pass
    plt.suptitle("Average Spectra for "+str(i*10+5)+str("th")+" to "+str(i*10+15)+str("th")+ " percentile in HD-EW")
    fig.subplots_adjust(left=0.05, right=0.95, wspace = 0.15, hspace=0)
    fig.set_size_inches(18.0,12.0)
    
    plt.savefig("spectra_avg_fine_"+str(i)+".png")'''

        
#########################################
###################################################
'''index=[1,10,100,1000,6041,6441,6841,7241,7400]

fig, axs = plt.subplots(nrows=8, ncols=1, sharex=True)

for i in range(1,9):
    n=index[i-1]
    j=810+i
    ax=plt.subplot(j)
    wls, fluxes, sn2s, sigmas, badpoints = getdata(n) #or succs
    plot(n)
    if i in range(1,8):
        plt.setp(ax.get_xticklabels(), visible=False)
    else:
        plt.setp(ax.get_xticklabels(), visible=True)
    
    
    if i==8:
        plt.xlabel("Wavelengths, Ang")
    else:
        pass
fig.subplots_adjust(left=0.05, right=0.95, wspace = 0.15, hspace=0)
plt.show()'''
