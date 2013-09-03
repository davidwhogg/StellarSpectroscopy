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

class datas:
    def __init__(self, hdew, gi, mjds, plates, fibers, extinction, conts, conterrs, fluxs, fluxerrs):
        self.hdew=hdew
        self.gi=gi
        self.mjds=mjds
        self.plates=plates
        self.fibers=fibers
        self.extinction=extinction
        self.conts=conts
        self.conterrs=conterrs
        self.fluxs=fluxs
        self.fluxerrs=fluxerrs
    def sort(self):
        ind=argsort(self.hdew)
        sorthdew=zeros(len(self.gi))
        sortgi=zeros(len(self.gi))
        sortplate=zeros(len(self.gi),int)
        sortmjd=zeros(len(self.gi),int)
        sortfiber=zeros(len(self.gi),int)
        sortext=zeros(len(self.gi),float)
        sortconts=zeros(len(self.gi),float)
        sortconterrs=zeros(len(self.gi),float)
        sortfluxs=zeros(len(self.gi),float)
        sortfluxerrs=zeros(len(self.gi),float)
        for i in range(len(self.gi)):
            sorthdew[i]=self.hdew[ind[i]]
            sortgi[i]=self.gi[ind[i]]
            sortmjd[i]=int(self.mjds[ind[i]])
            sortplate[i]=int(self.plates[ind[i]])
            sortfiber[i]=int(self.fibers[ind[i]])
            sortext[i]=self.extinction[ind[i]]
            sortconts[i]=self.conts[ind[i]]
            sortconterrs[i]=self.conterrs[ind[i]]
            sortfluxs[i]=self.fluxs[ind[i]]
            sortfluxerrs[i]=self.fluxerrs[ind[i]]
        plates=sortplate
        mjds=sortmjd
        fibers=sortfiber
        gi=sortgi
        hdew=sorthdew
        ext=sortext
        cont=sortconts
        conterr=sortconterrs
        flux=sortfluxs
        fluxerr=sortfluxerrs
        array=[plates,mjds,fibers,hdew,ext,cont,conterr,flux,fluxerr]
        alldata=[[1],[1],[1],[1],[1],[1],[1],[1],[1]]
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
            alldata[i]=array2 #alldata collates each sorted quantile.
        return alldata, array
    def getdata(self,plate,mjd,fiber):
            
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
    def deredshift(self, wls, fluxes, zm, badpoints):
        lamcorb=wls[0]/(1.0-zm)
        s = UnivariateSpline(wls[0], fluxes[0], k=3, s=0)
        xs=wls[0]
        ys=s(lamcorb)
        bady=s(badpoints)
        badx=badpoints
        return xs, ys, badx, bady
        

### Loading Data
f=open("datanewdr8bb","rb")
data=pickle.load(f)
f.close()

f=open("sorted","rb")
array1=pickle.load(f)
f.close()
extinction=np.array(array1[1])
magg=np.array(array1[15])
magr=np.array(array1[16])
magi=np.array(array1[17])
plates=np.array(array1[9])
mjds=np.array(array1[11])
fibers=np.array(array1[10])
gi=magg-magi
hdew=[]
conts=[]
conterrs=[]
fluxs=[]
fluxerrs=[]
ext2=[] #second set of extinctions, just to prove that the two opened files align
for i in range(len(magg)):
    ext2.append(data[0][i][5])
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



### Begin sorting: collect relevant plate-fiber-mjd and ext information
a=datas(hdew, gi, mjds, plates, fibers, extinction, conts, conterrs, fluxs, fluxerrs)
alldata, array = a.sort()
#alldata [i] has each EW quantile's data
superset=[]
for k in range(9): #within alldata
    intermediate=[]
    for l in range(9): #split quantiles
        indstart=37+l*75
        zone2=[]#this will contain the relevant data for each Extinction quantile \
                 #to be later summoned for making averages of the 75 spectra
        for m in range(5): #i.e. only appends plate/fiber/mjd/ext/HDEW
            zone2.append(alldata[k][m][indstart:indstart+75])
        intermediate.append(zone2)
    superset.append(intermediate)

### Using sorted plate-fiber-mjd, get spectra for averaging
averages_supercollection=[]
avg_avg_supercollection=[]
for n in range(9):
    averages_collection=[]
    avg_fluxavg=np.zeros(3120)
    for o in range(9):
        wlscomp=[]
        fluxcomp=[]
        zmcomp=[]
        badcomp=[]
        d=0
        wls_sum=np.zeros(3120)
        flux_sum=np.zeros(3120)
        for p in range(len(superset[0][0][0])): #should be 75
            plate=superset[n][o][0][p]
            mjd=superset[n][o][1][p]
            fiber=superset[n][o][2][p]
            wls, fluxes, sn2s, ivars, badpoints, zm = a.getdata(plate,mjd,fiber)
            b=(wls[0]>3900)
            c=(wls[0]<8000)
            wls_used = np.trim_zeros(c*b*wls[0])
            xs, ys, badx, bady = a.deredshift(wls, fluxes, zm, badpoints)
           
            flux_used = np.trim_zeros(c*b*ys)
            if len(wls_used)==3120:
                wls_sum += wls_used
                flux_sum += flux_used
                d+=1
                zmcomp.append(zm)
                badcomp.append(badpoints)
            else:
                pass
        print d
        fluxes_average=flux_sum/float(d)
        averages_collection.append(fluxes_average)
        avg_fluxavg+=fluxes_average
    avg_avg = avg_fluxavg/9.0
    avg_avg_supercollection.append(avg_avg)
    averages_supercollection.append(averages_collection)

### The above data was pickled and stored for convenience
'''
import pickle
f=open("avgavg","rb")
g=open("averages","rb")
h=open("superset","rb")
e=open("wls","rb")
wls_used=pickle.load(e)
superset=pickle.load(h)
averages_supercollection=pickle.load(g)
avg_avg_supercollection=pickle.load(f)
f.close()
g.close()
h.close()
e.close()'''

### Plot ratios
for l in range(9):
    fig, axs = plt.subplots(nrows=9, ncols=1, sharex=True)

    for m in range(9):
        j=911+m
        ax=plt.subplot(j)


        plt.step(wls_used, averages_supercollection[l][m], 'r')
        plt.step(wls_used, averages_supercollection[l][m]/avg_avg_supercollection[l], 'b')
        plt.xlim(3900,8000)
        plt.ylim(0.9*np.median(averages_supercollection[l][m]/avg_avg_supercollection[l]),1.1*np.median(averages_supercollection[l][m]/avg_avg_supercollection[l]))
        plt.tight_layout()
        plt.grid(True)

        plt.title("Ext range: "+str(superset[l][m][3][0])+" to "+str(superset[l][m][3][-1]))
        if m in range(0,8):
            plt.setp(ax.get_xticklabels(), visible=False)
        else:
            plt.setp(ax.get_xticklabels(), visible=True)
        if m==8:
            plt.xlabel("Wavelengths, Ang")
        else:
            pass
    plt.suptitle("Ratio of fluxes to average fluxes ranked by extinction, with EWs between the "+str(l*10+5)+"th to "+str(l*10+15)+"th percentile")
    fig.subplots_adjust(left=0.05, right=0.95, wspace = 0.15, hspace=0)
    fig.set_size_inches(18.0,12.0)
    plt.savefig("spectra_avg_ratio_"+str(l)+".png")


