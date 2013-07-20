import numpy as np
import os
from pylab import *
from scipy.interpolate import UnivariateSpline

import commands
import pyfits
import sys
os.chdir("/Users/admin/Desktop/Maser files")
directory=commands.getoutput("pwd")
sys.path.append(directory)
import sqlcl
import pickle

f=open("newlinesnewer","rb")
data=pickle.load(f)
f.close()
#Get the data via SQL query (we really only need plate-mjd-fiber for now)
def sort_data():
    objid=[]
    extinction=[]
    hbflux=[]
    hbreqw=[]
    hbcont=[]
    hdflux=[]
    hdreqw=[]
    hdcont=[]
    objs=[]
    plate=[]
    fiber=[]
    mjd=[]
    n2flux=[] #n2 is now h-a, but both are not used
    n2reqw=[]
    n2cont=[]
    magg=[]
    hgflux=[]
    hgreqw=[]
    hgcont=[]

    array1=[objid,extinction,hbflux,hbreqw,hbcont,hdflux,hdreqw,hdcont,objs,plate,fiber,mjd,\
            n2flux,n2reqw,n2cont,magg,hgflux,hgreqw,hgcont]
    #and p.extinction_g between "+str(exta)+" AND "+str(extb)+ ##
    query = "SELECT p.objID, \
    p.extinction_g, g.h_beta_reqw_err, g.h_beta_reqw, g.h_beta_cont, g.h_delta_reqw_err, g.h_delta_reqw, g.h_delta_cont, \
    p.obj, s.plate, s.fiberID, s.mjd, g.h_alpha_flux, g.h_alpha_reqw, g.h_alpha_cont, \
    p.type, g.h_gamma_flux, g.h_gamma_reqw, g.h_gamma_cont \
    FROM PhotoObj AS p \
    JOIN SpecObj as s ON s.specobjID=p.specobjID \
    JOIN galSpecLine as g on s.specobjID=g.specobjID \
    WHERE psfMag_r BETWEEN 15.0 and 19.0 \
    and p.type=6  \
     and  dbo.fPhotoStatus('PRIMARY')>0 and dbo.fPhotoFlags('STATIONARY')>0 \
    and calibStatus_r=1 \
    and s.elodieTEff!=0 and s.elodieFeH!=0 and s.elodieLogG!=0 \
    and ((flags&dbo.fPhotoFlags('BLENDED')) \
    +(flags&dbo.fPhotoFlags('DEBLEND_TOO_MANY_PEAKS')) + \
    (flags&dbo.fPhotoFlags('SATURATED')) \
    +(flags&dbo.fPhotoFlags('BADSKY'))+ \
    (flags&dbo.fPhotoFlags('COSMIC_RAY')) \
    +(flags&dbo.fPhotoFlags('PEAKS_TOO_CLOSE'))+ \
    (flags&dbo.fPhotoFlags('NOTCHECKED_CENTER')) \
    +(flags&dbo.fPhotoFlags('SATUR_CENTER'))+ \
    (flags&dbo.fPhotoFlags('INTERP_CENTER')) \
    +(flags&dbo.fPhotoFlags('INTERP'))+ \
    (flags&dbo.fPhotoFlags('PSF_FLUX_INTERP')))=0 \
    AND  (psfMag_u-psfmag_g) between 0.82-0.08 and 0.82+0.08 \
    AND (psfMag_g-psfmag_r) between 0.3-0.08 and 0.30+0.08 \
    AND (psfMag_r-psfmag_i) between 0.09-0.08 and 0.09+0.08 \
    AND (psfMag_i-psfmag_z) between 0.02-0.08 and 0.02+0.08 \
    ORDER BY extinction_g DESC"
    alldata=sqlcl.query(query).read()
    interim=alldata.replace("\n",",")
    nent=19 #(number of query columns)
    compiled=interim.split(",")

    #sort all values into arrays:
    for j in range(nent): 
        for i in range(j+nent,len(compiled)-1,nent):
            array1[j].append(float(compiled[i])) #note, plate/mjd/fiber should be int

    #some of the values should be type(int):
    for i in range(len(plate)): 
            plate[i]=int(plate[i])
            mjd[i]=int(mjd[i])
            fiber[i]=int(fiber[i])

    return plate, mjd, fiber, extinction, objid

#plate, mjd, fiber, extinction, objid = sort_data()
#print len(plate) #check number of stars surveyed

########take flux, wl data from fits files##########
def getdata(i):
        
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

    plateid=data[0][i][6]
    mjdid=data[0][i][7]
    fiberid=data[0][i][8]
    
    tab = pyfits.open(commands.getoutput("pwd")+'/spec-'+str(plateid).zfill(4)+'-'+str(mjdid)+'-'+str(fiberid).zfill(4)+'.fits')
    print "tab success", i
    tabs[0]=tab
    j=0 #-2700 ###########
    if type(tabs[0][2].data.field(63)[0])==float32: #distinguish SDSS,BOSS
        zm= tabs[0][2].data.field(63)[0] #redshift
        print i, "63"            
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
    for k in range(len(lam)): #redshift correction
        lamcor[k]=float(lam[k])/float((1+zm)) 
        
    wls[0]=lamcor
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
    return wls, fluxes, sn2s, ivars, badpoints

#wls, fluxes, sn2s, sigmas, badpoints = getdata(8)


########### Plot spectra ######## (need grouped_data) - can only plot 
def plot(n):
    #if data[0][n][7]!=mjd[n]:
        #print "mjd mismatch!"
        #sys.exit()
    s = UnivariateSpline(wls[0], fluxes[0], k=3, s=0)
    xs=linspace(min(wls[0]),max(wls[0]),len(wls[0])*10)
    # peaks are [hd, hg, hb,  TiO,  Na,  He-I,   K,    H ] [6306 O-I]
    peaks=[4102, 4340, 4861, 5582, 5896, 3890, 3934, 3970] 
    for w in range(5,6): #for each peak location..

        cont1=data[w][n][0]
        cont_err=data[w][n][1]
        flux=data[w][n][2]
        flux_err=data[w][n][3]
        eqw=data[w][n][4]
        peak_loc=peaks[w]
        ext=data[w][n][5]

        flux_domain = [x for x in xs if x<peak_loc+10 and x>peak_loc-10] #neighborhood of 20A
        ys = s(flux_domain)
        plt.fill_between(flux_domain, ys, cont1, color='gray', alpha=0.5)
        plt.axvline(x=peak_loc+10, color='b')
        plt.axvline(x=peak_loc-10, color='b')
        plt.axvline(x=peak_loc, color='r')
        #plt.axvline(x=peak_loc+100, color='k')
        #plt.axvline(x=peak_loc-100, color='k')
        #plt.axvline(x=peak_loc-20, color='k')
        #plt.axvline(x=peak_loc+20, color='k')


        plt.plot(np.array([peak_loc-100,peak_loc-20]),np.array([cont1]*2), color='k')
        plt.plot(np.array([peak_loc+20,peak_loc+100]),np.array([cont1]*2), color='k')

        #plt.plot(np.array([3850,3880]),np.array([cont1]*2), color='k')
        #plt.plot(np.array([3900,3920]),np.array([cont1]*2), color='k')

        #countdown = [x for x in badpoints if x<peak_loc and x>peaks[w]-100]
        #countup = [x for x in badpoints if x>peak_loc and x<peaks[w]+100]
        #print peaks[w], "number of non-zero ivar pixels in LHS: ", len(countdown)
        #print peaks[w], "number of non-zero ivar pixels in RHS: ", len(countup)
    k=5 #[hd, hg, hb,  TiO,  Na,  He-I,   K,    H ]
    s1="Cont: "+str(round(data[k][n][0],1)) #h-delta
    s2="Err: "+str(round(data[k][n][1],3))
    s3="Flux: "+str(round(data[k][n][2],1))
    s4="Err: "+str(round(data[k][n][3],1))
    s5="EW: "+str(round(data[k][n][4],2)) 
    plt.text(0.05,0.2, s1, fontsize=11, transform = ax.transAxes)
    plt.text(0.05,0.1, s2, fontsize=11, transform = ax.transAxes)
    plt.text(0.95,0.24, s3, fontsize=11, ha='right', transform = ax.transAxes)
    plt.text(0.95,0.16, s4, fontsize=11, ha='right', transform = ax.transAxes)
    plt.text(0.95,0.08, s5, fontsize=11, ha='right', transform = ax.transAxes)
    
    plt.step(wls[0],fluxes[0]+sigmas[0], 'g', linewidth=0.4, alpha=1)
    plt.step(wls[0],fluxes[0]-sigmas[0],  'g', linewidth=0.4, alpha=1)
    plt.step(xs,s(xs),'b', linewidth=0.5, alpha=1)
    #plot ivar=0 points
    plt.scatter(np.array(badpoints), np.array(s(badpoints)), c='r', marker='o')    
        
    #plt.xlabel("wavelengths (A)")
    #plt.ylabel("flux (E-17 ergs/s/cm^2/A)")

    plt.tight_layout()

    plateid=data[0][n][6]
    mjdid=data[0][n][7]
    fiberid=data[0][n][8]
    plt.title(str(plateid)+"-"+str(mjdid)+"-"+str(fiberid)+".fits")
    plt.grid(True)
    
    plt.xlim(peak_loc-139,peak_loc+139)
    return

'''plt.subplots(nrows=3, ncols=3)
plt.xlabel("Wavelengths, Ang")
plt.ylabel("flux (E-17 ergs/s/cm^2/A)")
plt.title("Examples of failures")'''


succ=[]
fail=[]

for i in range(5200):
    if data[5][i][3]>500:
        fail.append(i)
print len(fail), "fail"            
for i in range(5200):
    if data[5][i][4]<-2 and data[0][i][3]<100:
        succ.append(i)
print len(succ), "succ"


# peaks are [hd, hg, hb,  TiO,  Na,  He-I,   K,    H ] [6306 O-I]
peaks=[4102, 4340, 4861, 5582, 5896, 3890, 3934, 3970]

fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True)
for i in range(1,10):
    n=succs[i-1]#+240#mum[i-1]
    j=330+i
    ax=plt.subplot(j)
    wls, fluxes, sn2s, sigmas, badpoints = getdata(n) #or succs
    plot(n)
    if i in range(1,7):
        plt.setp(ax.get_xticklabels(), visible=False)
    else:
        plt.setp(ax.get_xticklabels(), visible=True)
    
    
    if i==4:
        plt.ylabel("flux (E-17 ergs/s/cm^2/A)")
    elif i==8:
        plt.xlabel("Wavelengths, Ang")
    else:
        pass
    
 

fig.suptitle('Examples with He-I', size='large')
fig.subplots_adjust(left=0.05, right=0.95, wspace = 0.15, hspace=0.15)

#plt.setp([a.get_xticklabels() for a in fig.axes[:-3]], visible=False)

plt.show()


print "The SQL Query has already been assembled and the variables plate, mjd, fiber are ready to be called"
print "Pick any number i between 0 and ", len(plate), " and try the functions: "
print "plot(i)"
