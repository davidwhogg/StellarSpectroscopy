import numpy as np
import pyfits
from pylab import *
from scipy.interpolate import UnivariateSpline

import commands
import os
import sys
os.chdir("/Users/admin/Desktop/Maser files")
directory=commands.getoutput("pwd")
sys.path.append(directory)
import sqlcl
import pickle

# Get Data

f=open("datas2","rb")
data=pickle.load(f)


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

plate, mjd, fiber, extinction, objid = sort_data()
print len(plate) #check number of stars surveyed

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
    if extinction[i]!=data[0][i][4]:
        print "mismatched extinction values!"
        sys.exit()
        
    ### use this by default
    #for q in range(1): #cf 132, 175, 198
    plateid=plate[i]
    mjdid=mjd[i]
    fiberid=fiber[i]
    
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
            badpoints.append(wls[0][v])
            
    tab.close()
    return wls, fluxes, sn2s, ivars, badpoints

#wls, fluxes, sn2s, ivars, badpoints = getdata()




##### Calculate cont, eqw, flux values ########### Don't bother with this 
def calc():
    betadata=[]
    gammadata=[]
    deltadata=[]
    
    grouped_data=[deltadata,gammadata,betadata] 
    for z in range(len(plate)-2700): #len(plate)-2700
        s = UnivariateSpline(wls[z], fluxes[z], k=3, s=0)
        xs=linspace(min(wls[z]),max(wls[z]),len(wls[z])*10)
        
        peaks=[4102, 4340, 4861] # Only using h-b, h-g, h-d

        for w in range(len(peaks)): #for each peak location..
            cont_zone=[x for x in xs if x>peaks[w]-100 and x<peaks[w]+100]
            cont_est=s(cont_zone)
            cont1=np.median(cont_est) #Cont value

            peak_loc_finder=[x for x in xs if x>peaks[w]-5 and x<peaks[w]+5]
            peak_locs=s(peak_loc_finder)
            for q in range(len(peak_locs)):
                if peak_locs[q]==min(peak_locs):
                    peak_loc=peak_loc_finder[q] #peak location found.
            flux_domain = [x for x in xs if x<peak_loc+10 and x>peak_loc-10] #neighborhood of 20A
            ys = s(flux_domain)
            ## integral:
            ys_corr = ys-cont1 #ready for integration
            #len flux_domain = len ys_corr
            flux = 0.5*(flux_domain[1]-flux_domain[0])*(2*sum(ys_corr)-ys_corr[0]-ys_corr[len(ys_corr)-1]) #trap rule
            eqw = flux/cont1
            data[w].append([cont1, flux, eqw, peak_loc, extinction[z+2700], plate[z+2700], mjd[z+2700], fiber[z+2700], str(int(objid[z+2700]))]) #grouped_data
            print "ok done", z, w

    return data #grouped_data and also 198
'''
data = calc() #len(tabs)?  ##grouped_data
f2=open("datas2","wb")    
pickle.dump(data,f2) #grouped_data
f2.close()'''

########### Plot spectra ######## (need grouped_data) - can only plot 
def plot(n):

    s = UnivariateSpline(wls[0], fluxes[0], k=3, s=0)
    xs=linspace(min(wls[0]),max(wls[0]),len(wls[0])*10)
    peaks=[4102, 4340, 4861] # Only using h-b, h-g, h-d
    for w in range(len(peaks)): #for each peak location..
        cont_zone=[x for x in wls[0] if x>peaks[w]-100 and x<peaks[w]+100 and x not in badpoints]
        
        cont_est=s(cont_zone)
        cont1=np.median(cont_est) #Cont value
        peak_loc_finder=[x for x in xs if x>peaks[w]-5 and x<peaks[w]+5]
        peak_locs=s(peak_loc_finder)
        
        for q in range(len(peak_locs)):
            if peak_locs[q]==min(peak_locs):
                peak_loc=peak_loc_finder[q] #peak location found.
        flux_domain = [x for x in xs if x<peak_loc+10 and x>peak_loc-10] #neighborhood of 20A
        ys = s(flux_domain)
        plt.fill_between(flux_domain, ys, cont1, color='gray', alpha=0.5)
        plt.axvline(x=peak_loc+10, color='k')
        plt.axvline(x=peak_loc-10, color='k')
        plt.axvline(x=peak_loc, color='r')
        plt.axvline(x=peak_loc+100, color='k')
        plt.axvline(x=peak_loc-100, color='k')
        plt.plot(np.array([peak_loc-100,peak_loc+100]),np.array([cont1]*2), color='k')
        countdown = [x for x in badpoints if x<peak_loc and x>peaks[w]-100]
        countup = [x for x in badpoints if x>peak_loc and x<peaks[w]+100]
        print peaks[w], "number of non-zero ivar pixels in LHS: ", len(countdown)
        print peaks[w], "number of non-zero ivar pixels in RHS: ", len(countup)
    plt.xlim(4002,5002)
    plt.step(wls[0],fluxes[0]+ivars[0], 'g', linewidth=0.4, alpha=1)
    plt.step(wls[0],fluxes[0]-ivars[0],  'g', linewidth=0.4, alpha=1)
    plt.step(xs,s(xs),'b', linewidth=0.5, alpha=1)
    #plot ivar=0 points
    plt.scatter(np.array(badpoints), np.array(s(badpoints)), c='r', marker='o')    
    
    plt.xlabel("wavelengths (A)")
    plt.ylabel("flux (E-17 ergs/s/cm^2/A)")

    plt.tight_layout()
    plt.title(str(plate[n])+"-"+str(mjd[n])+"-"+str(fiber[n])+".fits")
    plt.grid(True)
    plt.show()
    return

#fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True)

'''plt.subplots(nrows=3, ncols=3)
plt.xlabel("Wavelengths, Ang")
plt.ylabel("flux (E-17 ergs/s/cm^2/A)")
plt.title("Examples of failures")'''

'''ax=plt.subplot(331)
wls, fluxes, sn2s, ivars = getdata(8)
plot(8)

ax=plt.subplot(332)
wls, fluxes, sn2s, ivars = getdata(106)
plot(106)

ax=plt.subplot(333)
wls, fluxes, sn2s, ivars = getdata(736)
plot(736)

ax=plt.subplot(334)
wls, fluxes, sn2s, ivars = getdata(216)
plot(216)
plt.ylabel("flux (E-17 ergs/s/cm^2/A)")

ax=plt.subplot(335)
wls, fluxes, sn2s, ivars = getdata(351)
plot(351)

ax=plt.subplot(336)
wls, fluxes, sn2s, ivars = getdata(1363)
plot(1363)

ax=plt.subplot(337)
wls, fluxes, sn2s, ivars = getdata(1358)
plot(1358)

ax=plt.subplot(338)
wls, fluxes, sn2s, ivars = getdata(1465)
plot(1465)
plt.xlabel("Wavelengths, Ang")

ax=plt.subplot(339)
wls, fluxes, sn2s, ivars = getdata(1663)
plot(1663)

fig.suptitle('Failed examples')
fig.subplots_adjust(hspace=0.2)
#plt.setp([a.get_xticklabels() for a in fig.axes[:-3]], visible=False)
plt.show()
'''

print "The SQL Query has already been assembled and the variables plate, mjd, fiber are ready to be called"
print "Pick any number i between 0 and ", len(plate), " and try the functions: "
print "#wls, fluxes, sn2s, ivars, badpoints = getdata(i)"
print "and plot(i)"
