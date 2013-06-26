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
def eqw_estimation(exta,extb,n):
        
    objids=[] #this one is left empty, ignore

    objids=[]
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

    array1=[objids,extinction,hbflux,hbreqw,hbcont,hdflux,hdreqw,hdcont,objs,plate,fiber,mjd,\
            n2flux,n2reqw,n2cont,magg,hgflux,hgreqw,hgcont]

    query = "SELECT top "+str(n)+" p.objID, \
    p.extinction_g, g.h_beta_reqw_err, g.h_beta_reqw, g.h_beta_cont, g.h_delta_reqw_err, g.h_delta_reqw, g.h_delta_cont, \
    p.obj, s.plate, s.fiberID, s.mjd, g.h_alpha_flux, g.h_alpha_reqw, g.h_alpha_cont, \
    p.type, g.h_gamma_flux, g.h_gamma_reqw, g.h_gamma_cont \
    FROM PhotoObj AS p \
    JOIN SpecObj as s ON s.specobjID=p.specobjID \
    JOIN galSpecLine as g on s.specobjID=g.specobjID \
    JOIN galSpecInfo as info on s.specObjID=info.specobjID \
    WHERE psfMag_r BETWEEN 15.0 and 19.0 \
    and p.type=6 and p.extinction_g between "+str(exta)+" AND "+str(extb)+ \
    " and  dbo.fPhotoStatus('PRIMARY')>0 and dbo.fPhotoFlags('STATIONARY')>0 \
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
    AND (psfMag_i-psfmag_z) between 0.02-0.08 and 0.02+0.08"
    alldata=sqlcl.query(query).read()
    interim=alldata.replace("\n",",")
    compiled=interim.split(",")
    for i in range(20,len(compiled)-1,19):
        extinction.append(float(compiled[i]))
    for i in range(21,len(compiled)-1,19):
        hbflux.append(float(compiled[i]))
    for i in range(22,len(compiled)-1,19):
        hbreqw.append(float(compiled[i]))
    for i in range(23,len(compiled)-1,19):
        hbcont.append((compiled[i]))
    for i in range(24,len(compiled)-1,19):
        hdflux.append(float(compiled[i]))
    for i in range(25,len(compiled)-1,19):
        hdreqw.append(float(compiled[i]))
    for i in range(26,len(compiled)-1,19):
        hdcont.append(float(compiled[i]))
    for i in range(27,len(compiled)-1,19):
        objs.append(float(compiled[i]))
    for i in range(28,len(compiled)-1,19):
        plate.append(int(compiled[i]))
    for i in range(29,len(compiled)-1,19):
        fiber.append(int(compiled[i]))
    for i in range(30,len(compiled)-1,19):
        mjd.append(int(compiled[i]))
    for i in range(31,len(compiled)-1,19):
        n2flux.append(float(compiled[i]))
    for i in range(32,len(compiled)-1,19):
        n2reqw.append(float(compiled[i]))
    for i in range(33,len(compiled)-1,19):
        n2cont.append(float(compiled[i]))
    for i in range(34,len(compiled)-1,19):
        magg.append((compiled[i]))
    for i in range(35,len(compiled)-1,19):
        hgflux.append(float(compiled[i]))
    for i in range(36,len(compiled)-1,19):
        hgreqw.append(float(compiled[i]))
    for i in range(37,len(compiled)-1,19):
        hgcont.append(float(compiled[i]))

    tabs=[] #this will contain each fits file in one super-array
    fluxes=[]
    sn2s=[]
    errormags=[]
    wls=[]
    array2=[fluxes,wls,sn2s,errormags]
    for i in range(1) :
        plateid=plate[i]
        mjdid=mjd[i]
        fiberid=fiber[i]
        commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate='+str(plateid)+'&fiber='+str(fiberid)+'&mjd='+str(mjdid)+'"')
        
        tab = pyfits.open(commands.getoutput("pwd")+'/spec-'+str(plateid).zfill(4)+'-'+str(mjdid)+'-'+str(fiberid).zfill(4)+'.fits')
        tabs.append(tab)
        if type(tabs[i][2].data.field(63)[0])==float32: #distinguish SDSS,BOSS
            zm= tabs[i][2].data.field(63)[0] #redshift
            print i, "63"            
        elif type(tabs[i][2].data.field(37)[0])==float32:
            zm= taba[i][2].data.field(37)[0]
            print i, "37"
        else:
            print "error"
            
        flux=tabs[i][1].data.field(0) #flux
        fluxes.append(flux)
        loglam=tabs[i][1].data.field(1)
        loglam=np.array(loglam)
        lam=10**loglam
        lamcor=zeros(len(lam))
        for k in range(len(lam)): #redshift correction
            lamcor[k]=float(lam[k])/float((1+zm)) 
            
        wls.append(lamcor)
        #wls.append(lam)
        print lam-lamcor
        
        sn2=tabs[i][2].data.field(6)[0]+tabs[i][2].data.field(7)[0] #one of these entries is 0 always
        sn2s.append(sn2)
        errormag=1/sn2
        errormags.append(errormag)

    #below section graphs the spectrum
    betadata=[]
    gammadata=[]
    deltadata=[]
    grouped_data=[deltadata,gammadata,betadata]
    for z in range(1):
        s = UnivariateSpline(wls[z], fluxes[z], k=3, s=0)
        xs=linspace(min(wls[z]),max(wls[z]),len(wls[z])*10)
        
        peaks=[4102, 4340, 4861] # Only using h-b, h-g, h-d
        beta= [x for x in xs if x>4761 and x<4962]
        gamma= [x for x in xs if x>4240 and x<4441]
        delta= [x for x in xs if x>4002 and x<4202]
        windows=[delta, gamma, beta]


        for i in range(len(peaks)):
        
            for k in windows[i]:
                if s(k) == min(s(windows[i])):
                    peak_loc = k #location of minima, wl
                    print peak_loc, i
            flux_domain = [x for x in xs if x>peak_loc-10 and x<peak_loc+10]
            ys=s(flux_domain) #flux values
            ## ESTIMATING CONT: using median
            cont1 = np.median(s(windows[i]))
            ## integral:
            ys_corr = ys-cont1 #ready for integration
            #len flux_domain = len ys_corr
            flux = 0.5*(flux_domain[1]-flux_domain[0])*(2*sum(ys_corr)-ys_corr[0]-ys_corr[len(ys_corr)-1]) #trap rule
            eqw = flux/cont1
            grouped_data[i].append([cont1, flux, eqw, peak_loc,])
            plt.fill_between(flux_domain, ys, cont1, color='gray', alpha=0.5)

        plt.xlim(4000,5000)
        plt.step(xs,s(xs),'b', linewidth=0.5, alpha=1) 
        plt.xlabel("wavelengths (A)")
        plt.ylabel("flux (E-17 ergs/s/cm^2/A)")
        for j in range(len(peaks)): #integral neighborhoods
            plt.axvline(x=grouped_data[j][z][3]+10, color='g')
            plt.axvline(x=grouped_data[j][z][3]-10, color='g')
            plt.axvline(x=grouped_data[j][z][3], color='b')
        plt.axvline(x=4761) #cont neighborhoods
        plt.axvline(x=4962)
        plt.axvline(x=4240)
        plt.axvline(x=4441)
        plt.axvline(x=4002)
        plt.axvline(x=4202)
        ##################
        plt.plot(np.array([4761,4961]),np.array([grouped_data[2][z][0]]*2), color='k')
        plt.plot(np.array([4002,4202]),np.array([grouped_data[0][z][0]]*2), color='k')
        plt.plot(np.array([4240,4440]),np.array([grouped_data[1][z][0]]*2), color='k')

        plt.title("Reported continuum values for extinction ="+ str(extinction[z]))
        print z, "blue flux, cont, reqw @ h-b 4861 is", hbflux[z], hbcont[z], hbreqw[z]
        print z, "blue flux, cont, reqw @ h-d 4102 is", hdflux[z], hdcont[z], hdreqw[z]
        print z, "blue flux, cont, reqw @ h-a 6565 is", n2flux[z], n2cont[z], n2reqw[z]
        print z, "blue flux, cont, reqw @ h-gamma 4340 is", hgflux[z], hgcont[z], hgreqw[z]
        plt.grid(True)
        plt.show()
    return grouped_data
## grouped_data is separated into 3 columns, one for each peak location
## grouped_data[i] is further split into n entries for the n value chosen 
## grouped_data[i][j] is split into [0]=cont, [1]=flux, [2]=eqw, [3]=peak location
