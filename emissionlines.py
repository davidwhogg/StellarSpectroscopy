
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

diff= 0.08 #one of the four should always be commented out
diff1= 0.08
diff2= 0.08
diff3= 0.08
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
n2flux=[]
n2reqw=[]
n2cont=[]
magg=[]
hgflux=[]
hgreqw=[]
hgcont=[]

array1=[objids,extinction,hbflux,hbreqw,hbcont,hdflux,hdreqw,hdcont,objs,plate,fiber,mjd,\
        n2flux,n2reqw,n2cont,magg,hgflux,hgreqw,hgcont]


query = "SELECT top 10 p.objID, \
p.extinction_g, g.h_beta_reqw_err, g.h_beta_reqw, g.h_beta_cont, g.h_delta_reqw_err, g.h_delta_reqw, g.h_delta_cont, \
p.obj, s.plate, s.fiberID, s.mjd, g.h_alpha_flux, g.h_alpha_reqw, g.h_alpha_cont, \
p.type, g.h_gamma_flux, g.h_gamma_reqw, g.h_gamma_cont \
FROM PhotoObj AS p \
JOIN SpecObj as s ON s.specobjID=p.specobjID \
JOIN galSpecLine as g on s.specobjID=g.specobjID \
JOIN galSpecInfo as info on s.specObjID=info.specobjID \
WHERE psfMag_r BETWEEN 15.0 and 19.0 \
and p.type=6 and p.extinction_g BETWEEN 0.0215 and 0.226 \
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
AND  (psfMag_u-psfmag_g) between 0.82-" +str(diff)+ " and 0.82+"+str(diff)+ " \
AND (psfMag_g-psfmag_r) between 0.3-" +str(diff1)+ " and 0.30+"+str(diff1) + " \
AND (psfMag_r-psfmag_i) between 0.09-" +str(diff2)+  " and 0.09+" +str(diff2) + " \
AND (psfMag_i-psfmag_z) between 0.02-" +str(diff3)+ " and 0.02+" +str(diff3)
alldata=sqlcl.query(query).read()

#  r=0.53183
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
for i in range(len(hgcont)) :
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
'''

objids2=[] #this one is left empty, ignore
objids2=[]
extinction2=[]
hbflux2=[]
hbreqw2=[]
hbcont2=[]
hdflux2=[]
hdreqw2=[]
hdcont2=[]
objs2=[]
plate2=[]
fiber2=[]
mjd2=[]
n2flux2=[]
n2reqw2=[]
n2cont2=[]
magg2=[]
hgflux2=[]
hgreqw2=[]
hgcont2=[]
array12=[objids2,extinction2,hbflux2,hbreqw2,hbcont2,hdflux2,hdreqw2,hdcont2,objs2,plate2,fiber2,mjd2,\
        n2flux2,n2reqw2,n2cont2,magg2,hgflux2,hgreqw2,hgcont2]

query2 = "SELECT top 9 p.objID, \
p.extinction_g, g.h_beta_reqw_err, g.h_beta_reqw, g.h_beta_cont, g.h_delta_reqw_err, g.h_delta_reqw, g.h_delta_cont, \
p.obj, s.plate, s.fiberID, s.mjd, g.h_alpha_flux, g.h_alpha_reqw, g.h_alpha_cont,\
p.type, g.h_gamma_flux, g.h_gamma_reqw, g.h_gamma_cont \
FROM PhotoObj AS p \
JOIN SpecObj as s ON s.specobjID=p.specobjID \
JOIN galSpecLine as g on s.specobjID=g.specobjID \
JOIN galSpecInfo as info on s.specObjID=info.specobjID \
WHERE psfMag_r BETWEEN 15.0 and 19.0 \
and p.type=6 and p.extinction_g>0.35 \
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
AND  (psfMag_u-psfmag_g) between 0.82-" +str(diff)+ " and 0.82+"+str(diff)+ " \
AND (psfMag_g-psfmag_r) between 0.3-" +str(diff1)+ " and 0.30+"+str(diff1) + " \
AND (psfMag_r-psfmag_i) between 0.09-" +str(diff2)+  " and 0.09+" +str(diff2) + " \
AND (psfMag_i-psfmag_z) between 0.02-" +str(diff3)+ " and 0.02+" +str(diff3)
alldata2=sqlcl.query(query2).read()
#  r=0.53183
interim2=alldata2.replace("\n",",")
compiled2=interim2.split(",")

for i in range(20,len(compiled2)-1,19):
    extinction2.append(float(compiled2[i]))
for i in range(21,len(compiled2)-1,19):
    hbflux2.append(float(compiled2[i]))
for i in range(22,len(compiled2)-1,19):
    hbreqw2.append(float(compiled2[i]))
for i in range(23,len(compiled2)-1,19):
    hbcont2.append((compiled2[i]))
for i in range(24,len(compiled2)-1,19):
    hdflux2.append(float(compiled2[i]))
for i in range(25,len(compiled2)-1,19):
    hdreqw2.append(float(compiled2[i]))
for i in range(26,len(compiled2)-1,19):
    hdcont2.append(float(compiled2[i]))
for i in range(27,len(compiled2)-1,19):
    objs2.append(float(compiled2[i]))
for i in range(28,len(compiled2)-1,19):
    plate2.append(int(compiled2[i]))
for i in range(29,len(compiled2)-1,19):
    fiber2.append(int(compiled2[i]))
for i in range(30,len(compiled2)-1,19):
    mjd2.append(int(compiled2[i]))
for i in range(31,len(compiled2)-1,19):
    n2flux2.append(float(compiled2[i]))
for i in range(32,len(compiled2)-1,19):
    n2reqw2.append(float(compiled2[i]))
for i in range(33,len(compiled2)-1,19):
    n2cont2.append(float(compiled2[i]))
for i in range(34,len(compiled2)-1,19):
    magg2.append((compiled2[i]))
for i in range(35,len(compiled2)-1,19):
    hgflux2.append(float(compiled2[i]))
for i in range(36,len(compiled2)-1,19):
    hgreqw2.append(float(compiled2[i]))
for i in range(37,len(compiled2)-1,19):
    hgcont2.append(float(compiled2[i]))
    

print len(hgcont), len(hgcont2)

tabs2=[] #this will contain each fits file in one super-array
fluxes2=[]
sn2s2=[]
errormags2=[]
wls2=[]
array22=[fluxes2,wls2,sn2s2,errormags2]
for i in range(len(hgcont2)) :
    plateid2=plate2[i]
    mjdid2=mjd2[i]
    fiberid2=fiber2[i]
    commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate='+str(plateid2)+'&fiber='+str(fiberid2)+'&mjd='+str(mjdid2)+'"')
    
    tab2 = pyfits.open(commands.getoutput("pwd")+'/spec-'+str(plateid2).zfill(4)+'-'+str(mjdid2)+'-'+str(fiberid2).zfill(4)+'.fits')
    tabs2.append(tab2)

    if type(tabs2[i][2].data.field(63)[0])==float32: #distinguish SDSS,BOSS
        zm2= tabs2[i][2].data.field(63)[0] #redshift
        print i, "63"            
    elif type(tabs2[i][2].data.field(37)[0])==float32:
        zm2= taba2[i][2].data.field(37)[0]
        print i, "37"
    else:
        print "error"
        
    flux2=tabs2[i][1].data.field(0) #flux
    fluxes2.append(flux2)
    loglam2=tabs2[i][1].data.field(1)
    loglam2=np.array(loglam2)
    lam2=10**loglam2
    lamcor2=zeros(len(lam2))
    for k in range(len(lam2)): #redshift correction
        lamcor2[k]=float(lam2[k])/float((1+zm2))

    #wls2.append(lam2)
    wls2.append(lamcor2)
    sn22=tabs2[i][2].data.field(6)[0]+tabs2[i][2].data.field(7)[0] #one of these entries is 0 always
    sn2s2.append(sn22)
    errormag2=1/sn22
    errormags2.append(errormag2)

'''

##factors
#factor=s2(8000)/s(8000)

#below section graphs the spectrum
for z in range(10):
        
    s = UnivariateSpline(wls[z], fluxes[z], k=3, s=0)
    xs=linspace(min(wls[z]),max(wls[z]),len(wls[z])*10)
    ys=s(xs)

    #s2 = UnivariateSpline(wls2[z], fluxes2[z], k=3, s=0)
    #xs2=linspace(min(wls2[z]),max(wls2[z]),len(wls2[z])*10)
    #ys2=s2(xs2)
    plt.xlim(4000,5000)
    #plt.step(xs2,ys2,'r', linewidth=0.5, alpha=1)
    plt.step(xs,ys,'b', linewidth=0.2, alpha=1) 
    plt.xlabel("wavelengths (A)")
    plt.ylabel("flux (E-17 ergs/s/cm^2/A)")
    plt.axhline(y=hbcont[z], xmin=0.77, xmax=0.98, color='r', label="@ 4861A")
    plt.axhline(y=hdcont[z], xmin=0, xmax=0.2, color='k', label="@ 4102A")
    plt.axhline(y=hgcont[z], xmin=0.24, xmax=0.44, color='g', label="@ 4340A")
    plt.legend(loc=8)
    #plt.title("spectrum for high (red), low (blue) extinction")
    plt.title("Reported continuum values for extinction ="+ str(extinction[z]))
    #print z, "red flux, cont, reqw @ h-b 4861 is", hbflux2[z], hbcont2[z], hbreqw2[z]
    #print z, "red flux, cont, reqw @ h-d 4102 is", hdflux2[z], hdcont2[z], hdreqw2[z]
    #print z, "red flux, cont, reqw @ h-a 6565 is", n2flux2[z], n2cont2[z], n2reqw2[z]
    #print z, "red flux, cont, reqw @ h-gamma 4340 is", hgflux2[z], hgcont2[z], hgreqw2[z]
    print z, "blue flux, cont, reqw @ h-b 4861 is", hbflux[z], hbcont[z], hbreqw[z]
    print z, "blue flux, cont, reqw @ h-d 4102 is", hdflux[z], hdcont[z], hdreqw[z]
    print z, "blue flux, cont, reqw @ h-a 6565 is", n2flux[z], n2cont[z], n2reqw[z]
    print z, "blue flux, cont, reqw @ h-gamma 4340 is", hgflux[z], hgcont[z], hgreqw[z]
    plt.show()
