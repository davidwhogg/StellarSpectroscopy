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

plate, mjd, fiber, extinction, objid = sort_data()
print len(plate) #check number of stars surveyed




#######Download fits files##########

def downloadfits():
    for i in range(len(plate)):
        plateid=plate[i]
        mjdid=mjd[i]
        fiberid=fiber[i]
        
        if os.path.isfile('spec-'+str(plateid).zfill(4)+'-'+str(mjdid)+'-'+str(fiberid).zfill(4)+'.fits')==True:
            print "already exists", i
        else:
            commands.getoutput('wget --content-disposition "http://api.sdss3.org/spectrum?plate='+str(plateid)+'&fiber='+str(fiberid)+'&mjd='+str(mjdid)+'"')
            print "command success", i
    return

#downloadfits() #can be commented out if already downloaded

########take flux, wl data from fits files##########
def getdata():
        
    tabs=[] #this will contain each fits file in one super-array
    fluxes=[]
    sn2s=[]
    errormags=[]
    wls=[]
    ivars=[]
    badpoints=[]

    ### use this by default
    for i in range(2700):#2700,len(plate)):  #cf 125, 178, 222
        plateid=plate[i]
        mjdid=mjd[i]
        fiberid=fiber[i]
        
        tab = pyfits.open(commands.getoutput("pwd")+'/spec-'+str(plateid).zfill(4)+'-'+str(mjdid)+'-'+str(fiberid).zfill(4)+'.fits')
        print "tab success", i
        tabs.append(tab)
        j=i#-2700 #-2700 ###########
        if type(tabs[j][2].data.field(63)[0])==float32: #distinguish SDSS,BOSS
            zm= tabs[j][2].data.field(63)[0] #redshift
            print i, "63"            
        elif type(tabs[j][2].data.field(37)[0])==float32:
            zm= tabs[j][2].data.field(37)[0]
            print i, "37"
        else:
            print "error"
            
        flux=tabs[j][1].data.field(0) #flux
        fluxes.append(flux)
        loglam=tabs[j][1].data.field(1)
        loglam=np.array(loglam)
        lam=10**loglam
        ivar=tabs[j][1].data.field(2)
        ivars.append(ivar)
        lamcor=zeros(len(lam))
        for k in range(len(lam)): #redshift correction
            lamcor[k]=float(lam[k])/float((1+zm)) 
            
        wls.append(lamcor)
        sn2=tabs[j][2].data.field(6)[0]+tabs[j][2].data.field(7)[0] #one of these entries is 0 always
        sn2s.append(sn2)
        errormag=1/sn2
        errormags.append(errormag)
        bad=[]
        for v in range(len(ivars[j])):    
            if ivars[j][v]!=0:
                ivars[j][v]=1/sqrt(ivars[j][v])
            elif ivars[j][v]==0:
                bad.append(wls[j][v])
                ivars[j][v]=np.inf
        badpoints.append(bad)
        tab.close()
    return wls, fluxes, sn2s, ivars, badpoints

wls, fluxes, sn2s, ivars, badpoints = getdata()




##### Calculate cont, eqw, flux values ###########
def calc():

    #f=open("tightnew", "rb")
    #data=pickle.load(f)
    #f.close()
    data = [[],[],[],[],[],[],[],[]]
    for z in range(2700):#len(plate)-2700): #len(plate)-2700
        s = UnivariateSpline(wls[z], fluxes[z], k=3, s=0)
        xs=linspace(min(wls[z]),max(wls[z]),len(wls[z])*10)
        # peaks are [hd, hg, hb, He-I,   K,    H,   TiO,  O-I]
        peaks=[4102, 4340, 4861, 3890, 3934, 3970, 5582, 6306] # Only using h-delta, h-gamma, h-beta
                

        for w in range(len(peaks)): #for each peak location..
            good = (ivars[z]!=np.inf)
            cont_indx = good*(wls[z]>peaks[w]-100)*(wls[z]<peaks[w]+100)
            flux_indx = (wls[z]>peaks[w]-10)*(wls[z]<peaks[w]+10)

            cont_zone=wls[z][(cont_indx*(True-flux_indx))]

            
            #cont_zone=[x for x in wls[z] if x>peaks[w]-100 and x<peaks[w]+100 and x not in badpoints[z]]
            cont_est=s(cont_zone)
            cont1=np.median(cont_est) #Cont value


            bad = True-good
            fail_indx = bad*(wls[z]>peaks[w]-100)*(wls[z]<peaks[w]+100)
            fail_flag=wls[z][fail_indx]
            #fail_flag=[x for x in wls[z] if x>peaks[w]-100 and x<peaks[w]+100 and x in badpoints[z]]



            flux_zone=wls[z][flux_indx]
            ys = s(flux_zone)

            ###take ivars for the cont range
            errors = np.array(ivars[z][cont_indx])
                               
            cont_err = np.sqrt(sum(errors**2))/len(errors)
            
            ## integral:
            ys_corr = ys-cont1 #ready for integration
            #len flux_domain = len ys_corr
            flux = 0.5*(flux_zone[1]-flux_zone[0])*(2*sum(ys_corr)-ys_corr[0]-ys_corr[len(ys_corr)-1]) #trap rule
    
            #calculate flux error: weighted sum with weights=stepsize
            flux_errors = np.array(ivars[z][flux_indx])
            
            f_errors_squared = sum(flux_errors[1:-1]**2) #sum squares, not first or last value
            f_error_w = 4*(f_errors_squared)+(flux_errors[0]**2+flux_errors[len(flux_errors)-1]**2) #apply weights; 1 + 4 + ... + 4 + 1
            flux_err = np.sqrt(f_error_w)*(wls[z][1]-wls[z][0]) #weight with step-size            

            eqw = flux/cont1
            zz=z#+2700 #or z+2700
            data[w].append([cont1, cont_err, flux, flux_err, eqw, extinction[zz], plate[zz], mjd[zz], fiber[zz], str(int(objid[zz])), len(fail_flag)]) #grouped_data
            print "ok done", z, w
    return data #grouped_data and also 220

data = calc() #save data to file
f2=open("newlines","wb")    
pickle.dump(data,f2) 
f2.close()

########### Plot spectra ######## (need grouped_data) - can only plot 
def plot():
    f=open("datas2","rb")
    #grouped_data=pickle.load(f)
    f.close()
    for z in range(len(grouped_data[0])):# i.e. z=0
        s = UnivariateSpline(wls[z], fluxes[z], k=3, s=0)
        xs=linspace(min(wls[z]),max(wls[z]),len(wls[z])*10)
        peaks=[4102, 4340, 4861] # Only using h-b, h-g, h-d
        
        for w in range(len(peaks)): #for each peak location..
            cont_zone=[x for x in xs if x>peaks[w]-100 and x<peaks[w]+100]
            cont_est=s(cont_zone)
            cont1=np.median(cont_est) #Cont value
            peak_loc_finder=[x for x in xs if x>peaks[w]-2 and x<peaks[w]+2]
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
        plt.xlim(4000,5000)
        plt.step(xs,s(xs),'b', linewidth=0.5, alpha=1) 
        plt.xlabel("wavelengths (A)")
        plt.ylabel("flux (E-17 ergs/s/cm^2/A)")
    
    ################## Plot continuums
        plt.plot(np.array([4762,4962]),np.array([grouped_data[2][n][0]]*2), color='k')
        plt.plot(np.array([4002,4202]),np.array([grouped_data[0][n][0]]*2), color='k')
        plt.plot(np.array([4240,4440]),np.array([grouped_data[1][n][0]]*2), color='k')

        plt.title("Spectrum for plate-mjd-fiber = "+ str(plate[n])+"-"+str(mjd[n])+"-"+str(fiber[n]))
        plt.grid(True)
        plt.show()
    return

    
    #return
#calc()'''
    
#return grouped_data
## grouped_data is separated into 3 columns, one for each peak location
## grouped_data[i] is further split into n entries for the number of stars 
## grouped_data[i][j] is split into [0]=cont, [1]=flux, [2]=eqw, [3]=peak location
