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

#plate, mjd, fiber, extinction, objid = sort_data()
#print len(plate) #check number of stars surveyed

f=open("newlinesnew","rb")
data=pickle.load(f)
f.close()
plate=[]
mjd=[]
fiber=[]
objid=[]
extinction=[]

for i in range(5204):
    plate.append(data[0][i][6])
    mjd.append(data[0][i][7])
    fiber.append(data[0][i][8])
    extinction.append(data[0][i][5])
    objid.append(data[0][i][9])


    


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
def getdata(a,b):
        
    tabs=[] #this will contain each fits file in one super-array
    fluxes=[]
    sn2s=[]
    errormags=[]
    wls=[]
    ivars=[]
    badpoints=[]
    sigmas=[]


    ### use this by default
    for i in range(a,b):  #cf 125, 178, 222
        plateid=plate[i]
        mjdid=mjd[i]
        fiberid=fiber[i]
        
        tab = pyfits.open(commands.getoutput("pwd")+'/spec-'+str(plateid).zfill(4)+'-'+str(mjdid)+'-'+str(fiberid).zfill(4)+'.fits')
        tabs.append(tab)
        j=i-a ###########
        if type(tabs[j][2].data.field(63)[0])==float32: #distinguish SDSS,BOSS
            zm= tabs[j][2].data.field(63)[0] #redshift
        elif type(tabs[j][2].data.field(37)[0])==float32:
            zm= tabs[j][2].data.field(37)[0]
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
        sigma=np.zeros(len(ivars[j]))

        for v in range(len(ivars[j])):
            if ivars[j][v]!=0:
                sigma[v]=1/sqrt(ivars[j][v])
            elif ivars[j][v]==0:
                bad.append(wls[j][v])
                sigma[v]=np.inf
        sigmas.append(sigma)
        badpoints.append(bad)
        tab.close()
    return wls, fluxes, sn2s, sigmas, badpoints

wls, fluxes, sn2s, sigmas, badpoints = getdata(2700,5204)


lines=[\
    ["H-delta",4102,[4002,4082,4122,4202]],\
    ["H-gamma",4340,[4240,4320,4360,4440]],\
    ["H-beta",4861,[4761,4841,4881,4961]],\
    ["TiO",5582,[5482,5462,5602,5682]],\
    ["Na",5896,[5796,5876,5916,5996]],\
    ["He-I",3890,[3850,3880,3900,3920]],\
    ["K",3934,[3850,3880,3900,3920]],\
    ["H",3970,[3850,3880,3900,3920]]\
    ]
##### Calculate cont, eqw, flux values ###########
# lines=[line.....]
# line = ["name", peakloc, cont region]. line[2][0:3] has cont region. line[1] is peakloc.
def calc(a,b,lines): #the a and be should be same as the getdata(a,b)
    f=open("newlinesnewer", "rb")
    data=pickle.load(f)
    f.close()
    #data = [[],[],[],[],[],[],[],[]]
    for z in range(b-a): 
        s = UnivariateSpline(wls[z], fluxes[z], k=3, s=0)
        xs=linspace(min(wls[z]),max(wls[z]),len(wls[z])*10)
                

        for w in range(5): #for each peak location up to Na-5896
            good = (sigmas[z]!=np.inf)
            cont_prim = (wls[z]>lines[w][2][0])*(wls[z]<lines[w][2][1])
            cont_sec = (wls[z]>lines[w][2][2])*(wls[z]<lines[w][2][3])
            cont_indx = good*(cont_sec+cont_prim) #cont ignores +-20A from peak
            
            flux_indx = (wls[z]>lines[w][1]-10)*(wls[z]<lines[w][1]+10) #flux is +-10A from peak

            zz=z+a #or z+2700  ##need this for indexing subsequent parts


            #skip peaks if there are no wavelength-flux data for that peak
            if sum(flux_indx)==0 or sum(cont_indx)==0:
                data[w].append([0,0,0,0,0,extinction[zz],plate[zz],mjd[zz],fiber[zz],str(int(objid[zz])),0])
                print "skipped ", z, w

            else:

                ##Calculate Continuum
                cont_zone=wls[z][(cont_indx*(np.logical_not(flux_indx)))]
                cont_est=s(cont_zone) #apply spline to wavelengths
                cont1=np.median(cont_est) #take the median

                #Find bad points
                bad = np.logical_not(good)
                fail_indx = bad*(wls[z]>lines[w][1]-100)*(wls[z]<lines[w][1]+100)
                fail_flag=wls[z][fail_indx]

                #Calculate cont error
                if float(sum(good*cont_prim))/sum(cont_prim)<0.25 or float(sum(good*cont_sec))/sum(cont_sec)<0.25:
                    cont_err=np.inf
                else:
                    errors = np.array(sigmas[z][cont_indx])                   
                    cont_err = np.sqrt(sum(errors**2))/len(errors)
                
                
                #Calculate Flux
                flux_zone=wls[z][flux_indx]
                ys = s(flux_zone)
                ys_corr = ys-cont1 #subtract continuum
                flux = 0.5*(flux_zone[1]-flux_zone[0])*(2*sum(ys_corr)-ys_corr[0]-ys_corr[len(ys_corr)-1]) #trap rule
        
                #calculate flux error: weighted sum with weights=stepsize
                flux_errors = np.array(sigmas[z][flux_indx])
                f_errors_squared = sum(flux_errors[1:-1]**2) #sum squares, not first or last value
                f_error_w = 4*(f_errors_squared)+(flux_errors[0]**2+flux_errors[len(flux_errors)-1]**2) #apply weights; 1 + 4 + ... + 4 + 1
                flux_err = np.sqrt(f_error_w)*(wls[z][1]-wls[z][0]) #weight with step-size            

                #calculate EW
                eqw = flux/cont1


                data[w].append([cont1, cont_err, flux, flux_err, eqw, extinction[zz], plate[zz], mjd[zz], fiber[zz], str(int(objid[zz])), len(fail_flag)]) #grouped_data

        ######### TREAT He-I, H, K lines separately.
        for w in range(5,8): #this can also be if lines[w][1]==3890, etc
            good = (sigmas[z]!=np.inf)
            cont_prim = (wls[z]>3850)*(wls[z]<3880)
            cont_sec = (wls[z]<3920)*(wls[z]>3900)
            cont_indx = good*(cont_sec+cont_prim) #assume all three lines have same cont
            
            flux_indx = (wls[z]>lines[w][1]-10)*(wls[z]<lines[w][1]+10) #flux is +-10A from peak

            zz=z+a #or z+2700  ##need this for indexing subsequent parts


            #skip peaks if there are no wavelength-flux data for that peak
            if sum(flux_indx)==0 or sum(cont_indx)==0:
                data[w].append([0,0,0,0,0,extinction[zz],plate[zz],mjd[zz],fiber[zz],str(int(objid[zz])),0])
                print "skipped ", z, w

            else:

                ##Calculate Continuum
                cont_zone=wls[z][(cont_indx*(np.logical_not(flux_indx)))]
                cont_est=s(cont_zone) #apply spline to wavelengths
                cont1=np.median(cont_est) #take the median

                #Find bad points
                bad = np.logical_not(good)
                fail_indx = bad*(wls[z]>lines[w][1]-100)*(wls[z]<lines[w][1]+100)
                fail_flag=wls[z][fail_indx]

                #Calculate cont error
                if float(sum(good*cont_prim))/sum(cont_prim)<0.25 or float(sum(good*cont_sec))/sum(cont_sec)<0.25:
                    cont_err=np.inf
                else:
                    errors = np.array(sigmas[z][cont_indx])                   
                    cont_err = np.sqrt(sum(errors**2))/len(errors)
                
                #Canewlinesnewerlculate Flux
                flux_zone=wls[z][flux_indx]
                ys = s(flux_zone)
                ys_corr = ys-cont1 #subtract continuum
                flux = 0.5*(flux_zone[1]-flux_zone[0])*(2*sum(ys_corr)-ys_corr[0]-ys_corr[len(ys_corr)-1]) #trap rule
        
                #calculate flux error: weighted sum with weights=stepsize
                flux_errors = np.array(sigmas[z][flux_indx])
                f_errors_squared = sum(flux_errors[1:-1]**2) #sum squares, not first or last value
                f_error_w = 4*(f_errors_squared)+(flux_errors[0]**2+flux_errors[len(flux_errors)-1]**2) #apply weights; 1 + 4 + ... + 4 + 1
                flux_err = np.sqrt(f_error_w)*(wls[z][1]-wls[z][0]) #weight with step-size            

                #calculate EW
                eqw = flux/cont1


                data[w].append([cont1, cont_err, flux, flux_err, eqw, extinction[zz], plate[zz], mjd[zz], fiber[zz], str(int(objid[zz])), len(fail_flag)]) #grouped_data
        print "ok done", z
    return data #grouped_data and(sum(cont_prim)+0.01) also 220

data = calc(2700,5204,lines) #save data to file
f2=open("newlinesnewer","wb")    
pickle.dump(data,f2) 
f2.close()

    
#return grouped_data
## grouped_data is separated into 3 columns, one for each peak location
## grouped_data[i] is further split into n entries for the number of stars 
## grouped_data[i][j] is split into [0]=cont, [1]=flux, [2]=eqw, [3]=peak location

