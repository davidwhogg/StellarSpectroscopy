
#import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import UnivariateSpline


#rom pylab import *
#import commands
#import sys
#os.chdir("/Users/admin/Desktop/Maser files")
#directory=commands.getoutput("pwd")
#sys.path.append(directory)
#import sqlcl
import sqldata_notmain
array10=sqldata_notmain.sqldata(10) #the number is not important

wls=array10[1][1][0] #array of wavelengths for one particular star
fluxes=array10[1][0][0] #array of fluxes for one particular star

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

plt.subplot(111) #211 if making two plots
#plt.scatter(np.array(ras),np.array(decs),s=sa, c=ca,marker='o')
plt.scatter(maggr,magri,s=sa, c=ca, marker="o", alpha=0.1, edgecolors = 'none')

#fig = plt.gcf()
#fig.set_size_inches(6,3)
plt.xlabel("g-r")
plt.ylabel("r-i")

#below section graphs the spectrum
'''
plt.subplot(212)
plt.step(wls,fluxes)
plt.xlabel("wavelengths (A)")
plt.ylabel("flux (E-17 ergs/s/cm^2/A)")
fig = plt.gcf()
fig.set_size_inches(6,3)
'''


##### SPLINE
#s1 = UnivariateSpline(lamcor1, flux1, k=3) #these names are not defined
#plot wl, s1

 
#####

plt.show()


## below section plots several graphs together
'''

s=[1,2,3,8,12,13,14,15,16,17,18]
names=["extinction (mag)","TEff (K)","FeH","obj ID", "RA", "Dec", "mag_u","mag_g","mag_r","mag_i","mag_z"]

for i in range(len(s)):
	for j in range(len(s)):
		if s[i]!=s[j] and s[i]<s[j]:
			plt.scatter(np.array(array10[s[i]]),np.array(array10[s[j]]))
			plt.xlabel(names[i])
			plt.ylabel(names[j])
			plt.show()
			#plt.savefig(str(str(i)+"-",str(j)+".png"))
			
'''
