#import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os

#rom pylab import *
#import commands
#import sys
#os.chdir("/Users/admin/Desktop/Maser files")
#directory=commands.getoutput("pwd")
#sys.path.append(directory)
#import sqlcl
import sqldata_notmain
array10=sqldata_notmain.sqldata(10)

wls=array10[1][1][0] #just one spectrum
fluxes=array10[1][0][0]

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
c=['b']*len(objs)
c[12]='r'
s=[5]*len(objs)
s[12]=100
sa=np.array(s)
ca=np.array(c)
maggr=np.array(magg)-np.array(magr)
magug=np.array(magu)-np.array(magg)
magri=np.array(magr)-np.array(magi)

plt.subplot(211)
plt.scatter(np.array(ras),np.array(decs),s=sa, c=ca,marker='o')
fig = plt.gcf()
fig.set_size_inches(6,3)
plt.xlabel("ra")
plt.ylabel("dec")

plt.subplot(212)
plt.step(wls,fluxes)
plt.xlabel("wavelengths (A)")
plt.ylabel("flux (E-17 ergs/s/cm^2/A)")
fig = plt.gcf()
fig.set_size_inches(6,3)

plt.show()



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
