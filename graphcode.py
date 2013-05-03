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

objs=array10[8]
teff=array10[2]
feh=array10[3]
extinction=array10[1]
ras=array10[12]
decs=array10[13]
magu=array10[14]
magg=array10[15]
magr=array10[16]
magi=array10[17]
magz=array10[18]

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
			
