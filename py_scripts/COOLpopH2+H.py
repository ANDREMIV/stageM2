import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from numpy.random import randn
import os
#In each pop file, you need to manually add number of levels, collisionners in y description, level file name, nbofplots

os.chdir('../')
FFile= "levelsh2-h"
file = FFile+".txt"
f = open(file, "r")



S=f.read()


A=S.splitlines()
L=len(A)

for i in range(0,L,1):
    B=A[i].split('\t')
S=int((L-1)/2)
X=randn(S)
YP=randn(S)
YN=randn(S)
Z=randn(S)
WP=randn(S)
WN=randn(S)
i=0
    

for i in range(2,L,1):
    B=A[i].split('\t')
    I=int((i-2)/2)
    if i%2==0:
        X[I]=float(B[0])
        YP[I]=float(B[3])*1e6*1e7 #erg/s/cm^3
        if YP[I]<0:
            YN[I]=-1*YP[I]
            YP[I]=0
        else:
            YN[I]=0
        Z[I]=float(B[2])
        WP[I]=float(B[4])
        if WP[I]<0:
            WN[I]=-1*WP[I]
            WP[I]=0
        else:
            WN[I]=0
   

ax = plt.subplot(111)
plt.ylim(1e-42,1e-19)
plt.xlim(10,1000)
plt.xscale("log")
plt.yscale("log")
plt.ylabel('Collisionnal heating-cooling [erg/s]')
plt.xlabel('cosmological doppler shift Z')
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.plot(X,YP,label='>0')
plt.plot(X,YN,label='<0')
leg = plt.legend(loc='best', shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.savefig("coolingZ"+FFile+".png")
plt.show()



ax = plt.subplot(111)
plt.ylim(1e-42,1e-19)
plt.xlim(2,3000)
plt.xscale("log")
plt.yscale("log")
plt.ylabel('Collisionnal heating-cooling [erg/s]')
plt.xlabel('Kinetic Temperature (z) in K')
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.plot(Z,YP,label='>0')
plt.plot(Z,YN,label='<0')
leg = plt.legend(loc='best', shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)

plt.savefig("coolingT"+FFile+".png")
plt.show()




ax = plt.subplot(111)
plt.ylim(1e-24,1e-5)
plt.xlim(10,1000)
plt.xscale("log")
plt.yscale("log")
plt.ylabel('Temperature derivative [K/s]')
plt.xlabel('cosmological doppler shift Z')

plt.plot(X,WP,label='>0')
plt.plot(X,WN,label='<0')

leg = plt.legend(loc='best', shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)

plt.savefig("coolingdT"+FFile+".png")
plt.show()

f.close()