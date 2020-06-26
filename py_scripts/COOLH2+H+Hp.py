import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from numpy.random import randn
import os
#In each pop file, you need to manually add number of levels, collisionners in y description, level file name, nbofplots

os.chdir('../')



file = "PD2.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)

for i in range(0,L,1):
    B=A[i].split('\t')
    if B[1] == "1.#INF00e+000":
        B[1] = "3.017931e+255"
        break
L=i
X=randn(L)
Y=randn(L)

for i in range(0,L,1):
    B=A[i].split('\t')
    X[i]=float(B[0])
    Y[i]=float(B[1])

ax = plt.subplot(111)
plt.ylim(1e-46,1e-31)
plt.xlim(10,200)
plt.xscale("log")
ax.invert_xaxis()
plt.yscale("log")
plt.ylabel('Collisionnal heating-cooling [ $erg.s^{-1}.cm^{-3}$ ]')
plt.xlabel('cosmological doppler shift Z')
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.plot(X,Y,"r",label="Flower and Pineau des Forêt (2000)")
f.close()


FFile= "levelsh2-h-H+-rot"
file = FFile+".txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)

for i in range(0,L,1):
    B=A[i].split('\t')
S=int((L-1)/2)
X=randn(S)
V=randn(S)
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
        YP[I]=float(B[3])*1e7 #convert from J/s/cm^3 to erg/s/cm^3
        if YP[I]<0:
            YN[I]=-1*YP[I]
            YP[I]=0
        else:
            YN[I]=0
        Z[I]=float(B[2])
        WP[I]=float(B[4])
        V[I]=float(B[5])
        if WP[I]<0:
            WN[I]=-1*WP[I]
            WP[I]=0
        else:
            WN[I]=0
   



plt.plot(X,YP,label='RADEX')
#plt.plot(X,YN,label='<0')

file = "expansion2.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)
CST=3
S=int((L-CST))
XA=randn(S)
YA=randn(S)
ZA=randn(S)
CST2=10

for i in range(CST,L,1):
    B=A[i].split('\t')
    I=int(i-CST)
    XA[I]=float(B[4]) 
    YA[I]=float(B[CST2-1])*1e7 #convert from J/s/cm^3 to erg/s/cm^3
    ZA[I]=

plt.plot(XA,YA,label="lsoda")
f.close()


leg = plt.legend(loc='best', shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)

plt.savefig("coolingZ"+FFile+".png")
plt.show()



ax = plt.subplot(111)
plt.ylim(1e-46,1e-31)
plt.xlim(2,500)
plt.xscale("log")
plt.yscale("log")
plt.ylabel('Collisionnal heating-cooling [ $erg.s^{-1}.cm^{-3}$ ]')
plt.xlabel('Kinetic Temperature (z) in K')

plt.plot(Z,YP,label='our model')
#plt.plot(Z,YN,label='<0')
leg = plt.legend(loc='best', shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)

plt.savefig("coolingT"+FFile+".png")
plt.show()




ax = plt.subplot(111)
plt.ylim(1e-26,1e-15)
plt.xlim(10,200)
plt.xscale("log")
plt.yscale("log")
plt.ylabel('Temperature derivative [K/s]')
plt.xlabel('cosmological doppler shift Z')

plt.plot(X,WP)#,label='>0')
#plt.plot(X,WN,label='<0')

#leg = plt.legend(loc='best', shadow=True, fancybox=True)
#leg.get_frame().set_alpha(0.5)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
ax.minorticks_on()
ax.grid(b=True, which='minor', color='#999999',axis='both', linestyle='-', alpha=0.40)

plt.savefig("coolingdT"+FFile+".png")
plt.show()




f.close()

