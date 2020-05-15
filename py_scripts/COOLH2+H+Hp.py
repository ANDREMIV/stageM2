import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from numpy.random import randn
import os
#In each pop file, you need to manually add number of levels, collisionners in y description, level file name, nbofplots

os.chdir('../')
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
nbl=9
Y=randn(nbl,S)

i=0
    

for i in range(2,L,1):
    B=A[i].split('\t')
    I=int((i-2)/2)
    if i%2==0:
        X[I]=float(B[0])
        Y[I]=float(B[3])
   

ax = plt.subplot(111)
plt.ylim(1e-53,1e-33)
plt.xlim(10,3402)
plt.xscale("log")
plt.yscale("log")
ax.invert_xaxis()
plt.ylabel('cooling J/s/cm^3')
plt.xlabel('cosmological doppler shift Z')
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.plot(X,Y)
plt.savefig("cooling.png")
plt.show()

f.close()