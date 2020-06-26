import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from numpy.random import randn
import os
#In each pop file, you need to manually add number of levels, collisionners in y description, level file name, nbofplots

os.chdir('../')
FFile= "expansion2"
file = FFile+".txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)
CST=3
S=int((L-CST))
X=randn(S)
nbl=54
Y=randn(nbl,S)

i=2
QNs=A[i].split('\t')
    
CST2=10
for i in range(CST,L,1):
    B=A[i].split('\t')
    I=int(i-CST)
    X[I]=float(B[4])   
    for j in range(0,nbl,1):
        Y[j][I]=float(B[j+CST2])
        

nbseparation=6
dd=nbl//nbseparation
for k in range(0,nbseparation,1):
    

    for j in range(k*dd,(k+1)*dd,1):
        if j>=nbl:
            break
        else:
            plt.plot(X,Y[j],label=QNs[j])
    ax = plt.subplot(111)
    plt.ylim(1e-4,1)
    plt.xlim(10,3402)
    ax.invert_xaxis()
    plt.xscale("log")
    plt.yscale("log")
    plt.ylabel('populations of H2--> <--h out of statistical equilibrium')
    plt.xlabel('cosmological doppler shift Z')
    leg = plt.legend(loc='best', shadow=True, fancybox=True)
    leg.get_frame().set_alpha(0.5)
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
    s="%d" %(k)
    plt.savefig("NEQ"+FFile+"."+s+".png")
    plt.show()

f.close()