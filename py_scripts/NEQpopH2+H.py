import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from numpy.random import randn
import os
#In each pop file, you need to manually add number of levels, collisionners in y description, level file name, nbofplots

os.chdir('../')
FFile= "STATEQlevels"
file = FFile+".txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)

S=int((L-4))
X=randn(S)
nbl=54
Y=randn(nbl,S)

i=2
QNs=A[i].split('\t')
    

for i in range(4,L,1):
    B=A[i].split('\t')
    I=int(i-4)
    X[I]=float(B[0])   
    for j in range(0,nbl,1):
        Y[j][I]=float(B[j+1])
        

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
    plt.xlim(1e-17,1.1e-10)
    #plt.xscale("log")
    plt.yscale("log")
    plt.ylabel('populations of H2--> <--h')
    plt.xlabel('time t')
    leg = plt.legend(loc='best', shadow=True, fancybox=True)
    leg.get_frame().set_alpha(0.5)
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
    s="%d" %(k)
    plt.savefig("NEQ"+FFile+"."+s+".png")
    plt.show()

f.close()