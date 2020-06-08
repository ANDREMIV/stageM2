import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from numpy.random import randn
import os
os.chdir('../')

file = "GAH+.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)
X=randn(L)
Y=randn(L)

for i in range(0,L,1):
    B=A[i].split('\t')
    X[i]=float(B[0])
    Y[i]=float(B[1])



ax = plt.subplot(111)
plt.ylim(1e-28,1e-21)
plt.xlim(1e2,1e4)
plt.xscale("log")
plt.yscale("log")
#ax.invert_xaxis()
plt.ylabel('Cooling per molecule per collisionner [ $erg.s^{-1}.cm^3$ ]')
plt.xlabel('Kinetic Temperature T [K]')
plt.plot(X,Y,",:",label="GAH+")
f.close()

file = "GAHWF07.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)
X=randn(L)
Y=randn(L)

for i in range(0,L,1):
    B=A[i].split('\t')
    X[i]=float(B[0])
    Y[i]=float(B[1])

plt.plot(X,Y,",-.",label="GAHWF07")
f.close()

file = "GAHGP98.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)
X=randn(L)
Y=randn(L)

for i in range(0,L,1):
    B=A[i].split('\t')
    X[i]=float(B[0])
    Y[i]=float(B[1])

#plt.plot(X,Y,label="GAHGP98")
f.close()

file = "SQUAREH2HCP.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)
X=randn(L-1)
Y=randn(L-1)

for i in range(1,L,1):
    B=A[i].split('|\t')  
    X[i-1]=float(B[0])
    Y[i-1]=1e7*float(B[1]) #J/s.cm^3 to erg/s.cm^3

plt.plot(X,Y,",-",label="H")
f.close()

file = "SQUAREH2H+CP.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)
X=randn(L-1)
Y=randn(L-1)

for i in range(1,L,1):
    B=A[i].split('|\t')  
    X[i-1]=float(B[0])
    Y[i-1]=1e7*float(B[1]) #J/s.cm^3 to erg/s.cm^3

plt.plot(X,Y,",-",label="H+")
f.close()

file = "SQUAREH2H'CP.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)
X=randn(L-1)
Y=randn(L-1)

for i in range(1,L,1):
    B=A[i].split('|\t')  
    X[i-1]=float(B[0])
    Y[i-1]=1e7*float(B[1]) #J/s.cm^3 to erg/s.cm^3

plt.plot(X,Y,'-.',label="H'")
f.close()





leg = plt.legend(loc='best', shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.savefig("Cooling_power.png")
plt.show()
f.close()