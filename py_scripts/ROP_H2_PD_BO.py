import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from numpy.random import randn
import os
os.chdir('../')

file = "ropT.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)

for i in range(0,L,1):
    B=A[i].split('\t')
    if B[3] == "1.#INF00e+000":
        B[3] = "3.017931e+255"
        break
L=i
X=randn(L)
Y=randn(L)
W=randn(L)
Z=randn(L)

for i in range(0,L,1):
    B=A[i].split('\t')
    X[i]=float(B[0])
    Y[i]=float(B[1])
    W[i]=float(B[2])
    Z[i]=float(B[3])


ax = plt.subplot(111)
plt.ylim(0,3.2)
plt.xlim(10,3000)
plt.xscale("log")
ax.invert_xaxis()
plt.ylabel('Ortho Para ratio')
plt.xlabel('cosmological doppler shift Z')
plt.plot(X,Y,label="Radex h2--> <-- H + p")
plt.plot(X,W,'g-.',label="Boltzman's Trad")
plt.plot(X,Z,'r--',label="Boltzman's Tbar")
f.close()

file = "PD.txt"
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
    

#plt.plot(X,Y,"r.",label="Flower et Pineau des Forêt (2000)")
f.close()

file = "expansion2.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)
CST=3
S=int((L-CST))
X=randn(S)
Y=randn(S)
CST2=10

for i in range(CST,L,1):
    B=A[i].split('\t')
    I=int(i-CST)
    X[I]=float(B[4]) 
    Y[I]=float(B[CST2-2])

f.close()




#plt.plot(X,Y,label="lsoda h2--> <-- H + p")
leg = plt.legend(loc='best', shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.savefig("rop.png")
plt.show()
f.close()