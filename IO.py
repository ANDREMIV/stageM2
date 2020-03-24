import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from numpy.random import randn

"""
Simple Parser
"""

file = "expansion.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)-2

for i in range(0,L,1):
    B=A[i+2].split('\t')
    if B[3] == "-1.#IND000000000" or B[3] == "1.#QNAN00000000" :
        break
L=i-2
X=randn(L)
Y=randn(L)
Z=randn(L)
for i in range(0,L,1):
    B=A[i+2].split('\t')
    X[i]=float(B[0])
    Y[i]=float(B[1])
   
"""   
f.close()

C=file.split('.')
file=""
for i in range(0,len(C)-1,1):
    file= file+ C[i]
file=file+"trimmed"+C[len(C)-1]

f = open(file, "w")
for i in range(0,nbcase,1):
    f.write("%d\n" %Z[i])
   


f.close()
"""
plt.xlim(-1,0)
plt.ylim(0,1)
plt.ylabel('Scale factor')
plt.xlabel('Normalized cosmic time')
plt.grid(True)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.plot(X,Y)
plt.savefig("Scale_factor.png")
plt.show()

for i in range(0,L,1):
    B=A[i+2].split('\t')
    X[i]=-1/float(B[1])+1
    Y[i]=float(B[2])
    Z[i]=float(B[3])

ax = plt.subplot(111)
plt.ylabel('Temperatures')
plt.xlabel('Cosmic doppler shift z')
#plt.xlim(-7500,0)
plt.ylim(0,19000)
plt.plot(X,Y,label="Trad")
plt.plot(X,Z,label="Tbar")
leg = plt.legend(loc='best', ncol=2, shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid(True)
plt.savefig("Temperatures.png")
plt.show()


file = "XE.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)

for i in range(0,L,1):
    B=A[i].split('\t')
    if B[1] == "-1.#IND000":
        break
L=i
X=randn(L)
Y=randn(L)

for i in range(0,L,1):
    B=A[i].split('\t')
    X[i]=float(B[0])
    Y[i]=float(B[1])


plt.xlim(0,10)
plt.ylim(0,1)
plt.ylabel('ionization')
plt.xlabel('Ti/Tbar')
plt.grid(True)
plt.plot(X,Y)
plt.savefig("XE.png")
plt.show()

file = "expansionDS.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)

for i in range(0,L,1):
    B=A[i].split('\t')
    if B[1] == "-1.#IND000":
        break
L=i
X=randn(L)
Y=randn(L)

for i in range(0,L,1):
    B=A[i].split('\t')
    X[i]=float(B[0])
    Y[i]=float(B[1])


plt.xlim(-5,0)
plt.ylim(0,1)
plt.ylabel('scale factor')
plt.xlabel('cosmic time')
plt.grid(True)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.plot(X,Y)
plt.savefig("DeSitter.png")
plt.show()

file = "expansionEDS.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)

for i in range(0,L,1):
    B=A[i].split('\t')
    if B[1] == "-1.#IND000000000":
        break
L=i
X=randn(L)
Y=randn(L)

for i in range(0,L,1):
    B=A[i].split('\t')
    X[i]=float(B[0])
    Y[i]=float(B[1])


plt.xlim(-1,0)
plt.ylim(0,1)
plt.ylabel('scale factor')
plt.xlabel('cosmic time')
#plt.grid(True)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.plot(X,Y)
# Show the major grid lines with dark grey lines

plt.savefig("EinsteinDeSitter.png")
plt.show()

file = "rop.txt"
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
plt.xlim(10,300)
plt.xscale("log")
ax.invert_xaxis()
plt.ylabel('Ortho Para ratio')
plt.xlabel('cosmological doppler shift Z')
plt.plot(X,Y,label="Radex")
plt.plot(X,W,'g.',label="Boltzman's Trad")
plt.plot(X,Z,'r--',label="Boltzman's Tbar")

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

plt.plot(X,Y,"r.",label="Pineau De Flower")
leg = plt.legend(loc='upper left', ncol=4, shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.savefig("rop.png")
plt.show()

"""
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


ax = plt.subplot(111)
plt.ylim(0,3.5)
plt.xlim(1,1100)
plt.xscale("log")
ax.invert_xaxis()
plt.ylabel('Ortho Para ratio')
plt.xlabel('cosmological doppler shift Z')
plt.plot(X,Y,"r.",label="PD")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.savefig("PD.png")
plt.show()

"""

