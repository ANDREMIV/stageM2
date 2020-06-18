import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from numpy.random import randn
import os
os.chdir('../')


#Temperatures
file = "expansion.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)-2
for i in range(0,L,1):
    B=A[i+2].split('\t')
    if B[3] == "-1.#IND000000000" or B[3] == "1.#QNAN00000000" or B[3] == "-1.#INF000000000":
        break
L=i
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
   
"""


plt.xlim(0,1)
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
f.close()


for i in range(0,L,1):
    B=A[i+2].split('\t')
    X[i]=1/float(B[1])
    Y[i]=float(B[2])
    Z[i]=float(B[3])
    

ax = plt.subplot(111)
plt.ylabel('Temperatures')
plt.xlabel('Cosmic doppler shift z+1')
plt.xlim(1,3402)
plt.xscale("log")
ax.invert_xaxis()
plt.ylim(1,10000)
plt.yscale("log")
plt.plot(X,Y,label="Trad")
#plt.plot(X,Z,label="Tbar")


file = "expansion1e0.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)-2
for i in range(0,L,1):
    B=A[i+2].split('\t')
    if B[3] == "-1.#IND000000000" or B[3] == "1.#QNAN00000000" or B[3] == "-1.#INF000000000":
        break
L=i
X=randn(L)
Y=randn(L)
Z=randn(L)
for i in range(0,L,1):
    B=A[i+2].split('\t')
    X[i]=1/float(B[1])
    Y[i]=float(B[2])
    Z[i]=float(B[3])

plt.plot(X,Z,label="Tb avec xe = 1e-0")

f.close()

file = "expansion1e-1.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)-2
for i in range(0,L,1):
    B=A[i+2].split('\t')
    if B[3] == "-1.#IND000000000" or B[3] == "1.#QNAN00000000" or B[3] == "-1.#INF000000000":
        break
L=i
X=randn(L)
Y=randn(L)
Z=randn(L)
for i in range(0,L,1):
    B=A[i+2].split('\t')
    X[i]=1/float(B[1])
    Y[i]=float(B[2])
    Z[i]=float(B[3])

plt.plot(X,Z,label="Tb avec xe = 1e-1")

f.close()

file = "expansion1e-2.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)-2
for i in range(0,L,1):
    B=A[i+2].split('\t')
    if B[3] == "-1.#IND000000000" or B[3] == "1.#QNAN00000000" or B[3] == "-1.#INF000000000":
        break
L=i
X=randn(L)
Y=randn(L)
Z=randn(L)
for i in range(0,L,1):
    B=A[i+2].split('\t')
    X[i]=1/float(B[1])
    Y[i]=float(B[2])
    Z[i]=float(B[3])

plt.plot(X,Z,label="Tb avec xe = 1e-2")

f.close()

file = "expansion1e-3.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)-2
for i in range(0,L,1):
    B=A[i+2].split('\t')
    if B[3] == "-1.#IND000000000" or B[3] == "1.#QNAN00000000" or B[3] == "-1.#INF000000000":
        break
L=i
X=randn(L)
Y=randn(L)
Z=randn(L)
for i in range(0,L,1):
    B=A[i+2].split('\t')
    X[i]=1/float(B[1])
    Y[i]=float(B[2])
    Z[i]=float(B[3])

plt.plot(X,Z,label="Tb avec xe = 1e-3")

f.close()

file = "expansion1e-4.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)-2
for i in range(0,L,1):
    B=A[i+2].split('\t')
    if B[3] == "-1.#IND000000000" or B[3] == "1.#QNAN00000000" or B[3] == "-1.#INF000000000":
        break
L=i
X=randn(L)
Y=randn(L)
Z=randn(L)
for i in range(0,L,1):
    B=A[i+2].split('\t')
    X[i]=1/float(B[1])
    Y[i]=float(B[2])
    Z[i]=float(B[3])

plt.plot(X,Z,label="Tb avec xe = 1e-4")

f.close()

file = "expansion1e-5.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)-2
for i in range(0,L,1):
    B=A[i+2].split('\t')
    if B[3] == "-1.#IND000000000" or B[3] == "1.#QNAN00000000" or B[3] == "-1.#INF000000000":
        break
L=i
X=randn(L)
Y=randn(L)
Z=randn(L)
for i in range(0,L,1):
    B=A[i+2].split('\t')
    X[i]=1/float(B[1])
    Y[i]=float(B[2])
    Z[i]=float(B[3])

plt.plot(X,Z,label="Tb avec xe = 1e-5")

f.close()

file = "expansion_xe.txt"
f = open(file, "r")

S=f.read()


A=S.splitlines()
L=len(A)-2
for i in range(0,L,1):
    B=A[i+2].split('\t')
    if B[3] == "-1.#IND000000000" or B[3] == "1.#QNAN00000000" or B[3] == "-1.#INF000000000":
        break
L=i
X=randn(L)
Y=randn(L)
Z=randn(L)
for i in range(0,L,1):
    B=A[i+2].split('\t')
    X[i]=1/float(B[1])
    Y[i]=float(B[2])
    Z[i]=float(B[3])

plt.plot(X,Z,label="Tb avec xe = 0")

f.close()




leg = plt.legend(loc='best', shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid(True)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.savefig("Temperatures.png")
plt.show()


