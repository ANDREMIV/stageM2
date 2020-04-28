import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from numpy.random import randn
import os
os.chdir('../')

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
plt.ylabel('facteur d\'échelle a')
plt.xlabel('temps cosmique normalisé t')
#plt.grid(True)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.40)
plt.plot(X,Y,"b^",markevery=200,label="Solution intégrateur")
Y=np.power(3/2*(X+2/3),2.0/3.0)
plt.plot(X,Y,"rx",markevery=200,label="Solution analytique")
leg = plt.legend(loc='best', shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.savefig("EinsteinDeSitter.png")
plt.show()
f.close()