import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

indata=pd.read_table('densmap.dat', sep='\s+', header=None)
print(indata)
densdata=indata.values

plt.imshow(densdata, cmap='coolwarm', interpolation='nearest')
cb=plt.colorbar()
plt.show()
"""
indensdat=pd.read_table('heatmap', sep='\s+', header=None)
npdensdat=indensdat.values

for i in range(len(npdensdat)):
    if npdensdat[i+1,2] != npdensdat[i,2]:
        N=i
        break
print(N)
print(npdensdat[0:N,(1,3,4)])

x = npdensdat[0:N,3]
y = npdensdat[0:N,4]

X, Y = np.meshgrid(x,y)

print(X, Y)
plt.contourf(X,Y,npdensdat[0:N,1])
plt.show()

"""
