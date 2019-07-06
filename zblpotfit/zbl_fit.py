import math
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import os

#ZBL constants
b=[0.18175,0.50986,0.28022,0.02817]
c=[3.1998,0.94229,0.40290,0.20162]

#Buckingham potential, r in Angstroms
def buck(A,rho,C,r):
    bpot=A*math.exp(-r/rho)-C/r**6
    return bpot

#ZBL potential
def zbl(ZA,ZB,r):
    rdash=(0.88534*0.52917721067)/(ZA**(0.23)+ZB**(0.23)) #bohr radius in Angs
    zblsum=0.0
    for i in range(len(b)):
        zblsum=zblsum+b[i]*math.exp(-c[i]*r/rdash)
    zblpot= ((ZA*ZB)*1389354.835/r)*zblsum #Found in dlpoly constants r4pie0*10 to make Angs J/mol
    zblpot=zblpot*1.04E-5 #convert to eV
    return zblpot

#Switching funtion
def switchf(r,rm,eta):
    if r<rm:
        return (1-(math.exp((r-rm)/eta))/2)
    if r>=rm:
        return (math.exp((rm-r)/eta)/2)

#dUbuck/dr
def dBuck(A,rho,C,r):
    bpot=-(A/rho)*math.exp(-r/rho)+6*C/r**7
    return bpot

#dUZBL/dr
def dzbl(ZA,ZB,r):
    rdash=(0.88534*0.52917721067)/(ZA**(0.23)+ZB**(0.23)) #bohr radius in Angs
    dzblsum=0.0
    for i in range(len(b)):
        dzblsum=dzblsum-(b[i]*c[i]/rdash)*math.exp(-c[i]*r/rdash)
    dzbl= ((ZA*ZB)*1389354.835/r)*dzblsum #Found in dlpoly constants r4pie0*10 to make Angs J/mol
    dzbl= dzbl*1.04E-5 #convert to eV
    dzbl = (-zbl(ZA,ZB,r)/r) + dzbl
    return dzbl

#Switching function df/dr
def df(r,rm,eta):
    if r<rm:
        return (-math.exp((r-rm)/eta)/(2*eta))
    if r>=rm:
        return (-math.exp((rm-r)/eta)/(2*eta))

#dUtot/dr
def dU(rvalues,Z1,Z2,A,rho,C,rm,eta):
    dUtot=[]
    for r in rvalues:
        dUtot.append( df(r,rm,eta)*zbl(Z1,Z2,r)+switchf(r,rm,eta)*dzbl(Z1,Z2,r)-df(r,rm,eta)*buck(A,rho,C,r)+(1-switchf(r,rm,eta))*dBuck(A,rho,C,r) )
    return dUtot

#Calculate array of ZBL values from input r list
def zblcalc(rvalues,Z1,Z2):
  zblvalues=[]
  for r in rvalues:
    zblvalues.append(zbl(Z1,Z2,r))
  return zblvalues

#Calculate array of Buckhingham values from input r list
def buckcalc(rvalues,A,rho,C):
  buckvalues=[]
  for r in rvalues:
    buckvalues.append(buck(A,rho,C,r))
  return buckvalues

#Calculate convolution
def mixcalc(rvalues,Z1,Z2,A,rho,C,rm,eta):
  mixvalues=[]
  for r in rvalues:
    mixvalues.append(switchf(r,rm,eta)*zbl(Z1,Z2,r)+(1-switchf(r,rm,eta))*buck(A,rho,C,r))
  return mixvalues

#Calculate rough differential of any curve
def diffcalc(potvalues,rvalues):
  diffpot=np.diff(potvalues)
  diffr=np.diff(rvalues)
  dpotdr=diffpot/diffr
  newr=[]
  for i in range(len(rvalues)-1):
    r=(rvalues[i+1]+rvalues[i])/2
    newr.append(r)
  return [dpotdr,newr]
"""    
rvalues=np.linspace(1.5,1.7,400)

zblvalues=[]
for r in rvalues:
    zblvalues.append(zbl(14,8,r))

#plt.plot(rvalues,zblvalues)

rbuckvalues=np.linspace(1.5,1.7,400)
buckvalues=[]
for r in rbuckvalues:
    buckvalues.append(buck(50306.1,46.2978,r,0.161)) #Si-O Wang potential

#plt.plot(rbuckvalues,buckvalues)

rm=1.05
eta=0.09

print(buckvalues)

mixedpot=[]
for r in rvalues:
    mixedpot.append(switchf(r,rm,eta)*zbl(14,8,r)+(1-switchf(r,rm,eta))*buck(50306.1,46.2978,r,0.161))
diffpot=[]
for i in range(len(rvalues)):
    diffpot.append(mixedpot[i]-buckvalues[i])

#gradmix=np.diff(mixedpot,1)
plt.plot(rvalues,diffpot)
#plt.plot(rvalues[1:],gradmix)
#plt.plot(rvalues,mixedpot)

plt.show()
"""
