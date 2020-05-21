import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Defines Qi = rho* integral 4*pi*r*G(r)*sin(Qr) dr or integral D(r)sin(Qr) dr from Martin's total scattering formalism
def QiQ(Q,G,rlist,rho):
    dr = rlist[1]-rlist[0]
    integral=0
    for i in range(len(rlist)):
        integral=integral+(dr*rlist[i]*Gr[i]*np.sin(Q*rlist[i]))
    integral=4*np.pi*rho*integral
    return integral

print("Did you remain to change rho (number density N/V)")

grdata=pd.read_csv("Grtot.dat",header=None,delim_whitespace=True)
gr=np.array(grdata.values)
gr=np.transpose(gr)
rlist=gr[0]
Gr=gr[1]

Qlist=np.arange(0.1,45.1,0.05)

QiQlist=[]

rho=4997.0/40.1294647120**3

for Q in Qlist:
    QiQlist.append(QiQ(Q,Gr,rlist,rho))

outfile=open('Qi_tot.dat','w')

for i in range(len(Qlist)):
    info="%16.8f   %16.12f\n" % (Qlist[i],QiQlist[i])
    outfile.write(info)

outfile.close()
