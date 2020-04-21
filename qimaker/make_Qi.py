import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Defines Q*i(Q) = integral D(r)*sin(Qr) dr from Martin's total scattering formalism
def F(Q,Dr,rlist):
    dr = rlist[1]-rlist[0]
    integral=0
    for i in range(len(rlist)):
        integral=integral+(dr*Dr[i]*np.sin(Q*rlist[i]))
    return integral

Drinput=pd.read_csv("Drtot.dat",header=None,delim_whitespace=True)
Drdata=np.array(Drinput.values)
Drdata=np.transpose(Drdata)
rlist=Drdata[0]
Dr=Drdata[1]

Qlist=np.arange(0.1,45.1,0.05)

Flist=[]

#rho=20000.0/63.164860269**3

for Q in Qlist:
    Flist.append(F(Q,Dr,rlist))

outfile=open('Qi_tot.dat','w')

for i in range(len(Qlist)):
    info="%16.8f   %16.12f\n" % (Qlist[i],Flist[i])
    outfile.write(info)

outfile.close()
