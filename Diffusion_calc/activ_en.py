import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

data=pd.read_table('diffconsts.dat', sep='\s+', header=None)
datamat=data.values

origdata=datamat.copy()

outfile=open("Ea_D0_OUT","w")

for i in range(len(datamat)):
    datamat[i][0]=1/float(datamat[i][0])
    for j in range(1,len(datamat[i])):
        datamat[i][j]=np.log(datamat[i][j])

print(datamat)

    
"""
outfile=open("diffconsts.dat","w")
    
tempdir=os.listdir(os.getcwd())
for temp in tempdir:
    if os.path.isdir(temp) == True and RepresentsInt(temp)==True:
        os.chdir(temp)

        data=pd.read_table('statis_sum.dat', sep='\s+', header=None)
        data.columns=['Time' if x==0 else x for x in data.columns]
        #sdata=data.sort_values(by=['Time'])

        print(data.columns.values)

        x=data['Time']

        diffs=[]

        for j in range(1,len(data.columns)):
            y=data[j]

            x1=[]
            y1=[]
            d1=[]

            for i in range(len(x)):
                if x[i]>150.0:
                    x1.append(x[i])
                    D=(y[i]/x[i])/6 # 10^-8 m^2/s
                    d1.append(D)
                    y1.append(y[i])

            #plt.plot(x1,y1)
            plt.plot(x1,d1)

            print(j, np.mean(d1))
            diffs.append(np.mean(d1))

        plt.figure(int(temp))

        plt.show()

        os.chdir("../")
        outfile.write('%10s   ' % temp)
        outfile.write('%10f'*len(diffs) % tuple(diffs))
        outfile.write("\n")

outfile.close()
"""
"""
    z1=np.polyfit(x1,y1,1)
    p1=np.poly1d(z1)
    #plt.plot(x1,p1(x1),c="r",linewidth=1.0)
    y1hat=p1(x1)
    y1bar=np.sum(y1)/len(y1)
    ssreg1=np.sum((y1hat-y1bar)**2)
    sstot1=np.sum((y1-y1bar)**2)

    rsq1=ssreg1/sstot1

    print("Line fit rsq = ",rsq1)
"""
    


