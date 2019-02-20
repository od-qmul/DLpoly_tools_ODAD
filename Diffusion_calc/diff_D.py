import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

outfile=open("diffconsts.dat","w")
    
tempdir=os.listdir(os.getcwd())
for temp in tempdir:
    if os.path.isdir(temp) == True and RepresentsInt(temp)==True:
        os.chdir(temp)
        data=pd.read_table('statis_sum.dat', sep='\s+', header=None)
        data.columns=['Time' if x==0 else x for x in data.columns]

        x=data['Time']

        diffs=[]

        for j in range(1,len(data.columns)):
            y=data[j]
            x1=[]
            y1=[]
            d1=[]
            for i in range(len(x)):
                if y[i]>4.0 and x[i]>50.0:
                    x1.append(x[i])
                    y1.append(y[i])
#Purely for visual inspection, not used for Dfinal calculation 
                    Dapprox=(y[i]/x[i])/6 # 10^-8 m^2/s
                    d1.append(Dapprox)
    
            plt.plot(x1,y1)

            if len(x1) < 50:
                D = 0.0
            else:
                z1=np.polyfit(x1,y1,1)
                p1=np.poly1d(z1)
                plt.plot(x1,p1(x1),c="r",linewidth=1.0)
                y1hat=p1(x1)
                y1bar=np.sum(y1)/len(y1)
                ssreg1=np.sum((y1hat-y1bar)**2)
                sstot1=np.sum((y1-y1bar)**2)
                rsq1=ssreg1/sstot1
                D=z1[0]/6 # 10^-8 m^2/s
                
            print(j, D)
            diffs.append(D)

        plt.title(int(temp))
        plt.show()

        os.chdir("../")
        outfile.write('%10s   ' % temp)
        outfile.write('%10f'*len(diffs) % tuple(diffs))
        outfile.write("\n")

outfile.close()



