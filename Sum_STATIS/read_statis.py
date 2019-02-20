#Python script to read STATIS files
#Oliver Dicks 20/08/2018

import sys
import numpy as np
import copy

#Class that reads STATIS files and saves all data columns as a function of time and timestep
class Read_statis:
    def __init__(self, input_file):
        #Open STATIS file
        self.statis_file = open(input_file, 'r')

        #extract lines of CONFIG file
        self.lines = self.statis_file.readlines()

        #strip linebreak from strings and split in lines
        self.lines = [line.rstrip('\n') for line in self.lines]
        self.lines = [line.split() for line in self.lines]

        #Find number of entries per step
        nentry=int(self.lines[2][2])

        newstep= True
        self.data=[]
        datcnt=0

        #Extract data, first the step and time data, then the number of entries

        for line in self.lines[2:]:
            if newstep==True:
                stepdata=[]
                stepdata.append(int(line[0]))
                stepdata.append(float(line[1]))
                newstep=False
            else:
                for value in line:
                    if datcnt<nentry:
                        stepdata.append(float(value))
                        datcnt=datcnt+1
                        if datcnt==nentry:
                            datcnt=0
                            self.data.append(stepdata)
                            newstep=True

        self.statis_file.close()
        
    def write_columns(self,xcol,ycols):
        outfile=open("statis_sum.dat",'w')
        for i in range(len(self.data)):
            t=()
            t=t + (self.data[i][xcol],)
            for col in ycols:
                t=t+(self.data[i][col],)
            outfile.write('%10f  '*len(t) % t)
            outfile.write('\n')

        outfile.close()

#Reads in system arguments in form x column, then y columns as integers
        
inlist=sys.argv
for i in range(1,len(inlist)):
    inlist[i]=int(inlist[i])

stat=Read_statis("STATIS")
stat.write_columns(inlist[1],inlist[2:])


    
