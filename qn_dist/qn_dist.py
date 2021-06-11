import sys
import copy
import numpy as np

class Read_field:
    def __init__(self, input_file):
        #Open FIELD file
        self.field_file = open(input_file, 'r')

        #extract lines of FIELD file
        self.lines = self.field_file.readlines()

        #strip linebreak from strings in lines
        self.lines = [line.rstrip('\n') for line in self.lines]
        self.lines = [line.split() for line in self.lines]

        #identify number of atom types
        for line in self.lines:
            for i in range(len(line)):
                if (line[i] == "molecular") and (line[i+1] == "types"):
                    self.numatomtypes = int(line[i+2])
        #identify element types and number of each element type
        self.atomlist=[]
        for i in range(len(self.lines)):
            if len(self.lines[i]) == 0:
                pass
            else:
                if (self.lines[i][0] == "nummols"):
                    atomprops=[]
                    atomprops.append(self.lines[i+2][0]) #append atom type
                    atomprops.append(int(self.lines[i][1])) #append nummols
                    self.atomlist.append(atomprops)
                
field=Read_field("FIELD")
print(field.numatomtypes)
print(field.atomlist)

totnum = 0

for atom in field.atomlist:
  totnum = totnum + atom[1]

print(totnum)

f=open("ICOORD", "r")

f.readline()
f.readline()

icoord = np.empty((totnum+1,30),dtype=object)
qn_list = np.zeros((totnum+1,1),dtype=int)
qn_dist = np.zeros(6,dtype=int)
for i in range(totnum):
  line=f.readline().split()
  gnum=int(line[0])
  icoord[gnum][0] = line[1]
  icoord[gnum][1] = int(line[2])
  for j in range(int(line[2])):
    icoord[gnum][2+j]=int(line[3+j])
for i in range(1,totnum+1):
  if (icoord[i][0] != 'Si'):
    continue
  qcount = 0  
  coord = icoord[i][1]
  for j in range(coord):
    bnum = icoord[i][2+j]
    if icoord[bnum][0] == 'O':
      for k in range(icoord[bnum][1]):
        cnum = icoord[bnum][2+k]
        if (cnum != i) and (icoord[cnum][0] == 'Si'):
          qcount=qcount+1
  qn_list[i]=qcount
  qn_dist[qcount]=qn_dist[qcount]+1 

print(icoord[1])
print(icoord[1540])
print(icoord[194])

print(qn_list[54])
print(qn_dist)

qnout = open("qn_dist.out", "w")

for i in range(len(qn_dist)):
    info="{:2d}   {:10d}\n".format(i,qn_dist[i])
    qnout.write(info)

##Open rdf_all.dat produced by bash script from RDFDAT
#rdf_file = open("rdf_all.dat",'r')
#lines=rdf_file.readlines()
#lines = [line.split() for line in lines]
#rdf_file.close()
#
#pair_atoms=[]
#curpair=['','']
#cnt=0
#for atom in lines[0]:
#    if cnt==0 and atom!="#":
#        curpair[0]=atom
#        cnt=1
#    elif cnt==1:
#        curpair[1]=atom
#        pair_atoms.append(curpair.copy())
#        cnt=0
#
#Natom = copy.deepcopy(field.atomlist)
##Define total number of atoms
#Ntot=0
#for i in range(len(Natom)):
#    Ntot=Ntot+Natom[i][1]
#
##Define c1c2b1b2 for each pair
#coln=[]
#
#for pair in pair_atoms:
#    for element in Natom:
#        if pair[0]==element[0]:
#            N1=element[1]
#        if pair[1]==element[0]:
#            N2=element[1]
#    if pair[0]==pair[1]:
#        A=1.0
#    if pair[0]!=pair[1]:
#        A=2.0
#    c1c2=A*float(N1*N2)/(float(Ntot)*float(Ntot))
#    for element in listb:
#        if pair[0]==element[0]:
#            b1=element[1]
#        if pair[1]==element[0]:
#            b2=element[1]
#    b1b2= b1*b2
#    coln.append(c1c2*b1b2)
#
##Define Gr as Sum_ij (c_i*c_j*b_i*b_j*(g_ij(r)-1))   
#
#TotGr=[]
#
#for i in range(1,len(lines)):
#    Gr=0.0
#    for j in range(1,len(lines[i])):
#        Gr=Gr+(float(lines[i][j])-1.0)*coln[(j-1)]
#    TotGr.append([float(lines[i][0]),Gr])
#
#outfile=open('Grtot.dat','w')
#
#for rpos in TotGr:
#    info="%16.8f   %16.12f\n" % (rpos[0],rpos[1])
#    outfile.write(info)
#
#outfile.close()
#
##Define normalisation term for Gdash where Gdash-1=Gr/(Sum_ijc_i*b_i)^2 from Keen JAC (2000)
#
#sumbici = 0.0
#
#for atomspec in Natom:
#  for element in listb:
#    if atomspec[0]==element[0]:
#      bici=(float(atomspec[1])/float(Ntot))*element[1]
#  sumbici = sumbici + bici
#  
#normfactor=sumbici**(-2)
#print("(Sum_i bi*ci)^-2 = %16.8f" % (normfactor))
#Grdash =  copy.deepcopy(TotGr)
#for i in range(len(TotGr)):
#  Grdash[i][1]=(TotGr[i][1]*normfactor)+1.0
#
#outfile2=open('Grdash.dat','w')
#
#for rpos in Grdash:
#    info="%16.8f   %16.12f\n" % (float(rpos[0]),rpos[1])
#    outfile2.write(info)
#
#outfile2.close()
#
##grdata=pd.read_csv("Grtot.dat",header=None,delim_whitespace=True)
#gr=np.array(TotGr)
#gr=np.transpose(gr)
#print(gr)
#rlist=gr[0]
#Gr=gr[1]
#
#Qlist=np.arange(0.1,45.1,0.05)
#
#QiQlist=[]
#FQlist=[]
#SQlist=[]
#
#for Q in Qlist:
#    QiQlist.append(QiQ(Q,Gr,rlist,rho))
#    FQlist.append(QiQ(Q,Gr,rlist,rho)/Q)
#    SQlist.append((normfactor*QiQ(Q,Gr,rlist,rho)/Q)+1.0)
#
#Qoutfile=open('Qi_tot.dat','w')
#FQoutfile=open('FQ_tot.dat','w')
#SQoutfile=open('SQ_tot.dat','w')
#
#
#for i in range(len(Qlist)):
#    info="%16.8f   %16.12f\n" % (Qlist[i],QiQlist[i])
#    Qoutfile.write(info)
#
#Qoutfile.close()
#
#for i in range(len(Qlist)):
#    info="%16.8f   %16.12f\n" % (Qlist[i],FQlist[i])
#    FQoutfile.write(info)
#
#FQoutfile.close()
#
#for i in range(len(Qlist)):
#    info="%16.8f   %16.12f\n" % (Qlist[i],SQlist[i])
#    SQoutfile.write(info)
#
#SQoutfile.close()
