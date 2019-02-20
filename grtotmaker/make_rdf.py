import sys
import copy

listb=[['B', 5.30], ['Na', 3.64], ['Al', 3.449], ['Ca', 4.7], ['Si', 4.1491], ['O', 5.803],['Zr',7.160]]

#Class that reads CONFIG files and outputs the lattice parameters, a raw data array and an array with the fractional coordinates
class Read_config:
    def __init__(self, input_file):
        #Open CONFIG file
        self.config_file = open(input_file, 'r')

        #extract lines of CONFIG file
        self.lines = self.config_file.readlines()

        #strip linebreak from strings in lines
        self.lines = [line.rstrip('\n') for line in self.lines]

        #extract lattice parameters

        self.lat_par=[[0,0,0],[0,0,0],[0,0,0]]
        i=0
        for line in self.lines[2:5]:
                line_latpar = line.split()
                lat_par_temp=line_latpar[0:3]
                lat_par_temp=[float(element) for element in lat_par_temp]
                self.lat_par[i]=lat_par_temp
                i=i+1

        #Read off atom type, number and position and print to list as [[Element,number,x,y,z]...]
        self.atom_data=[]
        self.atom_types=[]
        new_atom=False
        for line in self.lines[5:]:
            line=line.split()
            if (type(line[0]) is str) and (len(line[0]) <= 2):
                exists = False
                for element in self.atom_types:
                    if line[0] == element:
                        exists = True
                if exists == False:
                    self.atom_types.append(line[0])
                this_atom=[]
                this_atom.append(line[0])
                this_atom.append(int(line[1]))
                new_atom=True
            else:
                if new_atom==True:
                    atom_coordinates=line
                    for i in range(len(atom_coordinates)):
                        atom_coordinates[i]=float(atom_coordinates[i])
                    this_atom.append(atom_coordinates)
                    self.atom_data.append(this_atom)
                    new_atom=False


config=Read_config('REVCON')
nelement=[]
for i in range(1,len(config.atom_data)):
    if config.atom_data[i][0] != config.atom_data[i-1][0]:
        nelement.append(config.atom_data[i-1][0:2])
nelement.append(config.atom_data[len(config.atom_data)-1][0:2])

Natom=copy.deepcopy(nelement)

for i in range(1,len(Natom)):
    Natom[i][1]=Natom[i][1]-nelement[i-1][1]

Ntot=0
for x in Natom:
    Ntot=Ntot+x[1]


rdf_file = open("rdf_all.dat",'r')

lines=rdf_file.readlines()
lines = [line.split() for line in lines]

rdf_file.close()


pair_atoms=[]
curpair=['','']
cnt=0
for atom in lines[0]:
    if cnt==0 and atom!="#":
        curpair[0]=atom
        cnt=1
    elif cnt==1:
        curpair[1]=atom
        pair_atoms.append(curpair.copy())
        cnt=0
        
print(pair_atoms)
coln=[]
for pair in pair_atoms:
    for element in Natom:
        if pair[0]==element[0]:
            N1=element[1]
        if pair[1]==element[0]:
            N2=element[1]
    if pair[0]==pair[1]:
        A=1.0
    if pair[0]!=pair[1]:
        A=2.0
    c1c2=A*float(N1*N2)/(float(Ntot)*float(Ntot))
#    print(c1c2)
    for element in listb:
        if pair[0]==element[0]:
            b1=element[1]
        if pair[1]==element[0]:
            b2=element[1]
    b1b2= b1*b2
#    print(b1b2)
    coln.append(c1c2*b1b2)
#    print(coln)

sumcoln=0
for term in coln:
    sumcoln=sumcoln+term

for i in range(len(coln)):
    coln[i]=coln[i]/sumcoln
    
print(coln)
Totgr=[]
print(lines[14])

for i in range(1,len(lines)):
    grtot=0.0
    for j in range(1,len(lines[i])):
        grtot=grtot+float(lines[i][j])*coln[(j-1)]
    Totgr.append([float(lines[i][0]),grtot])

outfile=open('grtot.dat','w')

for bit in Totgr:
    info="%16.8f   %16.12f\n" % (float(bit[0]),bit[1])
    outfile.write(info)

outfile.close()
    


