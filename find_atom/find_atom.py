#Python script to read CONFIG files and find atoms within R of a coordinate
#Oliver Dicks 26/04/2019

import sys
import numpy as np
import copy

#Transform cartesian coordinates into fractional coordinates
def frac_coord(cartpos,latpar):
    r=np.array(cartpos)
    latpar=np.array(latpar)
    latpar=latpar.transpose()
    rfrac=np.linalg.tensorsolve(latpar,r)
    return rfrac

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

        #Output as fractional coordinates
        self.frac_data=[]
        for i in range(len(self.atom_data)):
            frac_now=[self.atom_data[i][0:2],frac_coord(self.atom_data[i][2],self.lat_par)]
            self.frac_data.append(frac_now)

#Outputs a config file from an input of file name, atom_data and lattice parameters   
class write_config:
    def __init__(self,output_file,atdat,latpar):
        #Open file
        outfile=open(output_file,'w')

        outfile.write("CONFIG OUTPUT FILE\n")
        info = "  0   3  %i  \n" % (len(atdat))        
        outfile.write(info)
        for line in latpar:
            latstring = " %19.9f %19.9f %19.9f \n" % tuple(line)
            outfile.write(latstring)
        for atom in atdat:
            name_string = "%s %15s\n" % tuple(atom[0])
            outfile.write(name_string)
            coord_string = " %16.9f %16.9f %16.9f \n" % tuple(atom[1])
            outfile.write(coord_string)

        outfile.close()
        
inlist=sys.argv
for i in range(1,len(inlist)):
    inlist[i]=float(inlist[i])
R=inlist[1:4]
rcutoff=inlist[4]
rcutsq=rcutoff**2
config=Read_config("CONFIG")
for atom in config.atom_data:
  distsq=(atom[2][0]-R[0])**2+(atom[2][1]-R[1])**2+(atom[2][2]-R[2])**2
  if distsq<rcutsq:  
    print(atom)

