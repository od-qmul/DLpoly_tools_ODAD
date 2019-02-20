#Python script to read CONFIG files and resize by percentage or lattice parameter
#Oliver Dicks 05/02/2018

import sys
import numpy as np
import copy

amu=1.660539E-24 #g/a.u

def rand_coord():
    rand_coord=np.random.uniform(0,1.0,size=3)
    return rand_coord

def get_new_coords(fracin,new_latpar):
    frac=copy.deepcopy(fracin)
    final_atom_data = []
    lat=np.array(new_latpar)
    for i in range(len(frac)):
        temp_frac=np.array(frac[i][2])
        temp_newcoord = np.inner(lat,temp_frac)
        temp_line = frac[i]
        temp_line[2]=temp_newcoord
        final_atom_data.append(temp_line)
    return final_atom_data 

 
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
            name_string = "%s %15s\n" % tuple(atom[0:2])
            outfile.write(name_string)
            coord_string = " %16.9f %16.9f %16.9f \n" % tuple(atom[2])
            outfile.write(coord_string)

        outfile.close()

class get_input:
    def __init__(self):
        #Open file inrand
        in_file=open("inrand",'r')
        #extract lines of inrand
        lines = in_file.readlines()

        #strip linebreak from strings in lines
        lines = [line.rstrip('\n') for line in lines]

        dencnt = 0

        #read off density
        for line in lines:
            if len(line)>0:
                aline=line.split()
                if aline[0].lower() == "density" and dencnt==0:
                    self.density = float(aline [1]) #density in gcm-3
                    dencnt=dencnt+1
                elif aline[0].lower() == "density" and dencnt>0:
                    print("Too many densities specified")
                    exit()

        self.atom_dat=[]
        found_atoms = False
        #read off atom data
        for line in lines:
            if len(line)>0:
                bline=line.split()
                if bline[0].lower() == "endatoms":
                    break
                if found_atoms == True:
                    thisatom=bline
                    self.atom_dat.append(bline)
                if bline[0].lower() == "atoms":
                    found_atoms = True

        #define lattice parameters
        mass=0.0
        for atom_type in self.atom_dat:
            mass = mass + float(atom_type[1])*float(atom_type[2])*amu
        volume = mass/self.density
        lat_length = np.cbrt(volume)*1E-2*1E10
        self.lat_par=[[lat_length,0.0,0.0],[0.0,lat_length,0.0],[0.0,0.0,lat_length]]
          

indat=get_input()

print(indat.density)
print(indat.atom_dat)
print(indat.lat_par)

frac_data=[]
j=0

for element in indat.atom_dat:
    for i in range(int(element[2])):
        j=j+1
        this_atom=[]
        this_atom.append(element[0])
        this_atom.append(j)
        this_atom.append(rand_coord())
        frac_data.append(this_atom)

final_data=get_new_coords(frac_data,indat.lat_par) 

#print(final_data[0:10])

write_config("CONFIG",final_data,indat.lat_par)

