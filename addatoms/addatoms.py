import sys
import numpy as np
import copy


def rand_coord():
    rand_coord=np.random.uniform(-0.5,0.5,size=3)
    return rand_coord


#Transform cartesian coordinates into fractional coordinates
def frac_coord(cartpos,latpar):
    r=np.array(cartpos)
    latpar=np.array(latpar)
    latpar=latpar.transpose()
    rfrac=np.linalg.tensorsolve(latpar,r)
    return rfrac

#class that reads CONFIG files and outputs the lattice parameters, a raw data array and an array with the fractional coordinates
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

# read off atom type, number and position and print to list as [[Element,number,x,y,z]...]
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
            frac_now=[self.atom_data[i][0],self.atom_data[i][1],frac_coord(self.atom_data[i][2],self.lat_par)]
            self.frac_data.append(frac_now)

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
        in_file=open("input",'r')
        #extract lines of inrand
        lines = in_file.readlines()

        #strip linebreak from strings in lines
        lines = [line.rstrip('\n') for line in lines]
        
        self.numatoms=0
        for line in lines:
            if len(line)>0:
                aline=line.split()
                if aline[0].lower() == "num":
                    self.numatoms = int(aline [1])                 
                     
                    
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

          

indat=get_input()

config=Read_config("CONFIG")


#print(config.frac_data)

j=0
frac_data=[]
for element in indat.atom_dat:
    for i in range(int(element[1])):
        j=j+1
        this_atom=[]
        this_atom.append(element[0])
        this_atom.append(j+indat.numatoms)
        this_atom.append(rand_coord())
        config.frac_data.append(this_atom)
        
#print(config.lat_par)
#print(config.frac_data)
final_data=get_new_coords(config.frac_data,config.lat_par)
#print(final_data)
write_config("newCONFIG",final_data,config.lat_par)

