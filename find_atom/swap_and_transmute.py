#!/usr/bin/env python3
import collections
import six
import numpy as np
import argparse as cli
from operator import itemgetter


# python 3.8+ compatibility
try:
    collectionsAbc = collections.abc
except:
    collectionsAbc = collections

tokens = dict(config="CONFIG",outconfig="newconfig")

def update(d, u):
  for k, v in six.iteritems(u):
    dv = d.get(k, {})
    if not isinstance(dv, collectionsAbc.Mapping):
        d[k] = v
    elif isinstance(v, collectionsAbc.Mapping):
        d[k] = update(dv, v)
    else:
        d[k] = v
  return d

def checktokens(x):
    for key, val in x.items():
        if key in tokens:
            assert type(val) == type(tokens[key]), key + " has the wrong type"
        else:
           print("key {:s} is unknown check the manual".format(key))

class dlpoly(object):
  def __init__(self,**d):
    self.d={}
    update(self.d,d)
    checktokens(self.d)

    self.symbols = []
    self.positions = []
    self.velocities = []
    self.forces = []
    self.pbc=0
    self.cell = np.zeros((3, 3))

  def read_config(self):
    try:
        f=open(self.d['config'],'r')
    except IOError:
        return []

    title = f.readline()
    line = f.readline().split()
    levcfg = int(line[0])
    imcon = int(line[1])
    cell = np.zeros((3, 3))
    if imcon>0:
        for j in range(3):
            line = f.readline().split()
            for i in range(3):
                try:
                    cell[j, i] = float(line[i])
                except ValueError:
                    raise RuntimeError("error reading cell")
    symbols = []
    positions = []
    velocities = []
    forces = []
    line = f.readline()
    while line:
        symbol = line.split()[0]
        symbols.append(symbol)
        x, y, z = f.readline().split()[:3]
        positions.append([float(x), float(y), float(z)])
        if levcfg > 0:
            vx, vy, vz = f.readline().split()[:3]
            velocities.append([float(vx), float(vy), float(vz)])
        if levcfg > 1:
            fx, fy, fz = f.readline().split()[:3]
            forces.append([float(fx), float(fy), float(fz)])
        line = f.readline()

    self.symbols=symbols
    self.cell=cell
    self.positions=positions
    self.velocities=velocities
    self.forces=forces
    self.pbc=imcon
    print("read {0:s}".format(title))

  def write_config(self,title="no title",levcfg=0):

    with open(self.d['outconfig'],"w") as f: 
        f.write('{0:72s}\n'.format(title))
        natoms = len(self.symbols)
        f.write('{0:10d}{1:10d}{2:10d}\n'.format(levcfg, self.pbc, natoms))
        if self.pbc > 0:
            cell = self.cell
            for j in range(3):
                f.write('{0:20.10f}{1:20.10f}{2:20.10f}\n'.format(
                    cell[j, 0], cell[j, 1], cell[j, 2]))
                vels = []
        forces = []
        if levcfg > 0:
            vels = self.velocities
        if levcfg > 1:
            forces = self.forces

        labels = self.symbols

        for i, r in enumerate(self.positions):
            f.write("{0:8s}{1:10d}\n{2:20.10f}{3:20.10f}{4:20.10f}\n".format(
                labels[i], i+1, r[0], r[1], r[2]))
            if levcfg > 0:
                if len(vels) == 0:
                    f.write("{0:20.10f}{1:20.10f}{2:20.10f}\n".format(
                        0.0, 0.0, 0.0))
                else:
                    f.write("{0:20.10f}{1:20.10f}{2:20.10f}\n".format(
                        vels[i][0], vels[i][1], vels[i][2]))
            if levcfg > 1:
                if len(forces) == 0:
                    f.write("{0:20.10f}{1:20.10f}{2:20.10f}\n".format(
                                0.0, 0.0, 0.0))
                else:
                    f.write("{0:20.10f}{1:20.10f}{2:20.10f}\n".format(
                        forces[i][0], forces[i][1], forces[i][2]))
  def find_first(self,atom):
      k = -1
      for i,s in enumerate(self.symbols):
        if (s == atom):
          k = i
          break

      return k
      

  def find_atom(self,at,r0,around):
      candidates=[]
      for i, r in enumerate(self.positions):
          if ((r0[0]-around<r[0]<r0[0]+around) and 
              (r0[1]-around<r[1]<r0[1]+around) and
              (r0[2]-around<r[2]<r0[2]+around) and self.symbols[i] == at):
              candidates.append(i)
      if len(candidates) == 1:
          return candidates[0]
      elif len(candidates) == 0:
          raise Exception("no atoms of the type {} indicated found around position {}".format(at,r0))
      else:
          return choose_closest(candidates,self.positions,r0)

  def swap_atoms(self,i,j,newatom):

      r = self.positions[j]
      self.positions[j] = self.positions[i]
      self.positions[i] = r
      
      if len(self.velocities) > 0: 
        v = self.velocities[j]
        self.velocities[j] = self.velocities[i]
        self.velocities[i] = v
      
      if len(self.forces) > 0: 
        f = self.forces[j]
        self.forces[j] = self.forces[i]
        self.forces[i] = f

      self.symbols[i] = newatom

def choose_closest(c,r,r0):

    dists = [ (r0[0]-r[i][0])**2+(r0[1]-r[i][1])**2+(r0[2]-r[i][2])**2 for i in c ]
    return c[min(enumerate(dists), key=itemgetter(1))[0]]

def set_cli():

  parser = cli.ArgumentParser(description='simple code to swap two atoms in a config file')
  parser.add_argument('--config', help='config file to read from (default: %(default)s)',default="CONFIG")
  parser.add_argument('--newconfig', help='config file name to write (default: %(default)s)',default="newconfig")
  parser.add_argument('--position',nargs=3,type = float, help='position around which to check for atom (default: %(default)s)',
        metavar=('x','y','z'),default=[0.0,0.0,0.0])
  parser.add_argument('--atom',help='atom type to transmute (default: %(default)s)',default="Si")
  parser.add_argument('--swap',type=int,help='atom position to swap (default: %(default)s)',default=0)
  parser.add_argument('--newatom',help='atom type to which we transmute (default: %(default)s)',default="U")
  parser.add_argument('--level',type=int, choices=[0,1,2],help='level of detail to write the new config with (default: %(default)s)',default=2)
  parser.add_argument('--first',type=bool, choices=[True,False],help='swap first occurence of --atom to --newatom (default: %(default)s)',default=True)
  parser.add_argument('--around',type=float,help='defines the cube in which to search for atoms of type --atom (default: %(default)s)',default=2.0)

  return parser.parse_args()

if __name__ == '__main__':

    cli = set_cli()

    dlp=dlpoly(config=cli.config,outconfig=cli.newconfig)
    dlp.read_config()
    a = dlp.find_atom(cli.atom,cli.position,cli.around)
    print("Note: indeces start at 0")
    print("Found closest {} at index {}".format(cli.atom,a))
    sw=cli.swap
    if (cli.first):
        sw = dlp.find_first(cli.atom)
    print("Swap atoms {} with {} and transmute {} to {}".format(sw, a,sw,cli.newatom))
    dlp.swap_atoms(sw,a,cli.newatom)
    print("write the config to {}".format(cli.newconfig))
    dlp.write_config(title="new file with U",levcfg=cli.level)
