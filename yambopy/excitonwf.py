# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from itertools import product
import matplotlib.pyplot as plt
import json
import numpy as np

def red_car(red,lat): return np.array(map( lambda coord: coord[0]*lat[0]+coord[1]*lat[1]+coord[2]*lat[2], red))

def jump_to(f,tag):
    """ Jump to a line in file
    """
    while True:
        line = f.readline()
        if tag in line:
            break
def v2str(v):
    return ("%12.8lf "*len(v))%tuple(v)

class ExcitonWaveFunction():
    """ Class to read excitonic wavefunctions from yambo in the 3D xsf format
    """
    def __init__(self):
        self.cell = []
        self.atoms = []
        self.lattice = []
        self.initialized = False

    def read_file(self,filename):
        f = open(filename)
        jump_to(f,"PRIMVEC")
        self.lattice.append( map(float,f.readline().strip().split()) )
        self.lattice.append( map(float,f.readline().strip().split()) )
        self.lattice.append( map(float,f.readline().strip().split()) )

        jump_to(f,"PRIMCOORD")
        self.natoms = int(f.readline().split()[0])-1

        #read the hole position
        self.hole = map(float,f.readline().strip().split())

        #read the atoms positions
        self.atoms = []
        for i in range(self.natoms):
            self.atoms.append( map(float,f.readline().strip().split()) )

        #get atypes
        self.atypes = np.unique([a[0] for a in self.atoms]).tolist()
        atypes_dict = dict([(a,n) for n,a in enumerate(self.atypes)])
        self.atoms = [ [atypes_dict[a[0]]]+a[1:] for a in self.atoms]

        jump_to(f,"BEGIN_DATAGRID_3D")
        self.nx, self.ny, self.nz = map(int, f.readline().strip().split())
        f.readline() #ignore

        #read cell
        self.cell.append( map(float,f.readline().strip().split()) )
        self.cell.append( map(float,f.readline().strip().split()) )
        self.cell.append( map(float,f.readline().strip().split()) )

        #read data
        self.datagrid = np.zeros([self.nz,self.ny,self.nx])
        for k,j,i in product(range(self.nz),range(self.ny),range(self.nx)):
            self.datagrid[k,j,i] = float(f.readline())
        self.initialized = True

    def plot_slice_x(self,n):
        """ plot a slice of the 3d grid
        """
        plt.imshow(self.datagrid[:,:,n])
        plt.show()

    def plot_slice_z(self,n):
        """ plot a slice of the 3d grid
        """
        plt.imshow(self.datagrid[n,:,:])
        plt.show()

    def write_xsf(self,filename):
        f = open(filename,'w')
        #structure
        f.write('CRYSTAL\n')
        f.write('PRIMVEC\n')
        for vlat in self.lattice:
            f.write(v2str(vlat)+'\n')
        f.write('PRIMCOORD\n')
        f.write('%d 1\n'%len(self.atoms))
        f.write(v2str(self.hole)+'\n')
        for atom in self.atoms:
            f.write("%d "%self.atypes[atom[0]]+v2str(atom[1:])+'\n')
        #datagrid
        f.write('BEGIN_BLOCK_DATAGRID_3D\n')
        f.write('excitonwf\n')
        f.write('BEGIN_DATAGRID_3D\n')
        f.write('%d %d %d\n'%(self.nx,self.ny,self.nz))
        f.write('0.00 0.00 0.00\n')
        for vlat in self.lattice:
            f.write(v2str(vlat)+'\n')
        for x in self.datagrid.flatten():
            f.write('%lf\n'%x)
        f.write('END_DATAGRID_3D\n')
        f.write('END_BLOCK_DATAGRID_3D\n')
        f.close()
    
    def get_data(self):
        return { "datagrid": self.datagrid.flatten().tolist(),
                 "lattice": self.lattice,
                 "atoms": self.atoms,
                 "atypes": self.atypes,
                 "hole": hole,
                 "nx": self.nx,
                 "ny": self.ny,
                 "nz": self.nz }

    def write_json(self):
        """ Write as a json file
        """
        f = open("datagrid.json","w")
        json.dump(self.get_data(),f)
        f.close()

    def read_json_file(self,filename):
        f = open(filename,"r")
        data = json.load(f)
        f.close()

        self.read_json(data)        

    def read_json(self,data):
        """ Write as a json file
        """
        self.datagrid = data["datagrid"] 
        self.lattice = data["lattice"]
        self.atoms = data["atoms"]
        self.atypes = data["atypes"]
        self.nx = data["nx"]
        self.ny = data["ny"]
        self.nz = data["nz"]
        self.natoms = len(self.atoms)

        self.initialized = True

    def __str__(self):
        s = ""
        s += "lattice:\n"
        for i in range(3):
            s += ("%12.8lf "*3)%tuple(self.lattice[i])+"\n"
        s += "natoms:\n"
        s += "atoms:\n"
        for i in range(self.natoms):
            s += ("%3d "+"%12.8lf "*3)%tuple(self.atoms[i])+"\n"
        s += "atypes:\n"
        for n,a in enumerate(self.atypes):
            s += "%3d %3d\n"%(n,a)
        s += "nx: %d\n"%self.nx
        s += "ny: %d\n"%self.ny
        s += "nz: %d\n"%self.nz
        return s
