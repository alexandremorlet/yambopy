# Copyright (c) 2015, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from ase import Atoms
#we try to use matplotlib, if not present we won't use it
try:
    from matplotlib import pyplot as plt
except ImportError:
    _has_matplotlib = False
else:
    _has_matplotlib = True
import json
import numpy as np
import os

class YamboBSEAbsorptionSpectra():
    def __init__(self,job_string,threshold=0.2):
        """ From a database of excitons
        """
        self.job_string = job_string
        self.threshold = threshold
        self.data = {"excitons":[]}
        self.atoms = None

        #use YamboOut to read the absorption spectra
        y = YamboOut('.')
        print y
        for key,value in y.data.items():
            if "eps" in key: self.data["eps"] = value.tolist()
            if "eel" in key: self.data["eel"] = value.tolist()

    def get_excitons(self):
        """ Obtain the excitons using ypp
        """
        filename = "o-%s.exc_I_sorted"%self.job_string
        if not os.path.isfile(filename):
            os.system("ypp -e s -J %s"%self.job_string)
        self.excitons = np.loadtxt(filename)
        return self.excitons[self.excitons[:,1]>self.threshold]

    def get_wavefunctions(self):
        """ Collect all the wavefuncitons with an intensity larger than self.threshold
        """
        self.filtered_excitons = self.excitons[self.excitons[:,1]>self.threshold]
        self.filtered_excitons

        #read the ypp file using YamboIn
        ypp = YamboIn()
        ypp.read_file("ypp.in")

        keywords = ["lattice", "atoms", "atypes", "nx", "ny", "nz"]
        for exciton in self.filtered_excitons:
            #get info
            e,intensity,i = exciton

            #create ypp input file and run
            ypp["States"] = "%d - %d"%(i,i)
            ypp.write("ypp_%d.in"%i)
            filename = "o-%s.exc_3d_%d.xsf"%(self.job_string,i)
            print filename
            if not os.path.isfile(filename):
                os.system("ypp -F ypp_%d.in -J %s"%(i,self.job_string))

            #read the excitonic wavefunction
            ew = YamboExcitonWaveFunction()
            ew.read_file(filename)
            data = ew.get_data()
            for word in keywords:
                self.data[word] = data[word]
            self.data["excitons"].append({"energy": e,
                                          "intensity": intensity,
                                          "hole": data["hole"],
                                          "index": i,
                                          "datagrid": np.array(data["datagrid"])})

    def get_atoms(self):
        """ Get a ase atoms class
        """
        if "lattice" in self.data.keys():
            self.atypes = self.data["atypes"]
            self.lat = self.data["lattice"]
            self.atom_types = [self.atypes[a[0]] for a in self.data["atoms"]]
            self.pos = [a[1:] for a in self.data["atoms"]]
            self.atoms = Atoms(self.atom_types, self.pos, pbc=[1,1,1])
            self.atoms.set_cell(self.lat)

    def write_json(self):
        """ Write a jsonfile with the absorption spectra and the wavefunctions of certain excitons
        """
        print "writing json file...",
        JsonDumper(self.data,"absorptionspectra.json")
        print "done!"
