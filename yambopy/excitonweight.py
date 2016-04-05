# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
import numpy as np
from numpy.linalg import norm, inv
from netCDF4 import Dataset
from math import sqrt
#we try to use matplotlib, if not present we won't use it
try:
    from matplotlib import pyplot as plt
    from matplotlib.mlab import griddata
except ImportError:
    _has_matplotlib = False
else:
    _has_matplotlib = True
from itertools import product

class YamboExcitonWeight():
    def __init__(self,filename,save='SAVE'):
        #read excitons file
        self.excitons = np.loadtxt(filename)

        #read database
        self.nc_db    = Dataset('%s/ns.db1'%save)
        self.sym_car  = self.nc_db.variables['SYMMETRY'][:]
        self.kpts_iku = self.nc_db.variables['K-POINTS'][:].T
        self.lat      = self.nc_db.variables['LATTICE_VECTORS'][:].T
        self.alat     = self.nc_db.variables['LATTICE_PARAMETER'][:].T

        #caclulate the reciprocal lattice
        self.rlat  = rec_lat(self.lat)
        self.nsym  = len(self.sym_car)

        #convert form internal yambo units to cartesian lattice units
        self.kpts_car = np.array([ k/self.alat for k in self.kpts_iku ])

        #convert cartesian transformations to reduced transformations
        self.sym_rlu = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            a = np.dot(s.T,inv(self.rlat))
            self.sym_rlu[n] = np.dot(inv(self.lat.T),a)

        #convert cartesian transformations to reciprocal transformations
        self.sym_rec = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            self.sym_rec[n] = inv(s).T

        self.weights = None

    def write_irr(self,filename="irr.dat"):
        """ write the list of kpoints from the irreducible brillouin zone
        """
        f = open("irr.dat",'w')
        for k in self.kpts_car:
            f.write(("%12.8lf"*3)%(k[1],k[0],k[2])+"\n")
        f.close()

    def write_full(self,filename="full.dat"):
        """ write the list of kpoints in the full brillouin zone
        """
        #generate all the possible points
        kpoints = self.kpts_car
        kpts = []
        for k in self.kpts_car:
            for sym in self.sym_car:
                kpts.append(np.dot(sym,k))

        f = open("full.dat",'w')
        for q in kpts:
            f.write(("%12.8lf "*3)%tuple(q)+"\n")
        f.close()

    def calc_kpts_weights(self):
        """ calculate the weights and kpoints of the excitons
        """
        self.weights = dict()

        #first run set everything to zero
        for line in self.excitons:
            v,c,k,sym,w,e = line
            self.weights[(int(k),int(sym))] = 0

        #add weights
        for line in self.excitons:
            v,c,k,sym,w,e = line
            self.weights[(int(k),int(sym))] += w

        #rename symmetries and kpoints
        sym = self.sym_car
        kpoints = self.kpts_car

        qpts = []
        weights = []
        for r in product(xrange(3),repeat=3):
            for k,s in self.weights.keys():
                w   = self.weights[(k,s)]
                weights.append( w )
                qpt = np.dot(sym[s-1],kpoints[k-1])+red_car([r],self.rlat)[0]
                qpts.append( qpt )
        return np.array(qpts), np.array(weights)

    def plot_contour(self,resX=500,resY=500):
        """ plot a contour
            idea taken from http://stackoverflow.com/questions/18764814/make-contour-of-scatter
        """
        kpts, z = self.calc_kpts_weights()
        x,y = kpts[:,0],kpts[:,1]
        xi = np.linspace(min(x), max(x), resX)
        yi = np.linspace(min(y), max(y), resY)
        Z = griddata(x, y, z, xi, yi, interp='cubic')
        X, Y = np.meshgrid(xi, yi)

        plt.contourf(X, Y, Z, cmap='gist_heat_r')
        plt.show()

    def plot_weights(self):
        """ Plot the weights in a scatter plot of this exciton
        """
        cmap = plt.get_cmap("gist_heat_r")

        kpts, weights = self.calc_kpts_weights()
        plt.scatter(kpts[:,0], kpts[:,1], s=50, marker='H', color=[cmap(sqrt(c)) for c in weights])
        plt.axes().set_aspect('equal', 'datalim')
        plt.show()

    def __str__(self):
        s = ""
        s += "reciprocal lattice:\n"
        s += "\n".join([("%12.8lf "*3)%tuple(r) for r in self.rlat])+"\n"
        s += "lattice:\n"
        s += "\n".join([("%12.8lf "*3)%tuple(r) for r in self.lat])+"\n"
        s += "alat:\n"
        s += ("%12.8lf "*3)%tuple(self.alat)+"\n"
        return s

if __name__ == "__main__":
    ye = YamboExciton('o-yambo.exc_weights_at_1_02')
    print ye
    ye.write_irr()
    ye.write_full()
    #ye.plot_contour()
    ye.plot_weights()
