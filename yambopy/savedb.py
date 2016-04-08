# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yambopy.netcdf import *
from yambopy.plot import *
from itertools import product

ha2ev = 27.211396132

def expand_kpts_val(kpts,syms,val):
    """ Take a list of qpoints and symmetry operations and return the full brillouin zone
    with the corresponding index in the irreducible brillouin zone
    """
    full_kpts = []
    full_val  = []
    print "nkpoints:", len(kpts)
    for nk,k in enumerate(kpts):
        for sym in syms:
            full_kpts.append((nk,np.dot(sym,k)))
            full_val.append(val[nk])

    return full_kpts, full_val

def expand_kpts(kpts,syms):
    """ Take a list of qpoints and symmetry operations and return the full brillouin zone
    with the corresponding index in the irreducible brillouin zone
    """
    full_kpts = []
    print "nkpoints:", len(kpts)
    for nk,k in enumerate(kpts):
        for sym in syms:
            full_kpts.append((nk,np.dot(sym,k)))

    return full_kpts

class YamboSaveDB():
    """ Read information from the SAVE database in Yambo
    """
    def __init__(self,save='SAVE'):
        #read database
        try:
            database = '%s/ns.db1'%save
            self.nc_db    = Dataset(database)
        except:
            print "Error reading %d database"%database
            exit()
        self.atomic_numbers   = self.nc_db.variables['atomic_numbers'][:]
        self.atomic_positions = self.nc_db.variables['ATOM_POS'][:]
        self.eigenvalues      = self.nc_db.variables['EIGENVALUES'][:]
        self.sym_car          = self.nc_db.variables['SYMMETRY'][:]
        self.kpts_iku         = self.nc_db.variables['K-POINTS'][:].T
        self.lat              = self.nc_db.variables['LATTICE_VECTORS'][:].T
        self.alat             = self.nc_db.variables['LATTICE_PARAMETER'][:].T

        #caclulate the reciprocal lattice
        self.rlat  = rec_lat(self.lat)
        self.nsym  = len(self.sym_car)

        #convert form internal yambo units to cartesian lattice units
        self.kpts_car = np.array([ k/self.alat for k in self.kpts_iku ])

        #convert cartesian transformations to reduced transformations
        inv = np.linalg.inv
        self.sym_rlu = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            a = np.dot(s.T,inv(self.rlat))
            self.sym_rlu[n] = np.dot(inv(self.lat.T),a)

        #convert cartesian transformations to reciprocal transformations
        self.sym_rec = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            self.sym_rec[n] = inv(s).T

    def expand_kpts(self,repx=range(3),repy=range(3),repz=range(3)):
        """ Take a list of qpoints and symmetry operations and return the full brillouin zone
        with the corresponding index in the irreducible brillouin zone
        """
        kpoints_indexes = []
        kpoints_full    = []

        #expand using symmetries
        for nk,k in enumerate(self.kpts_car):
            for r in product(repx,repy,repz):
                d = red_car([r],self.rlat)[0]
                for sym in self.sym_car:
                    kpoints_full.append(np.dot(sym,k)+d)
                    kpoints_indexes.append(nk)

        print kpoints_full[:5]

        return np.array(kpoints_full), np.array(kpoints_indexes)

    def plot_bs(self,size=20,bandc=1,bandv=None,expand=True,repx=range(3),repy=range(3),repz=range(3)):
        """ Plot the weights in a scatter plot of this exciton
        """
        if bandv is None: bandv = self.nbandsv

        cmap = plt.get_cmap("gist_heat_r")

        eigenvalues = self.eigenvalues
        print "tansitions %d -> %d"%(bandv,bandc+self.nbandsv)
        weights = (eigenvalues[0,:,self.nbandsv+bandc-1]-eigenvalues[0,:,bandv-1])*ha2ev
        print "min:", min(weights)
        print "max:", max(weights)
        weights = weights/max(weights)

        if expand:
            kpts, nks = self.expand_kpts(repx=repx,repy=repy,repz=repz)
            weights = weights[nks]
        else:
            kpts = self.kpts_car

        fig = plt.figure(figsize=(10,10))
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
        #for r in product(range(3),repeat=3):
        #        d = red_car([r],self.rlat)[0]
        plt.scatter(kpts[:,0], kpts[:,1], s=size, marker='H', color=[cmap(c) for c in weights])
        plt.axes().set_aspect('equal', 'datalim')
        plt.show()
