# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yambopy.netcdf import *
from yambopy.plot import *
from math import sqrt

I = complex(0,1)

class DipolesDB(YamboSaveDB):
    def __init__(self,save='SAVE',filename='ndb.dip_iR_and_P'):
        """ Open all the dipoles databases
        """
        YamboSaveDB.__init__(self,save=save)
        self.filename = filename

        #read dipoles
        self.database = Dataset(filename, 'r')
        self.nk_ibz, self.nk_bz, _, _ = self.database.variables['HEAD_R_LATT']
        self.spin,_ = self.database.variables['SPIN_VARS'][:]
        spin, self.nbands, self.nbandsv = self.database.variables['PARS'][:3]
        self.nbandsc = self.nbands - self.nbandsv

    def get_databases(self):
        self.databases = [Dataset("%s_fragment_%d"%(self.filename,nk+1)) for nk in range(self.nk_ibz)]

        self.dipoles = np.zeros([self.nk_ibz,3,self.spin,self.nbandsc,self.nbandsv],dtype=complex)
        for nk,db in enumerate(self.databases):
            for i in xrange(3):
                for s in xrange(self.spin):
                    #dip = db.variables['DIP_iR_k_%04d_xyz_%04d_spin_%04d'%(nk+1,i+1,s+1)][:]
                    dip = db.variables['DIP_P_k_%04d_xyz_%04d_spin_%04d'%(nk+1,i+1,s+1)][:].T
                    self.dipoles[nk,i,s] = dip[:,:,0]+I*dip[:,:,1]
        return self.dipoles

    def plot(self,size=20,dir=0,spin=0,bandc=1,bandv=None,expand=True):
        """ Plot the weights in a scatter plot of this exciton
        """
        cmap = plt.get_cmap("gist_heat_r")

        #check if bandc is an integer
        if bandc is int: bandc = (bandc,)
        #check if bandv in an integer
        if bandv is None: bandv = (self.nbandsv,)

        dipoles = self.get_databases()

        print "tansitions %d -> %d"%(bandv,bandc+self.nbandsv)
        weights = np.zeros([self.nk_bz])
        for c,v in bandc,bandv:
            weights += dipoles[:,dir,spin,c-1,v-1]
        weights = np.array([abs(x) for x in weights])
        weights = weights/max(weights)
        kpts = self.kpts_car

        if expand:
            kpts, nks = self.expand_kpts()
            weights = weights[nks]
        else:
            kpts = self.kpts_car

        fig = plt.figure(figsize=(10,10))
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
        plt.scatter(kpts[:,0], kpts[:,1], s=size, marker='H', color=[cmap(sqrt(c)) for c in weights])
        plt.axes().set_aspect('equal', 'datalim')
        plt.show()

    def __str__(self):
        s = ""
        s += "nk_ibz    : %d\n"%self.nk_ibz
        s += "nk_bz     : %d\n"%self.nk_bz
        s += "nbandsc   : %d\n"%self.nbandsc
        s += "nbandsv   : %d\n"%self.nbandsv
        return s

if __name__ == "__main__":
    ddb = DipolesDB()
    ddb.get_databases()
    print ddb
