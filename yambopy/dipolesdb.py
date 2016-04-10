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

class YamboDipolesDB(YamboSaveDB):
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
        self.dipoles = None

    def get_databases(self):
        self.dipoles = np.zeros([self.nk_ibz,3,self.spin,self.nbandsc,self.nbandsv],dtype=complex)
        for nk in range(self.nk_ibz):

            #open database
            db = Dataset("%s_fragment_%d"%(self.filename,nk+1))

            for i in xrange(3):
                for s in xrange(self.spin):
                    #dip = db.variables['DIP_iR_k_%04d_xyz_%04d_spin_%04d'%(nk+1,i+1,s+1)][:]
                    dip = db.variables['DIP_P_k_%04d_xyz_%04d_spin_%04d'%(nk+1,i+1,s+1)][:].T
                    self.dipoles[nk,i,s] = dip[:,:,0]+I*dip[:,:,1]

            #close database
            db.close()

        return self.dipoles

    def get_dipoles(self,dir=0,expand=True):
        #read databases
        if self.dipoles is None:
            self.get_databases()

        #check if we need to expand the dipoles to the full BZ
        if expand:
            field_dir = [[1,0,0],[0,1,0],[0,0,1]][dir]
            kpts, nks, nss = self.expand_kpts(repx=range(1),repy=range(1),repz=range(1))

            #rotate filed according to the different symmetries
            field_dir_rot = [np.dot(s,field_dir) for s in self.sym_car]

            #generate the new dipoles
            dipoles = np.zeros([len(nks),3,self.spin,self.nbandsc,self.nbandsv],dtype=complex)
            for nk_full,nk,ns in zip(xrange(len(nks)),nks,nss):
                for s,c,v in product(xrange(self.spin),xrange(self.nbandsc),xrange(self.nbandsv)):
                    #rotate dipoles
                    dipoles[nk_full,:,s,c,v] = np.dot(self.sym_car[ns],self.dipoles[nk,:,s,c,v])
        else:
            kpts = self.kpts_car
            dipoles = self.dipoles

        return dipoles, kpts

    def plot(self,size=20,dir=0,spin=0,bandc=1,bandv=1,expand=True):
        """ Plot the weights in a scatter plot of this exciton
        """
        cmap = plt.get_cmap("gist_heat_r")

        dipoles, kpts = self.get_dipoles(expand=expand)

        #check if bandc is an integer
        if type(bandc) is int: bandc = (bandc,)
        #check if bandv in an integer
        if type(bandv) is int: bandv = (bandv,)

        print "tansitions %s -> %s"%(str(bandv),str(bandc))
        #weights = np.zeros([self.nk_bz])
        weights = np.zeros([len(kpts)])
        for c,v in product(bandc,bandv):
            weights += np.absolute(dipoles[:,dir,spin,c-1,v-1])
        weights = np.array([x.real for x in weights])
        weights = weights

        fig = plt.figure(figsize=(10,10))
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
        plt.scatter(kpts[:,0], kpts[:,1], s=size, marker='H', lw=0, cmap=cmap, c=weights)
        plt.axes().set_aspect('equal', 'datalim')
        plt.colorbar()
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
