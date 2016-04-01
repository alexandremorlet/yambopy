from netCDF4 import Dataset
from math import sqrt
import numpy as np

I = complex(0,1)

class DipolesDB():
    def __init__(self,filename='ndb.dip_iR_and_P'):
        """ Open all the dipoles databases
        """
        self.filename = filename

        #read dipoles
        self.database = Dataset(filename, 'r')
        self.nk_ibz, self.nk_bz, _, _ = self.database.variables['HEAD_R_LATT']
        self.spin,_ = self.database.variables['SPIN_VARS'][:]
        spin, self.nbands, self.nbandsv = self.database.variables['PARS'][:3]
        self.nbandsc = self.nbands - self.nbandsv

    def get_databases(self):
        self.databases = [Dataset("%s_fragment_%d"%(self.filename,nk+1)) for nk in range(self.nk_ibz)]

        self.dipoles = np.zeros([self.nk_ibz,3,self.spin,self.nbandsv,self.nbandsc],dtype=complex)
        for nk,db in enumerate(self.databases):
            for i in xrange(3):
                for s in xrange(self.spin):
                    dip = db.variables['DIP_iR_k_%04d_xyz_%04d_spin_%04d'%(nk+1,i+1,s+1)][:]
                    self.dipoles[nk,i,s] = dip[0,:,:]+I*dip[1,:,:]
    
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
