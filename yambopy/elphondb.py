# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from netCDF4 import Dataset
from math import sqrt
import numpy as np

I = complex(0,1)

class ElectronPhononDB():
    """ Python class to read the electron-phonon matrix elements generated with yambo
    """
    def __init__(self,filename='ndb.elph_gkkp_expanded'):
        self.filename = filename

        self.database = Dataset(filename)

        kgrid = self.database.variables['PH_K'][:].T
        qgrid = self.database.variables['PH_Q'][:].T

        self.nk,_ = kgrid.shape
        self.nq,_ = qgrid.shape

        database = Dataset('%s_fragment_1'%filename)
        
        nq,_,self.nmodes,nbands2 = database.variables['ELPH_GKKP_Q1'][:].shape
        nmodes = database.variables['PH_FREQS1'][:].shape[0]
        self.nbands = sqrt(nbands2)
        database.close()       
 
        #some checks
        if nmodes != self.nmodes: print 'variable nmodes has different dimensions %d != %d'%(nmodes,self.nmodes)
        if nq     != self.nq:     print 'variable nmodes has different dimensions'

    def get_databases(self):
        """ Load all the gkkp databases to memory
        """
        self.databases = [ Dataset('%s_fragment_%d'%(self.filename,nq+1)) for nq in range(self.nq) ]
        
        gkkp = np.array([ db['ELPH_GKKP_Q%d'%(n+1)] for n,db in enumerate(self.databases) ])
        self.gkkp = gkkp[:,:,0,:,:] + I*gkkp[:,:,1,:,:]

    def __str__(self):
        s = ''
        s+= 'nkpoints: %d\n'%self.nk
        s+= 'nqpoints: %d\n'%self.nq
        s+= 'nmodes: %d\n'%self.nmodes
        s+= 'nbands: %d\n'%self.nbands
        return s

if __name__ == '__main__':
    elph = ElectronPhononDB()
    print elph
    elph.get_databases()
