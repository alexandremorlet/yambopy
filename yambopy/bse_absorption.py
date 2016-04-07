# Copyright (c) 2015, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yambopy.plot  import *
import os

class YamboBSEAbsorptionSpectra(YamboSaveDB):
    """ Create a file with information about the excitons from Yambo files
    """
    def __init__(self,job_string,save='SAVE'):
        YamboSaveDB.__init__(self,save=save)
        self.job_string = job_string
        self.data = {"excitons":[],
                     "lattice": self.lat,
                     "atypes": self.atomic_numbers,
                     "atoms": self.atomic_positions}

        self.atoms = None
        self.excitons = None

        #use YamboOut to read the absorption spectra
        y = YamboOut('.')
        # we obtain all the bse spectra
        absorptionspectra = y.get_data(('eps','diago'))
        #we just use one of them
        key = list(absorptionspectra)[0]
        for key,value in absorptionspectra[key].items():
            self.data[key] = value

    def get_excitons(self,min_intensity=0.1,max_energy=4,Degen_Step=0.0100):
        """ Obtain the excitons using ypp
        """
        filename = "o-%s.exc_E_sorted"%self.job_string
        if not os.path.isfile(filename):
            os.system("ypp -e s -J %s"%self.job_string)
        self.excitons = np.loadtxt(filename)

        #filter with degen
        if Degen_Step:
            new_excitons = []
            prev_exc = 0
            for exc in self.excitons:
                e,i,index = exc
                if abs(e-prev_exc)>Degen_Step:
                    new_excitons.append(exc)
                prev_exc = e
            self.excitons = np.array(new_excitons)

        #filter with energy
        self.excitons = self.excitons[self.excitons[:,0]<max_energy]

        #filter with intensity
        self.excitons = self.excitons[self.excitons[:,1]>min_intensity]

        return self.excitons

    def get_wavefunctions(self,FFTGvecs=30,Degen_Step=0.0100,repx=range(3),repy=range(3),repz=range(3),wf=False):
        """ Collect all the wavefuncitons with an intensity larger than self.threshold
        """
        if self.excitons is None:
            print "Excitons not present use get_exitons() first"
            exit()

        #create a ypp file using YamboIn for reading the wavefunction
        yppwf = YamboIn('ypp -e w',filename='ypp.in')
        yppwf['Format'] = "x"
        yppwf['Direction'] = "123"
        yppwf['FFTGvecs'] = FFTGvecs
        yppwf['Degen_Step'] = Degen_Step

        #create a ypp file using YamboIn for reading the excitonic weights
        yppew = YamboIn('ypp -e a',filename='ypp.in')
        yppew['MinWeight'] = 1e-6
        yppew['Degen_Step'] = Degen_Step

        keywords = ["lattice", "atoms", "atypes", "nx", "ny", "nz"]
        for exciton in self.excitons:
            #get info
            e,intensity,i = exciton

            if wf:
                ##############################################################
                # Excitonic Wavefunction
                ##############################################################
                #create ypp input for the wavefunction file and run
                yppwf["States"] = "%d - %d"%(i,i)
                yppwf.write("yppwf_%d.in"%i)

                filename = "o-%s.exc_3d_%d.xsf"%(self.job_string,i)
                if not os.path.isfile(filename):
                    os.system("ypp -F yppwf_%d.in -J %s"%(i,self.job_string))

                #read the excitonic wavefunction
                ewf = YamboExcitonWaveFunction()
                ewf.read_file(filename)
                data = ewf.get_data()
                for word in keywords:
                    self.data[word] = data[word]

            ##############################################################
            # Excitonic Amplitudes
            ##############################################################
            #create ypp input for the amplitudes file and run
            yppew["States"] = "%d - %d"%(i,i)
            yppew.write("yppew_%d.in"%i)

            filename = "o-%s.exc_weights_at_%d"%(self.job_string,i)
            if not os.path.isfile(filename):
                os.system("ypp -F yppew_%d.in -J %s"%(i,self.job_string))

            #read the excitonic wavefunction
            ew = YamboExcitonWeight(filename)
            qpts, weights = ew.calc_kpts_weights(repx=repx,repy=repy,repz=repz)

            ############
            # Save data
            ############
            exciton = {"energy": e,
                       "intensity": intensity,
                       "weights": weights,
                       "qpts": qpts,
                       "index": i}
            if wf:
                exciton["hole"] = data["hole"]
                exciton["datagrid"] = np.array(data["datagrid"])

            self.data["excitons"].append(exciton)

    def write_json(self):
        """ Write a jsonfile with the absorption spectra and the wavefunctions of certain excitons
        """
        print "writing json file...",
        JsonDumper(self.data,"absorptionspectra.json")
        print "done!"
