#
# Author: Alejandro Molina-Sanchez
# Run real-time simulations with yambo
# We put this script in the same folder than the SAVE directory 
#
#
#from __future__ import print_function
from yambopy.inputfile import *
from pwpy.inputfile import *
from pwpy.outputxml import *

yambo    = 'yambo'
yambo_rt = 'yambo_rt'
ypp_rt   = 'ypp_rt'

'''
# check if the database is present
if not os.path.isdir('database'):
    os.mkdir('database')

#check if the nscf cycle is present
if os.path.isdir('nscf/si.save'):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit() 

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    log_p2y   = shell.run(['p2y'],  cwd="nscf/si.save")
    print(log_p2y.output,file=open("p2y.log","w"))
    log_yambo = shell.run(['yambo'],cwd="nscf/si.save")
    print(log_yambo.output,file=open("yambo.log","w"))
    shell.run('mv nscf/si.save/SAVE database'.split())
'''

def realtime(calculation,name_folder):
  if calculation == 'collision':
    return YamboIn('%s -r -e -v c -V all'%yambo_rt)    
  if calculation == 'tdsex':
    return YamboIn('%s -q p -v c -V all'%yambo_rt)    
  if calculation == 'pumping':
    return YamboIn('%s -q p -v c -V all'%yambo_rt)    
  if calculation == 'dissipation':
    return YamboIn('%s -p c -g n -V all'%yambo_rt)    

calculation = 'collision'
#calculation = 'tdsex'
#calculation = 'pumping'
#calculation = 'diss-elph'
#calculation = 'diss-elel'
#calculation = 'diss-both'

y = realtime('collision','.')

# System Common variables
y['FFTGvecs']  = [15,'Ha']
y['EXXRLvcs']  = [15,'Ha']
y['SCBands']   = [[25,28],'']

# Collision
y['NGsBlkXp']  = [1200,'mHa']
y['BndsRnXs' ] = [[1,70],'']
y['CUTGeo']    = 'box z'
y['CUTBox']    = [[0,0,38], '']
y.write('03_COLLISION')
os.system('rm yambo.in')

# Time-dependent
y = realtime('tdsex','.')
# System
y['GfnQP_Wv']  = [ [0.04,0.00,0.00],'' ]
y['GfnQP_Wc']  = [ [0.04,0.00,0.00],'' ]
y['GfnQPdb']   = 'none' 
# Time-propagation 
y['RTstep']     = [0.01,'fs']
y['NETime']     = [  80,'fs']
y['Integrator'] = "HEUN EXP"
y['IOtime']     = [ [0.05, 0.05, 0.01], 'fs' ] 
 # Pump Pulse
y['Probe_Int']  = [ 1E5 , 'kWLm2']
y['Probe_kind'] = "DELTA"
y['Probe_Damp'] = [0,'fs']
y['Probe_Freq'] = [[0.0,0.0],'eV']
y.write('04_TDSEX')
os.system('rm yambo.in')


exit()


if calculation == 'collision':
  print('running yambo-collision')
  #log_yambo = shell.run(('%s -F 03_COLLISION -J '%yambo,'COLLISION').split(),cwd='gw')
  print ('%s -F yambo_run.in -J '%yambo + ' COLLISION'  )
