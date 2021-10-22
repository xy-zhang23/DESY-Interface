# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 17:19:37 2020

@author: lixiangk
"""

from interface import *

P0 = 17 # MeV/c
gamma0 = M2G(P0)

lambda0 = 100e-6
lam_u = 30e-3
delz = lam_u/2

setup = Namelist('setup', 
                 rootname = 'test',
                 lattice = 'gen4lat.lat',
                 beamline = 'tunnel',
                 gamma0 = gamma0,
                 lambda0 = lambda0,
                 delz = delz,
                 seed = 1,
                 npart = 8192,
                 nbins = 4,
                 one4one = False,
                 shotnoise = True)

time = Namelist('time',
                s0 = 0,
                slen = 5e-3,
                sample = 1,
                time = True)

profile_gauss = Namelist('profile_gauss', 
                         label = 'beampro',
                         c0 = 150.0,
                         s0 = 2.5e-3,
                         sig = 1e-3)

profile_step = Namelist('profile_step', 
                         label = 'beampro',
                         c0 = 150.0,
                         s_start = 0e-3,
                         s_end = 5e-3)

beam = Namelist('beam',
                gamma = gamma0,
                delgam = gamma0*0.5e-2,
                current = '@beampro',
                ex = 4e-6,
                ey = 4e-6,
                betax = 19.99072,
                betay = 0.74887,
                alphax = 10.93,
                alphay = 3.25,
                bunch = 0,
                emod = 0)

importbeam = Namelist('importbeam',
                      file = 'test.0.par.h5')

efield = Namelist('efield',
                  rmax = 0,
                  nz = 1,
                  nphi = 2,
                  ngrid = 170)

field = Namelist('field',
                 power = 0,
                 phase = 0,
                 dgrid = 20e-3,
                 ngrid = 201)

track = Namelist('track',
                 output_step = 1,
                 field_dump_step = 0,
                 beam_dump_step = 0)

g4 = Genesis4(setup, time, profile_step, importbeam, field, track)
#print(g4.output)

g4.output = g4.output.lower()
g4.write('test.in')
#g4.qsub('test.sh', 'test.in')

#%% Astra format to hdf5 format used by Genesis 1.3 V4


#%% slice-wise parameters from txt to hdf format
#fname = get_file('.*')
import h5py
def slice2hdf5(inputName = None, ext = '.h5'):
    try:
        import h5py
    except Exception as err:
        print(err)
        
    if inputName is None:
        inputName = get_file('.*')
        
    file = open(inputName, 'r')
    line = file.readline()
    i = 1
    while True:
        if not line:
            break
        match = re.search(r'COLUMNS' , line.upper())
        if match:
            keys = line.split()[2:]
            break
        line = file.readline()
        i += 1
    file.close()
    
    data = np.loadtxt(inputName, skiprows = i)  
    
    path, base, _ = fileparts(inputName)
    outputName = os.path.join(path, base+ext); print(outputName)
    
    #outputName = 'test.in.h5'
    with h5py.File(outputName, 'w') as file:
        for i, key in enumerate(keys):
            dset = file.create_dataset(key, data = data[:,i])
    
    return

slice2hdf5('slice@beam_und.ini.001')        

#file = h5py.File(outputName, 'r')
#plt.plot(file.get('ZPOS')[:], file.get('CURPEAK')[:])
#file.close()