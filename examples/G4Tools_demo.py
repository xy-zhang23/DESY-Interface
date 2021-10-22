# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 15:58:57 2021

@author: lixiangk
"""

from interface import *

workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2021\Genesis-demo'
os.chdir(workdir)

sliceLength = 300e-6
sliceSpacing = sliceLength
mpart = 4096*2
mslice = 120 # slices ahead of the bunch, only to save the radiation fields

# Load the astra distribution
fname = 'run2.0451.002'
dist = pd_loadtxt(fname)
dist[0,2] = 0
dist[1:,5] += dist[0,5]

# Convert from ev/C to unitless momentum for px and py and 
# from ev/C to gamma for pz
dist[:,3:6] /= (g_mec2*1e6)
gamma = np.sqrt(1+dist[:,3]**2+dist[:,4]**2+dist[:,5]**2)
dist[:,5] = gamma
gamma0 = np.mean(gamma)

# Find the range of distribution
zmin = np.min(dist[:,2])
dist[:,2] -= zmin; zmin = 0 # set bunch tail to 0
zmax = np.max(dist[:,2])
nslice = int((zmax-zmin)/sliceLength) # slices for the bunch itself

# Create a hdf5 file for writing
fout = h5py.File('test.par.h5', 'w')
# General setting for the slices
keys = ['beamletsize', 'one4one', 'refposition',
        'slicecount', 'slicelength', 'slicespacing']
_ = fout.create_dataset('beamletsize',
                        data = np.array([4], dtype = int))
_ = fout.create_dataset('one4one', 
                        data = np.array([0], dtype = int))
_ = fout.create_dataset('refposition', 
                        data = np.array([0], dtype = float))

islice = 0
nmin = 1
ztemp = dist[:,2]
stat = []
for i in np.arange(nslice):
    # select particles within a slice
    select = (ztemp>=i*sliceLength) * (ztemp<(i+1)*sliceLength)
    npar = np.sum(select)
    if npar<nmin: # too small partile number for resampling
        continue
    else:
        inputSliceDist = dist[select,0:6]
        
        inputSliceDist[:,2] -= i*sliceLength
        inputSliceDist[:,2] *= (2*np.pi/sliceLength)
        
        # Calculate charge and current of the slice
        charge = np.sum(dist[select,-3])
        current = -charge*1e-9/sliceLength*g_c
        
        islice += 1; print(islice, '/', nslice+mslice, ':', current, ' A')
        
        sliceDist = resampleParticles(inputSliceDist, mpart)
        #xx, yy, theta, px, py, gamma = sliceDist.T[:]
        theta = sliceDist.T[2]
        bf0 = np.mean(np.exp(-theta[:]*1j))
        bf2, bf2_mean = np.abs(bf0)**2, 1.0/(-charge*1e-9/g_qe)
        
        stat.append([islice, charge, current, bf2, bf2_mean])
        #print('bunching factor: ', 1/bf2, 1/bf2_mean, bf2/bf2_mean)
        
        # add to hdf5
        sliceName = str.format('slice%06d/' % islice)
        
        _ = fout.create_dataset(sliceName+'current',
                                data = np.array([current], dtype = float))
        
        coors = ['x', 'y', 'theta', 'px', 'py', 'gamma']
        for i, coor in enumerate(coors):
            xx = sliceDist.T[i]
            _ = fout.create_dataset(sliceName+coor,
                                    data = np.array(xx, dtype = float))

# add more slices ahead of the bunch
current = 0
sliceDist = np.zeros((mpart, 6))
sliceDist[:,2] = np.random.rand(mpart)*np.pi*2
sliceDist[:,5] = gamma0

for i in np.arange(mslice):
    
    # add to hdf5
    islice += 1; print(islice, '/', nslice+mslice, ':', current, ' A') 
    sliceName = str.format('slice%06d/' % islice)
    
    stat.append([islice, 0, 0, 0, 0])
    
    _ = fout.create_dataset(sliceName+'current',
                            data = np.array([current], dtype = float))
    
    coors = ['x', 'y', 'theta', 'px', 'py', 'gamma']
    for i, coor in enumerate(coors):
        xx = sliceDist.T[i]
        _ = fout.create_dataset(sliceName+coor,
                                data = np.array(xx, dtype = float))
        
_ = fout.create_dataset('slicecount', 
                        data = np.array([islice], dtype = int))
_ = fout.create_dataset('slicelength', 
                        data = np.array([sliceLength], dtype = float))
_ = fout.create_dataset('slicespacing', 
                        data = np.array([sliceSpacing], dtype = float))

fout.close()
stat = np.array(stat)

# Plot
fig, ax = plt.subplots(nrows = 2, sharex = True, figsize = (6, 4))
ax[0].plot(stat[:,0]*sliceLength*1e3, stat[:,2], '-')
ax[0].set_ylabel(r'Current (A)')
ax[0].legend(['Resampling'])

ax[1].plot(stat[:,0]*sliceLength*1e3, stat[:,3], '-')
ax[1].plot(stat[:,0]*sliceLength*1e3, [1.0/mpart for _ in np.arange(len(stat))], '-')
ax[1].plot(stat[:,0]*sliceLength*1e3, stat[:,4], '-')

ax[1].set_yscale('log')
#ax[1].set_ylim(1e-9, 1)

ax[1].set_xlabel(r'# of slices')
ax[1].set_xlabel(r'$z$ (mm)')
ax[1].set_ylabel(r'Bunching factor')
ax[1].legend(['Resampling', '1/mpart', '1/ne'], ncol = 2)

for axis in ax:
    axis.grid()

#%% Calculate the bunching factors from Genesis1.3 V4 particle output
import h5py
import numpy.random as rd

f = h5py.File('test1.0.par.h5', 'r')

# General setting for the slices
keys = ['beamletsize', 'one4one', 'refposition',
        'slicecount', 'slicelength', 'slicespacing']
for key in keys:
    a = f.get(key)
    print(a[:])

sliceLength = f.get('slicelength')[0]

stat = []
for i in np.arange(nslice):
    name = 'slice%6.6d' % (i+1)
    
    x = f[name+'/x'][:]
    y = f[name+'/y'][:]
    theta = f[name+'/theta'][:]
    
    px = f[name+'/px'][:]
    py = f[name+'/py'][:]
    gamma = f[name+'/gamma'][:]
    
    current = f[name+'/current'][0]
    charge = current*sliceLength/g_c
    
    if charge>0:
        bf0 = np.mean(np.exp(-theta*1j))
        bf2, bf2_mean = np.abs(bf0)**2, 1.0/(charge/g_qe)
    else:
        bf2, bf2_mean = 0, 0
    
    stat.append([i, charge, current, bf2, bf2_mean])
        
f.close()
stat = np.array(stat)
