# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 15:15:22 2023

@author: lixiangk
"""
from interface import *

def astra2gen(inputName = None, inputDist = None, outputName = None):
    '''
    Astra input particle format to hdf5 format used in Genesis1.3 Version 4

    Parameters
    ----------
    inputName : string, optional
        Name of the input particle file. If defined, prefered than inputDist. 
        None by default. 
    inputDist : 6D array, optional
        Particle distribution from Astra, with absolute z and Pz. None by default. 
    outputName: string, optional
        Name of the output particle file.
    Returns
    -------
    None.
    
    Examples
    --------
    ```
    fname = 'beam_modulated.ini'
    fout = 'beam_modulated.h5'
    kwargs = dict(inputName = fname,
                 outputName = fout)
    astra2hdf5(**kwargs)
    ```
    '''
    
    header = r'? VERSION = 1.0\n'+\
    '? COLUMNS X PX Y PY T GAMMA\n'

    if inputName != None:   
        data = pd_loadtxt(inputName)
        data[1:,2] += data[0,2]
        data[1:,5] += data[0,5]
    elif inputDist != None:
        data = inputDist[:,0:6]
    else:
        print('No input file or data!')
        return
        
    if outputName == None:
        outputName = 'temp.dat'
    print('The distribution is saved to '+outputName)
    
    x, y, z = data[:,0:3].T[:]
    t = z/g_c
    
    tminps, tmaxps = t.min()*1e12, t.max()*1e12
    print('Bunch length in ps: ', tmaxps-tminps)
    
    Px, Py, Pz = data[:,3:6].T[:]
    xp = Px/Pz
    yp = Py/Pz
    P = np.sqrt(1+(Px**2+Py**2+Pz**2)/g_mec2/1e6/g_mec2/1e6) # gamma
    
    data[:,0] = x
    data[:,1] = Px/g_mec2/1e6
    data[:,2] = y
    data[:,3] = Py/g_mec2/1e6
    data[:,4] = t
    data[:,5] = P/g_mec2/1e6
    
    np.savetxt(outputName, data[:,0:6], header = header, fmt = '%-15.6e')
    
    # # saving to hdf5
    # f = h5py.File(outputName, 'w')
    
    # dx = f.create_dataset("x", data = x)
    # dy = f.create_dataset("y", data = y)
    # dt = f.create_dataset("t", data = t)
    
    # dxp = f.create_dataset("xp", data = xp)
    # dyp = f.create_dataset("yp", data = yp)
    # dp = f.create_dataset("p", data = p)
    
    # f.close()
    return


#%%

def Halton(i, j):
    '''
    Hamseley sequence reproduced from Mithra2.

    Parameters
    ----------
    i : TYPE
        DESCRIPTION.
    j : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    if (i > 20):
        print(" dimension can not be larger than 20. ")
        return;
    
    prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
    #int p0, p, k, k0, a;
    
    x = 0.0;
    k0 = j + 1;
    p = prime[i];

    p0 = p;
    k  = k0;
    x  = 0.0;
    while (k > 0):
        a = k%p;
        x += a/float(p0);
        k = int(k/p);
        p0 *= p;
      
    return 1.0 - x;
halton = Halton;

# def HaltonNorm(j, mu = 0, sig = 1):
#     r1 = Halton(8, j)
#     r2 = Halton(9, j)   
#     return sig*np.sqrt(-2*np.log(r1))*np.cos(2*np.pi*r2)+mu

def HaltonNorm(j, mu = 0, sig = 1, c_sig = 0):
    r1 = Halton(8, j)
    if c_sig == 0:
        rr = sp.special.erfinv(2*r1-1)*np.sqrt(2)*sig+mu
    else:
        ccc = sp.special.erf(c_sig/np.sqrt(2))
        rr = sp.special.erfinv(2*r1*ccc-ccc)*np.sqrt(2)*sig+mu
    return rr


erf = sp.special.erf
erfinv = sp.special.erfinv

def HaltonNorm2(j, mu = 0, sig = 1, a = -10, b = 10):
    # F(x, a) = (erf(x/np.sqrt(2))-erf(a/np.sqrt(2)))/2/ccc -> 0, 1
    Fa = erf(a/np.sqrt(2))
    Fb = erf(b/np.sqrt(2))
    r1 = Halton(8, j)
    r2 = erfinv(2*r1*(Fb-Fa)/2+Fa)*np.sqrt(2)
    return r2

# rr = np.array([HaltonNorm2(j, a = -3, b = 0) for j in np.arange(10000)])
# fig, ax = plt.subplots(figsize = (5, 4))
# ax.hist(rr, bins = 20, histtype = r'step')
#%%
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2021\Genesis-demo\Minimal'
os.chdir(workdir)

import h5py

fname = 'test.0.par.h5'
f = h5py.File(fname, 'r')

Qb = 2e-9
sigma_z = 1.85e-3
sigma_t = sigma_z/g_c
I = Qb/np.sqrt(2*np.pi)/sigma_t

lambda0 = 100e-6
k0 = 2*np.pi/lambda0

Nm_slice = len(theta)

rr = np.array([HaltonNorm2(j, a = a, b = b) for j in np.arange(Nm_slice)])

fG = lambda x, mu = 0, sig = sigma_z:1/np.sqrt(2*np.pi)/sig*np.exp(-x**2/2/sig**2)
Ns_b = int(sigma_z*3/lambda0)//2*2+2

#fig, ax = plt.subplots()

Ns = f['slicecount'][0]

Nm_slices = 0
for islice in np.arange(0, Ns_b):
    name = str.format('slice%06d' % (islice+2+6))
    
    current = f[name+'/current'][0]
    if current > 0:
        
        x1 = -3*sigma_z+islice*lambda0
        a = (x1-0.5*lambda0)/sigma_z
        b = (x1+0.5*lambda0)/sigma_z
        
        theta = f[name+'/theta'][:]; Nm_slice = len(theta)
        
        rr = np.array([HaltonNorm2(Nm_slice+j, a = a, b = a+1) for j in np.arange(Nm_slice)])
        
        theta1 = 2*np.pi*rr
        current1 = fG(x1)
        
        Nm_slices += Nm_slice
        
        print(islice, current1)
    
        #ax.hist(theta/k0+x1, bins = 200, histtype = r'step')
        #ax.hist(theta1/k0+x1, bins = 200, histtype = r'step')
    
#%%
for 
current = f['slice000010/current']
gamma = f['slice000010/gamma']
theta = f['slice000010/theta']







fig, ax = plt.subplots(figsize = (5, 4))
ax.hist(theta, bins = 20, histtype = r'step')



#%%

workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2022\THzComm2\2nC\SASE'
os.chdir(workdir)

inputName = '368A.2809.003'
astra2hdf5(inputName)