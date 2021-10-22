# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 11:36:54 2021

@author: lixiangk
"""

from .Numeric import *

from scipy import fft
import h5py

#%% Convert 6D particle distributions in ascii format (used in Astra) to hdf5
#   format for Genesis1.3 V4
def astra2hdf5(inputName = None, inputDist = None, outputName = None, ext = '.h5'):
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

    '''
    
    if inputName != None:   
        data = pd_loadtxt(inputName)
        data[1:,2] += data[0,2]
        data[1:,5] += data[0,5]
        if outputName == None:
            path, base, _ = fileparts(inputName)
            outputName = path+os.sep+base+ext 
    elif inputDist != None:
        data = inputDist[:,0:6]
        if outputName == None:
            outputName = 'temp'+ext
    else:
        print('No input file or data!')
        return
    print('The distribution is saved to '+outputName)
    
    x, y, z = data[:,0:3].T[:]
    t = z/g_c
    
    Px, Py, Pz = data[:,3:6].T[:]
    xp = Px/Pz
    yp = Py/Pz
    p = np.sqrt(Px**2+Py**2+Pz**2)/1e6/g_mec2
    
    f = h5py.File(outputName, 'w')
    
    dx = f.create_dataset("x", data = x)
    dy = f.create_dataset("y", data = y)
    dt = f.create_dataset("t", data = t)
    
    dxp = f.create_dataset("xp", data = xp)
    dyp = f.create_dataset("yp", data = yp)
    dp = f.create_dataset("p", data = p)
    
    f.close()
    
#%% Convert slice-wise parameters from ascii (used in Genesis1.3 V2) to hdf5
#   format for V4
def slice2hdf5(inputName = None, outputName = None, ext = '.h5'):
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
    
    if outputName == None:
        path, base, _ = fileparts(inputName)
        outputName = os.path.join(path, base+ext); print(outputName)
    
    #outputName = 'test.in.h5'
    with h5py.File(outputName, 'w') as file:
        for i, key in enumerate(keys):
            dset = file.create_dataset(key, data = data[:,i])
    
    return

#%% Resample 6D particle distribution in such a way that every slice has the 
#   same number of macro particles, as input for Genesis1.3 V4.
#   This is a remake of the Genesis1.3 V4 C++ codes
def resampleParticles(inputBeamDist, mpart = 8192):
    '''
    Resampling to add samples to a known set, rewriting the C++ code in Genesis 1.3
    V4 into Python.
    Add particles to the input distribution `beamdist0` and return a new distribution
    of size `mpart`.

    Parameters
    ----------
    inputBeamDist : 2D array of size (ndist0, 6)
        6D particle distribution, with columns storing coordinates of x, y, theta, 
        px, py, gamma, respectively.
    mpart : int
        Number of particles of the new distribution to return.

    Returns
    -------
    beamDist : 2D array of size (mpart, 6)
        6D particle distribution, with columns storing coordinates of x, y, theta,
        px, py, gamma, respectively.
        
    Examples
    -------
    ```
    mpart = 4096
    beamdist0 = np.random.rand(4000, 6) # x, y, theta, px, py, gamma
    beamdist = resampleParticles(beamdist0, mpart)
    ```
    '''
    
    nlen, ndim = inputBeamDist.shape
    if ndim<6:
        print('Wrong size, at least 6D!')
        return None
    inputBeamDist = inputBeamDist[:,0:6]
    
    if nlen>=mpart:
        # No loop
        while nlen>mpart:
            n1 = np.floor(np.random.rand(nlen-mpart)*nlen).astype(int)
            #import pdb; pdb.set_trace()
            inputBeamDist = np.delete(inputBeamDist, n1, axis = 0)
            nlen = len(inputBeamDist)
        # Use loop as done in Genesis
        # while nlen>mpart:
        #     n1 = int(np.random.rand()*nlen)
        #     inputBeamDist = np.delete(inputBeamDist, n1, axis = 0)
        #     nlen = nlen-1
        beamDist = inputBeamDist
    else:
        beamDist = np.zeros((mpart, 6))
        beamDist[:nlen] = inputBeamDist  
    
        # Step 1 - calculate the center and the rms beam size
        avg = np.mean(inputBeamDist, axis = 0)
        std = np.std( inputBeamDist, axis = 0)
        
        # Step 2 - invert the beam size for normalization and check for "cold" dimensions, e.g. zero energy spread
        std_inv = np.where(std == 0, 1, 1/std)
        
        # Step 3 - normalize distribution so that it is aligned to the origin and has an rms size of unity in all dimensions
        beamDist[:nlen] = (beamDist[:nlen]-avg)*std_inv
        
        # Step 4 - add particles
        for i in np.arange(nlen, mpart):
            n1 = int(np.random.rand()*nlen)
            rmin = 1e9
            
            ## no loop
            selected = beamDist[n1]
            distance = np.sum((np.array([selected]*nlen)-beamDist[:nlen])**2, axis = 1)
            n2 = np.argmin(distance)
            # # Use loop as done in Genesis
            # for j in np.arange(nlen):
            #     distance = np.sqrt(np.sum((beamDist[n1]-beamDist[j])**2))
            #     if distance<rmin and j != n1:
            #         n2 = j
            #         rmin = distance
            
            #print(i, n1, n2)      
            beamDist[i] = 0.5*(beamDist[n1]+beamDist[n2])+(2*np.random.rand()-1)*(beamDist[n1]-beamDist[n2])
            
        # Step 5 - scale back
        beamDist = beamDist/std_inv+avg
    
    return beamDist

#%% Add modulation to a distribution
def modulation1D(x, funcMod, args = (), centered = True, xc = None):
    '''
    Modulate the input dataset with the predefined modulation function.

    Parameters
    ----------
    x : 1D array-like
        The dataset to be modulated.
    funcMod : 1D function
        The probability density function of the modulated profile.
    args: tuple, optional
        Extra arguments passed to the function.
    centered: boolean, optional
        True by default. If True, center the distribution, that is to shift it
        by its geometric average before applying the modulation.
    xc: float, optional  
        If not None, then center the distribution at xc. To be implemented.
    Returns
    -------
    xnew: 1D array-like
        The modulated dataset.

    Examples:
    ```
    x = np.random.rand(100000)
    funcMod = lambda x, k:np.sin(k*x)**2
    args = [ 2*np.pi/0.5 ]
    
    xnew = modulation1D(x, funcMod, *args)
    plt.figure()
    
    _ = plt.hist(x, bins = 100, histtype = r'step')
    _ = plt.hist(xnew, bins = 100, histtype = r'step')
    ```
    '''
    
    xmin = x.min()
    xmax = x.max()
    
    if centered:
        xavg = np.mean(x)
        if xmax-xavg>xavg-xmin:
            xmin = xavg-(xmax-xavg)
        else:
            xmax = xavg+(xavg-xmin)
    
    # normalize the dataset
    xnorm = (x-xmin)/(xmax-xmin)
    
    xbin = np.linspace(0, 1, 1001)
    freq = funcMod(xbin*(xmax-xmin)+xmin, *args)
    freq = freq/np.sum(freq)
    for i in np.arange(1, len(freq)):
        freq[i] += freq[i-1]
    
    from scipy.interpolate import interp1d
    funcInterp = interp1d(freq, xbin, bounds_error = False, fill_value = 0)
    
    xnew = funcInterp(xnorm)
    # scale back
    xnew = xnew*(xmax-xmin)+xmin
    
    return xnew

def modulation2D(x, y, funcMod, args = (), centered = True, xc = None):
    '''
    Modulate the input dataset with the predefined modulation function.

    Parameters
    ----------
    x : 1D array-like
        The first dataset to be modulated.
    y : 1D array-like
        The second dataset to be modulated.
    funcMod : modulation function
        The probability density function of the modulated profile.
    args: tuple, optional
        Extra arguments passed to the function.
    centered: boolean, optional
        True by default. If True, center the distribution, that is to shift it
        by its geometric average before applying the modulation.
    xc: float, optional  
        If not None, then center the distribution at xc. To be implemented.
    Returns
    -------
    xnew: 1D array-like
        The modulated dataset.

    Examples:
    ```
    
    ```
    '''
    # to be implemented
    return None
