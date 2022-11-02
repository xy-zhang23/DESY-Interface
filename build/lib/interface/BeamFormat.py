# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 21:22:18 2020

@author: lixiangk
"""
from .Numeric import *
from .BeamDiagnostics import *

from scipy import fft
import h5py

fmt = '%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%4d%4d'
astra_format = fmt

#%% Scale the input distribution according to the given match parameters
#   This is a remake of the method in Genesis-1.3-Version4
def Matching(inputName = None, inputDist = None, outputName = None,
             betax = 0, betay = 0, alphax = 0, alphay = 0, flipZ = False):
    '''
    Astra input particle format to hdf5 format used in Genesis1.3 Version 4

    Parameters
    ----------
    inputName : string, optional
        Name of the input particle file. If defined, prefered than inputDist. 
        None by default. 
    inputDist : 10D array, optional
        Particle distribution from Astra, with absolute z and Pz. None by default. 
    outputName: string, optional
        Name of the output particle file.
    betax, betay, alphax, alphay: Twiss parameters to be matched
    Returns
    -------
    dist : 10D array 
        Particle distribution in Astra format
    
    --------------
    Examples
    --------
    ```
    fname = 'ast.2529.003'
    fout = 'ast_.2529.013'
    kwargs = {}
    kwargs.update(
        inputName = fname,
        outputName = fout,
        betax = 5.115,
        betay = 1.634,
        alphax =  7.326,
        alphay = -0.897)
    Matching(**kwargs)
    ```
    '''
    # Load the astra distribution
    if inputName != None:   
        dist = pd_loadtxt(inputName)
        dist[1:,2] += dist[0,2]
        dist[1:,5] += dist[0,5]
    elif inputDist != None:
        dist = inputDist
    else:
        print('No input file or data!')
        return
    
    if outputName == None:
        outputName = 'temp.dat'
    print('The distribution is saved to '+outputName)
    
    diag = BeamDiagnostics(dist = dist)
    print('Before matching:')
    for key in ['nemit_x', 'nemit_y', 'beta_x', 'beta_y', 'alpha_x', 'alpha_y', 'cor_Ekin']:
        i = diag.keyIndex.get(key)
        unit = diag.keyUnit.get(key)
        print(str('%17s: %15.6E %s' % (key, diag.x[i], unit)))
    #diag.demo()
    
    # Define the Twiss parameters of the matched distribution
    match = 1
    gammax, gammay = (1+alphax**2)/betax, (1+alphay**2)/betay
    
    if match:
        
        bx, by, ax, ay = diag.beta_x, diag.beta_y, diag.alpha_x, diag.alpha_y
        px = dist[:,3]/dist[:,5] # px->xp
        py = dist[:,4]/dist[:,5] # py->yp
        
        px += (ax/bx)*dist[:,0]
        py += (ay/by)*dist[:,1]
        
        dist[:,0] *= np.sqrt(betax/bx)
        dist[:,1] *= np.sqrt(betay/by)
        
        px *= np.sqrt(bx/betax)
        py *= np.sqrt(by/betay)
        
        px -= (alphax/betax)*dist[:,0]
        py -= (alphay/betay)*dist[:,1]
        
        dist[:,3] = px*dist[:,5]
        dist[:,4] = py*dist[:,5]
       
    diag = BeamDiagnostics(dist = dist)
    print('After matching:')
    for key in ['nemit_x', 'nemit_y', 'beta_x', 'beta_y', 'alpha_x', 'alpha_y', 'cor_Ekin']:
        i = diag.keyIndex.get(key)
        unit = diag.keyUnit.get(key)
        print(str('%17s: %15.6E %s' % (key, diag.x[i], unit)))

    # Use relative coordinates to reference particle in z and Pz
    dist[1:,2] -= dist[0,2]
    dist[1:,5] -= dist[0,5]

    if flipZ:
        dist[1:,2] *= -1

    np.savetxt(outputName, dist, fmt = astra_format)
    return dist

#%% Convert 6D particle distributions in ascii format (used in Astra) to hdf5
#   format for Genesis1.3 V4
def astra2hdf5(inputName = None, inputDist = None, outputName = None):
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
        outputName = 'temp.h5'
    print('The distribution is saved to '+outputName)
    
    x, y, z = data[:,0:3].T[:]
    t = z/g_c
    
    tminps, tmaxps = t.min()*1e12, t.max()*1e12
    print('Bunch length in ps: ', tmaxps-tminps)
    
    Px, Py, Pz = data[:,3:6].T[:]
    xp = Px/Pz
    yp = Py/Pz
    p = np.sqrt(1+(Px**2+Py**2+Pz**2)/1e12/g_mec2/g_mec2) # gamma
    
    
    # saving to hdf5
    f = h5py.File(outputName, 'w')
    
    dx = f.create_dataset("x", data = x)
    dy = f.create_dataset("y", data = y)
    dt = f.create_dataset("t", data = t)
    
    dxp = f.create_dataset("xp", data = xp)
    dyp = f.create_dataset("yp", data = yp)
    dp = f.create_dataset("p", data = p)
    
    f.close()
    return

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
#%% Convert Astra distribution to slice-wise parameters for Genesis1.3-Version2
def astra2slice(inputName, outputName = None, nslice = 100, nc = 9):
    '''
    The output file follows the format required for Genesis 1.3 simulation.
    From left the right, the column are:
        ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS
    Parameters
    ----------
    inputName: string
        filename of Astra output
    outputName: string
        filename of the output
    nslice: int, optional
        number of slices
    nc: int, optional
        number of slices to discard on both sides
    
    Returns
    -------
    slicePara: 2D array
        2D array of the slice-wise parameters.
    
    Examples
    --------
    ```
    # Convert astra distribution to slice-wise parameters with 100 slices
    r = astra2slice('ast.0528.001', nslice = 100, nc = 2)
    ```
    '''
    
    header = '? VERSION = 1.0\n'+\
    '? COLUMNS ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS'
    
    data = np.loadtxt(inputName)
    data[1:,2] += data[0,2]
    data[1:,5] += data[0,5]
    
    nop = len(data[:,0])
    #select = (data[:,0]>=-3e-3)*(data[:,0]<=3e-3)*(data[:,1]>=-1e-3)*(data[:,1]<=1e-3)
    select = (data[:,9]>0); data = data[select]
    nop = len(data[:,0])
    
    nc = nc
    counts, edges = np.histogram(data[:,2], bins = nslice)
    edges = (edges[:-1]+edges[1:])/2.0
    #print(counts, edges, nc)
    
    zmin, zmax = edges[0], edges[-1]
    #print(zmin, zmax)
    
    for i in np.arange(len(counts[0:nslice//2])):
        if counts[i] < nc:
            zmin = edges[i+1]; print(zmin)
    for i in np.arange(len(counts[nslice//2:])):
        if counts[nslice-i-1] < nc:
            zmax = edges[nslice-i-2]
    #print(zmin, zmax)
    
    
    dz = (zmax-zmin)/nslice
    
    r = []; zpos = 0; counts = 0
    for i in np.arange(nslice):
        zmin_i, zmax_i = zmin+i*dz, zmin+i*dz+dz
        select = (data[:,2]>=zmin_i)*(data[:,2]<=zmax_i)
        beam_i = data[select]; #print(len(beam_i))
        
        diag = BeamDiagnostics(dist = beam_i);
        Ipeak = -diag.Q_b/dz*g_c
        
        if i == 0:
            tmp = [zpos, kinetic2gamma(diag.Ekin), diag.std_Ekin*1e-3/g_mec2, \
                   diag.nemit_x*1e-6, diag.nemit_y*1e-6, diag.beta_x, diag.beta_y, \
                   0, 0, 0, 0, diag.alpha_x, diag.alpha_y, 0, 0]
            r.append(tmp); #counts += 1; print(counts)
            zpos += dz/2.0
            
            
        tmp = [zpos, kinetic2gamma(diag.Ekin), diag.std_Ekin*1e-3/g_mec2, \
             diag.nemit_x*1e-6, diag.nemit_y*1e-6, diag.beta_x, diag.beta_y, \
             0, 0, 0, 0, diag.alpha_x, diag.alpha_y, Ipeak, 0]
        r.append(tmp); #counts += 1; print(counts)
        
        zpos += dz
        
        if i == nslice-1:
            zpos -= dz/2.0
            tmp = [zpos, kinetic2gamma(diag.Ekin), diag.std_Ekin*1e-3/g_mec2, \
                   diag.nemit_x*1e-6, diag.nemit_y*1e-6, diag.beta_x, diag.beta_y, \
                   0, 0, 0, 0, diag.alpha_x, diag.alpha_y, 0, 0]
            r.append(tmp); #counts += 1; print(counts)
            
            
    if outputName == None:
        outputName = 'slice@'+inputName+'.dat'
    print('The distribution is saved to '+outputName)
    np.savetxt(outputName, np.array(r), fmt = '%-15.6e', \
               header = header, comments='')
    return np.array(r)

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
def resampleParticles(inputBeamDist, mpart = 8192, force = False, closest = None):
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
    force : boolean
        If true, always returned the resampled distribution
    Returns
    -------
    beamDist : 2D array of size (mpart, 6)
        6D particle distribution, with columns storing coordinates of x, y, theta,
        px, py, gamma, respectively.
        
    Examples
    --------
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
    inputBeamDist = np.copy(inputBeamDist[:,0:6])
    
    if nlen>=mpart and force == False:
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
        if force == False:
            beamDist[:nlen] = np.copy(inputBeamDist)
            start = nlen
        else:
            start = 0
            
        # Step 1 - calculate the center and the rms beam size
        avg = np.mean(inputBeamDist, axis = 0)
        std = np.std( inputBeamDist, axis = 0)
        
        # Step 2 - invert the beam size for normalization and check for "cold" dimensions, e.g. zero energy spread
        std_inv = np.where(std == 0, 1, 1/std)
        
        # Step 3 - normalize distribution so that it is aligned to the origin and has an rms size of unity in all dimensions
        inputBeamDist[:nlen] = (inputBeamDist[:nlen]-avg)*std_inv
        
        # sort by z
        #temp = np.sum(inputBeamDist[:,2:3]*inputBeamDist[:,2:3], axis = 1)
        temp = inputBeamDist[:,2]
        inputBeamDist = inputBeamDist[np.argsort(temp)]
        
        if force == False:
            beamDist[:nlen] = (beamDist[:nlen]-avg)*std_inv
        
        # Step 4.0 - find the most closest particle of each particle
        if np.any(closest) == None:
            hw = 1024
            closest = np.zeros(nlen)
            for i in np.arange(nlen):
                l1 = np.max([0, i-hw])
                l2 = np.min([nlen, i+hw])
                ln = l2-l1
                
                n1 = i
                selected = inputBeamDist[n1]
                
                v1 = np.array([selected]*ln)
                v2 = inputBeamDist[l1:l2]
                
                distance = np.sum((v1-v2)**2, axis = 1)
                distance[n1-l1] = 1e20; #print(distance)
                
                n2 = np.argmin(distance); #print(n1, n2, l1, l2)
                closest[n1] = l1+n2
        
        # Step 4.1 - add particles
        #nn = int(np.random.rand(mpart-start)*nlen)
        
        # Step 4 - add particles
        for i in np.arange(start, mpart):
            n1 = int(np.random.rand()*nlen)
            
            ## no loop
            # selected = inputBeamDist[n1]
            # distance = np.sum((np.array([selected]*nlen)-inputBeamDist[:nlen])**2, axis = 1)
            # distance[n1] = 1e20
            # n2 = np.argmin(distance)
            
            n2 = int(closest[n1])
            
            ## Use loop as done in Genesis
            # rmin = 1e9
            # for j in np.arange(nlen):
            #     distance = np.sqrt(np.sum((beamDist[n1]-beamDist[j])**2))
            #     if distance<rmin and j != n1:
            #         n2 = j
            #         rmin = distance
            
            if(i%10000 == 0): print(i, n1, n2)      
            beamDist[i] = 0.5*(inputBeamDist[n1]+inputBeamDist[n2])+\
                          (2*np.random.rand()-1)*(inputBeamDist[n1]-inputBeamDist[n2])
            
        # Step 5 - scale back
        beamDist = beamDist/std_inv+avg
    
    return beamDist, closest

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

    Examples
    --------
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
    
    #
    xbin = np.linspace(0, 1, 64*120+1)
    freq = funcMod(xbin*(xmax-xmin)+xmin, *args)
    freq = freq/np.sum(freq)
    for i in np.arange(1, len(freq)):
        freq[i] += freq[i-1]
    #
    # from scipy.integrate import quad
    # xbin = np.linspace(0, 1, 10001)
    # freq = np.zeros(len(xbin))
    # for i in np.arange(1, len(freq)):
    #     y, err = quad(funcMod, xbin[i-1], xbin[i], *args)
    #     freq[i] = freq[i-1]+y
    # print(freq[-1])
    
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

    Examples
    --------
    ```
    
    ```
    '''
    # to be implemented
    return None


#%% Astra and Warp formats conversion
def astra2warp(fname, fout = None, Q_coef = -1.0, High_res = True):
    '''
    The output file follows the format required for Warp simulation.
    From left the right, the column are:
      X Y Z UX UY UZ W, where Ux Uy Uz are dimentionless momentum, W is macro particle charge
    Parameters
      fname: filename of Astra output
      fout: filename of the output
      Q_coef: an coefficient to scale the bunch charge, default is -1.0
    '''
    data = np.loadtxt(fname)
    data[1:,2] += data[0,2]
    data[1:,5] += data[0,5]
    
    select = (data[:,9]>0); data = data[select]
    data[:,3] = data[:,3]/g_mec2/1e6
    data[:,4] = data[:,4]/g_mec2/1e6
    data[:,5] = data[:,5]/g_mec2/1e6
    
    data[:,6] = data[:,7]*1e-9/g_qe*Q_coef # number of electrons for each macro particle
    
    if fout is None:
        fout = fname+'.warp'
    print('The distribution is saved to '+fout)
    if High_res:
        fmt = '%20.12E'
    else:
        fmt = '%14.6E'
    np.savetxt(fout, data[:,0:7], fmt = fmt)
    return data[:,0:7]

def warp2astra(fname, fout = None, Run = 1, Q_coef = -1.0, ratio = 100):
    '''
    The output file follows the format required for Astra simulation.
    Parameters
      fname: filename of Warp output, which includes columns of X Y Z UX UY UZ W, where Ux Uy Uz are dimentionless momentum, W is macro particle charge
      fout: filename of the output
      Q_coef: an coefficient to scale the bunch charge, default is -1.0
      ratio: a factor used to scale the position of electron bunch to be used in Astra file name
    '''
    
    data = np.loadtxt(fname)
    data[:,3:6] *= g_mec2*1e6
    
    z0 = weighted_mean(data[:,2], data[:,-1])
    pz0 = weighted_mean(data[:,5], data[:,-1])
    
    data[:,2] -= z0
    data[:,5] -= pz0; #print pz0; return
    data[:,6] *= (-g_qe*1e9) # convert to nC
    
    nop = len(data)
    d1 = np.zeros((nop+1, 10))
    d1[0, 2] = z0
    d1[0, 5] = pz0
    
    d1[1:,0:6] = data[:,0:6]
    d1[1:,7] = data[:,6]
    d1[:, -2] = 1
    d1[:, -1] = 5
    
    if fout is None:
        print('The current bunch center is at ', z0, ' meters')
        fid = z0*ratio
        while round(fid)<1:
            fid *= 10
        fid = round(fid)
        fout = 'ast.%04d.%03d' % (fid, Run)
    print('The distribution is saved to '+fout)
    
    np.savetxt(fout, d1, fmt = '%12.4E%12.4E%12.4E%12.4E%12.4E%12.4E%12.4E%12.4E%4d%4d')
    return d1
