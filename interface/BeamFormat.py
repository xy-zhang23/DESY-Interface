# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 21:22:18 2020

@author: lixiangk
"""
from .Numeric import *

def astra2slice(fname, fout = None, nslice = 100, nc = 9):
    '''
    The output file follows the format required for Genesis 1.3 simulation.
    From left the right, the column are:
      ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS
    Parameters
      fname: filename of Astra output
      fout: filename of the output
      nslice: number of slices
      nc: number of slices to discard on the sides
    '''
    
    header = '? VERSION = 1.0\n'+\
    '? COLUMNS ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS'
    
    data = np.loadtxt(fname)
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
            
            
    if fout is None:
        fout = 'slice@'+fname+'.dat'
    print('The distribution is saved to '+fout)
    np.savetxt(fout, np.array(r), fmt = '%-15.6e', \
               header = header, comments='')
    return np.array(r)

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
