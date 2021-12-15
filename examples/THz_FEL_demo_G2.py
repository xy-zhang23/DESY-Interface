# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 17:19:37 2020

@author: lixiangk
"""

from interface import *

P0 = 17 # MeV/c
#P0 = 2.054595E+01

Q = 4000e-12 # Coulumb

lam_s = 100e-6
Nu, lam_u, B = 120, 3.0e-2, 1.2799839

K = 3.49

FWHM = 20e-12 # second
sigma_z = FWHM*g_c
curpeak = -Q/FWHM; print ('peak current: ', curpeak)

Ek = momentum2kinetic(P0) # MeV
gamma = kinetic2gamma(Ek); # unitless
delgam = gamma*0.5e-2

freq_s = g_c/lam_s
print ('resonant wavelength from kinetic energy: ', lam_s*1e6, ' um')

# define the # of slices, usually the slices in the electron bunch + the number of undulator period
nslice = Nu+sigma_z/lam_s; print ('# of slice: ', nslice-Nu)
nslice = 150
ntail = 0

# define electron parameters, which are dominated by the input particle file, thougth
emit_x, emit_y = 4e-6, 4e-6
sigma_x, sigma_y = 1.58e-3, 0.18543e-3
alpha_x, alpha_y = 7.18, 1.72

delz = 0.5
ipseed = 1
fname = 'pithz.%d' % ipseed

beamfile = '../pithz.2745.019'
maginfile = '../LCLS-I.lat'

gen = Genesis2()

# set undulator parameters
gen.set(aw0   = K/np.sqrt(2.),
        awd   = K/np.sqrt(2.),
        nwig  = Nu,
        xlamd = lam_u,
        iertyp= 0,
        delaw = 0,
        iseed = -1)

# set electron beam parameters
gen.set(gamma0  = gamma,
        delgam  = delgam,
        curpeak = curpeak,
        rxbeam  = sigma_x,
        rybeam  = sigma_y,
        emitx   = emit_x,
        emity   = emit_y,
        alphax  = alpha_x,
        alphay  = alpha_y,
        npart   = 8192*2)

# set particle-loading parameters
gen.set(ipseed = 1,
        nbins  = 16)

# set mesh paremeters
gen.set(nptr = 170,
        dgrid= 0.02, 
        nscr = 2, 
        nscz = 1)

# set time-dependent parameters/grids
gen.set(itdp   = 1,
        nslice = nslice,
        ntail  = ntail,
        iotail = 1,
        curlen = sigma_z)

# set radiation parameters
gen.set(xlamds = lam_s, 
        zrayl  = 3.5251E-02,
        ncar = 151,
        prad0 = 0) 

# set simulation control parameters
gen.set(delz  = delz,
        zstop = Nu*lam_u)

# set input and ouput control parameters
gen.set(ffspec = 1,
        ipradi = 0,
        beamfile   = beamfile,
        outputfile = fname+'.out',
        maginfile = maginfile,
        lout = [1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1])

#gen.write(fname+'.in')
#gen.qsub(fname+'.sh')
print(gen.output)

#%% Batch generating inputs for parameter scan
var1 = np.arange(1, 21)
var2 = [0, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]

combi = np.array([[v1, v2] for v1 in var1 for v2 in var2])

for x in combi:
    
    # undulator error 
    # iseed = np.random.randint(5000, 10000)
    # Kerr = 0 # K/np.sqrt(2.)*(0.01+0.01*rr)
    # gen.set(iseed = iseed, delaw = Kerr, iertyp = 0)
    
    # shot noise
    ipseed = int(x[0])
    
    # bunching factor
    bunch = x[1]
    
    # output name
    fname = 'pithz.%d.%.6f' % (ipseed, bunch)
    
    gen.set(ipseed = ipseed)
    gen.set(bunch = bunch)
    
    gen.set(ippart = 30, ispart = 10)
    
    gen.set(outputfile = fname+'.out')

    gen.write(fname+'.in')
    gen.qsub(fname+'.sh', fname+'.in')

exit()
