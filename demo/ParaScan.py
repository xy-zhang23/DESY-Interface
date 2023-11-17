
from ObjectFunc import *

from timeit import default_timer
import numpy as np

import time

from scipy.interpolate import interp1d, interp2d

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

#%% Post processing of Astra parameter scan
def post_ParaScan(x, *args):
    '''
    Photocathode to the EMSY1 at 5.277 m
    Parameters:
      x: an array or list of the variables to be optimized
    Returns:
      energy: the quantity that is envisoned as the "energy" of the sample
    '''
    
    sigma_x = sigma_y = x[1]/4.
    phi_gun, phi_booster = x[3], x[5]
    Imain = x[6]
    
    phi_gun, phi_booster = x[3], x[5]
    
    MaxE_gun = get_MaxE_gun(phi_gun, x[2])
    #MaxE_gun = x[2]
    
    MaxE_booster = get_MaxE_booster(MaxE_gun, phi_gun, phi_booster, x[4])
    #MaxE_booster = x[4]
    
    MaxB = I2B(Imain)
    
    Q_total = x[0]/1e3
    #Ipart = int(Q_total*50e3)

    direc = str.format('Q-%.2fpC-D-%.2fmm-E1-%.2fMV_m-phi1-%.2fdeg-E2-%.2fMV_m-phi2-%.2fdeg-I-%.2fA' %\
                       (Q_total*1e3, x[1], MaxE_gun, phi_gun, MaxE_booster, phi_booster, Imain))
    print(direc)
    
    
    if 1:

        os.chdir(direc)
                
        try:
            fname = 'ast.0528.001'
            diag = BeamDiagnostics(fname = fname, energy = False)
        except:
            fname = 'ast.0527.001'
            diag = BeamDiagnostics(fname = fname, energy = False)

        
        res = list(x) + list(diag.x)
        
        os.chdir('..'); #os.system('rm -r '+direc)

        with open('ParaScan_250pC@5.28m.dat','a') as f_handle:
            np.savetxt(f_handle,np.atleast_2d(res),fmt='%14.6E')
            
    return 0


#%% Simple parameter scan

# Bunch charge
var0 = [250] 

# Diameter of BSA
var1 = np.linspace(0.9, 1.5, 7) 
var1 = [1.2]

# Gun phase w.r.t. MMMG phase
var2 = [0] 

# Booster phase w.r.t. MMMG phase
var3 = [0]

# Solenoid current
var4 = np.linspace(350, 370, 11)
var4 = [368]

# Define scan parameters, from left to right
# By default: charge, BSA, momentum after gun, gun phase, momentum after booster, booster phase, solenoid current
# charge, BSA, gun gradient, gun phase, booster gradient, booster phase, solenoid current
combi = np.array([[v0, v1, 6.3, v2, 17, v3, v4] for v0 in var0 
                  for v1 in var1 for v2 in var2 for v3 in var3 for v4 in var4])
                                  
for x in combi:
    # Produce Astra inputs
    obj_ParaScan(x)
    
    # Run post processor
    #post_MBIScan(x)
    pass

#%%
exit()

#%% Analyze results from post-processing
import pandas as pd

data = np.loadtxt('ParaScan_250pC@5.28m.dat')
scanned = ['charge', 'BSA', 'E1', 'phi1', 'E2', 'phi2', 'Imain']
diaged = list(BeamDiagnostics().keyIndex.keys())
df = pd.DataFrame(data, columns = scanned+diaged)

#%% nemit-vs-Imain-vs-
flag = 'BSA'
var = np.unique(df[flag])
var = [1.2]

fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(df[flag]-v0)<0.01; print(select)
    ds = df[select]
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds['Imain'], np.sqrt(ds['nemit_x']*ds['nemit_y']))

ax.grid()
ax.set_xlabel(r'Imain (A)')
ax.set_ylabel(r'Norm. emittance ($\mu$m)')

ax.legend(['%s = %.1f mm' % (flag, v0) for v0 in var])
#ax.legend(['%s = %.0f deg' % (flag, v0) for v0 in var])

fig.savefig('nemit-vs-Imain-vs-BSA.png')

#%% nemit-vs-Imain-vs-BSA
#select = np.abs(df['phi2']+5)<0.01; print(select)
#ds0 = df[select]

flag = 'BSA'
var = [0.6, 0.8, 1.0, 1.2]
fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(df[flag]-v0)<0.01; print(select)
    ds = df[select]
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds['Imain'], np.sqrt(ds['nemit_x']*ds['nemit_y']))

ax.grid()
ax.set_xlabel(r'Imain (A)')
ax.set_ylabel(r'Norm. emittance ($\mu$m)')

ax.set_xlim(350, 370)
ax.set_ylim(0, 1)

ax.legend(['%s = %.1f mm' % (flag, v0) for v0 in var])
#ax.legend(['%s = %.0f deg' % (flag, v0) for v0 in var])

fig.savefig('nemit-vs-Imain-vs-BSA2.png')

#%% rms-vs-Imain-vs-
fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(1-df[flag]/v0)<0.05
    ds = df[select]
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds['Imain'], np.sqrt(ds['std_x']*ds['std_y']))

ax.grid()
ax.set_xlabel(r'Imain (A)')
ax.set_ylabel(r'RMS size (mm)')

ax.set_xlim(350, 370)
ax.set_ylim(0, 1)

ax.legend(['%s = %.1f mm' % (flag, v0) for v0 in var])

fig.savefig('rms-vs-Imain-vs-BSA.png')

#%% cor_Ek-vs-phi2-vs-

#flag = 'Imain'; var = [384, 388, 392]

fig, ax = plt.subplots()
for v0 in var:
    #select = np.abs(1-df[flag]/v0)<0.01
    #ds = df[select]
    ds = df
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds['phi2'], ds['cor_Ekin'])
    ax.plot(ds['phi2'], ds['std_Ekin'])

ax.grid()
ax.set_xlabel(r'$\phi_{boo}$ (degree)')
ax.set_ylabel(r'Energy spread (keV)')

ax.legend([r'Correlated', 'RMS'])

fig.savefig('cor_Ekin-vs-phi2-BSA1.8mm3.png')