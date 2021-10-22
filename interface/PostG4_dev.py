# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 14:57:04 2020

@author: lixiangk
"""

from universal import *
from interface.postG4 import *

#%% Load the output of version 2
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2020\Genesis-demo\version2'
os.chdir(workdir)

fname = 'pithz.1.out'
pg = PostGenesis(fname, debug = 1)

#%% OR load the output of version 4
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2020\Genesis-demo\feedback2'
os.chdir(workdir)

fname = 'test.out.h5'
pg = PostGenesis(fname, debug = 1)
#%% Animate the radiation profile along the undulator
#from IPython import display

zplot = pg.zplot
temp = pg.get_fielddata('power')/1e6

ppower = np.max(temp, axis = 1)
ymax = np.max(temp)/1e6

fig = plt.figure()
ax1 = fig.add_subplot(5, 1, (1, 3))
ax2 = fig.add_subplot(3, 1, 3)

step = 5
for i in np.arange(0, temp.shape[0], step):
    ax1.plot(temp[i], '-')
    
    ax1.set_title('z = : %.3f m' % zplot[i])
    ax1.set_xlabel('# of slice')
    ax1.set_ylabel('Power (MW)')
    #ax.set_ylim(0, ymax)
    ax1.grid()
    
    ax2.plot(zplot[:i], ppower[:i], '-')
    ax2.set_xlabel(r'$z$ (m)')
    ax2.set_ylabel(r'Peak power (MW)')
    ax2.grid()
    
    #display.display(plt.gcf())
    #display.clear_output(wait=True)
    
    plt.pause(0.1)
    if i<temp.shape[0]-step:
        ax1.cla()
        ax2.cla()

#%% Animate the bunching factor along the undulator
zplot = pg.zplot

temp = pg.get_fielddata('power')/1e6
ppower = np.max(temp, axis = 1)

temp = pg.get_beamdata('bunching')

ymax = np.max(temp)

fig = plt.figure()
ax1 = fig.add_subplot(5, 1, (1, 3))
ax2 = fig.add_subplot(3, 1, 3)

step = 5
for i in np.arange(0, temp.shape[0], step):
    ax1.plot(temp[i], '-')
    
    ax1.set_title('z = : %.3f m' % zplot[i])
    ax1.set_xlabel('# of slice')
    ax1.set_ylabel('Bunching')
    #ax.set_ylim(0, ymax)
    ax1.grid()
    
    ax2.plot(zplot[:i], ppower[:i]/1e6, '-')
    ax2.set_xlabel(r'$z$ (m)')
    ax2.set_ylabel(r'Total power (MW)')
    ax2.grid()
    
    plt.pause(0.1)
    if i<temp.shape[0]-step:
        ax1.cla()
        ax2.cla()
        
#%% Plot the bunching factor of signle slice along the undulator

current = pg.current
zplot = pg.zplot

# print(pg.outputs) # check the name of the property of interest
temp = pg.get_beamdata('bunching')

fig, [ax1, ax2] = plt.subplots(nrows = 2)

ax1.plot(current, 'k-')
ax1.set_xlabel(r'# of slices')
ax1.set_ylabel(r'Current (A)')
ax1.grid()

slices = [20, 30, 40]
for i in slices:
    ax2.plot(zplot, temp[:,i], '-')

ax2.set_xlabel(r'$z$ (m)')
ax2.set_ylabel(r'Bunching factor')
ax2.legend(['# %d' % i for i in slices])
ax2.grid()

fig.tight_layout()
#fig.savefig('bunching-factor.eps')

#%% Plot size of the whole radiation along the undulator

zplot = pg.zplot
current = pg.current

name = 'r_size'
name = 'xsize' # or 'ysize'
r_size = (pg.get_fielddata(name) @ current)/np.sum(current)
                      
fig, ax = plt.subplots(figsize = (5, 3))

ax.plot(zplot, r_size*1e3, 'r-')
ax.grid()

ax.set_xlabel(r'$z$ (m)')
ax.set_ylabel(r'Radiation trans. size (mm)')

#%% Plot beam size of the whole electron bunch along the undulator

x = pg.zplot
current = pg.current

# version 2
#xrms = (pg.get_fielddata('xrms') @ current)/np.sum(current)
#yrms = (pg.get_fielddata('yrms') @ current)/np.sum(current)

# version 4
xrms = (pg.get_data('Beam', 'xsize') @ current)/np.sum(current)
yrms = (pg.get_data('Beam', 'ysize') @ current)/np.sum(current)
                   
fig, ax = plt.subplots(figsize = (5, 3))

ax.plot(x, xrms*1e3, 'r-')
ax.plot(x, yrms*1e3, 'b-')
ax.grid()

ax.set_xlabel(r'$z$ (m)')
ax.set_ylabel(r'RMS size (mm)')
ax.legend([r'$x$', r'$y$'])

#%% Batch process outputs of parameter scan, e.g., the random seed
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2019\LCLS-I-matching-100um\_n200k-sig_x-1.55mm-sig_y-0.30mm-alp_x-10.93-alp_y-3.25\SASE-002'
os.chdir(workdir)

#fname = 'pithz.1.out'
#pg = PostGenesis(fname, debug = 1)

ext = '_test.txt'
for i in np.arange(1, 11):
    try:
        ipseed = i+0; print(ipseed, end = ' ')
        fname = 'pithz.%d.out' % ipseed
        # fname = 'pithz.%dA.%.1fps.%03d.out' % (Ipeak, FWHM, ipseed)
        
        pg = PostGenesis(fname)

        zplot = pg.zplot
        zpower = pg.zpower
        zenergy = pg.zenergy
        
        ppower = np.max(pg.get_fielddata('power'), axis = 1)
        
        tt = pg.zbunch/g_c
        tspec = pg.get_fielddata('power', at = 3.6)
        ww, fspec = pg.get_spectrum(at = 3.6)

        if i == 1:
            with open('./power-z-ipseed'+ext, 'w') as f_handle:
                np.savetxt(f_handle,np.atleast_2d(zplot),fmt='%15.6E')
            with open('./energy-z-ipseed'+ext, 'w') as f_handle:
                np.savetxt(f_handle,np.atleast_2d(zplot),fmt='%15.6E')
            with open('./peak-power-z-ipseed'+ext, 'w') as f_handle:
                np.savetxt(f_handle,np.atleast_2d(zplot),fmt='%15.6E')
                
            with open('./spectrum-lamds-ipseed'+ext, 'w') as f_handle:
                np.savetxt(f_handle,np.atleast_2d(ww),fmt='%15.6E')
            with open('./power-t-ipseed'+ext, 'w') as f_handle:
                np.savetxt(f_handle,np.atleast_2d(tt),fmt='%15.6E')
            
        with open('./power-z-ipseed'+ext, 'a') as f_handle:
            np.savetxt(f_handle,np.atleast_2d(zpower),fmt='%15.6E')
        with open('./energy-z-ipseed'+ext, 'a') as f_handle:
            np.savetxt(f_handle,np.atleast_2d(zenergy),fmt='%15.6E')
        with open('./peak-power-z-ipseed'+ext, 'a') as f_handle:
            np.savetxt(f_handle,np.atleast_2d(ppower),fmt='%15.6E')
            
        with open('./spectrum-lamds-ipseed'+ext, 'a') as f_handle:
            np.savetxt(f_handle,np.atleast_2d(fspec),fmt='%15.6E')
        with open('./power-t-ipseed'+ext, 'a') as f_handle:
            np.savetxt(f_handle,np.atleast_2d(tspec),fmt='%15.6E')
    except Exception as err:
        print(err)
        print('Error got at: ', i, )
        pass
print('\n')

#%% Plot: energy and power along z
ext = '_test'
pp = np.loadtxt('peak-power-z-ipseed'+ext+'.txt'); pp = pp.T; print (pp.shape)
EE = np.loadtxt('energy-z-ipseed'+ext+'.txt'); EE = EE.T; print (EE.shape)

end = pp.shape[1]

output = f'P =: {np.mean(pp[-1,1:end])/1e6:.2f} +/- {np.std(pp[-1,1:end])/1e6:.2f} MW\n'
output += f'E =: {np.mean(EE[-1,1:end])*1e6:.2f} +/- {np.std(EE[-1,1:end])*1e6:.2f} uJ\n'
print(output)

with open('stat'+ext+'.dat', 'w') as f_handle:
    f_handle.write(output)

fig, ax1 = plt.subplots(figsize = (4, 4))
for i in np.arange(1, end):
    ax1.plot(pp[:,0], pp[:,i]/1e6, '-', color='grey') # energy: MW
ax1.plot(pp[:,0], np.mean(pp[:,1:], 1)/1e6, 'k-')

ax1.set_xlabel(r'$z$ (m)')
ax1.set_ylabel(r'power (MW)')
ax1.set_yscale('log')
ax1.grid()
fig.savefig('ppower-vs-z-ipseed'+ext+'.png')

fig, ax2 = plt.subplots(figsize = (4, 4))
for i in np.arange(1, end):
    print(i, end = ' ')
    ax2.plot(EE[:,0], EE[:,i]*1e6, '-', color='grey') # energy: uJ
ax2.plot(EE[:,0], np.mean(EE[:,1:], 1)*1e6, 'k-')
#res.append(EE[:,0])
#res.append(np.mean(EE[:,1:], 1))

ax2.set_xlabel(r'$z$ (m)', fontsize = 12)
ax2.set_ylabel(r'$E$ ($\mu$J)', fontsize = 12)
ax2.set_yscale('log')
#ax2.set_ylim(0, 150)
#ax2.set_yticks([1e-7, 1e-6, 1e-6, 1e])
ax2.grid()

#locmin = mpl.ticker.LogLocator(base = 10.0, subs = (0.2, 0.4, 0.6, 0.8), numticks = 12)
#ax2.yaxis.set_minor_locator(locmin)
#ax2.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

#ax2.text(0.925, 0.95, '($b$)', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)

fig.savefig('energy-vs-z-ipseed'+ext+'.png', density = 300)

#%% Plot: spectrum and waveform profile
ext = '_test'
ss = np.loadtxt('spectrum-lamds-ipseed'+ext+'.txt'); ss = ss.T; print (ss.shape)
pp = np.loadtxt('power-t-ipseed'+ext+'.txt'); pp = pp.T; print (pp.shape)

end = ss.shape[1]

xlamds = []; width = []
fig, ax = plt.subplots(figsize = (4, 4))
for i in np.arange(1, end):
    ax.plot(ss[:,0]*1e6, ss[:,i]/1e8, '-', color='grey') # energy: MW
    xlamds.append([weighted_mean(ss[:,0], ss[:,i])*1e6])
    width.append([weighted_std(ss[:,0], ss[:,i])*1e6])
ax.plot(ss[:,0]*1e6, np.mean(ss[:,1:], 1)/1e8, 'k-')

output =  f'Centre wavelength: {np.mean(xlamds):.2f} +/- {np.std(xlamds):.2f} um\n'
output += f'Spectrum width: {np.mean(width):.2f} +/- {np.std(width):.2f} um\n'

ax.set_xlabel(r'$\lambda_s$ ($\mu$m)', fontsize = 12)
ax.set_ylabel(r'intensity (arb. unit)', fontsize = 12)
#ax.set_yscale('log')
ax.set_xlim(80, 120)
#ax.set_ylim(-0.05, 1.1)
ax.grid()

#ax.text(0.925, 0.95, '($b$)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
fig.savefig('spectrum-vs-lamds-ipseed'+ext+'.png')


jitter = []
fig, ax = plt.subplots()
for i in np.arange(1, end):
    ax.plot(pp[:,0]*1e12-20, pp[:,i]/1e8, '-', color='grey') # energy: MW
    jitter.append([weighted_mean(pp[:,0], pp[:,i])/g_c*1e12-20])
ax.plot(pp[:,0]*1e12-20, np.mean(pp[:,1:], 1)/1e8, 'k-')

output += f'Arrival time jitter: {np.mean(jitter):.2f} +/- {np.std(jitter):.2f} ps\n'
print(output)

with open('stat'+ext+'.dat', 'a') as f_handle:
    f_handle.write(output)

ax.set_xlabel(r'$t$ (ps)')
ax.set_ylabel(r'intensity (arb. unit)')
#ax.set_yscale('log')
ax.set_xlim(-20, 20)
#ax.set_ylim(-0.005, 0.6)
ax.grid()
fig.savefig('power-vs-t-ipseed'+ext+'.png')

#%%
f = h5py.File('test.out.h5', 'r')
f.keys()
#%%
f = h5py.File('test.200.par.h5', 'r')

slice = f['slice000001']

current = slice.get('current')[0]
print('The current is ', current)

x = slice.get('x')[:]
y = slice.get('y')[:]
px = slice.get('px')[:]
py = slice.get('py')[:]
theta = slice.get('theta')[:]
gamma = slice.get('gamma')[:]

# fig, [ax, ax2] = plt.subplots(ncols = 2, figsize = (8 ,4))
# ax.plot(x, px, '.')
# ax.plot(y, py, '.')
# ax.grid()

ax2.plot(theta, gamma, '.')
ax2.grid()

#%%

plt.figure()
plt.xlabel('theta')
plt.ylabel('counts')

for step in np.linspace(0, 240, 240//5+1):
    f = h5py.File('test.%d.par.h5' % step, 'r')

    slice = f['slice000001']
    theta = slice.get('theta')[:]

    plt.hist(np.mod(theta, 2*np.pi), 100, histtype = r'step')
    plt.xlim(-30, 30)
    
    plt.pause(0.1)
    if step<240:
        plt.cla()


#%%
#f = h5py.File('../Benchmark.out.h5', 'r')
f = h5py.File('test.out.h5', 'r')

beam = f.get('Beam')

xsize = beam.get('xsize')[:]
ysize = beam.get('ysize')[:]
current = beam.get('current')[:]
bunching = beam.get('bunching')[:]

lattice = f.get('Lattice')
z = lattice.get('z')[:]
aw = lattice.get('aw')[:]

fig, [ax1, ax2, ax3] = plt.subplots(ncols = 3, figsize = (12, 4))

ax1.plot(z, aw, '-o')
ax1.set_xlabel(r'$z$ (m)')
ax1.set_ylabel(r'$a_W$')
ax1.grid()

ax2.plot(xsize, '-')
ax2.plot(ysize, '-')
ax2.set_xlabel(r'$z$ (m)')
ax2.set_ylabel(r'RMS size (mm)')
ax2.legend(['x', 'y'])
ax2.grid()

#ax3.plot(current, '-')
ax3.plot(bunching, '-')
ax3.set_xlabel(r'$z$ (m)')
ax3.set_ylabel(r'bunching')
ax3.grid()

fig.tight_layout()
fig.savefig('test.png')


#%% Field
f = h5py.File('test.out.h5', 'r')

Field = f.get('Field')
power = Field.get('power')[:]

lam_s = 100e-6
energy = np.sum(power, axis = 1)*lam_s/g_c

zpos = np.linspace(0, 3.6, power.shape[0])

fig, ax = plt.subplots()
ax.plot(zpos, energy*1e6, '-o')
ax.set_xlabel(r'$z$ (m)')
ax.set_ylabel(r'Energy ($\mu$J)')
ax.set_yscale('log')

#%%
tmp = f.get('Beam')
for key in list(tmp.keys()):
    try:
        print(key, ': ', tmp.get(key)[:].T)
    except Exception as err:
        print(key, ': ')
        pass
    
#%%
betax = 20.772
betay = 0.270
alphax = 7.180000
alphay = 1.720000

print(233.52666069/alphax, 57.17006118/alphay, 20.30348682/betax, 0.26965361/betay)


#%%
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2020\Warp'
os.chdir(workdir)

beam1 = pd_loadtxt('beam_und_250k_4_shot.ini')
