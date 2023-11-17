import sys
if sys.platform == "linux" or sys.platform == "linux2":
    sys.path.append(r'/afs/ifh.de/group/pitz/data/lixiangk/work/sync/python')
    rootdir = r'//afs/ifh.de/group/pitz/data/lixiangk/work'
elif sys.platform == "win32":
    sys.path.append(r'\\afs\ifh.de\group\pitz\data\lixiangk\work\sync\python')
    rootdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\work'
elif sys.platform == "darwin":
    print("OS X")
else:
    print("Unknown platform!")

from interface import *

from timeit import default_timer
import time

import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import curve_fit

field_maps = os.path.join(rootdir, 'sync', 'field-maps')
gun_profile = 'gun51cavity.txt'

I2B = lambda I: -(0.0000372+0.000588*I)
B2I = lambda B: (-B-0.0000372)/0.000588

# Prepare the interpolation function for the get_MaxE_booster()
#data_gun = np.loadtxt(field_maps+os.sep+'phiGun42_scan.dat')
#fEG_gun = interp1d(data_gun[:,0], data_gun[:,1])

#data_gun = np.loadtxt(field_maps+os.sep+'phiGun51_scan.dat', usecols = ())
#fEG_gun = interp1d(data_gun[:,0], data_gun[:,1])

### 2D gun gradient and phase scan and interpolation, updated 27.02.2023
data_gun = np.loadtxt(field_maps+os.sep+'gun51_scan2d.dat')
shape = (51, 91)
E_gun   = np.reshape(data_gun[:,0], shape)
phi_gun = np.reshape(data_gun[:,1], shape)
EG_gun  = np.reshape(data_gun[:,2], shape)

# fEG2d_gun = interp2d(E_gun, phi_gun, EG_gun, bounds_error = False)
from scipy.interpolate import RegularGridInterpolator
fEG2d_gun_t = RegularGridInterpolator((np.unique(E_gun), np.unique(phi_gun)), EG_gun,
                                bounds_error = False)
fEG2d_gun = lambda E, phi: fEG2d_gun_t(([E], [phi]))[0]

def get_MaxE_gun(phi_gun = 0, EG = 6.3):
    
    EG1 = np.array([fEG2d_gun(E0, phi_gun) for E0 in np.unique(E_gun)])
    fMaxE = interp1d(EG1, np.unique(E_gun),
                     bounds_error = False, fill_value = 'extrapolate')
    
    #EG1 = fEG2d_gun(phi_gun, EG)
    #MaxE_gun = fEgun2d(phi_gun, EG)
    MaxE_gun = fMaxE(EG)
    return MaxE_gun.item()


###

### 2D booster gradient and phase scan and interpolation, updated 01.03.2023
# Somehow, the three columns in booster 2d scan data is phi, grad, and momentum,
# which is different from that for the gun
data_booster = np.loadtxt(field_maps+os.sep+'phi2_scan.dat')
shape = (71, 121)
phi_booster = np.reshape(data_booster[:,0], shape)
E_booster   = np.reshape(data_booster[:,1], shape)
EG_booster  = np.reshape(data_booster[:,2], shape)

# fEG2d_booster = interp2d(E_booster, phi_booster, EG_booster)
fEG2d_booster_t = RegularGridInterpolator((np.unique(E_booster), np.unique(phi_booster)), EG_booster.T,
                                          bounds_error = False)
fEG2d_booster = lambda E, phi: fEG2d_booster_t(([E], [phi]))[0]

def get_MaxE_booster(MaxE_gun = 60, phi_gun = 0, phi_booster = 0, Ef = 17.05):
    
    Eb = np.linspace(10, 22, 121*5)
    EG = fEG2d_booster(Eb, phi_booster)
    fEb_EG = interp1d(EG, Eb)
    
    E1 = fEG2d_gun(MaxE_gun, phi_gun)
    
    MaxE_booster = fEb_EG(Ef-E1)
    return MaxE_booster.item()

###

astra_to_sco = 1.590444459638447
TDS_1MV_ratio = 0.07345637 # TO use this ratio, the input Efield will be the TDS voltage in the unit of MV
ast_fmt = '%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%4d%4d'

Lquad = 0.0675
quads = {'HIGH1.Q1':4.79, 'HIGH1.Q2':5.005, 'HIGH1.Q3':5.6025, 'HIGH1.Q4':5.8525, 'HIGH1.Q5':6.6475,
        'HIGH1.Q6':6.8925, 'HIGH1.Q7':8.18, 'HIGH1.Q8':8.655, 'HIGH1.Q9':10.208, 'HIGH1.Q10':10.388,
        'PST.QM1':12.088, 'PST.QM2':12.468, 'PST.QM3':12.848, 'PST.QT1':13.228, 'PST.QT2':13.608,
        'PST.QT3':13.988, 'PST.QT4':14.368, 'PST.QT5':14.748, 'PST.QT6':15.128, 'HIGH2.Q1':16.635,
        'HIGH2.Q2':16.735, 'HIGH2.Q3':19.587, 'HIGH2.Q4':21.600, 'HIGH2.Q5':23.013, 'HIGH3.Q1':26.350,
        'HIGH3.Q2':26.750, 'HIGH3.Q3':27.150} # 2019.10.05

quads.update({'HIGH2.Q3':19.587, 'HIGH2.Q4':21.600, 'HIGH2.Q5':23.013,
              'HIGH3.Q1':26.350, 'HIGH3.Q2':26.750, 'HIGH3.Q3':27.150})

quads.update({'BIO.Q1':28.10, 'BIO.Q2':28.25, 'BIO.Q3':28.40, 
              'V.Q1':0.23, 'V.Q2':1.5175, 'V.Q3':2.5575, 
              'TEMP.Q1':5.25, 'TEMP.Q2':5.50, 'TEMP.Q3':5.75}) # 2020

quads.update({'HIGH2.Q3':19.587, 'HIGH2.Q4':22.1935, 'HIGH2.Q5':23.0785, 
              'HIGH3.Q1':27.0185, 'HIGH3.Q2':27.4185, 'HIGH3.Q3':27.8185,
              'BACK3.Q1':0.2685, 'BACK3.Q2':0.2685+0.4, 'BACK3.Q3':0.2685+0.4*2}) # 2021.04.21

quads.update({'HIGH2.Q3':19.250, 'HIGH2.Q4':21.460, 'HIGH2.Q5':23.220, 
              'HIGH3.Q1':27.108, 'HIGH3.Q2':27.338, 'HIGH3.Q3':27.778}) # 2021.08.17

# 2022.01.04
quads['BACK.Q1'] = 29.887-1.8-quads['HIGH3.Q3']
quads['BACK.Q2'] = quads['HIGH3.Q3']-quads['HIGH3.Q2']+quads['BACK.Q1']
quads['BACK.Q3'] = quads['HIGH3.Q2']-quads['HIGH3.Q1']+quads['BACK.Q2']

# 2022.02.02
zmatch = 25.293
quads['BACK.QM1'] = zmatch - quads['HIGH2.Q5']
quads['BACK.QM2'] = zmatch - quads['HIGH2.Q4']
quads['BACK.QM3'] = zmatch - quads['HIGH2.Q3']
quads['BACK.QM4'] = zmatch - quads['HIGH2.Q2']

# 2021.10.27
quads['BIO.QM0'] = 0.1
quads['BIO.QM1'] = 0.2
quads['BIO.QM2'] = 0.5
quads['BIO.QM3'] = 0.8

quads['BIO.QM4'] = 2.05
quads['BIO.QM5'] = 2.35
quads['BIO.QM6'] = 2.65

# 2023.04.27
quads['BIO.Q1'] = 0.340
quads['BIO.Q2'] = 0.655
quads['BIO.Q3'] = 1.105
quads['BIO.Q4'] = 1.420


scrns = {'LOW.SCR1':0.8030, 'LOW.SCR2':1.3790, 'LOW.SCR3':1.7080, 'HIGH1.SCR1':5.2770, 'HIGH1.SCR2':6.2500,\
        'HIGH1.SCR3':7.1250, 'HIGH1.SCR4':8.4100, 'HIGH1.SCR5':8.9200, 'PST.SCR1':12.2780, 'PST.SCR2':13.0380,\
        'PST.SCR3':13.7980, 'PST.SCR4':14.5580, 'PST.SCR5':15.3180, 'HIGH2.SCR1':16.3030, 'HIGH2.SCR2':18.2620,\
        'HIGH2.SCR3':22.86, 'HIGH3.SCR1':26.950, 'BACK2.SCR1':0.46, 'BACK3.SCR1':5.227}
    
scrns = {'LOW.SCR1':0.8030, 'LOW.SCR2':1.3790, 'LOW.SCR3':1.7080, 'HIGH1.SCR1':5.2770, 'HIGH1.SCR2':6.2500,\
        'HIGH1.SCR3':7.1250, 'HIGH1.SCR4':8.4100, 'HIGH1.SCR5':8.9200, 'PST.SCR1':12.2780, 'PST.SCR2':13.0380,\
        'PST.SCR3':13.7980, 'PST.SCR4':14.5580, 'PST.SCR5':15.3180, 'HIGH2.SCR1':16.3030, 'HIGH2.SCR2':18.2620,\
        'HIGH2.SCR3':23.450, 'HIGH3.SCR1':27.627, 'HIGH3.UND':28.087, 'BACK2.SCR1':0.46, 'BACK3.SCR1':5.227}

# 2022.07.14
scrns['HIGH3.SCR1'] = 27.558
scrns['HIGH3.SCR2'] = 32.040
scrns['HIGH3.SCR3'] = 33.040

# 2021.10.27
scrns['BIO.SCR1'] = 1.425
scrns['BIO.WINDOW'] = 4.275

# 2021.11.07
scrns['HIGH3.D1_ENT'] = 25.793
scrns['HIGH3.D2_EXIT'] = 27.4269

# Imaging quads after exit window
Ldrift = 0.12
quads['BIO.QI1'] = scrns['BIO.WINDOW']+Ldrift*1+Lquad/2.0
quads['BIO.QI2'] = scrns['BIO.WINDOW']+Ldrift*2+Lquad*1.0+Lquad/2.0
quads['BIO.QI3'] = scrns['BIO.WINDOW']+Ldrift*4+Lquad*2.0+Lquad/2.0
quads['BIO.QI4'] = scrns['BIO.WINDOW']+Ldrift*5+Lquad*3.0+Lquad/2.0

scrns['BIO.WATER'] = scrns['BIO.WINDOW']+Ldrift*6+Lquad*4

scrns['BACK.SCR1'] = quads['BACK.Q3']+(quads['HIGH3.Q1']-scrns['HIGH3.SCR1'])
scrns['BACK.SCR2'] = quads['BACK.Q3']+(quads['HIGH3.Q1']-scrns['HIGH2.SCR3'])

# quads from the exit of High3.D1, the first dipole of the dogleg
l1, l2, l3 = 0.45, 0.379286, 0.15
quads['DOGLEG.Q1'] = l1+Lquad/2
quads['DOGLEG.Q2'] = quads['DOGLEG.Q1']+l2+Lquad
quads['DOGLEG.Q3'] = quads['DOGLEG.Q2']+l3*2+Lquad
quads['DOGLEG.Q4'] = quads['DOGLEG.Q3']+l2+Lquad
scrns['DOGLEG.LEND'] = 2*(l1+Lquad+l2+Lquad+l3)

# quads backwards from around window
quads['BACK.QW1'] = 0.05
quads['BACK.QW2'] = quads['BACK.QW1']+0.1
quads['BACK.QW3'] = quads['BACK.QW2']+0.1


steerers = {'LOW.ST1':0.4920, 'LOW.ST2':0.9630, 'LOW.ST3':1.2700, 
            'LOW.ST4':1.4630, 'LOW.ST5-A':1.9860,
            'HIGH1.ST1':4.8950, 'HIGH1.STA1':5.4270, 'HIGH1.ST2':7.0400, 'HIGH1.STA2':7.2975, 
            'HIGH1.ST3':8.1280, 'HIGH1.ST4':9.8230, 
            'PST.ST1':11.600, 'PST.ST2 ':12.390, 'PST.ST3 ':13.150, 
            'PST.ST4':13.910, 'PST.ST5':14.670, 'PST.ST6':14.981, 
            'HIGH2.ST1':16.453, 'HIGH2.ST2':16.832, 'HIGH2.ST3':19.137, 
            'HIGH2.ST4':21.400, 'HIGH2.ST5':21.972, 
            'HIGH3.ST1':26.250, 'HIGH3.ST2':26.550}

try: # Redefine as case insensitive dictionary
    from requests.structures import CaseInsensitiveDict
    pitz = CaseInsensitiveDict(**scrns, **quads, **steerers)
    scrns = CaseInsensitiveDict(**scrns)
    quads = CaseInsensitiveDict(**quads)
    steerers = CaseInsensitiveDict(**steerers)
except Exception as err:
    print(err)
    pitz = {}
    pitz.update(**scrns, **quads, **steerers)


# Create a folder, mainly for photoinjector optimization
def CreateDir(x, flag = 'injector', **kwargs):
    Ipart = 250000
    if len(kwargs)>0:
        if 'Ipart' in kwargs.keys():
            Ipart = kwargs['Ipart']
            
    if flag.upper() in ['INJECTOR', 'PI']:
        
        sigma_x = sigma_y = x[1]/4.
        phi_gun, phi_booster = x[3], x[5]
        Imain = x[6]

        MaxE_gun = x[2]
        phi_gun, phi_booster = x[3], x[5]
        MaxE_gun = get_MaxE_gun(phi_gun, 6.3)
    
        MaxE_booster = get_MaxE_booster(MaxE_gun, phi_gun, phi_booster, 17)
        
        Q_total = x[0]/1e3
        #Ipart = int(Q_total*50e3)

        direc = str.format('Q-%.2fpC-D-%.2fmm-E1-%.2fMV_m-phi1-%.2fdeg-E2-%.2fMV_m-phi2-%.2fdeg-I-%.2fA' %\
                           (Q_total*1e3, x[1], MaxE_gun, phi_gun, MaxE_booster, phi_booster, Imain))
    
    elif flag.upper() in ['QUADS', 'TRANSPORT']:
        direc = ''
        for i in np.arange(len(x[:])):
            direc += str.format('Q%d-%.2fT_m-' %  (i+1, x[i]))
        direc = direc[:-1]
    elif flag.upper() in ['TWISS', 'MATCHING']:        
        direc = str.format('n%.0fk-sig_x-%.2fmm-sig_y-%.2fmm-alp_x-%.2f-alp_y-%.2f-nemit_x-%.2fum-nemit_y-%.2fum' % (Ipart/1000., *x))
    return direc

def obj_ParaScan(x, *args):
    
    '''
    Created on Nov 17, 2019
    Simulation from photocathode via EMSY1 at 5.277 m to 20 m
    Parameters:
      x: an array or list of the variables to be optimized
    Returns:
      inputs for running Astra simulations on Farm
    '''
    
    Ipart = 200000
    
    Q_total = -x[0]/1e3
    BSA = x[1]
    FWHM = 7e-3 # for Gaussian
    #FWHM = 20e-3 # for Flattop
    
    sigma_x = sigma_y = BSA/4.
    C_sigma_x = C_sigma_y = 0
    
    #sigma_x = sigma_y = 3.0/2.355 # Gaussian truncated
    #C_sigma_x = C_sigma_y = BSA/2./sigma_x
    
    sigma_x = sigma_y = 0.75 #0.96 # Gaussian truncated
    C_sigma_x = C_sigma_y = BSA/2./sigma_x
    
    phi_gun, phi_booster = x[3], x[5]
    
    MaxE_gun = get_MaxE_gun(phi_gun, x[2])
    #MaxE_gun = x[2]
    
    MaxE_booster = get_MaxE_booster(MaxE_gun, phi_gun, phi_booster, x[4])
    #MaxE_booster = x[4]
    
    Imain = x[6]
    MaxB = I2B(Imain)
    
    # path where the field maps are saved
    field_maps = rootdir+os.sep+'sync'+os.sep+'field-maps'
    
    Distribution = 'beam.ini'
    #Distribution = 'ast.0528.011'
    generator = Generator1(FNAME = Distribution, IPart = Ipart, Species = 'electrons', Q_total = Q_total,
                           Ref_Ekin = 0.0e-6, LE = 0.55e-3*1, dist_pz = 'i',
                           Dist_z = 'g', sig_clock = FWHM/2.355, Cathode = True,
                           Dist_x = '2', sig_x = sigma_x, Dist_px = 'g', Nemit_x = 0,
                           Dist_y = '2', sig_y = sigma_y, Dist_py = 'g', Nemit_y = 0,
                           C_sig_x = C_sigma_x, C_sig_y = C_sigma_y) # Gaussian by default
    
    generator.set(Dist_x = 'r', sig_x = BSA/4.0, Dist_y = 'r', sig_y = BSA/4.0) # Transverse uniform
    #generator.set(Dist_z = 'p', Lt = FWHM, Rt = 1e-3) # Temporal flattop
    
    #generator.set(Ref_Ekin = 1e-6, LE = 0, dist_pz = 'g', Nemit_x = 0.17, Nemit_y = 0.18)
    
    Run = 1
    newrun = Module('Newrun', Run = Run, Head = 'PITZ beam line simulation', Distribution = Distribution, CathodeS = True,
                    Auto_Phase = True, Track_All = True, check_ref_part = False, Lprompt = False, Max_step=200000)
    #newrun.set(Xoff = 0.5, Yoff = 0.5)
    #newrun.set(Qbunch = Q_total)
    #newrun.set(Run = 1, XYrms = sigma_x)
    
    charge = Module('Charge', LSPCH = True, Lmirror = True, Nrad = 40, Nlong_in = 50, N_min = 50, Max_scale = 0.05, Max_count = 20)
    charge.set(N_min= 50)
    #charge.set(L2D_3D = True, Z_TRANS = 5, NXF = 32, NYF = 32, NZF = 32, min_grid_trans = 0.03e-3)
    #charge.set(LSPCH3D = True, NXF = 32, NYF = 32, NZF = 32)
    
    cavity = Module('Cavity', LEfield = True, File_Efield = 
                    [field_maps+os.sep+gun_profile, field_maps+os.sep+'CDS14_15mm.txt'],
                    MaxE = [MaxE_gun, MaxE_booster], C_pos = [0., 2.675], Nue = [1.3, 1.3], Phi = [phi_gun, phi_booster])
    
    soleno = Module('Solenoid', LBfield = True, File_Bfield = [field_maps+os.sep+'gunsolenoidsPITZ.txt'],
                    MaxB = [-MaxB], S_pos = [0.], S_xrot = [0*1e-3], S_yrot = [0])
    
    Zstart, Zstop = 0, 5.28
    Zemit = int((Zstop-Zstart)/0.01)
    output = Module('Output', Zstart = Zstart, Zstop = Zstop, Zemit = Zemit, Zphase = 1, RefS = True, EmitS = True,
                    PhaseS = True, TrackS = False, LandFS = True, C_EmitS = True, LPROJECT_EMIT = True,
                    LOCAL_EMIT = False, Screen = [2.455, 5.28])
    #output.set(Zstop = 0.5, Zemit = 50) # emission
    
    apertu = Module('Aperture', LApert = True, File_Aperture = [field_maps+os.sep+'app.txt'])
    
    #Qsol = 0 #-0.0034*x[-2]/365.0
    #Qsol_zrot = 0 # -0.4363 #x[-1]/180.0*np.pi
    #quadru = Module('Quadrupole', LQuad = True, Q_type = ['../Qsol.data'], Q_pos = [0], Q_grad = [Qsol], Q_zrot = [Qsol_zrot])
    
    astra = Astra()
    astra.add_modules([newrun, charge, cavity, soleno, output])
    
    
    direc = str.format('Q-%.2fpC-D-%.2fmm-E1-%.2fMV_m-phi1-%.2fdeg-E2-%.2fMV_m-phi2-%.2fdeg-I-%.2fA' %\
                       (np.abs(Q_total)*1e3, BSA, MaxE_gun, phi_gun, MaxE_booster, phi_booster, Imain))
    
    
    #os.system('mkdir -p '+direc)
    try:
        os.mkdir(direc)
    except Exception as err:
        print(err.errno == 17)
        if err.errno == 17:
            pass
        else:
            return
    
    job_name = direc+'-%03d' % Run
    gen_name = 'gen' #+`Run`
    ast_name = 'ast' #+`Run`
    
    generator.write(direc+os.sep+gen_name+'.in')
    astra.write(direc+os.sep+ast_name+'.in')
    
    astra.submit(job_name, ast_name, gen_name, direc)
    
    return
