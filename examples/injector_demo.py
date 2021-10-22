# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:45:16 2020

@author: lixiangk
"""

from interface import *
from scipy.interpolate import interp1d, interp2d

field_maps = os.path.join('.', 'field-maps')

I2B = lambda I: -(0.0000372+0.000588*I)
B2I = lambda B: (-B-0.0000372)/0.000588

###### Prepare the interpolation function for the get_MaxE_booster()
data_gun = np.loadtxt(field_maps+os.sep+'phiGun42_scan.dat')
fEG_gun = interp1d(data_gun[:,0], data_gun[:,1])

data_booster = np.loadtxt(field_maps+os.sep+'phiBooster_scan.dat')
EG_booster = np.reshape(data_booster[:,2], (71, 121)) #-6.2678

phi_booster = np.linspace(-40, 30, 71)
E_booster   = np.linspace(10, 22, 121)

fEG_booster = interp2d(E_booster, phi_booster, EG_booster)

def get_MaxE_booster(MaxE_gun = 60, phi_gun = 0, phi_booster = 0, Ef = 17.05):

    Eb = np.linspace(10, 22, 121*5)
    EG = fEG_booster(Eb, phi_booster)
    fEb_EG = interp1d(EG, Eb)

    MaxE_booster = fEb_EG(Ef-fEG_gun(phi_gun))
    return np.asscalar(MaxE_booster)
######

field_maps = os.path.join('..', 'field-maps')
def obj_THz4nC_MOGA(x):
    '''
    Created on November 13, 2019
    Simulation from photocathode via EMSY1 at 5.277 m to 20 m
    Goal function is a combination of average emittance and correlated energy spread exterpolated at undulator center, 29.25 m
    Laser spot fixed at 4 mm, gun phase at MMMG
    Variables to be optimized: gun and booster phases, solenoid current
    Parameters:
      x: an array or list of the variables to be optimized
    Returns:
      energy: the quantity that is envisoned as the "energy" of the sample
    '''

    Ipart = 1000
    BSA = x[0]
    sigma_x = sigma_y = BSA/4.
    C_sigma_x = C_sigma_y = 0
    
    #sigma_x = sigma_y = 3.0/2.355 # Gaussian truncated
    #C_sigma_x = C_sigma_y = BSA/2./sigma_x
    
    phi_gun, phi_booster = x[1], x[2]
    Imain = x[3]
    
    MaxE_gun = 60
    MaxE_booster = get_MaxE_booster(MaxE_gun, phi_gun, phi_booster, 17.0)
    MaxB = I2B(Imain)

    Q_total = -4.0
    
    generator = Generator1(FNAME = 'beam.ini', IPart = Ipart, Species = 'electrons',
                           Q_total = Q_total, Cathode = True,
                           Ref_Ekin = 0.0e-6, LE = 0.55e-3, dist_pz = 'i',
                           Dist_z = 'p', Lt = 21.5e-3, rt = 2e-3, 
                           Dist_x = 'r', sig_x = sigma_x, C_sig_x = C_sigma_x,
                           Dist_px = 'g', Nemit_x = 0,
                           Dist_y = 'r', sig_y = sigma_y, C_sig_y = C_sigma_y,
                           Dist_py = 'g', Nemit_y = 0)
    
    newrun = Module('Newrun', Run = 1, Head = 'PITZ beamline simulation',
                    Distribution = 'beam.ini', Qbunch = Q_total,
                    Auto_Phase = True, Track_All = True, check_ref_part = False,
                    Lprompt = False, Max_step=200000)
    newrun.set(Run = 1)
    
    charge = Module('Charge', LSPCH = True, Lmirror = True, 
                    Nrad = 40, Nlong_in = 50, N_min = 50, 
                    Max_scale = 0.05, Max_count = 20)
    charge.set(LSPCH = False, Lmirror = False)
    
    File_Efield = [field_maps+os.sep+fp for fp in ['gun42cavity.txt', 'CDS14_15mm.txt']]
    cavity = Module('Cavity', LEfield = True, File_Efield =  File_Efield,
                    MaxE = [MaxE_gun, MaxE_booster], C_pos = [0., 2.675],
                    Nue = [1.3, 1.3], Phi = [phi_gun, phi_booster])
    
    soleno = Module('Solenoid', LBfield = True, 
                    File_Bfield = [field_maps+os.sep+'gunsolenoidsPITZ.txt'], 
                    MaxB = MaxB, S_pos = [0.])
    
    output = Module('Output', Zstart = 0, Zstop = 20, Zemit = 400, Zphase = 1,
                    RefS = True, EmitS = True, PhaseS = True, TrackS = False, 
                    LandFS = True, C_EmitS = True, LPROJECT_EMIT = True,
                    LOCAL_EMIT = True, Screen = [5.28])
    apertu = Module('Aperture', LApert = True,
                    File_Aperture = [field_maps+os.sep+'app.txt'])
    
    
    #astra = Astra(newrun, charge, cavity, soleno, output)
    astra = Astra()
    astra.add_modules([newrun, charge, cavity, soleno, output])
    
    direc = str.format('Q%.1fpC-D-%.2fmm-E1-%.2fMV_m-phi1-%.2fdeg-E2-%.2fMV_m-phi2-%.2fdeg-I-%.2fA' %
                       (Q_total*1e3, BSA, MaxE_gun, phi_gun, MaxE_booster, phi_booster, Imain))
    #os.system('mkdir -p '+direc)
    try: 
        os.mkdir(direc) 
    except OSError as error: 
        print(error)  
        
    os.chdir(direc)
    
    astName = 'ast.in'
    genName = 'gen.in'
    jobName = 'myjob.%s.sh' % direc
    
    generator.write(genName)
    astra.write(astName)
    
    #### Either submit to Farm
    #os.chdir('..')
    #astra.qsub(job_basename, ast_basename, gen_basename, direc, submit = True)
    #return 0
    ####
    
    ### or run directly
    os.system('generator %s 2>&1 | tee gen.log' % genName)
    os.system('astra %s 2>&1 | tee ast.log' % astName)
    
    #print(x)
    # Postprocessing and calculate the goal functions
    try:
        
        fname = 'ast.Xemit.001'
        xemit = np.loadtxt(fname)
        select = (xemit[:,0]>5.0)
        xemit = xemit[select]
        nemit_x_std = np.std(xemit[:,5])
        nemit_x_avg = np.mean(xemit[:,5])
        
        
        fname = 'ast.Zemit.001'
        zemit = np.loadtxt(fname)       
        select = (zemit[:,0]>5.0)
        zemit = zemit[select]
        
        f = lambda x, a, b: a+b*x
        popt, pcov = curve_fit(f, zemit[:,0], zemit[:,6])
        cov = f(29.25, *popt)
        
        fname = 'ast.0528.001'
        diag1 = BeamDiagnostics(fname = fname)
        
        fname = 'ast.2000.001'
        diag2 = BeamDiagnostics(fname = fname)
        I1 = diag2.I1
        
        if diag1.loss_cathode<=500 and diag1.loss_aperture<=500: 
            # if there is little beam loss to the cathode or to the aperture
            obj = [4*nemit_x_avg*nemit_x_std, (np.fabs(cov+50.0)/50.0)**2, 200.0/I1]
        else:
            # if there is too much beam loss, then set the objective to 1000
            obj = [1e3, 1e3, 1e3]   
        
        # Save intermediate diagnosing results
        res = [Ipart, BSA, MaxE_gun, phi_gun, MaxE_booster, phi_booster, Imain]+list(diag1.x)
        with open('..'+os.sep+'diag@5.28m.dat','a') as f_handle:
            np.savetxt(f_handle, np.atleast_2d(res), fmt='%14.6E')

        res = [Ipart, BSA, MaxE_gun, phi_gun, MaxE_booster, phi_booster, Imain]+list(diag2.x)
        with open('..'+os.sep+'diag@20m.dat','a') as f_handle:
            np.savetxt(f_handle, np.atleast_2d(res), fmt='%14.6E')
        
    except:
        # if no output beam file, then set the objective to 10000
        nemit_x, cov, I1 = 999, 999, 999
        obj = [nemit_x, cov, I1]
        print('Error: not simulated well!')
    
    os.chdir('..'); 
    
    #removing intermediate results
    #os.system('rm -r '+direc)
    
    # Save parameters and goal functions in another file
    res = [Ipart, BSA, MaxE_gun, phi_gun, MaxE_booster, phi_booster, Imain, nemit_x_avg, nemit_x_std, cov, I1]+obj
    with open('results.dat','a') as f_handle:
        np.savetxt(f_handle, np.atleast_2d(res), fmt='%14.6E')
    
    print(obj)
    return obj

#obj_THz4nC_MOGA[4, 0, -20, 380])