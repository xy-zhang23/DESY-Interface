# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:45:16 2020

@author: lixiangk
"""

from interface import *

def RotationMatrix(theta):
    cc = np.cos(theta)
    ss = np.sin(theta)
    return np.array([[cc, -ss], [ss, cc]])

def RotatedBy(X, theta, Xc = [0, 0]):
    M = RotationMatrix(theta)
    X = np.reshape(X, (len(X), 1))
    Xc = np.reshape(Xc, (len(Xc), 1))
    return M @ (X-Xc)+Xc

def PrepareDipoleCoordinates(rho, theta, poleWidth = 0, X0 = [0, 0],
                             theta0 = 0, beta1 = 0, beta2 = 0):
    '''
    Return four pairs of coordinates for the four corners of the pole face, in the form of 
    (x, z) in unit of meter.
    |--
    |                 
    A1             |
    |            B1  
    A          |
    |        B
    A2     | 
    |   B2
    | |
    |
    '''
    if poleWidth == 0:
        poleWidth = rho
    width = poleWidth/2
    sg = np.sign(theta)
    
    X0 = np.reshape(X0, (2, 1))
    Xa = np.reshape([0, 0], (2, 1))
    Xc = np.reshape([0, rho*sg], (2, 1))
    
    Xa1 = np.reshape([0,  width], (2, 1))+Xa
    Xa2 = np.reshape([0, -width], (2, 1))+Xa
    
    Xb = RotatedBy(Xa-Xc, theta)+Xc
    Xb1 = RotatedBy(Xa1-Xc, theta)+Xc
    Xb2 = RotatedBy(Xa2-Xc, theta)+Xc
    
    Xa1 = RotatedBy(Xa1-Xa, beta1)+Xa
    Xa2 = RotatedBy(Xa2-Xa, beta1)+Xa
    
    Xb1 = RotatedBy(Xb1-Xb, beta2)+Xb
    Xb2 = RotatedBy(Xb2-Xb, beta2)+Xb
    
    #
    V = np.concatenate((Xa1, Xa2, Xb1, Xb2), axis=1) # 2x4
    M = RotationMatrix(theta0) # 2x2
    V = (M@V).T # (2x4).T -> 4x2
    
    V[:,0] += X0[1]
    V[:,1] += X0[0]
    
    return [vi[::-1] for vi in V]

def obj_dogleg(x, *args, **kwargs):
    '''
    Photocathode to the EMSY1 at 5.277 m
    Parameters:
      x: an array or list of the variables to be optimized, here the gradients of quads
      *args: the names of quads
      **kwargs: key-value pairs of Astra module property, e.g., 'Run', 1
    Returns:
      direc: the directory where Astra simulation runs
    '''
    
    z0 = 26
    lq = 0.0675
    
    L, rho, theta = 1.5, 0.4, np.pi/3
    
    poleWidth = 0.36
    gap = 50e-3
    Bdip = 0.1077634*1.5
    
    l1, l2, l3, g1, g2, rr = x
    
    Run = 10
    Distribution = '..'+os.sep+'ast.2550.007'
    Zstop = 28
    Screen = [27.5, 27.8]
    if len(kwargs)>0:
        if 'Run' in list(kwargs.keys()):
            Run = kwargs['Run']

        if 'Distribution' in list(kwargs.keys()):
            Distribution = kwargs['Distribution']

        if 'Zstop' in list(kwargs.keys()):
            Zstop = kwargs['Zstop']

        if 'Screen' in list(kwargs.keys()):
            Screen = kwargs['Screen']

    newrun = Module('Newrun', Run = Run, Head = 'PITZ beam line simulation', Distribution = Distribution,\
                    Auto_Phase = True, Track_All = True, check_ref_part = False, Lprompt = False, Max_step=200000)
    #newrun.set(Run = Run)

    charge = Module('Charge', LSPCH = True, Lmirror = True, Nrad = 30, Nlong_in = 50, N_min = 10,\
                    NXF = 32, NYF = 32, NZF = 32, L2D_3D = True, Z_TRANS = 5.0, Max_scale = 0.05, Max_count = 20)
    #charge.set(LSPCH = False, L2D_3D = False)
    #charge.set(NXF = 16, NYF = 16, NZF = 16)
    
    output = Module('Output', Zstart = 0, Zstop = Zstop, Zemit = int(Zstop*50), Zphase = 1, RefS = True, EmitS = True,\
                    PhaseS = True, TrackS = False, LandFS = True, Screen = Screen)
    
    # Quadrupole
    Xb = np.array([rho*(np.cos(theta)-1), z0+rho*np.sin(theta)]) # (x, z)
    L0 = (L+Xb[0]*2)/np.sin(theta);
    
    L1 = np.array([l1+lq/2, l1+lq+l2+lq/2, L0/2+l3+lq/2, L0-l1-lq/2]); print(L1, L0, L)
    Qs = np.array([Xb[0]-L1*np.sin(theta), Xb[1]+L1*np.cos(theta)]) # (x, z)
    
    quadru = Module('Quadrupole', LQuad = True, 
                    Q_type = ['../Q3.data' for _ in np.arange(4)],
                    Q_pos = Qs[1], 
                    Q_grad = [g1*rr, g2, g2, g1*rr],
                    Q_xoff = Qs[0],
                    Q_xrot = [-theta for _ in np.arange(4)])
    
    # Dipole
    X11, X12, X13, X14 = PrepareDipoleCoordinates(rho, -theta, poleWidth = poleWidth, theta0 = 0, X0 = [0, z0]) # (x, z)
    
    X0 = [Xb[0]-L0*np.sin(theta), Xb[1]+L0*np.cos(theta)] # (x, z)
    X21, X22, X23, X24 = PrepareDipoleCoordinates(rho, theta, poleWidth = poleWidth, theta0 = -theta, X0 = X0)
    
    D_type = ['hor', 'hor']; D_strength = [-Bdip, Bdip]; D_radius = [rho, -rho]
    D1 = [X11, X21]; D2 = [X12, X22]; D3 = [X13, X23]; D4 = [X14, X24]; #print(D1)
    
    D_Gap = [[gap, gap], [gap, gap]]
    
    dipole = Module('Dipole', LDipole = True, D_type = D_type, D1 = D1, D2 = D2, D3 = D3, D4 = D4,
                    D_Gap = D_Gap, D_strength = D_strength, D_radius = D_radius)
                    
    
    astra = Astra()
    astra.add_modules([newrun, charge, quadru, dipole, output])
    
    direc = ''
    for i in np.arange(len(x)):
        direc += str.format('%.2f-' %  (x[i]))
    direc = direc[:-1]
    
    ###
    
    #os.system('mkdir -p '+direc)
    try: 
        os.mkdir(direc) 
    except OSError as error: 
        print(error)   
    
    #os.chdir(direc)
    
    job_name = 'myjob_'+direc
    #gen_name = 'gen' #+`Run`
    ast_name = 'ast' #+`Run`
    
    #generator.write(direc+os.sep+gen_name+'.in')
    astra.write(direc+os.sep+ast_name+'.in'); print(direc+os.sep+ast_name+'.in')
    astra.qsub(job_name, ast_name, None, direc)
    
    return direc

obj_dogleg([0.5, 0.3, 0.1, 1, -0.5, 1])