from .Plots import *

def nemixrms(x, xp, w = None):
    '''
    Parameters
      x, xp: phase space coordinates. xp could be either divergence or normalized momentum
      w: weighting factor or charge or the particles
    Returns
      emittance: phase space emittance
    '''
    if w is None:
        w = np.ones(len(x))
    x2 = weighted_mean(x*x, w)
    xp2 = weighted_mean(xp*xp, w)
    xxp = weighted_mean(x*xp, w)
    emit_x = np.sqrt(x2*xp2-xxp*xxp)
    return emit_x
nemiyrms = nemixrms

class BeamDiagnostics:
    '''
    Postprocessing class for Astra and/or Warp simulations
    '''
    def __init__(self, fname = None, dist = None, start_index = 0, plot = False, **kwargs):
        '''
        Parameters
          fname: full name of the data file
          dist: 2D arrays that stores the distributions
          start_index: 0 by default, the starting number to refer to the calculated beam parameters
          plot: False by default, if plot current profile
          **kwargs:
            key_value = None
            tracking = 'Astra'
            momentum = 'eV/c'
            em3d = None
        '''
        
        keys = ['avg_z', 'nemit_x', 'nemit_y', 'std_x', 'std_y', 'std_z', 
                'Ekin', 'std_Ekin', 'nemit_z', 'Q_b', 'avg_x', 'avg_y', 
                'alpha_x', 'beta_x', 'gamma_x', 'emit_x', 
                'alpha_y', 'beta_y', 'gamma_y', 'emit_y',
                'cor_Ekin', 'loss_cathode', 'loss_aperture', 
                'FWHM', 'NoP', 'I1', 'I0', 'cov_xxp', 'cov_yyp']
        
        indices = np.arange(len(keys))
        units = ['m', 'um', 'um', 'mm', 'mm', 'mm', 'MeV', 'keV', 'keV mm', 'C',
                 'm', 'm', ' ', 'm', 'm^-1', 'm', ' ', 'm', 'm^-1', 'm', 'keV', ' ',
                 ' ', 'mm', ' ', 'A', 'A', 'mmmrad', 'mmmrad']
        keyIndex = {}
        keyUnit = {}
        for i, key in enumerate(keys):
            keyIndex[key] = indices[i] + start_index
            keyUnit[key] = units[i]
            
        self.keyIndex = keyIndex
        self.keyUnit = keyUnit
        
        kwargs.update(plot = plot)
        if fname is not None:
            try:
                import pandas as pd
                dist = pd.read_csv(fname, delimiter = '\s+', header = None).values
            except Exception as err:
                print(err)
                dist = np.loadtxt(fname)
                
            if 'tracking' in kwargs.keys():
                if kwargs['tracking'].upper() in ['ASTRA']:
                    dist[1:,2] += dist[0,2]; dist[1:,5] += dist[0,5]
                else:
                    pass
            else:
                dist[1:,2] += dist[0,2]; dist[1:,5] += dist[0,5]
            
            self.dist = dist
            self.diagnostics(**kwargs)
        elif dist is not None:
            self.dist = dist
            self.diagnostics(**kwargs)
        else:
            self.x = np.zeros(len(self.keyIndex))

    def __getitem__(self, key):
        return self.get(key)
    
    def __getattr__(self, key):
        return self.get(key)
    
    @property
    def x(self):
        return self.__x

    @x.setter
    def x(self, x):
        self.__x = x
        
    def get(self, key):
        '''
        Return the corresponding parameter if key is one of `self.kv.keys()`, or return the index of the key 
        if requesting by "idx_"+key and key is one of `self.kv.keys()`
        '''
        if key in self.keyIndex.keys():
            index = self.keyIndex[key]
            return self.x[index]
        elif key in ['idx_'+ i for i in self.keyIndex.keys()]:
            index = self.keyIndex[key[4:]]
            return index
        else:
            return -999
    
    def get_index(self, key):
        index = self.keyIndex[key]
        return index
        
    def diagnostics(self, key_value = None, tracking = 'Astra',
                    momentum = 'eV/c', energy = True,
                    em3d = None, plot = False, fig_ext = '.eps'):
        '''
        Parameters
          x y z Px Py Pz w: 6D phase space of the particles
        Outputs
          [avg_z nemit_x nemit_y std_x std_y std_z Ekin std_Ekin nemit_z Q_b\
          avg_x avg_y alpha_x beta_x gamma_x emit_x alpha_y beta_y gamma_y emit_y\
          cor_zEk loss_cathode loss_aperture FWHM NoP I1 I0]: array
        '''
        
        tracking = tracking.upper()
        if tracking in ['ASTRA', 'A']:
            ix, iy, iz, iux, iuy, iuz, iw = 0, 1, 2, 3, 4, 5, 7
            weight_to_charge = 1e-9
        elif tracking in ['WARP', 'W']:
            ix, iy, iz, iux, iuy, iuz, iw = 0, 1, 2, 3, 4, 5, 6
            weight_to_charge = g_qe
            default_momentum = 'Dimentionless'
        else:
            ix, iy, iz, iux, iuy, iuz, iw = 0, 1, 2, 3, 4, 5, 6
            weight_to_charge = 1.0
        
        if key_value != None:
            keys = key_value.keys()
            if 'x' in keys: ix = key_value['x']
            if 'y' in keys: iy = key_value['y']
            if 'z' in keys: iz = key_value['z']
            if 'ux' in keys: iux = key_value['ux']
            if 'uy' in keys: iuy = key_value['uy']
            if 'uz' in keys: iuz = key_value['uz']
            if 'w' in keys: iw = key_value['w']
        
        dist = np.copy(self.dist)
        NoP = len(dist[:,1])
        
        if tracking in ['ASTRA', 'A']:
            select = (dist[:,9]==-89); n_lost_cathode = np.sum(select)
            select = (dist[:,9]==-15); n_lost_aperture = np.sum(select)

            select = (dist[:,9]>0) # active particles
            if np.sum(select) == 0:
                select = (dist[:,9]==-1); 
                if np.sum(select) == 0:
                    self.x = np.zeros(len(self.keyIndex))
                    return
                else:
                    dist[:,iz] = dist[:,6]*1e-9*g_c

            dist = dist[select]
        else:
            n_lost_cathode = -1
            n_lost_aperture = -1
        
        
        x, y, z = dist[:,ix], dist[:,iy] ,dist[:,iz] # positions in meter
        
        momentum = momentum.upper()
        if momentum == 'EV/C':
            bgx, bgy, bgz = dist[:,iux]/g_mec2/1e6, dist[:,iuy]/g_mec2/1e6, dist[:,iuz]/g_mec2/1e6 # dimentionless momentum
        elif momentum == 'DIMENTIONLESS':
            bgx, bgy, bgz = dist[:,iux], dist[:,iuy], dist[:,iuz]
        elif momentum == 'M/S':
            bgx, bgy, bgz = dist[:,iux]/g_c, dist[:,iuy]/g_c, dist[:,iuz]/g_c
            gamma = 1.0/np.sqrt(1.0-bgx**2-bgy**2-bgz**2)
            bgx, bgy, bgz = bgx*gamma, bgy*gamma, bgz*gamma
        
        
        w = dist[:,iw]*weight_to_charge # C
        
        if energy:
            Ek = np.sqrt(1+bgx**2+bgy**2+bgz**2)*g_mec2-g_mec2 # kinetic energy in MeV
        else:
            Ek = np.sqrt(bgx**2+bgy**2+bgz**2)*g_mec2 # P*c in MeV/c
            
        Qtot = np.sum(w) # C
        
        x2, y2, z2 = np.apply_along_axis(np.cov, arr = dist[:,0:3], axis = 0, aweights = np.abs(w), bias = True)
        std_x, std_y, std_z = np.sqrt([x2, y2, z2])
        
        xc, yc, zmean = np.apply_along_axis(weighted_mean, arr = dist[:,0:3], axis = 0, weights = w)
        
        xp = bgx/bgz; yp = bgy/bgz
        xp2 = weighted_cov(xp, xp, w); yp2 = weighted_cov(yp, yp, w)
        xxp = weighted_cov(x, xp, w);  yyp = weighted_cov(y, yp, w)
        
        emit_x = np.sqrt(x2*xp2-xxp*xxp); emit_y = np.sqrt(y2*yp2-yyp*yyp)
        alpha_x = -xxp/emit_x; beta_x = x2/emit_x; gamma_x = xp2/emit_x
        alpha_y = -yyp/emit_y; beta_y = y2/emit_y; gamma_y = yp2/emit_y
        
        fwhm = get_FWHM(z)
        [I1, I0] = zdist(z, weights = w, plot = plot, fig_ext = fig_ext)
        
        select = (z-zmean>-1.*std_z)*(z-zmean<1.*std_z); Ek1 = Ek[select]; z1 = z[select]; w1 = w[select]
        Ek1mean = weighted_mean(Ek1, w1); z1mean = weighted_mean(z1, w1); std_z1 = weighted_std(z1, w1);
        cor_zEkin = weighted_mean(Ek1*z1, w1)-z1mean*Ek1mean; # MeV/c*m
        cor_Ekin = cor_zEkin/std_z1**2; # MeV/c/m
        
        Ek2 = Ek-Ek1mean-cor_Ekin*(z-z1mean) # higher order energy spread
        std_Ek2 = weighted_std(Ek2, w)*1e3
        #import pdb; pdb.set_trace()
        
        if plot:
            fig, [ax1, ax2] = plt.subplots(ncols = 2, figsize = (6, 3))
            ax1.plot((z[::979]-zmean)*1e3, Ek[::979], '.')
            ztmp = np.linspace(-5, 5)*1e-3
            Ektmp = Ek1mean+cor_Ekin*ztmp
            ax1.plot(ztmp*1e3, Ektmp, '-')
        
            ax2.plot((z[::979]-zmean)*1e3, Ek2[::979], '.')
            fig.tight_layout()
            fig.savefig('Ek-z.png')
        
        
        if em3d is not None:
            Bz = em3d.get_field(x, y, z)[:,5]
            bgx += 0.5*Bz*y*g_c/g_mec2/1e6
            bgy -= 0.5*Bz*x*g_c/g_mec2/1e6
            
        bgx2 = weighted_cov(bgx, bgx, w); bgy2 = weighted_cov(bgy, bgy, w)
        xbgx = weighted_cov(x, bgx, w);   ybgy = weighted_cov(y, bgy, w)
        
        nemit_x = np.sqrt(x2*bgx2-xbgx**2); nemit_y = np.sqrt(y2*bgy2-ybgy**2)
        
        x = np.array([zmean, nemit_x*1e6, nemit_y*1e6, std_x*1e3, std_y*1e3, std_z*1e3, 
                      weighted_mean(Ek, w), weighted_std(Ek, w)*1e3, nemixrms(z, Ek, w)*1e6, 
                      Qtot, xc, yc,
                      alpha_x, beta_x, gamma_x, emit_x, alpha_y, beta_y, gamma_y, emit_y,
                      cor_Ekin, n_lost_cathode, n_lost_aperture, fwhm*1e3, NoP, I1, I0, xxp*1e6, yyp*1e6, std_Ek2])
        
        self.x = x
    
    def demo(self, basename = None):
        if basename is not None:
            filename = 'diag@'+basename+'.txt'
        else:
            filename = 'diag@'+'.txt'
        
        a = sorted(self.keyIndex.items(), key=lambda x: x[1])
        ss = ''
        for i in np.arange(len(a)):
            ss += str('%17s: %15.6E %s\n' % (a[i][0], self.x[i], self.keyUnit.get(a[i][0])))
        print(ss)

        ff = open(filename, 'w')
        ff.write(ss)
        ff.close()
        return 

def astra_demo(fname):
    pass

def warp_post(fname = None, dist = None, **kwargs):
    '''
    Parameters
      x y z bgx bgy bgz w: 6D phase space of the particles
    Outputs
      see class `BeamDiagnostics`
    Examples
    
      data = pd_loadtxt('ast.xxxx.001')
      data[1:,2] += data[0,2]
      data[1:,5] += data[0,5]
      astra_post(dist = data)
      
      or 
      
      astra_post('ast.xxxx.001')
    '''
    
    kwargs.update(tracking = 'Warp')
    kwargs.update(momentum = 'Dimentionless')
    
    diag = BeamDiagnostics(fname = fname, dist = dist, **kwargs)
    return diag.x

def astra_post(fname = None, dist = None, **kwargs):
    '''
    Parameters
      x y z Px Py Pz t w: 6D phase space of the particles
    Outputs
      see class `BeamDiagnostics`
    Examples
    
      data = pd_loadtxt('ast.xxxx.001')
      warp_post(dist = data)
      
      or 
      
      warp_post('ast.xxxx.001')
    '''
    
    kwargs.update(tracking = 'Astra')
    kwargs.update(momentum = 'eV/c')
    
    diag = BeamDiagnostics(fname = fname, dist = dist, **kwargs)
    return diag.x

### Functions to access the default statistics from Astra simulation
def get_nemixrms(fname = 'ast.Xemit.001'):
    xemit = np.loadtxt(fname)
    return xemit[-1,5]
def get_nemixrms_at(z, fname = 'ast.Xemit.001'):
    xemit = np.loadtxt(fname)
    fxemit = interp1d(xemit[:,0], xemit[:,5])
    return fxemit(z)

def get_nemiyrms(fname = 'ast.Yemit.001'):
    yemit = np.loadtxt(fname)
    return yemit[-1,5]
def get_nemiyrms_at(z, fname = 'ast.Yemit.001'):
    yemit = np.loadtxt(fname)
    fyemit = interp1d(yemit[:,0], yemit[:,5])
    return fyemit(z)

def get_ref_momentum(fname = 'ast.ref.001'):
    ref = np.loadtxt(fname)
    return ref[-1,2]

### Functions to calculate statistics using particle distributions from Astra simulation 
def cal_momentum(fname = 'beam.ini'):
    beam = np.loadtxt(fname)
    beam[1:,5] += beam[0,5]
    return np.mean(np.sqrt(beam[:,3]**2+beam[:,4]**2+beam[:,5]**2))
def cal_kinetic(fname = 'beam.ini'):
    momentum = cal_momentum(fname)
    return momentum2kinetic(momentum)

def cal_nemixrms(fname = 'beam.ini'):
    beam = np.loadtxt(fname)
    return nemixrms(beam[:,0], beam[:,3])/g_mec2/1e6*1e6
def cal_nemiyrms(fname = 'beam.ini'):
    beam = np.loadtxt(fname)
    return nemixrms(beam[:,1], beam[:,4])/g_mec2/1e6*1e6