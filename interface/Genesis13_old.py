# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 18:45:15 2020

@author: lixiangk
"""
from universal import *
from numpy.fft import *

class Genesis:
    '''
    Class for Genesis 1.3, obsolete
      1. Producing Genesis 1.3 input for simulation
      2. Read Genesis 1.3 output (raw format) into memories and
         1) get the radiation powers/energies along the undulator;
         2) get the temporal/spectrum distribution at given positions
    '''
    
    def __init__(self, **kwargs):
        self.name = 'newrun'
        
        undulator = {'aw0':1.0, 'xlamd':3.0e-2, 'nwig':120, 'nsec':1, 'xkx':0, 'xky':1, 'iwityp':0, 'awd':1.0}
        electron  = {'npart':8192, 'gamma0':20, 'delgam':0.5e-2, 'rxbeam':1e-3, 'rybeam':1e-3, 'emitx':2e-6,\
                     'emity':2e-6, 'alphax':0, 'alphay':0, 'curpeak':1e2, 'bunch':0}
        radiation = {'xlamds':100e-6, 'prad0':0, 'zrayl':0.1, 'zwaist':0}
        grid      = {'ncar':201,'rmax0':9, 'dgrid':0, 'nscz':1, 'nscr':1, 'nptr':40}
        control   = {'delz':1, 'zstop':0, 'nbins':4, 'iall':0, 'ipseed':-1, 'shotnoise':1}
        time      = {'itdp':0, 'curlen':1e-3, 'zsep':1, 'nslice':408, 'ntail':-253, 'iotail':0}
        io        = {'iphsty':1, 'ishsty':1, 'ippart':0, 'ispart':0, 'ipradi':0, 'isradi':0, 'idump':0, 'idmpfld':0,\
                     'idmppar':0, 'ffspec':0}
        lout      = {'lout':[1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1]}
        focusing  = {'quadf':0}
        
        modules = [undulator, electron, radiation, grid, control, time, io, lout, focusing]
        self.kv = {}
        
        for module in modules:
            self.kv.update(module)
        
        for k, v in kwargs.items():
            self.kv.update({k:v})
        self.build()
        return
    def set(self,**kwargs):
        '''
        Set Genesis 1.3 input, `**kwargs` should be in the list of input parameters defined in the manual
        '''
        for k, v in kwargs.items():
            self.kv.update({k:v})
        self.build()
        return
    def reset(self,**kwargs):
        '''
        Reset Genesis 1.3 input, `**kwargs` should be in the list of input parameters defined in the manual
        '''
        for k, v in kwargs.items():
            self.kv.update({k:v})
        self.build()
        return
    def build(self):
        '''
        Prepare for output by iterating through the `**kwargs`
        '''
        output = ' $'+self.name.upper()+' \n'
        for k, v in self.kv.items():
            if isinstance(v, (int, np.int8, np.int16, np.int32, np.int64)) and not isinstance(v, bool):
                output += str.format(' %-9s = %d\n' % (k.upper(), v))
            elif isinstance(v, (float, np.float32, np.float64)):
                output += str.format(' %-9s = %f\n' % (k.upper(), v))
            elif isinstance(v, str):
                output += str.format(' %-9s = \'%s\'\n' % (k.upper(), v))
            elif isinstance(v, bool):
                output += str.format(' %-9s = %s\n' % (k.upper(), str(v))) 
            elif isinstance(v, (list, np.ndarray)):
                output += str.format(' %-9s = ' % k.upper())
                for v1 in v:
                    output += str.format('%s ' % (str(v1)))
                output += '\n'
        output += ' $END\n'
        self.output = output
        return
    def readdata(self, fname):
        '''
        Read the standard output (in raw format and stored in `fname`) into memories for further analysis
        '''
        data=open(fname, 'r')
        line=data.readline()

        kv = {}

        while True:
            if not line:
                break

            match = re.search(r'^[^=]*=[^=]*$' , line)
            if match:
                # print line
                key, value = line.split('=')
                match = re.search(r'file', key)
                if match:
                    kv.update({key.strip():value})
                    line=data.readline()
                    continue

                if len(value.split()) == 1:
                    value = np.float(value.replace("D", "E"))
                else:
                    value = [np.float(v.replace("D", "E")) for v in value.split()]
                kv.update({key.strip():value})          

            match = re.search(r'^[ ]*z\[m\][ ]*aw[ ]*qfld[ ]*$', line)
            if match:
                # print line
                nz = np.int(kv['zstop']/kv['delz']/kv['xlamd'])+1
                field = np.zeros((nz, 3))
                for i in np.arange(nz):
                    line = data.readline()
                    field[i] = np.array([np.float(v.replace("D", "E")) for v in line.split()])
                break
            line=data.readline()

        nt = np.int(kv['nslice'])
        nc = np.int(np.sum(kv['lout'])); #print nc

        current = np.zeros((nt, 2))
        results = np.zeros((nt, nz, nc)); #print results.shape

        islice = 0
        while True:
            if not line:
                break
            match = re.search(r'[ E]* current', line)
            if match:
                icurrent, tmp = line.split()
                icurrent = np.float(icurrent)
                current[islice, 0] = islice+1; current[islice, 1] = icurrent

        
                line = data.readline()
                line = data.readline()
                line = data.readline()

                outputs = line.split()

                for j in np.arange(nz):
                    line = data.readline(); #print line
                    line = line.replace('NaN', '  0')
                    results[islice, j] = np.array([np.float(v.replace("D", "E")) for v in line.split()])
                #break
                islice += 1
            line = data.readline()
        self.kv = kv
        self.nt = nt; self.nz = nz; self.nc = nc
        self.field = field
        self.current = current
        self.data = results
        return
    
    def get_slice(self, i):
        return self.data[i,:,:]
    
    def plot_current(self):
        tmp = self.current
        
        fig, ax = plt.subplots()
        ax.plot(tmp[:,0], tmp[:,1], '-*')
        ax.set_xlabel(r'# of slice')
        ax.set_ylabel(r'current (A)')
        
        
    def get_power_at(self, z = None):
        '''
        Get the radiation power distribution at position z
        Parameters
          z: longitudinal position along the undulator
        Returns
          a 2D array that stores the radition temporal distribution at position z
        '''
    
        if z == None:
            z = self.kv['zstop']
        col = np.int(z/self.kv['delz']/self.kv['xlamd'])
        
        nt = self.nt
        dist = np.zeros((nt, 2))
        
        period_s = self.kv['xlamds']/g_c
        dist[:,0] = np.arange(nt)*period_s*g_c # meter
        dist[:,1] = self.data[:, col, 0]
        
        return dist
    def get_power_density_at(self, z = None):
        '''
        Get the radiation power density at position z
        Parameters
          z: longitudinal position along the undulator
        Returns
          a 2D array that stores the radition power density temporal distribution at position z
        '''
    
        if z == None:
            z = self.kv['zstop']
        col = np.int(z/self.kv['delz']/self.kv['xlamd'])
        
        nt = self.nt
        dist = np.zeros((nt, 2))
        
        period_s = self.kv['xlamds']/g_c
        dist[:,0] = np.arange(nt)*period_s*g_c # meter
        dist[:,1] = self.data[:, col, 2]
        
        return dist
    def get_radiation_size_at(self, z = None):
        '''
        Get the radiation power density at position z
        Parameters
          z: longitudinal position along the undulator
        Returns
          a 2D array that stores the radition power density temporal distribution at position z
        '''
    
        if z == None:
            z = self.kv['zstop']
        col = np.int(z/self.kv['delz']/self.kv['xlamd'])
        
        nt = self.nt
        dist = np.zeros((nt, 2))
        
        period_s = self.kv['xlamds']/g_c
        dist[:,0] = np.arange(nt)*period_s*g_c # meter
        dist[:,1] = self.data[:, col, 4]
        
        return dist
    def get_spectrum_at(self, z = None):
        '''
        Get the radiation spectrum at position z
        Parameters
          z: longitudinal position along the undulator
        Returns
          a 2D array that stores the radition spectrum at position z
        '''
        
        if z == None:
            z = self.kv['zstop']
        col = np.int(z/self.kv['delz']/self.kv['xlamd'])
        
        nt = self.nt
        dist = np.zeros((nt, 2))
        
        amp = self.data[:, col, 2]
        phase = self.data[:, col, 3]
        
        signal = np.zeros(nt, dtype = complex)
        
        for i in np.arange(nt):
            signal[i] = np.sqrt(amp[i])*np.complex(np.cos(phase[i]), np.sin(phase[i]))

        spectrum = np.abs(fftshift(fft(signal, nt)))
        spectrum = spectrum*spectrum

        period_s = self.kv['xlamds']/g_c
        freq_s = 1./period_s
        F = 1.0*np.arange(-nt/2, nt/2,)/nt*freq_s+freq_s # Frequency
        F = g_c/F # Frequency to wavelength

        dist[:,0] = F; dist[:,1] = spectrum
        return dist
    def get_power(self):
        '''
        Get the radiation power as a function of position z
        Parameters
          None
        Returns
          a 2D array that stores the positions and powers, respectively
        '''
        
        nz = self.nz; # print nz
        dist = np.zeros((nz, 2))
        dist[:,0] = self.field[:,0]
        
        dist[:,1] = self.data[:,:,0].sum(axis=0)
        return dist
    def get_peak_power(self):
        '''
        Get the radiation power as a function of position z
        Parameters
          None
        Returns
          a 2D array that stores the positions and powers, respectively
        '''
        
        nz = self.nz; # print nz
        dist = np.zeros((nz, 2))
        dist[:,0] = self.field[:,0]
        
        dist[:,1] = self.data[:,:,0].max(axis=0)
        return dist
    def get_energy(self):
        '''
        Get the radiation energy as a function of position z
        Parameters
          None
        Returns
          a 2D array that stores the positions and energies, respectively
        '''
        
        nz = self.nz
        dist = np.zeros((nz, 2))
        dist[:,0] = self.field[:,0]
        
        period_s = self.kv['xlamds']/g_c
        dist[:,1] = self.data[:,:,0].sum(axis=0)*period_s
        return dist
    def get_rms_x(self):
        '''
        Get the rms beam size in x direction as a function of position z
        Parameters
          None
        Returns
          a 2D array that stores the positions and rms beam sizes in x direction, respectively
        '''
        
        nz = self.nz
        dist = np.zeros((nz, 2))
        dist[:,0] = self.field[:,0]
        
        period_s = self.kv['xlamds']/g_c
        
        xrms = self.data[:,:,7]
        current = self.current
        
        for i in np.arange(xrms.shape[1]):
            dist[i, 1] = np.sqrt(np.sum(xrms[:,i]*xrms[:,i]*current[:,1])/np.sum(current[:,1]))
    
        return dist
    def get_rms_y(self):
        '''
        Get the rms beam size in y direction as a function of position z
        Parameters
          None
        Returns
          a 2D array that stores the positions and rms beam sizes in y direction, respectively
        '''
        
        nz = self.nz
        dist = np.zeros((nz, 2))
        dist[:,0] = self.field[:,0]
        
        period_s = self.kv['xlamds']/g_c
        
        xrms = self.data[:,:,8]; print(xrms.shape)
        current = self.current
        
        for i in np.arange(xrms.shape[1]):
            dist[i, 1] = np.sqrt(np.sum(xrms[:,i]*xrms[:,i]*current[:,1])/np.sum(current[:,1]))
    
        return dist
    def write(self, filename = 'fel', suffix = None):
        '''
        Write the Genesis input file to the `filename`
        '''
        if suffix is None:
            suffix = '.in'
        ff = open(filename+suffix, 'w')
        ff.write(self.output)
        ff.close()
        return
    def qsub(self, filename):
        '''
        Write the batch file to the `filename`, which could be submitted to a server by "qsub filename"
        '''
        con = '''\
#!/bin/zsh
#
#$ -cwd
#$ -o '''+filename+'''.o
#$ -e '''+filename+'''.e
#$ -V
#$ -l h_cpu=12:00:00
#$ -l h_rss=2G
#$ -P pitz
##$ -pe multicore 32
##$ -R y
echo '''+filename+'''.in | genesis'''
        ff = open(filename+'.sh', 'w')
        ff.write(con)
        ff.close()
        return
