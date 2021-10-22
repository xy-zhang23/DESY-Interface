# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 14:57:04 2020

@author: lixiangk
"""

import h5py
from numpy.fft import fftshift,fft

from .Numeric import *

def calcSpectrum(amp, phase, lambda0 = 100e-6, freq0 = None):
    '''
    Calculate the spectrum from samples

    Parameters
    ----------
    amp : 1D or 2D array
        Amplitudes of samples of the signal. In the case of 2D, the first dimension 
        is along the slices and the second dimension is along the undulator.
    phase : 1D or 2D array
        Phases of samples of the signal. 
    lambda0 : double, optional
        Seperation of sampling (usually the wavelength) in meter. The default is 100e-6.
    freq0 : double, optional
        Sampling frequency. If defined, it dominates `lambda0`. The default is None.
        
    Returns
    -------
    wavelength : 1D array
        Wavelength of the signal transformed in frequency domain.
    spectrum : 1D or 2D array
        Spectra intensity of the signal transformed in frequency domain.

    '''
    
    nsample = len(amp) # number of samples
    
    signal = np.sqrt(amp)*(np.cos(phase)+np.sin(phase)*1j)
    
    axis = 0
    spectrum = np.abs(fftshift(fft(signal, nsample, axis), axis))
    spectrum = spectrum*spectrum
    
    if freq0 is None and lambda0 is not None:
        freq0 = 1./lambda0*g_c # sampling frequency
    
    F = 1.0*np.arange(-nsample/2, nsample/2,)/nsample*freq0+freq0 # Frequency
    wavelength = g_c/F # Frequency to wavelength
    
    return wavelength, spectrum
    
class PostGenesis13:
    version = 4 # default
    def __init__(self, fname = None, **kwargs):
        '''
        
        Parameters
        ----------
        fname : str
            Genesis V4 main output file.
        **kwargs : TYPE
            DESCRIPTION.
            
        Returns
        -------
        None.

        '''
        
        debug = 0
        if len(kwargs)>0:
            if 'debug' in kwargs.keys():
                debug = kwargs['debug']    
            if 'harmonic' in kwargs.keys():
                harmonic = kwargs['harmonic']
                self.fieldname = 'Field%d' % harmonic
            else:
                self.fieldname = 'Field'
                
        if fname is None:
            fname = get_file('.h5'); print(fname)
        
        _, _, ext = fileparts(fname); print(ext)
        
        if ext.upper() in ['.H5']:
            self.version = 4
            self.load4(fname)
        elif ext.upper() in ['.OUT']:
            self.version = 2
            self.load2(fname)    
        else:
            print('Unknown file extension!')
            return
        
        if debug:  
            if 'fig_ext' in kwargs.keys():
                fig_ext = kwargs['fig_ext']
            else:
                fig_ext = '.png'
            self.plot_current(fig_ext = fig_ext)
            self.plot_power(fig_ext = fig_ext)
            self.plot_energy(fig_ext = fig_ext)
            self.plot_spectrum(fig_ext = fig_ext)
        
    def print_all(self, prop = None, order = 2):
        if self.version == 2:
            print(self.outputs)
            return
        
        print('./')
        for key in self.file.keys():
            if prop is None:
                print('  - %s' % key)
                if order == 1:
                    continue
                for subkey in self.file.get(key).keys():
                    print('    - %s' % subkey)
            elif prop.upper() == key.upper():
                print('  - %s' % key)
                for subkey in self.file.get(key).keys():
                    print('    - %s' % subkey)
                break
            
    def load4(self, fname):
        file = h5py.File(fname, 'r')
        
        tmp = {}
        for k, v in file.get('Beam').items():
            tmp.update({k:v[:]})
        self.beam = tmp
        
        tmp = {}
        for k, v in file.get(self.fieldname).items():
            tmp.update({k:v[:]})
        self.field = tmp
        
        nstep, nslice = file.get(self.fieldname).get('power').shape
        
        meta = file.get('Meta').get('InputFile')[0].decode().split()
        
        kv = {}
        for _, a in enumerate(meta):
            tmp = a.split('=')
            if len(tmp)>1:
                kv.update({tmp[0]:tmp[1]})
                
        self.file = file
        self.nslice = nslice
        self.nstep = nstep
        self.kv = kv
        
        self.outputs = [temp for temp in self.file.get(self.fieldname).keys()] # Output properties for radiation fields
        
        self.sample = file.get('Global/sample')[0]
        self.lambdaref = file.get('Global/lambdaref')[0]
        self.gamma0 = file.get('Global/gamma0')[0]
        self.one4one = file.get('Global/one4one')[0]
        self.scan = file.get('Global/scan')[0]
        self.slen = file.get('Global/slen')[0]
        self.time = file.get('Global/time')[0]
        
        lslice = self.lambdaref*self.sample
        
        self.lslice = lslice # length of a slice
        #self.zbunch = np.linspace(lslice/2, self.slen-lslice/2, nslice)
        self.zbunch = np.arange(nslice)*self.lslice+self.lslice/2
        self.current = file.get('Beam/current')[:].flatten()
        
        #self.bunching = file.get('Beam/bunching')
        
        self.zplot = file.get('Lattice/zplot')[:] # z coordinates along the undulator
        self.zstep = self.zplot[1]-self.zplot[0]  
        
        self.zpower = np.sum(file.get(self.fieldname).get('power')[:], axis = 1) # power along the undulator
        self.zenergy = self.zpower*self.lslice/g_c # energy along the undulator
        
        self.power = np.sum(self.zpower)    # Total power
        self.energy = np.sum(self.zenergy)  # Total energy
        
        # spectrum at the end of the undulator
        self.wavelength, self.spectrum = self.get_spectrum(self.zplot[-1])
        return
    
    load = load4
    
    def load2(self, fname):
        
        '''
        Read the standard output (in ascii format and stored in `fname`) 
        into memories for further analysis
        '''
        
        file = open(fname, 'r')
        line = file.readline()

        # First, read the input namelist from the header of the file 
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
                    line = file.readline()
                    continue

                if len(value.split()) == 1:
                    value = np.float(value.replace("D", "E"))
                else:
                    value = [np.float(v.replace("D", "E")) for v in value.split()]
                kv.update({key.strip():value})          

            match = re.search(r'^[ ]*z\[m\][ ]*aw[ ]*qfld[ ]*$', line)
            if match:
                # print line
                nstep = np.int(kv['zstop']/kv['delz']/kv['xlamd'])+1
                field = np.zeros((nstep, 3))
                for i in np.arange(nstep):
                    line = file.readline()
                    field[i] = np.array([np.float(v.replace("D", "E")) for v in line.split()])
                break
            line = file.readline()
        
        nslice = np.int(kv['nslice'])
        nc = np.int(np.sum(kv['lout'])); #print nc # number of output items

        # Initialize the arrays for storing the data blocks
        current = np.zeros((nslice, 1))
        data = np.zeros((nslice, nstep, nc)); #print data.shape
        
        # Then read the blocks of data from the file
        islice = 0
        while True:
            if not line:
                break
            match = re.search(r'[ E]* current', line)
            if match:
                icurrent, tmp = line.split()
                icurrent = np.float(icurrent)
                current[islice, 0] = icurrent
        
                line = file.readline()
                line = file.readline()
                line = file.readline()

                outputs = line.split()

                for j in np.arange(nstep):
                    line = file.readline(); #print line
                    line = line.replace('NaN', '  0')
                    data[islice, j] = np.array([np.float(v.replace("D", "E")) for v in line.split()])
                #break
                islice += 1
            line = file.readline()
            
        file.close()
        
        self.kv = kv # dict with key-value pairs
        self.nslice = nslice; self.nstep = nstep; self.nc = nc
        self.field = field
        
        self.outputs = outputs # Output properties for radiation fields
        
        #self.data = data
        # to make the first two dimensions consistent with version 4
        # now nstep, nslice, nc
        self.data = np.transpose(data, [1, 0, 2]) 
        
        self.gamma0 = self.kv['gamma0']
        self.time = self.kv['itdp']
        self.lambdaref = self.kv['xlamds']
        self.sample = 1.0/self.kv['zsep']
        
        self.lslice = self.lambdaref/self.sample # length of a slice
        
        self.zbunch = np.arange(nslice)*self.lslice+self.lslice/2
        self.current = current
        
        self.zplot = field[:,0]; self.zstep = self.zplot[1]-self.zplot[0] # z coordinates along the undulator
        
        self.zpower = self.data[:,:,0].sum(axis = 1) # power along the undulator
        self.zenergy = self.zpower*self.lslice/g_c   # energy along the undulator
        
        self.power = np.sum(self.zpower)    # Total power
        self.energy = np.sum(self.zenergy)  # Total energy
        
        self.wavelength, self.spectrum = self.get_spectrum(self.zplot[-1]) # spectrum at the end of the undulator
        return
    
    def get_fielddata(self, name, at = None):
        '''
        Get the field data/property defined by `name` along the slices at position z = at
        Parameters
          name: string
              Name of a property such as 'power', 'increment', 'p_mid', and so on for version 2
              and 'power', 'intensity-farfield', 'phase-farfield' and so on for version 4
          at: float or None
              Longitudinal position along the undulator, used to calculate the nth step of 
              interest; if None, return all
        Returns
        -------
        data : 1D array
            The requested data along the slices
        '''
        if at is not None:
            if at < 0:
                col = -at
            else:
                #col = np.int(at/self.zstep)
                col = np.where(np.abs(self.zplot-at)<self.zstep/2)[0][0]
        else:
            col = slice(0, self.nstep)
        #print(col)
        
        props = [temp.upper() for temp in self.outputs]
        if name.upper() not in props: # check if `name` is included in the outputs
            print('Available properties are: ', self.outputs)
            return None
            
        if self.version<4:
            third = self.outputs.index(name)
            return self.data[col,:,third]
        else:
            return self.get_data(self.fieldname, name)[:][col]

    def get_beamdata(self, name, at = None):
        '''
        Get the beam data/property defined by `name` along the slices at position z = at
        Parameters
          name: string
              Name of a property such as 'power', 'increment', 'p_mid', and so on for version 2
              and 'power', 'intensity-farfield', 'phase-farfield' and so on for version 4
          at: float or None
              Longitudinal position along the undulator, used to calculate the nth step of 
              interest; if None, return all
        Returns
        -------
        data : 1D array
            The requested data along the slices
        '''
        if at is not None:
            if at < 0:
                col = -at
            else:
                #col = np.int(at/self.zstep)
                col = np.where(np.abs(self.zplot-at)<self.zstep/2)[0][0]
        else:
            col = slice(0, self.nstep)
        #print(col)
        
        if self.version < 4:
            props = [temp.upper() for temp in self.outputs]
        else:
            props = [temp.upper() for temp in self.file.get('Beam').keys()]
            
        if name.upper() not in props: # check if `name` is included in the outputs
            print('Available properties are: ', self.outputs)
            return None
            
        if self.version<4:
            third = self.outputs.index(name)
            return self.data[col,:,third]
        else:
            return self.get_data('Beam', name)[:][col]
    
    def get_spectrum(self, at = None):
        '''
        Get the radiation spectrum at position z = at
        Parameters
          at: float or None
              Longitudinal position along the undulator, used to calculate the
              nth step of interest; if None, return all
        Returns
        -------
        wavelength : 1D array
            Wavelength of the signal transformed in frequency domain.
        spectrum : 1D or 2D array
            Intensity of the signal transformed in frequency domain. if at == None, 
            then return a 2D array, with the first dimension along the slices and
            the second dimension along the undulator.
        '''
        
        if self.version<4:
            n1, n2 = 'p_mid', 'phi_mid'
        else:
            n1, n2 = 'intensity-nearfield', 'phase-nearfield'
            
        amp = self.get_fielddata(n1, at)
        phase = self.get_fielddata(n2, at)
        
        return calcSpectrum(amp, phase, self.lslice)
    
    def get_data(self, path, *args):
        '''
        Get data by the `path` from the hdf5 file, only valid for version 4

        Parameters
        ----------
        path : string
            Path of the variable of interest, e.g., 'Field/power'.
        *args: string
            subpath, e.g., get_data('Field', 'power') is equiverlent to get_data('Field/power')
            
        Returns
        -------
        data : array
            The dimension of the array depends on the variable type.

        '''
        if self.version == 2:
            print('Only work for version 4!')
            return
        
        if self.version == 4:
            if len(args)>0:
                for temp in args:
                    path = path+'/'+temp
            try:
                data = self.file[path][:]
            except Exception as err:
                print(err)
                data = None
        return data
        
    def plot_current(self, x = 'time', fig_ext = '.png'):
        
        fig, ax = plt.subplots()
        if x.upper() == 'TIME':
            ax.plot(self.zbunch/g_c*1e12, self.current, '-')
            ax.set_xlabel(r'Time (ps)')
        elif x.upper() == 'LENGTH':
            ax.plot(self.zbunch, self.current, '-')
            ax.set_xlabel(r'Position (mm)')
        elif x.upper() == 'SLICE':
            ax.plot(self.current, '-')
            ax.set_xlabel(r'# of slice')
            
        ax.set_ylabel(r'Current (A)')
        ax.grid()
        fig.savefig('current'+fig_ext)
        
    def plot_spectrum(self, fig_ext = '.png'):
        
        fig, ax = plt.subplots()
        ax.plot(self.wavelength*1e6, self.spectrum, '-')
        ax.set_xlabel(r'Wavelength ($\mu$m)')
        ax.set_ylabel(r'Intensity (arb. unit)')
        ax.grid()
        fig.savefig('spectrum'+fig_ext)
        
    def plot_power(self, fig_ext = '.png'):
        
        fig, ax = plt.subplots(ncols = 1, figsize = (4, 4))
        ax.plot(self.zplot, self.zpower/1e6, '-')
        ax.set_xlabel(r'$z$ (m)')
        ax.set_ylabel(r'Power (MW)')
        ax.set_yscale('log')
        ax.grid()
        fig.savefig('power-z'+fig_ext)
        
    def plot_energy(self, fig_ext = '.png'):
        
        fig, ax = plt.subplots(ncols = 1, figsize = (4, 4))
        ax.plot(self.zplot, self.zenergy*1e6, '-')
        ax.set_xlabel(r'$z$ (m)')
        ax.set_ylabel(r'Energy ($\mu$J)')
        ax.set_yscale('log')
        ax.grid()
        fig.savefig('energy-z'+fig_ext)

PostGenesis = PostGenesis13