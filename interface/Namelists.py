# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:16:04 2020
Interfaces for Generator, Astra and Genesis version 2 and 4, and so on
@author: lixiangk
"""

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import numpy as np
#import subprocess

from .Job import *

try:
    from requests.structures import CaseInsensitiveDict
    ADict = CaseInsensitiveDict
except Exception as err:
    print(err)
    
    class CaseInsensitiveDict(dict):
        @classmethod
        def _k(cls, key):
            return key.lower() if isinstance(key, str) else key

        def __init__(self, *args, **kwargs):
            super(CaseInsensitiveDict, self).__init__(*args, **kwargs)
            self._convert_keys()
        def __getitem__(self, key):
            return super(CaseInsensitiveDict, self).__getitem__(self.__class__._k(key))
        def __setitem__(self, key, value):
            super(CaseInsensitiveDict, self).__setitem__(self.__class__._k(key), value)
        def __delitem__(self, key):
            return super(CaseInsensitiveDict, self).__delitem__(self.__class__._k(key))
        def __contains__(self, key):
            return super(CaseInsensitiveDict, self).__contains__(self.__class__._k(key))
        def has_key(self, key):
            return super(CaseInsensitiveDict, self).has_key(self.__class__._k(key))
        def pop(self, key, *args, **kwargs):
            return super(CaseInsensitiveDict, self).pop(self.__class__._k(key), *args, **kwargs)
        def get(self, key, *args, **kwargs):
            return super(CaseInsensitiveDict, self).get(self.__class__._k(key), *args, **kwargs)
        def setdefault(self, key, *args, **kwargs):
            return super(CaseInsensitiveDict, self).setdefault(self.__class__._k(key), *args, **kwargs)
        def update(self, E={}, **F):
            super(CaseInsensitiveDict, self).update(self.__class__(E))
            super(CaseInsensitiveDict, self).update(self.__class__(**F))
        def _convert_keys(self):
            for k in list(self.keys()):
                v = super(CaseInsensitiveDict, self).pop(k)
                self.__setitem__(k, v)
    ADict = CaseInsensitiveDict
    #ADict = dict



def toStr(v):
    '''
    Convert the input argument to string.
    Parameters:
      v: any type
    Returns:
      str: string
    '''
    t = type(v)
    if t == str:
        v = v
    elif isinstance(v, (bool, np.bool_)):
        v = '{}'.format(v).lower()
    elif isinstance(v, (int, float, np.integer, np.floating)):
        v = '{:.7g}'.format(v)
    return v
  
class Namelist():
    '''
    Namelist used in Astra: newrun, output, scan, modules, error, charge, aperture,
                      wake, cavity, solenoid, quadrupole, dipole, laser
          or used in Genesis 1.3 V4: setup, lattice, field, beam
    Example:
      To add a `newrun` module, call
          newrun = Module('Newrun', Run = 1, Head = 'PITZ beamline simulation',
                           Distribution = 'beam.ini', Auto_Phase = True, 
                           Track_All = True, check_ref_part = False, 
                           Lprompt = False, Max_step=200000)
      The inputs for newrun could be changed by 
          newrun.set(Run = 2, Track_all = False)
    '''
    starter = '&'
    ending = '&end'
    case = ''
    
    def __init__(self, name = 'Input', **kwargs):
        '''
        Parameters
          name: name of the namelist, e.g., 'Input', 'Newrun', case insensitive
          **kwargs: key value pairs defining the properties of the `Namelist`
        '''
        self.kv = ADict()
        self.set(name, **kwargs)
        
    def set(self, *args, **kwargs):
        '''
        Set or update the properties of the namelist
        Parameters
        ----------
        *args : the name of the `Namelist` if not empty
        **kwargs : key value pairs defining the properties of the `Namelist`
        Returns
        -------
        None.
        '''
        if len(args)>0:
            self.name = args[0]
            
        if ADict == dict:
            for key, value in kwargs.items():
                self.kv.update({key.upper():value})
        else:
            self.kv.update(**kwargs)
        
        self.update()
        
    def delete(self, key, *args):
        '''
        Delete one or more properties
        Parameters
        ----------
        key : String
            Name of the property, e.g., `MaxE` of `Cavity` for `Astra`
        *args: String
            More key names 
        Returns
        -------
        None.

        '''
        key = key.upper()
        if key in self.kv.keys():
            self.kv.pop(key)
        if len(args)>0:
            for key in args:
                if key in self.kv.keys():
                    self.kv.pop(key)
        self.update()
        
    def update(self, quoting = True):
        '''
        Update output to be written to a file
        Returns
        -------
        None.

        '''
        output = self.starter+self.name+' \n'
        for k, v in self.kv.items():
            if self.case.upper() in ['UPPER', 'BIG', 'BIGGER']:
                k = k.upper()
            elif self.case.upper() in ['LOWER', 'SMALL', 'SMALLER']:
                k = k.lower()
            else:
                k = k
                
            #print(k)
            if isinstance(v, (list, tuple, np.ndarray)):
                dim = len(np.asarray(v).shape)
                if dim == 1:
                    for i, vi in enumerate(v):
                        if isinstance(vi, str):
                            if quoting:
                                vi = '\''+vi+'\''
                        else:
                            vi = toStr(vi)
                        output += ' {}({})={}\n'.format(k, i+1, vi)
                elif dim == 2:
                    for i, vi in enumerate(v):
                        if k.upper() in ['D_GAP', 'MODULE', 'AP_GR', 'Q_MULT_A', 'Q_MULT_B']: # only for `Astra`->Cavity
                            for j, vij in enumerate(vi):
                                if isinstance(vij, str):
                                    if quoting:
                                        vij = '\''+vij+'\''
                                else:
                                    vij = toStr(vij)
                                output += ' {}({},{})={}\n'.format(k, j+1, i+1, toStr(vij))
                        elif k.upper() in ['D1', 'D2', 'D3', 'D4']: # only for `Astra`->Dipole
                            output += ' {}({})=({},{})\n'.format(k, i+1, toStr(vi[0]), toStr(vi[1]))
                        else:
                            print('Unknown key!')
            else:
                if isinstance(v, str):
                    if quoting:
                        v = '\''+v+'\''
                else:
                    v = toStr(v)
                output += ' {}={}\n'.format(k, v)
        
        output += self.ending+'\n\n'
        self.output = output
    
    def write(self, inputName, case = None):
        
        if case == None:
            case = self.case
            
        if case.upper() in ['UPPER', 'BIG', 'BIGGER']:
            output = self.output.upper()
        elif case.upper() in ['LOWER', 'SMALL', 'SMALLER']:
            output = self.output.lower()
        else:
            output = self.output
            
        ff = open(inputName , 'w')
        ff.write(output)
        ff.close()
        

# Make an alias of the class `Namelist`
Module = Namelist


class Namelists:
    '''
    `Namelists` is a class of a `dict` type, with each of its value a `Namelist`
    instance
    This class is used to generate input file for `Astra`, since it usually 
    includes more than one `Namelist`.
    '''
    
    def __init__(self, namelists = None, *args):
        self.kv = ADict()
        #print(namelists, args)
        if namelists is not None:
            namelists = np.asarray([namelists]).flatten()
            #self.addNamelists(np.asarray([namelists]).flatten())    
            self.add(*namelists)
        if len(args)>0:
            self.add(*args)
        
    def add(self, namelist, *args, **kwargs):
        '''
        Add at least one or more `Namelist` instance to the class
        Parameters
        ----------
        namelist : Namelist
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        self.kv.update({namelist.name:namelist})
        if len(args)>0:
            for namelist in args:
               self.kv.update({namelist.name:namelist})
        self.update(**kwargs)
    
    def addNamelists(self, namelists):
        '''
        Add a series of `Namelist`

        Parameters
        ----------
        namelists : array or tuple or list
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        for namelist in namelists:
            self.add(namelist)
                
    # Create a alias, to be compatible with old versions
    add_namelists = addNamelists
    add_modules = add_namelists
    #update = add_namelists
    
    def update(self):
        output = ''
        for _, namelist in self.kv.items():
            #print('updated here 2: ', namelist.name)
            output += namelist.output
        self.output = output
        
    def delete(self, key, *args):
        '''
        Delete one or more properties
        Parameters
        ----------
        key : String
            Name of the `Namelist` to be deleted
        *args: String
            More names 
        Returns
        -------
        None.

        '''
        key = key.upper()
        if key in self.kv.keys():
            self.kv.pop(key)
        if len(args)>0:
            for key in args:
                if key in self.kv.keys():
                    self.kv.pop(key)
        self.update()
        
    def write(self, inputName, case = 'upper'):
        
        if case.upper() in ['UPPER', 'BIG', 'BIGGER']:
            output = self.output.upper()
        elif case.upper() in ['LOWER', 'SMALL', 'SMALLER']:
            output = self.output.lower()
        else:
            output = self.output
            
        ff = open(inputName, 'w')
        ff.write(output)
        ff.close()
        
# Make a alias of the class `Namelists`
Modules = Namelists