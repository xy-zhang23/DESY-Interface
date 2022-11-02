# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:16:04 2020
Interfaces for Generator, Astra and Genesis version 2 and 4, and so on
@author: lixiangk

"""

from .Namelists import *

# 'Genesis2` also has only one `Namelist`, like Generator in Astra
class Genesis2(Namelist):
    starter = ' $'
    ending = ' $end'
    case = ''
    
    def __init__(self, *args, **kwargs):
        
        super(Genesis2, self).__init__(name = 'newrun', **kwargs)
        
    def update(self, quoting = True):
        '''
        Update output string to write into a file
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
                        if i == 0:
                            output += ' {}={} '.format(k, vi)
                        else:
                            output += '{} '.format(vi)
                    output += '\n'
            else:
                if isinstance(v, str):
                    if quoting:
                        v = '\''+v+'\''
                else:
                    v = toStr(v)
                output += ' {}={}\n'.format(k, v)
        
        output += self.ending+'\n\n'
        self.output = output
    
    def write(self, inputName = 'gene.in', case = ''):   
        super(Genesis2, self).write(inputName, case)
        
    def qsub(self, jobName = None, inputName = 'gene.in', direc = '.',
             submit = False, command = 'genesis', **kwargs):
        job = QsubJob(command = command, echo = True)
        job.create(jobName, inputName, direc, submit, **kwargs)
        
        
# `Genesis4` has a list of `Namelist` as `Astra`, therefore just inherit from
# the `Astra` class
class Genesis4(Namelists):
    '''
    Interface class to generate input file for `Genesis` version 4. 
    Since its input includes usually more than one `Namelist`, the class is 
    defined as a list of `Namelist`.
    '''
    def __init__(self, namelists = None, *args):
        super(Genesis4, self).__init__(namelists, *args)
        #self.update() # It's already called in __init__
        
    def update(self):
        output = ''
        for _, namelist in self.kv.items():
            #print('updated here: ', namelist.name)
            
            namelist.starter = '&'
            namelist.ending = '&end'
            namelist.case = 'lower'

            namelist.name = namelist.name.lower()
            
            namelist.update(quoting = False)
            output += namelist.output
        self.output = output
        
    def write(self, inputName = 'gene.in', case = ''):
        super(Genesis4, self).write(inputName, case)
    
    def qsub(self, jobName = None, inputName = 'gene.in', direc = '.',
             submit = False, command = 'gencore', 
             cmd1 = 'module add gnu openmpi phdf5/1.10.6', **kwargs):
        
        job = QsubJob(command = command, echo = False)
        job.create(jobName, inputName, direc, submit, cmd1 = cmd1, **kwargs)
    
    def sbatch(self, jobName = None, inputName = 'gene.in', direc = '.',
             submit = False, command = 'gencore', 
             cmd1 = 'module add gnu openmpi phdf5/1.10.6', **kwargs):
        
        job = SbatchJob(command = command, echo = False)
        job.create(jobName, inputName, direc, submit, cmd1 = cmd1, **kwargs)