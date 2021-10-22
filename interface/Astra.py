# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:16:04 2020
Interfaces for Generator, Astra and Genesis version 2 and 4, and so on
@author: lixiangk
"""

from .Namelists import *

def _exist_(st, fname):
    '''
    Check if a string `st` exists in the file `fname`
    '''
    exist = False
    if os.path.exists(fname):
        log = open(fname)
        for i, line in enumerate(log):
            if st in line:
                exist = True
    return exist

class Generator1(Namelist):
    '''
    `Generator` has only one `Namelist`
    '''
    def __init__(self, *args, **kwargs):
        
        super(Generator1, self).__init__(name = 'Input', **kwargs)
        self.set(Lprompt = False)
    
    def write(self, inputName = 'gen.in'): 
        super(Generator1, self).write(inputName)
    
    def run(self, inputName = 'gen.in', force = True, command = 'generator'):
        '''
        If force is True, it will always start Astra simulation, otherwise,
        it will check if the simualtion has been already done before
        '''
        _, baseName, ext = fileparts(inputName)
        if ext == '':
            ext = '.in'
        inputName = baseName+ext
        
        logfile = baseName+'.log'        
        finished = _exist_('phase-space distribution saved', logfile)
        print('Is finished?', finished)
        
        r = -1
        if (not finished) or force:
            cmd = command+' '+inputName
            #os.systme(cmd+' 2>&1 | tee '+logfile)
            with open(logfile, 'w') as fout:
                r = subprocess.call(cmd, stdout = fout)     
        return r
        
    def qsub(self, jobName = None, inputName = 'gen.in', direc = '.',
             submit = False, command = 'generator', **kwargs):
        job = QsubJob(command = command, echo = False)
        job.create(jobName, inputName, direc, submit, **kwargs)
        

class Astra(Namelists):
    '''
    Interface class to generate input file for `Astra`. Since its input includes 
    usually more than one `Namelist`, the class is defined as a list of `Namelist`.
    '''
    
    def __init__(self, namelists = None, *args):
        super(Astra, self).__init__(namelists, *args)
        
    def write(self, inputName = 'ast.in'): 
        super(Astra, self).write(inputName)
        
    def run(self, inputName = 'ast.in', force = True, command = 'astra'):
        '''
        If force is True, it will always start Astra simulation, otherwise,
        it will check if the simualtion has been already done before
        '''
        _, baseName, ext = fileparts(inputName)
        if ext == '':
            ext = '.in'
        inputName = baseName+ext
        
        logfile = baseName+'.log'        
        finished = _exist_('finished simulation', logfile)
        print('Is finished?', finished)
        
        r = -1
        if (not finished) or force:
            cmd = command+' '+inputName
            #os.systme(cmd+' 2>&1 | tee '+logfile)
            with open(logfile, 'w') as fout:
                r = subprocess.call(cmd, stdout = fout)     
        return r
    
    def qsub(self, jobName = None, inputName = 'ast.in', genName = None, 
              direc = '.', submit = False, command = 'astra', **kwargs):
        '''
        Write the batch file to the `filename`, which could be submitted to 
        a server by "qsub filename". Only tested in Zeuthen site
        '''
        
        if genName is not None:
            _, baseName, ext = fileparts(genName)
            if ext == '':
                ext = '.in'
            
            cmd1 = '''generator '''+baseName+ext+''' 2>&1 | tee '''+baseName+'''.log
'''
        else:
            cmd1 = '''\n'''
        
        job = QsubJob(command = command, echo = False)
        job.create(jobName, inputName, direc, submit, cmd1 = cmd1, **kwargs)
        
