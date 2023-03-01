# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 16:40:22 2021

@author: lixiangk
"""

import os

def fileparts(fullname):
    '''
    Split the full name of a file into parts
    Parameters
      fullname: full name of the file
    Return
      [path, name, ext]: path, name and extension of the file
    '''
    
    [path, name] = os.path.split(fullname)
    [name, ext] = os.path.splitext(name)
    return [path, name, ext]

class CondorJob:
    '''
    Create a batch file to submit a job to the server by the `qsub` command.
    Only tested in Zeuthen site.
    '''
    def __init__(self, command = 'astra', echo = False):
        self.command = command
        self.echo = echo
            
    def create(self, jobName = None, inputName = 'ast.in', direc = '.',
               submit = False, **kwargs):
        '''
        Parameters
          jobName: name of the job
          inputName: input file name for e.g. `Astra` or `Generator`
          direc: directory name to open before running the program
          submit: if True, directly submit the job to the server
          **kwargs: more command to run before the program
        '''
        
        _, baseName, ext = fileparts(inputName)
        if ext == '':
            ext = '.in'
            
        if jobName == None:
            jobBaseName = 'myjob@'+baseName
        else:
            _, jobBaseName, _ = fileparts(jobName)
        
        con1 = '''\
#HTC condor job

executable = /bin/zsh
arguments = '''+jobBaseName+'''.sh

output = '''+jobBaseName+'''.o
error = '''+jobBaseName+'''.e

getenv = True
request_memory = 2G
#request_cpus = 32

queue 1
'''
        
        chdir = '''cd '''+direc+'''
'''
        con2 = chdir
        
        if len(kwargs)>0:
            for _, cmd in kwargs.items():
                con2 += cmd + '''
'''
        if self.echo:
            cmd = '''echo '''+baseName+ext+''' |'''+self.command
        else:
            cmd = self.command+''' '''+baseName+ext
        cmd += ''' 2>&1 | tee '''+baseName+'''.log
'''
        con2 += cmd
        
	
        subName = jobBaseName+'.submit'
        ff = open(subName, 'w')
        ff.write(con1)
        ff.close()
        
        jobName = jobBaseName+'.sh'
        ff = open(jobName, 'w')
        ff.write(con2)
        ff.close()
        
        if submit:
            os.system('condor_submit '+subName)
            
        return    

class QsubJob:
    '''
    Create a batch file to submit a job to the server by the `qsub` command.
    Only tested in Zeuthen site.
    '''
    def __init__(self, command = 'astra', echo = False):
        self.command = command
        self.echo = echo
            
    def create(self, jobName = None, inputName = 'ast.in', direc = '.',
               submit = False, **kwargs):
        '''
        Parameters
          jobName: name of the job
          inputName: input file name for e.g. `Astra` or `Generator`
          direc: directory name to open before running the program
          submit: if True, directly submit the job to the server
          **kwargs: more command to run before the program
        '''
        
        _, baseName, ext = fileparts(inputName)
        if ext == '':
            ext = '.in'
            
        if jobName == None:
            jobBaseName = 'myjob@'+baseName
        else:
            _, jobBaseName, _ = fileparts(jobName)
        
        con = '''\
#!/bin/zsh
#
#$ -cwd
#$ -o '''+jobBaseName+'''.o
#$ -e '''+jobBaseName+'''.e
#$ -V
#$ -l h_cpu=12:00:00
#$ -l h_rss=2G
#$ -P pitz
##$ -pe multicore 32
##$ -R y

'''
        
        chdir = '''cd '''+direc+'''
'''
        con += chdir
        
        if len(kwargs)>0:
            for _, cmd in kwargs.items():
                con += cmd + '''
'''
        if self.echo:
            cmd = '''echo '''+baseName+ext+''' |'''+self.command
        else:
            cmd = self.command+''' '''+baseName+ext
        cmd += ''' 2>&1 | tee '''+baseName+'''.log
'''
        con += cmd
        
        jobName = jobBaseName+'.sh'
        ff = open(jobName, 'w')
        ff.write(con)
        ff.close()
        
        if submit:
            os.system('qsub '+jobName)
            
        return        

class SbatchJob:
    '''
    Create a batch file to submit a job to the server by the `qsub` command.
    Only tested in Zeuthen site.
    '''
    def __init__(self, command = 'astra', echo = False):
        self.command = command
        self.echo = echo
            
    def create(self, jobName = None, inputName = 'ast.in', direc = '.',
               submit = False, **kwargs):
        '''
        Parameters
          jobName: name of the job
          inputName: input file name for e.g. `Astra` or `Generator`
          direc: directory name to open before running the program
          submit: if True, directly submit the job to the server
          **kwargs: more command to run before the program
        '''
        
        _, baseName, ext = fileparts(inputName)
        if ext == '':
            ext = '.in'
            
        if jobName == None:
            jobBaseName = 'myjob@'+baseName
        else:
            _, jobBaseName, _ = fileparts(jobName)
        
        con = '''\
#!/bin/bash

# one process for each physical core:
#SBATCH -c 2

# run on the haswell partition (pax11):
##SBATCH -p haswell
#SBATCH -p broadwell

# runtime 
#SBATCH --time=47:59:59

# number of tasks/cores, processes
#SBATCH --ntasks=1

# Job name:
#SBATCH --job-name=test
##SBATCH --error=test-%N-%j.err
##SBATCH --output=test-%N-%j.out
#SBATCH --error=test.err
#SBATCH --output=test.out

# copy environment variables from submit environment
#SBATCH --get-user-env

# send mail on all occasions:
##SBATCH --mail-type=ALL

'''
        
        chdir = '''cd '''+direc+'''
'''
        con += chdir
        
        if len(kwargs)>0:
            for _, cmd in kwargs.items():
                con += cmd + '''
'''
        if self.echo:
            cmd = '''echo '''+baseName+ext+''' |'''+self.command
        else:
            cmd = self.command+''' '''+baseName+ext
        cmd += ''' 2>&1 | tee '''+baseName+'''.log
'''
        con += cmd
        
        jobName = jobBaseName+'.sh'
        ff = open(jobName, 'w')
        ff.write(con)
        ff.close()
        
        if submit:
            os.system('sbatch '+jobName)
            
        return

