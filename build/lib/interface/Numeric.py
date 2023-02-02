# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 20:54:52 2020

@author: lixiangk
"""

import sys
    
# from IPython.display import Image
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import curve_fit

from scipy.special import jn, jn_zeros, j0, j1
from scipy.integrate import quad, ode

from scipy.constants import codata

import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.ticker import AutoMinorLocator

import os, re

from timeit import default_timer
import time

# import other packages needed by the users
# import h5py
# from opmd_viewer import OpenPMDTimeSeries

### Default plot setup
def plot_config():
    from cycler import cycler
    from matplotlib.ticker import AutoMinorLocator

    fsize = 12 # a quarter of the paper width: 20 pt; half of the paper width: 12
    font = {'size' : fsize, 'family' : 'serif'}
    color_cycle = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
    linestyle_cycle = ['-', '--', '-.', ':', (0, (5, 2, 5, 2)), (0, (10, 2, 5, 2, 2, 2)), (0, (12, 2, 2, 2))]
    marker_cycle = ['o', 'd', 'v', '^', '<', '>', '*']
    
    mpl.rc('font', **font)
    mpl.rc('xtick', labelsize = 10, direction = 'in', top   = True)
    mpl.rc('ytick', labelsize = 10, direction = 'in', right = True)
    mpl.rc('xtick.major', size = 5, width = 1)
    mpl.rc('ytick.major', size = 5, width = 1)
    mpl.rc('xtick.minor', size = 3, width = 0.7, visible = True)
    mpl.rc('ytick.minor', size = 3, width = 0.7, visible = True)
    
    mpl.rc('lines', linewidth=2, markersize=6, color='r')
    # mpl.rc('lines', linestyle = 'solid')
    mpl.rc('axes', labelpad = 0, prop_cycle=(cycler('color', color_cycle) + cycler('linestyle', linestyle_cycle) + cycler('marker', marker_cycle)))
    mpl.rc('legend', fontsize = 12, labelspacing = 0.05, handletextpad=0.4, frameon=False, handlelength=2.1)
    
    mpl.rc('figure', dpi = 100, figsize = (4, 4))
    mpl.rc('figure.subplot', bottom = 0.15, top = 0.9, left = 0.15, right = 0.9)
    
    mpl.rc('image', cmap = 'jet')
    
    return
plot_config()

def reset_margin(bottom = 0.15, top = 0.9, left = 0.15, right = 0.9):
    mpl.rc('figure.subplot', bottom = bottom, top = top, left = left, right = right)
    return
###

class Dict(dict):
    '''
    Postprocessing class for Astra and/or Warp simulations
    '''
    def __init__(self, **kwargs):
        '''
        Parameters
          **kwargs: key-value pairs
        '''
        
        self.x = dict(**kwargs)
        
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
        Return the corresponding parameter if key is one of `self.x.keys()`
        '''
        if key in self.x.keys():
            return self.x.get(key)
        
# commonly used scientific constants
class Const(Dict):
    
    h   =codata.value('Planck constant')
    hbar=codata.value('Planck constant over 2 pi')
    
    g   =codata.value('standard acceleration of gravity')
    c   =codata.value('speed of light in vacuum')
    
    qe  =codata.value('elementary charge')
    me  =codata.value('electron mass')
    re  =codata.value('classical electron radius')
    mec2=codata.value('electron mass energy equivalent in MeV')
    eV  =codata.value('electron volt')
    
    kB  =codata.value('Boltzmann constant')
    sigma =codata.value('Stefan-Boltzmann constant')
    
    eps0=codata.value('electric constant')
    mu0 =codata.value('mag. constant')
    
    def __init__(self):
        pass

class Unit:
    
    x     = r'$x$ (mm)'
    y     = r'$y$ (mm)'
        
    r     = r'$r$ (mm)'
    xi    = r'$\xi$ (mm)'
        
    z     = r'$z$ (m)'
        
    rms_x = r'$\sigma_x$ (mm)'
    rms_y = r'$\sigma_y$ (mm)'
    rms_z = r'$\sigma_z$ (mm)'
    
    rms_xy  = r'$\sigma_{x/y}$ (mm)'
    rms_xyz = r'$\sigma_{x/y/z}$ (mm)'
    
    emi_x = r'$\varepsilon_{n,x}$ (mm mrad)'
    emi_y = r'$\varepsilon_{n,y}$ (mm mrad)'
    emi_z = r'$\varepsilon_{n,y}$ (keV mm)'
    
    emi_xy = r'$\varepsilon_{n,x/y}$ (mm mrad)'
    
    momentum    = r'$P$ (MeV/c)'
    
    momentum_x  = r'$P_x$ (MeV/c)'
    momentum_y  = r'$P_y$ (MeV/c)'
    momentum_z  = r'$P_z$ (MeV/c)'
    
    kinetic     = r'$E_{\rm k}$ (MeV)'
    kinetic_rel = r'$\Delta E/E$ (%)'
    rel_kinetic = r'$\Delta E/E$ (%)'
    
    energy_spread = r'$\Delta E$ (keV)'
    kinetic_spread = r'$\Delta E$ (keV)'
    
    momentum_rel = r'$\Delta P/P$ (%)'
    rel_momentum = r'$\Delta P/P$ (%)'
    
    momentum_spread = r'$\Delta P$ (keV)/c'
    
    cs_alpha = r'$\alpha$'
    cs_beta  = r'$\beta$ (m)'
    cs_gamma = r'$\gamma$ (m$^{-1}$)'
    
    def __init__(self):
        pass
        
g_h   =codata.value('Planck constant')
g_hbar=codata.value('Planck constant over 2 pi')

g_g   =codata.value('standard acceleration of gravity')
g_c   =codata.value('speed of light in vacuum')

g_qe  =codata.value('elementary charge')
g_me  =codata.value('electron mass')
g_re  =codata.value('classical electron radius')
g_mec2=codata.value('electron mass energy equivalent in MeV')
g_eV  =codata.value('electron volt')

g_kB  =codata.value('Boltzmann constant')
g_sigma =codata.value('Stefan-Boltzmann constant')

g_eps0=codata.value('electric constant')
g_mu0 =codata.value('mag. constant')

u_x     = r'$x$ (mm)'
u_y     = r'$y$ (mm)'

u_r     = r'$r$ (mm)'
u_xi    = r'$\xi$ (mm)'

u_z     = r'$z$ (m)'

u_rms_x = r'$\sigma_x$ (mm)'
u_rms_y = r'$\sigma_y$ (mm)'
u_rms_z = r'$\sigma_z$ (mm)'

u_rms_xy  = r'$\sigma_{x/y}$ (mm)'
u_rms_xyz = r'$\sigma_{x/y/z}$ (mm)'

u_emi_x = r'$\varepsilon_{n,x}$ (mm mrad)'
u_emi_y = r'$\varepsilon_{n,y}$ (mm mrad)'
u_emi_z = r'$\varepsilon_{n,y}$ (keV mm)'

u_emi_xy = r'$\varepsilon_{n,x/y}$ (mm mrad)'

u_momentum    = r'$P$ (MeV/c)'

u_momentum_x  = r'$P_x$ (MeV/c)'
u_momentum_y  = r'$P_y$ (MeV/c)'
u_momentum_z  = r'$P_z$ (MeV/c)'

u_kinetic     = r'$E_{\rm k}$ (MeV)'
u_kinetic_rel = r'$\Delta E/E$ (%)'
u_rel_kinetic = r'$\Delta E/E$ (%)'

u_cs_alpha = r'$\alpha$'
u_cs_beta  = r'$\beta$ (m)'
u_cs_gamma = r'$\gamma$ (m$^{-1}$)'

try:
    import datetime
    def timestamp(option = True):   
    
        today = datetime.datetime.now()
    
        if option:
            stamp = '__%02d.%02d.%02d %02d:%02d:%02d__' % (today.day, today.month, today.year,
                                                           today.hour, today.minute, today.second)
        else:
            stamp = '__%02d.%02d.%02d__%02d_%02d_%02d__' % (today.day, today.month, today.year,
                                                                today.hour, today.minute, today.second)
        return stamp
except Exception as err:
    print(err)


try:
    import pandas as pd
    def pd_loadtxt(filename, **kwargs):
        if len(kwargs)>0:
            if 'delimiter' not in kwargs.keys() and 'sep' not in kwargs.keys():
                kwargs.update(delimiter = '\s+')
            if 'header' not in kwargs.keys():
                kwargs.update(header = None)
        else:
            kwargs.update(delimiter = '\s+')
            kwargs.update(header = None)
        
        data = pd.read_csv(filename, **kwargs).values
        return data
except Exception as err:
    print(err)


try:
    import tkinter as tk
    from tkinter import filedialog
    
    def get_path():
        '''
        Return a selected path by the user
        '''
        
        root = tk.Tk()
        #root.lift()
        root.attributes("-topmost", True)
        root.withdraw()
    
        path = filedialog.askdirectory(parent = root)
        
        if path == '':
            print('No path is selected!')
            path = '.'
    
        return path
    
    def get_file(filetype = '.txt', initialdir = "/"):
        '''
        Return a selected file name (including the path) by the user 
        '''
        print(initialdir)
        
        root = tk.Tk()
        #root.lift()
        root.attributes("-topmost", True)
        root.withdraw()
    
        if filetype != None:
            filetypes = (("user files", '*'+filetype), ("all files", "*.*"))
        else:
            filetypes = (("all files", "*.*"), ("text files", "*.txt"), ("jpeg files", "*.jpg"))
        filename =  filedialog.askopenfilename(parent = root, initialdir = initialdir, title = "Select file",\
                                                 filetypes = filetypes)
    
        return filename
    
    def get_file_to_save(filetype = None, initialdir = "/"):
        '''
        Return a selected file name (including the path) by the user 
        '''
        
        root = tk.Tk()
        #root.lift()
        root.attributes("-topmost", True)
        root.withdraw()
    
        if filetype != None:
            filetypes = (("user files", '*'+filetype), ("all files", "*.*"))
        else:
            filetypes = (("all files", "*.*"), ("text files", "*.txt"), ("jpeg files", "*.jpg"))
        filename =  filedialog.asksaveasfilename(parent = root, initialdir = initialdir, title = "Select file",\
                                                 filetypes = filetypes)
    
        return filename
    
except Exception as err:
    print(err)
    

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

def remove_files(file_list, direc = './'):
    '''
    Remove the files in `file_list` under director `direc`
    '''
    
    for f in file_list:
        if os.path.exists(direc+f):
            os.remove(direc+f)
    return

try:
    def convert(in_fig = None, out_format = None, **kwargs):
        '''
        Convert the figure to another format.
        Parameters
          **kwargs:
            quality: 90 in default
            dpi: (300, 300) in default
        '''
        
        from PIL import Image
        
        if in_fig is None:
            in_fig = get_file()
        
        if out_format is None:
            out_format = '.png'
        [path, name, ext] = fileparts(in_fig)
        print(r''+path)
        if path == '':
            out_fig = name+out_format
        else:
            out_fig = path+os.sep+name+out_format
        try:
            img = Image.open(in_fig)
            img.save(out_fig, **kwargs)
            print('Converted '+in_fig+' to '+out_fig)
        except:
            print('Error: cannot open the figure')
except Exception as err:
    print(err)
#convert("backup\2019\42\20191014\peak_current@PITZ__14.10.2019__12_03_36__.eps")
        
# # Relativistic equations

# In[5]:


def beta2gamma(beta):
    return 1.0/np.sqrt(1.0-beta**2)
def gamma2beta(gamma):
    return np.sqrt(1.0-1.0/gamma**2)
def bg2beta(bg):
    return bg/np.sqrt(1.0+bg**2)
def beta2bg(beta):
    return beta/np.sqrt(1.0-beta**2)
def bg2gamma(bg):
    return np.sqrt(1.0+bg**2)
def gamma2bg(gamma):
    return np.sqrt(gamma**2-1.0)

def kinetic2gamma(Ek, E0 = g_mec2):
    '''
    Parameters
      Ek: kinetic energy, in unit of MeV
    Returns
      gamma: relativistic factor
    '''   
    return 1.0+Ek/E0
K2G = kinetic2gamma

def kinetic2beta(Ek, E0 = g_mec2):
    '''
    Parameters
      Ek: kinetic energy, in unit of MeV
    Returns
      beta: ratio of velocity of particle to light
    ''' 
    return gamma2beta(kinetic2gamma(Ek, E0))
K2B = kinetic2beta

def momentum2gamma(pc, E0 = g_mec2):
    '''
    Parameters
      pc: momentum, MeV/c
    Returns
      gamma: relativistic factor
    '''
    return np.sqrt(pc**2/E0**2+1)
M2G = momentum2gamma

def momentum2kinetic(pc, E0 = g_mec2):
    '''
    Parameters
      pc: momentum, MeV/c
    Returns
      Ek: kinetic energy, MeV
    '''
    return E0*(np.sqrt(pc**2/E0**2+1)-1)
M2K = momentum2kinetic

def kinetic2momentum(Ek, E0 = g_mec2):
    '''
    Parameters
      Ek: kinetic energy, MeV
    Returns
      pc: momentum, MeV/c
    '''
    return np.sqrt((Ek+g_mec2)**2-g_mec2**2)
K2M = kinetic2momentum

# # Return the min or max of an array and the index
# ---
# To calculate the weighted std from a 1D array $w$, the first thing is to get the mean value. Taking the index of the array as the samples, then the mean value can be obtained from
# 
# $$ x_c = \frac{\sum_{i=0}^{n-1} i\cdot w[i]}{  \sum_{i=0}^{n-1} w[i]} $$
# 
# And the standard variation will be
# 
# $$\sigma = \sqrt{\frac{\sum_{i=0}^{n-1} (i-x_c)^2 w[i]}{\sum_{i=0}^{n-1} w[i]}}$$

def index_min(x):
    ''' 
    Parameters
      x: an 1-D array
    Returns
      [i, v]: the index of the minimum and itself
    '''
    i,v=0,x[0]
    for index,value in enumerate(x):
        if value<v:
            i,v=index,value
    return [i,v]
def index_max(x):
    ''' 
    Parameters
      x: an 1-D array
    Returns
      [i, v]: the index of the maximum and itself
    '''
    i,v=0,x[0]
    for index,value in enumerate(x):
        if value>v:
            i,v=index,value
    return [i,v]


# # Statistics with weighting factor

def weighted_sum(x, weights = None, axis = 0, returned = False):
    '''
    Parameters
      x: 1D array of the samples
      w: 1D array of the weighting factor
      returned: False by default, if True, also return the sum of weights
    Returns
      sum_x: sum of the samples
    '''
    if weights is None:
        sum_w = np.sum(x == x, axis = axis)
        #sum_w = sum_w[0]
        sum_x = np.sum(x, axis = axis)
    else:
        sum_w = np.sum(weights)
        sum_x = np.sum(x*weights, axis = axis)
    if returned:
        return sum_x, sum_w
    else:
        return sum_x
    
def weighted_mean(x, weights = None, axis = 0, returned = False):
    '''
    Parameters
      x: 1D array of the samples
      w: 1D array of the weighting factor
      returned: False by default, if True, also return the sum of weights
    Returns
      xc: mean value of the samples
    '''
    s, sum_w = weighted_sum(x, weights = weights, axis = axis, returned = True)
    xc = s/sum_w
    if returned:
        return xc, sum_w
    else:
        return xc

def weighted_std(x, weights = None, axis = 0, returned = False):
    '''
    Parameters
      x: 1D array of the samples
      w: 1D array of the weighting factor
    Returns
      sigma_x: standard viariation of the samples
    '''
    xc, sum_w = weighted_mean(x, weights = weights, axis = axis, returned = True)
    s = weighted_sum(x*x, weights = weights, axis = axis)
    
    std_x = np.sqrt(s/sum_w-xc*xc)
    if returned:
        return std_x, sum_w
    else:
        return std_x

def weighted_cov(x, y = None, weights = None, axis = 0, returned = False):
    '''
    Parameters
      x, y: 1D array of the samples
      w: 1D array of the weighting factor
    Returns
      <x*y>: covariance of x and y
    '''
    if y is None:
        y = x
    xc, sum_w = weighted_mean(x, weights = weights, axis = axis, returned = True)
    yc, sum_w = weighted_mean(y, weights = weights, axis = axis, returned = True)
    #print(xc, yc)
    s = weighted_sum(x*y, weights = weights, axis = axis)
    cov = s/sum_w-xc*yc
    #cov = np.sum((x-xc)*(y-yc)*weights, axis = axis)/np.sum(weights)
    if returned:
        return cov, sum_w
    else:
        return cov


# # Conversion between eV and nm, mm and fs

def eV2nm(eV):
    '''
    Parameters
      eV: energy of a photon in unit of electron volt
    Returns
      nm: wavelength of a photon in unit of nanometer
    '''
    return g_h*g_c/(eV*g_qe)*1e9
def nm2eV(nm):
    '''
    Parameters
      nm: wavelength of a photon in unit of nanometer
    Returns
      eV: energy of a photon in unit of electron volt
    '''
    return g_h*g_c/(nm*1e-9*g_qe)
def mm2fs(mm,beta=1):
    '''
    Parameters
      mm: bunch length in unit of minimeter
      beta: bunch velocity over light velocity in vacuum, default 1
    Returns
      fs: bunch length in unit of femtosecond
    '''
    return mm*1e-3/(beta*g_c)*1e15
def fs2mm(fs,beta=1):
    '''
    Parameters
      fs: bunch length in unit of femtosecond
      beta: bunch velocity over light velocity in vacuum, default 1
    Returns
      mm: bunch length in unit of minimeter
    '''
    return fs*1e-15*(beta*g_c)*1e3

def test():
    print (g_mec2,pc2E(50))
    print (eV2nm(1),nm2eV(800))
    print (mm2fs(1e-3),fs2mm(5000/2.355))
    gamma=kinetic2gamma(150.)
    print (gamma2bg(gamma))
# test()


# # Conversion between dB and W

def dBm2mW(x):
    ''' dBm->mW '''
    return 10**(x/10.)
def dBm2W(x):
    ''' dBm->W '''
    return 10**(x/10.)/1000.
def dB2P(x):
    ''' dB->linear for power'''
    return 10**(x/10.)
def dB2U(x):
    ''' dB->linear for voltage'''
    return 10**(x/20.)


# # Fit functions

f1 = lambda x,a,b:a*x+b
f2 = lambda x,a,b,c:a*x**2+b*x+c
f3 = lambda x,a,b,c,d:a*x**3+b*x**2+c*x+d

linear = lambda x, a, b: a+b*x
binomial= lambda x, a, b, c: a+b*x+c*x*x
gaussian = lambda x, sigma, mean, amp: amp*np.exp(-(x-mean)**2/2./sigma**2)

def lsq(x,y):
    ''' 
    Least square fitting
    Parameters
      x, y: 1-D array of coordinates to be fitted
    Returns
      [a, b, r]: a is slope, b the interseption and the closer r is to 1.0, the better the fitting is
    '''
    xc,yc=np.mean(x),np.mean(y)
    xyc=np.mean(x*y)
    x2c=np.mean(x*x)
    y2c=np.mean(y*y)
    
    a=(xyc-xc*yc)/(x2c-xc**2)
    b=yc-a*xc
    r=np.abs((xyc-xc*yc)/np.sqrt((x2c-xc**2)*(y2c-yc**2)))
    return np.array([a,b,r])


# # Linear Interpolation

def linear_interp(x, p1, p2):
    '''
    Linear interpolation. Given two points p1 = [x1, y1, ..], p2 = [x2, y2, ...], get the y-value at x
    Parameters
      x: the position of the point in the known dimension
      p1 = [x1, y1, ..]: the position of the first known point
      p2 = [x2, y2, ..]: the position of the second known point
    Returns
      [y, ...]: the position of the point in the dimension to be interpolated 
    '''
    return p2[1:]+(x-p2[0])/(p2[0]-p1[0])*(p2[1:]-p1[1:])


# # Calculate the FWHM of a distribution

from scipy.interpolate import splrep, sproot, splev

class MultiplePeaks(Exception): pass
class NoPeaksFound(Exception): pass

def cal_FWHM(x, y, k = 3):
    """
    Determine full-with-half-maximum of a peaked set of points, x and y.

    Assumes that there is only one peak present in the datasset.  The function
    uses a spline interpolation of order k.
    """
    
    xc = weighted_mean(x, y)
    sigma_x = weighted_std(x, y)
    select = (x>xc-sigma_x)*(x<xc+sigma_x)
    half_max = np.max(y[select])/2.; #print half_max
    
    #half_max = max(y)/2.0
    
    s = splrep(x, y - half_max, k = k)
    roots = sproot(s)

    if len(roots) > 2:
         r = -1
         #raise MultiplePeaks("The dataset appears to have multiple peaks, and "
         #                    "thus the FWHM can't be determined.")
    if len(roots) < 2:
        r = 0
        #raise NoPeaksFound("No proper peaks were found in the data set; likely "
        #                   "the dataset is flat (e.g. all zeros).")
    else:
        r = np.max(roots)-np.min(roots)
    return r

def get_FWHM(x, k = 3, bins = 50, weights = None):
    '''
    Parameters
      x: 1D array of the samples
      k: spline interpolation order
    Returns
      FWHM: full width at half maximum of the distribution of the samples
    '''
    try:
        counts, edges = np.histogram(x, bins = bins, weights = weights)
        centers = edges[1:]-0.5*(edges[1]-edges[0])
        counts = smooth_easy(counts, 5)
        r = cal_FWHM(centers, counts, k)
    except:
        r = -2
    return r

# # 包络与平滑

def envelope(s):
    '''
    Parameters
      s: 1-D array of source signal
    Returns
      q_u: 1-D array of upper envelope
      q_l: 1-D array of lower envelope
    '''
    # Prepend the first value of (s) to the interpolating values. 
    # This forces the model to use the same starting point for both the upper and lower envelope models.
    u_x = [0,]
    u_y = [s[0],]

    l_x = [0,]
    l_y = [s[0],]

    # Detect peaks and troughs and mark their location in u_x,u_y,l_x,l_y respectively.
    for k in np.arange(1,len(s)-1):
        if (np.sign(s[k]-s[k-1])==1) and (np.sign(s[k]-s[k+1])==1):
            u_x.append(k)
            u_y.append(s[k])

        if (np.sign(s[k]-s[k-1])==-1) and ((np.sign(s[k]-s[k+1]))==-1):
            l_x.append(k)
            l_y.append(s[k])

    # Append the last value of (s) to the interpolating values. 
    # This forces the model to use the same ending point for both the upper and lower envelope models.
    u_x.append(len(s)-1)
    u_y.append(s[-1])

    l_x.append(len(s)-1)
    l_y.append(s[-1])

    # Fit suitable models to the data. Here I am using cubic splines, similarly to the MATLAB example given in the question.
    u_fit = interp1d(u_x, u_y, kind = 'linear', bounds_error = False, fill_value=0.0)
    l_fit = interp1d(l_x, l_y, kind = 'linear', bounds_error = False, fill_value=0.0)

    # Evaluate each model over the domain of (s)
    u_env = u_fit(np.arange(len(s)))
    l_env = l_fit(np.arange(len(s)))

    return u_env, l_env

def smooth_easy(y, box_pts):
    '''
    Smooth the data points using a moving box
    Parameters
      y: 1-D array to be smoothed
      box_pts: number of points in the moving box
    Returns
      y_smooth: smoothed 1-D array
    '''
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def smooth(x,window_len=11,window='hanning'):
    """
    smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise (ValueError, "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise (ValueError, "Input vector needs to be bigger than window size.")


    if window_len < 3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise (ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s = np.r_[x[window_len-1:0:-1], x, x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(int(w/w.sum()), s, mode = 'valid')
    return y[int(window_len/2):-int(window_len/2)]

