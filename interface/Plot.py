from .Misc import *

### Generate 1D and 2D statistics from raw samples
def hist_centered(x, bins = None, range = None, weights = None, **kwargs):
    '''
    Returns
      counts, centers: 1d arrays including each bin's count and center, respectively
    '''
    if bins is None:
        bins = 100

    counts, edges = np.histogram(x, bins = bins, range = range, weights = weights, **kwargs)
    centers = edges[1:]-0.5*(edges[1]-edges[0])
    return counts, centers

def hist2d_scattered(x, y, bins = [100, 100], range = None, weights = None, **kwargs):
    '''
    Returns
      xs, ys, cs: 1d arrays, (xs, ys) are the coordinates and cs the number of particles at (xs, ys)
    '''
    counts, xedges, yedges = np.histogram2d(x, y, bins = bins, range = range, weights = weights, **kwargs)
    xs = np.array([xedges[i]+0.5*(xedges[1]-xedges[0]) for i in np.arange(len(xedges)-1) for j in np.arange(len(yedges)-1)])
    ys = np.array([yedges[j]+0.5*(yedges[1]-yedges[0]) for i in np.arange(len(xedges)-1) for j in np.arange(len(yedges)-1)])
    cs = np.array([counts[i,j] for i in np.arange(len(xedges)-1) for j in np.arange(len(yedges)-1)])
    select = (cs>0)
    return xs[select], ys[select], cs[select]

def hist2d_contourf(x, y, weights = None, bins = [100, 100], range = None):
    '''
    Returns
      xs, ys, cs: 1d arrays, (xs, ys) are the coordinates and cs the number of particles at (xs, ys)
    '''
    counts, xedges, yedges = np.histogram2d(x, y, bins = bins, range = range, weights = weights)
    xs = np.array([xedges[i]+0.5*(xedges[1]-xedges[0]) for i in np.arange(len(xedges)-1)])
    ys = np.array([yedges[j]+0.5*(yedges[1]-yedges[0]) for j in np.arange(len(yedges)-1)])
    ax, ay = np.meshgrid(xs, ys, indexing = 'ij')
    return ax, ay, counts


### plot 1D distribution
def zdist(x, bins = 100, range = None, weights = None, smooth = False, fig_ext = '.eps', plot = True, **fig_kw):
    '''
    Parameters
      x: longitudinal coordinates in m
    '''
    
    counts, centers = hist_centered(x, bins = bins, range = range, weights = weights)
    dz = centers[1]-centers[0]
    
    currents = np.abs(counts/dz*g_c); Ipeak0 = np.max(currents)
    if Ipeak0 > 1000:
        ratio = 1e-3
        ylabel = r'$I$ (kA)'
    else:
        ratio = 1
        ylabel = r'$I$ (A)'
        
    if plot:
        fig, ax = plt.subplots(**fig_kw)
        ax.plot(centers/g_c*1e12, currents*ratio, '-')
    
    if smooth:
        currents = smooth_easy(currents, 5); Ipeak1 = np.max(currents)
    else:
        Ipeak1 = Ipeak0
        
    if plot:
        if smooth:
            ax.plot(centers/g_c*1e12, currents*ratio, '--')
    
        ax.grid()
        #ax.set_xlabel(r'$\xi$ (mm)')
        ax.set_xlabel(r'$t_0$ (ps)')
        ax.set_ylabel(ylabel)
        #ax.set_xlim(-15, 15)
        #ax.set_ylim(0, 200)
        fig.savefig('dQ_dt-t'+fig_ext)
    return [Ipeak1, Ipeak0]

def xdist(x, bins = 50, range = None, weights = None, fig_ext = '.eps', **fig_kw):
    '''
    Parameters
      x: horizontal coordinates in mm
    '''
    fig, ax = plt.subplots(**fig_kw)
    h = ax.hist(x, bins = bins, range = range, weights = weights, histtype = 'step')
    #ax.grid()
    ax.set_xlabel(r'$x$ (mm)')
    if range == None:
        ax.set_ylabel(r'd$Q/$d$x$ (arb. unit)')
    else:
        ax.set_ylabel(r'd$Q/$d$x$ (pC/mm)')
    #ax.set_xlim(-15, 15)
    #ax.set_ylim(0, 200)
    fig.savefig('dQ_dx-x'+fig_ext)
    return

def ydist(y, bins = 50, range = None, weights = None, fig_ext = '.eps', **fig_kw):
    '''
    Parameters
      x: vertical coordinates in mm
    '''
    fig, ax = plt.subplots(**fig_kw)
    h = ax.hist(x, bins = bins, range = range, weights = weights, histtype = 'step')
    #ax.grid()
    ax.set_xlabel(r'$y$ (mm)')
    if range == None:
        ax.set_ylabel(r'd$Q/$d$y$ (arb. unit)')
    else:
        ax.set_ylabel(r'd$Q/$d$y$ (pC/mm)')
    #ax.set_xlim(-15, 15)
    #ax.set_ylim(0, 200)
    fig.savefig('dQ_dy-y'+fig_ext)
    return

def pdist(z, bins = 50, range = None, weights = None, fig_ext = '.eps', **fig_kw):
    '''
    Parameters
      p: momentum in MeV/c
    '''
    fig, ax = plt.subplots(**fig_kw)
    h = ax.hist(z, bins = bins, range = range, weights = weights, histtype = 'step')
    #ax.grid()
    ax.set_xlabel(r'$P$ (MeV/c)')
    if range == None:
        ax.set_ylabel(r'd$Q/$d$P$ (arb. unit)')
    else:
        ax.set_ylabel(r'd$Q/$d$P$ (pC/mm)')
    #ax.set_xlim(-15, 15)
    #ax.set_ylim(0, 200)
    fig.savefig('dQ_dP-P'+fig_ext)
    return

def tdist(t, bins = 50, range = None, weights = None, fig_ext = '.eps', **fig_kw):
    '''
    Parameters
      t: temperal coordinates in ps
    '''
    fig, ax = plt.subplots(**fig_kw)
    h = ax.hist(t, bins = bins, range = range, weights = weights, histtype = 'step')
    #ax.grid()
    ax.set_xlabel(r'$t$ (ps)')
    if range == None:
        ax.set_ylabel(r'd$Q/$d$t$ (arb. unit)')
    else:
        ax.set_ylabel(r'd$Q/$d$t$ (pC/ps)')
    #ax.set_xlim(-15, 15)
    #ax.set_ylim(0, 200)
    fig.savefig('dQ_dt-t'+fig_ext)
    return


### plot slice parameters
def plot_slice(fname, fig_ext = '.eps'):
    
    d1 = astra2slice(fname, 'slice@'+fname+'.dat', nc = 1)
    
    zz = (d1[:,0]-np.mean(d1[:,0]))*1e3
    II = d1[:,-2]; print('Peak current is ', II.max(), ' A')
    ee = d1[:,3]*1e6

    fig, ax = plt.subplots()
    ax.plot(zz, II, 'r-')
    ax.plot([-1], [-1], 'b--')
    ax.set_xlabel(r'$\xi$ (mm)')
    ax.set_ylabel(r'Current (A)', color = 'r')
    #ax.set_xlim(-5, 5)
    #ax.set_ylim(0, 200)
    ax.tick_params(axis='y', colors='r')
    ax.grid()
    ax.legend(['current', 'emittance'])

    ax1 = ax.twinx()
    ax1.plot(zz, ee, 'b--')
    ax1.set_ylabel(u_emi_x, color = 'b')
    #ax1.set_ylim(0., 10)
    ax1.tick_params(axis='y', colors='b')

    fig.savefig('slice@'+fname+fig_ext)
    

### plot 2D phase space distributions
def plot2d(x, y, weights = None, xlabel = r'$x$ (mm)', ylabel = r'$y$ (mm)', figname = 'xy2d',\
           fig_ext = '.eps', vmin = 0.001, vmax = 1., bins = None, extent = None, returned = False, **kwargs):
    '''
    Parameters
      bins: None or [xbins, ybins], set the sampling frequencies for x and y
      extent: None or [xmin, xmax, ymin, ymax], set the lower limit and upper limit of the x-axis and y-axis
    '''
    if bins == None:
        xbins = ybins = 100
    else:
        xbins, ybins = bins
    if extent == None:
        xs = np.sort(x[::5])
        xavg = np.mean(xs); dx = xs[-int(0.005*len(xs))]-xs[int(0.005*len(xs))]
        xavg = np.mean(x); dx = np.max(x)-np.min(x)
        xmin, xmax = xavg-dx*0.75, xavg+dx*0.75
        
        ys = np.sort(y[::5])
        yavg = np.mean(ys); dy = ys[-int(0.005*len(ys))]-ys[int(0.005*len(ys))]
        yavg = np.mean(y); dy = np.max(y)-np.min(y)
        ymin, ymax = yavg-dy*0.75, yavg+dy*0.75
    else:
        xmin, xmax, ymin, ymax = extent
    print('extent: ', xmin, xmax, ymin, ymax)
    
    fig = plt.figure(figsize=(4, 3.375))
    ax1 = fig.add_axes([0.18, 0.15, 0.675, 0.15]); ax1.patch.set_alpha(0)
    ax3 = fig.add_axes([0.18, 0.15, 0.15, 0.80]); ax3.patch.set_alpha(0)
    ax2 = fig.add_axes([0.18, 0.15, 0.75, 0.80]); ax2.patch.set_alpha(0)

    cnts, cens = hist_centered(x, bins=xbins, range=(xmin,xmax), weights=weights)
    ax1.plot(cens, smooth_easy(cnts, 8), 'b-')

    # vmin is the lower limit of the colorbar, and is free to change
    # xs, ys, cs = hist2d_scattered(x, y, bins=[xbins, ybins], range=[[xmin, xmax],[ymin, ymax]], weights=weights)
    # cax = ax2.scatter(xs, ys, s=2, c=cs/np.max(cs), edgecolor='', vmin=vmin, norm=mpl.colors.LogNorm())
   
    ax, ay, signal = hist2d_contourf(x, y, bins=[xbins, ybins], range=[[xmin, xmax],[ymin, ymax]], weights=weights)
    signal = np.abs(signal)
    signal = signal/np.max(signal)
    v = np.linspace(vmin, vmax, 99)
    cax = ax2.contourf(ax, ay, signal, v)
    
    cbar = fig.colorbar(cax, fraction=0.09, pad=0.01, format = '%.1f')
    cbar.set_ticks([vmin, 1.0])
    cbar.ax.tick_params(labelsize=11, pad=0)

    cnts, cens = hist_centered(y, bins=ybins, range=(ymin,ymax), weights=weights)
    ax3.plot(smooth_easy(cnts, 4), cens, 'b-')

    ax1.axis('off')
    ax1.set_xlim(xmin, xmax)

    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel, labelpad=-1)
    ax2.minorticks_on()
    ax2.grid()

    ax3.axis('off')
    ax3.set_ylim(ymin, ymax)

    fig.savefig(figname+fig_ext)
    
    if returned:
        return [ax, ay, signal]
    return

# In[41]:

def plot2dx2(x, y, x1, y1, weights = None, xlabel = r'$x$ (mm)', ylabel = r'$y$ (mm)', figname = 'xy2d',\
           fig_ext = '.eps', vmin = 0.01, vmax = 1.0, bins = None, extent = None, **kwargs):
    '''
    Parameters
      bins: None or [xbins, ybins], set the sampling frequencies for x and y
      extent: None or [xmin, xmax, ymin, ymax], set the lower limit and upper limit of the x-axis and y-axis
    '''
    if bins == None:
        xbins = ybins = 100
    else:
        xbins, ybins = bins
    if extent == None:
        xavg = np.mean(x); dx = np.max(x)-np.min(x)
        xmin, xmax = xavg-dx*0.75, xavg+dx*0.75
        yavg = np.mean(y); dy = np.max(y)-np.min(y)
        ymin, ymax = yavg-dy*0.75, yavg+dy*0.75
    else:
        xmin, xmax, ymin, ymax = extent
    
    fig = plt.figure(figsize=(4, 3.375))
    ax1 = fig.add_axes([0.18, 0.15, 0.675, 0.15]); ax1.patch.set_alpha(0)
    ax3 = fig.add_axes([0.18, 0.15, 0.15, 0.80]); ax3.patch.set_alpha(0)
    ax2 = fig.add_axes([0.18, 0.15, 0.75, 0.80]); ax2.patch.set_alpha(0)

    cnts, cens = hist_centered(x, bins=xbins, range=(xmin,xmax), weights=weights)
    #ax1.plot(cens, smooth_easy(cnts, 8), 'b-')

    #xs, ys, cs = hist2d_scattered(x, y, bins=[xbins, ybins], range=[[xmin, xmax],[ymin, ymax]], weights=weights)
    ## vmin is the lower limit of the colorbar, and is free to change
    #cax = ax2.scatter(xs, ys, s=2, c=cs/np.max(cs), edgecolor='', vmin=vmin, norm=mpl.colors.LogNorm())
    
    ax, ay, signal = hist2d_contourf(x, y, bins=[xbins, ybins], range=[[xmin, xmax],[ymin, ymax]], weights=weights)
    signal = np.abs(signal)
    signal = signal/np.max(signal)
    v = np.linspace(vmin, vmax, 99)
    cax = ax2.contourf(ax, ay, signal, v)
    
    #xs1, ys1, cs1 = hist2d_scattered(x1, y1, bins=[xbins, ybins], range=[[xmin, xmax],[ymin, ymax]], weights=weights)
    ## vmin is the lower limit of the colorbar, and is free to change
    #cax = ax2.scatter(xs1, ys1, s=2, c=cs1/np.max(cs1), edgecolor='', vmin=vmin, norm=mpl.colors.LogNorm())
    
    ax, ay, signal = hist2d_contourf(x1, y1, bins=[xbins, ybins], range=[[xmin, xmax],[ymin, ymax]], weights=weights)
    signal = np.abs(signal)
    signal = signal/np.max(signal)
    v = np.linspace(vmin, vmax, 99)
    cax = ax2.contourf(ax, ay, signal, v)
    
    cbar = fig.colorbar(cax, fraction=0.09, pad=0.01)
    cbar.set_ticks([vmin, 1.0])
    cbar.ax.tick_params(labelsize=11, pad=0)

    cnts, cens = hist_centered(y, bins=ybins, range=(ymin,ymax), weights=weights)
    #ax3.plot(smooth_easy(cnts, 4), cens, 'b-')

    ax1.axis('off')
    ax1.set_xlim(xmin, xmax)

    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel, labelpad=0)
    ax2.minorticks_on()
    ax2.grid()

    ax3.axis('off')
    ax3.set_ylim(ymin, ymax)

    fig.savefig(figname+fig_ext)
    return

def plot2d_xpx_ypy(fname = 'ast.0500.001', fig_ext = '.eps', **kwargs):
    '''
    Parameters:
      **kwargs: bins = (xbins, ybins), extent = [xmin, xmax, ymin, ymax]
    Returns:
      A saved figure
    '''

    beam = np.loadtxt(fname)
    beam[1:,2] += beam[0, 2]
    beam[1:,5] += beam[0, 5]
    
    select = (beam[:,9] == 5)|(beam[:,9] == -1)|(beam[:,9] == -3)
    beam = beam[select]
    
    x, y, w = beam[:,0]*1e3, beam[:,3]/1e3, -beam[:,7]*1e3
    x1, y1 = beam[:,1]*1e3, beam[:,4]/1e3
    plot2dx2(x = x, y = y, x1 = x1, y1 = y1, weights = w, xlabel = r'$x$ or $y$ (mm)',\
             ylabel = r'$p_x$ or $p_y$ (keV/$c$)', figname = 'xpx+ypy', fig_ext = fig_ext, **kwargs)
    
    return

def plot2d_all(fname = 'ast.0500.001', dist = None, fig_ext = '.eps', **kwargs):
    '''
    Parameters:
      **kwargs: bins = (xbins, ybins), extent = [xmin, xmax, ymin, ymax]
    Returns:
      A saved figure
    '''
    if dist is None:
        dist = np.loadtxt(fname)
        dist[1:,2] += dist[0, 2]
        dist[1:,5] += dist[0, 5]
    
    select = (dist[:,9] > 0)|(dist[:,9] == -1)|(dist[:,9] == -3)
    dist = dist[select]
    
    #tdist(t = beam[:,6]*1e3, bins = 60, weights = -beam[:,7]*1e3*2, range = (-15, 15))
    
    x, y, w = dist[:,0]*1e3, dist[:,1]*1e3, -dist[:,7]*1e3
    plot2d(x = x, y = y, weights = w, fig_ext = fig_ext, **kwargs)

    x, y, w = dist[:,0]*1e3, dist[:,3]/1e6, -dist[:,7]*1e3
    plot2d(x = x, y = y, weights = w, ylabel = r'$p_x$ (MeV/$c$)', figname = 'xpx', fig_ext = fig_ext, **kwargs)

    x, y, w = dist[:,2]*1e3, dist[:,5]/1e6, -dist[:,7]*1e3
    plot2d(x = x, y = y, weights = w, xlabel = r'$z$ (mm)', ylabel = r'$p_z$ (MeV/$c$)', figname = 'zpz', fig_ext = fig_ext, **kwargs)
    
    return

def plot2d_xy(fname = 'ast.0500.001', dist = None, fig_ext = '.eps', **kwargs):
    '''
    Parameters:
      **kwargs: bins = (xbins, ybins), extent = [xmin, xmax, ymin, ymax]
    Returns:
      A saved figure
    '''

    if dist is None:
        dist = np.loadtxt(fname)
        dist[1:,2] += dist[0, 2]
        dist[1:,5] += dist[0, 5]
    
    select = (dist[:,9] > 0)|(dist[:,9] == -1)|(dist[:,9] == -3)
    dist = dist[select]
    
    x, y, w = dist[:,0]*1e3, dist[:,1]*1e3, -dist[:,7]*1e3
    plot2d(x = x, y = y, weights = w, fig_ext = fig_ext, **kwargs)
   
    return

def plot2d_xpx(fname = 'ast.0500.001', dist = None, fig_ext = '.eps', **kwargs):
    '''
    Parameters:
      **kwargs: bins = (xbins, ybins), extent = [xmin, xmax, ymin, ymax]
    Returns:
      A saved figure
    '''

    if dist is None:
        dist = np.loadtxt(fname)
        dist[1:,2] += dist[0, 2]
        dist[1:,5] += dist[0, 5]
    
    select = (dist[:,9] > 0)|(dist[:,9] == -1)|(dist[:,9] == -3)
    dist = dist[select]
    
    x, y, w = dist[:,0]*1e3, dist[:,3]/1e3, -dist[:,7]*1e3
    plot2d(x = x, y = y, weights = w, ylabel = r'$p_x$ (keV/$c$)', figname = 'xpx', fig_ext = fig_ext, **kwargs)
    
    return

def plot2d_xxp(fname = 'ast.0500.001', dist = None, fig_ext = '.eps', **kwargs):
    '''
    Parameters:
      **kwargs: bins = (xbins, ybins), extent = [xmin, xmax, ymin, ymax]
    Returns:
      A saved figure
    '''

    if dist is None:
        dist = np.loadtxt(fname)
        dist[1:,2] += dist[0, 2]
        dist[1:,5] += dist[0, 5]
    
    select = (dist[:,9] > 0)|(dist[:,9] == -1)|(dist[:,9] == -3)
    dist = dist[select]
    
    x, y, w = dist[:,0]*1e3, dist[:,3]/dist[:,5]*1e3, -dist[:,7]*1e3
    plot2d(x = x, y = y, weights = w, ylabel = r'$x^{\prime}$ (mrad)', figname = 'xxp', fig_ext = fig_ext, **kwargs)
    
    return

def plot2d_ypy(fname = 'ast.0500.001', dist = None, fig_ext = '.eps', **kwargs):
    '''
    Parameters:
      **kwargs: bins = (xbins, ybins), extent = [xmin, xmax, ymin, ymax]
    Returns:
      A saved figure
    '''

    if dist is None:
        dist = np.loadtxt(fname)
        dist[1:,2] += dist[0, 2]
        dist[1:,5] += dist[0, 5]
    
    select = (dist[:,9] > 0)|(dist[:,9] == -1)|(dist[:,9] == -3)
    dist = dist[select]
    
    x, y, w = dist[:,1]*1e3, dist[:,4]/1e3, -dist[:,7]*1e3
    plot2d(x = x, y = y, weights = w, xlabel = r'$y$ (mm)', ylabel = r'$p_y$ (keV/$c$)', figname = 'ypy', fig_ext = fig_ext, **kwargs)
    
    return

def plot2d_yyp(fname = 'ast.0500.001', dist = None, fig_ext = '.eps', **kwargs):
    '''
    Parameters:
      **kwargs: bins = (xbins, ybins), extent = [xmin, xmax, ymin, ymax]
    Returns:
      A saved figure
    '''

    if dist is None:
        dist = np.loadtxt(fname)
        dist[1:,2] += dist[0, 2]
        dist[1:,5] += dist[0, 5]
    
    select = (dist[:,9] > 0)|(dist[:,9] == -1)|(dist[:,9] == -3)
    dist = dist[select]
    
    x, y, w = dist[:,1]*1e3, dist[:,4]/dist[:,5]*1e3, -dist[:,7]*1e3
    plot2d(x = x, y = y, weights = w, xlabel = r'$y$ (mm)', ylabel = r'$y^{\prime}$ (mrad)', figname = 'yyp', fig_ext = fig_ext, **kwargs)
    
    return

def plot2d_zpz(fname = 'ast.0500.001', dist = None, fig_ext = '.eps', **kwargs):
    '''
    Parameters:
      **kwargs: bins = (xbins, ybins), extent = [xmin, xmax, ymin, ymax]
    Returns:
      A saved figure
    '''

    if dist is None:
        dist = np.loadtxt(fname)
        #dist[1:,2] += dist[0, 2]
        dist[0, 2] = 0
        dist[1:,5] += dist[0, 5]
    
    select = (dist[:,9] > 0)
    if np.sum(select)>0:
        pass
    else:
        select = (dist[:,9] == -1)|(dist[:,9] == -3)
        if np.sum(select) > 0:
            dist[:,2] = dist[:,6]*1e-9*g_c
        else:
            return
    
    dist = dist[select]
    
    x, y, w = dist[:,2]*1e3, dist[:,5]/1e6, -dist[:,7]*1e3
    plot2d(x = x, y = y, weights = w, xlabel = r'$\xi$ (mm)', ylabel = r'$p_z$ (MeV/$c$)', figname = 'zpz', fig_ext = fig_ext, **kwargs)
    
    return

### plot evolutions of beam parameters along z
def plot_avg_xy(prefix = 'ast', suffix = '001', fig_ext = '.eps', extent = None, **fig_kw):
    xemit = np.loadtxt(prefix+'.Xemit.'+suffix)
    yemit = np.loadtxt(prefix+'.Yemit.'+suffix)
    
    fig, ax = plt.subplots(**fig_kw)
    ax.plot(xemit[:,0], xemit[:,2], '-')
    ax.plot(yemit[:,0], yemit[:,2], '-')
    ax.grid()
    ax.set_xlabel(r'$z$ (m)')
    ax.set_ylabel(r'$x, y$ (mm)')
    if extent is not None:
        xmin, xmax, ymin, ymax = extent
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    ax.legend(['$x$', '$y$'])
    fig.savefig('avg_xy-z'+fig_ext)
    return

def plot_rms_xy(prefix = 'ast', suffix = '001', fig_ext = '.eps', extent = None, **fig_kw):
    xemit = np.loadtxt(prefix+'.Xemit.'+suffix)
    yemit = np.loadtxt(prefix+'.Yemit.'+suffix)
    
    fig, ax = plt.subplots(**fig_kw)
    ax.plot(xemit[:,0], xemit[:,3], '-')
    ax.plot(yemit[:,0], yemit[:,3], '-')
    ax.grid()
    ax.set_xlabel(r'$z$ (m)')
    ax.set_ylabel(r'$\sigma_{x/y}$ (mm)')
    if extent is not None:
        xmin, xmax, ymin, ymax = extent
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    ax.legend(['$x$', '$y$'])
    fig.savefig('rms_xy-z'+fig_ext)
    return

def plot_rms_xyz(prefix = 'ast', suffix = '001', fig_ext = '.eps', extent = None, **fig_kw):
    xemit = np.loadtxt(prefix+'.Xemit.'+suffix)
    yemit = np.loadtxt(prefix+'.Yemit.'+suffix)
    zemit = np.loadtxt(prefix+'.Zemit.'+suffix)
    
    fig, ax = plt.subplots(**fig_kw)
    ax.plot(xemit[:,0], xemit[:,3], '-')
    ax.plot(yemit[:,0], yemit[:,3], '-')
    ax1 = ax.twinx()
    ax1.plot(zemit[:,0], zemit[:,3], 'g-')
    ax1.set_ylim()
    ax1.set_ylabel(u_rms_z)
    ax.grid()
    ax.set_xlabel(r'$z$ (m)')
    ax.set_ylabel(r'$\sigma_{x,y}$ (mm)')
    if extent is not None:
        xmin, xmax, ymin, ymax = extent
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    ax.legend(['$x$', '$y$', '$z$'])
    fig.savefig('rms_xyz-z'+fig_ext)
    return

def plot_emi_xy(prefix = 'ast', suffix = '001', fig_ext = '.eps', extent = None, **fig_kw):
    xemit = np.loadtxt(prefix+'.Xemit.'+suffix)
    yemit = np.loadtxt(prefix+'.Yemit.'+suffix)
    
    fig, ax = plt.subplots(**fig_kw)
    ax.plot(xemit[:,0], xemit[:,5], '-')
    ax.plot(yemit[:,0], yemit[:,5], '-')
    ax.grid()
    ax.set_xlabel(r'$z$ (m)')
    ax.set_ylabel(r'$\varepsilon_{n, x/y}$ (mm)')
    if extent is not None:
        xmin, xmax, ymin, ymax = extent
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    ax.legend(['$x$', '$y$'], loc = 'lower right')
    fig.savefig('emi_xy-z'+fig_ext)
    return

def plot_kin(prefix = 'ast', suffix = '001', fig_ext = '.eps', extent = None, **fig_kw):
    '''
    Plot Ek and Delta_Ek 
    '''
    zemit = np.loadtxt(prefix+'.Zemit.'+suffix)
    
    fig, ax = plt.subplots(**fig_kw)
    ax.plot(zemit[:,0], zemit[:,2], 'r-', label = r'$E_k$')
    ax2 = ax.twinx()
    ax2.plot(zemit[:,0], zemit[:,4]/1e3/zemit[:,2]*100, 'b-', label = r'$\sigma_E/E$')
    ax.grid()
    ax.set_xlabel(r'$z$ (m)')
    ax.set_ylabel(u_kinetic)
    ax2.set_ylabel(u_kinetic_rel)
    if extent is not None:
        xmin, xmax, ymin, ymax = extent
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    fig.savefig('kin-z'+fig_ext)
    return

def plot_kin2(prefix = 'ast', suffix = '001', fig_ext = '.eps', extent = None, **fig_kw):
    '''
    Plot Ek and Cor_zEk 
    '''
    zemit = np.loadtxt(prefix+'.Zemit.'+suffix)
    
    fig, ax = plt.subplots(**fig_kw)
    ax.plot(zemit[:,0], zemit[:,6], 'r-', label = r'$\langle z\cdot E_{\rm k}\rangle$')
    ax2 = ax.twinx()
    ax2.plot(zemit[:,0], zemit[:,4]/1e3/zemit[:,2]*100, 'b-', label = r'$\sigma_E/E$')
    ax.grid()
    ax.set_xlabel(r'$z$ (m)')
    ax.set_ylabel(r'$\langle z\cdot E_{\rm k}\rangle$')
    ax2.set_ylabel(u_kinetic_rel)
    if extent is not None:
        xmin, xmax, ymin, ymax = extent
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    fig.savefig('kin-z'+fig_ext)
    return

def plot_zEk(prefix = 'ast', suffix = '001', fig_ext = '.eps', extent = None, **fig_kw):
    '''
    Plot only Cor_zEk
    '''
    fig, ax = plt.subplots(**fig_kw)
    
    zemit = np.loadtxt(prefix+'.Zemit.'+suffix)
    
    ax.plot(zemit[:,0], zemit[:,6], 'r-', label = r'$\langle z\cdot E_{\rm k}\rangle$')
    ax.grid()
    ax.set_xlabel(r'$z$ (m)')
    ax.set_ylabel(r'$\langle z\cdot E_{\rm k}\rangle$ (keV)', color = 'r')
    
    ax.tick_params(axis='y', colors='r', which = 'both')
    # ax.spines['left'].set_color('blue')
    
    ax2 = ax.twinx()
    ax2.plot(zemit[:,0], zemit[:,4]/1e3/zemit[:,2]*100, 'b-', label = r'$\sigma_E/E$')
    ax2.set_ylabel(u_kinetic_rel, color = 'b')
    
    ax2.tick_params(axis='y', colors='b', which = 'both')
    ax2.spines['left'].set_color('red')
    ax2.spines['right'].set_color('blue')
    
    fig.savefig('zEkin-z'+fig_ext)
    return

def plot_all_z(prefix = 'ast', suffix = '001', fig_ext = '.eps', extent = None, **fig_kw):

    fig, [ax1, ax2, ax3] = plt.subplots(nrows = 3, **fig_kw)
    
    xemit = np.loadtxt(prefix+'.Xemit.'+suffix)
    yemit = np.loadtxt(prefix+'.Yemit.'+suffix)
    zemit = np.loadtxt(prefix+'.Zemit.'+suffix)
    
    ax1.plot(xemit[:,0], xemit[:,3], '-')
    ax1.plot(yemit[:,0], yemit[:,3], '-')
    ax1.grid()
    ax1.set_xlabel(r'$z$ (m)')
    ax1.set_ylabel(r'$\sigma_{x/y}$ (mm)')
    ax1.legend(['$x$', '$y$'])
    ax1.set_ylim(0, 10)
    
    ax12 = ax1.twinx()
    ax12.plot(zemit[:,0], zemit[:,3], '-')
    ax12.set_ylabel(r'$\sigma_{z}$ (mm)')
    ax12.set_ylim(0, 1)
    
    ax2.plot(xemit[:,0], xemit[:,5], '-')
    ax2.plot(yemit[:,0], yemit[:,5], '-')
    ax2.grid()
    ax2.set_xlabel(r'$z$ (m)')
    ax2.set_ylabel(r'$\varepsilon_{n, x/y}$ (mm)')
    ax2.legend(['$x$', '$y$'])
    
    ax3.plot(zemit[:,0], zemit[:,6], 'r-', label = r'$\langle z\cdot E_{\rm k}\rangle$')
    ax32 = ax3.twinx()
    ax32.plot(zemit[:,0], zemit[:,4]/1e3/zemit[:,2]*100, 'b-', label = r'$\sigma_E/E$')
    ax3.grid()
    ax3.set_xlabel(r'$z$ (m)')
    ax3.set_ylabel(r'$\langle z\cdot E_{\rm k}\rangle$ (keV)')
    ax32.set_ylabel(u_kinetic_rel)
    
    fig.savefig('all-along-z'+fig_ext)
    return