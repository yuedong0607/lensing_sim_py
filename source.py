import numpy as np
from math import pi

def get_source():

    # Randomly generate uniformly distributed points
    Nsource = 500
    xsmin = -1.     # Mpc/h
    xsmax = 1.      # Mpc/h
    ysmin = -1.     # Mpc/h
    ysmax = 1.      # Mpc/h

    xs = np.random.uniform(xsmin, xsmax, size = Nsource)
    ys = np.random.uniform(ysmin, ysmax, size = Nsource)
    zs = np.random.uniform(0.3, 0.7, size = Nsource)

    ind = (np.sqrt(xs**2 + ys**2) <= 1.5)
    xs = xs[ind]
    ys = ys[ind]
    zs = zs[ind]

    return xs, ys, zs

def get_source_around(xl, yl):

    # Randomly generate uniformaly distributed points around a particular lens
    Nsource = 500
    rc = 1. #Mpc/h
    phi = 2*pi*np.random.uniform(size = Nsource)
    r = np.sqrt(np.abs(np.random.normal(scale = 0.15, size = Nsource)))
    xs = (rc*r) * np.cos(phi) + xl
    ys = (rc*r) * np.sin(phi) + yl
    zs = np.random.uniform(0.3, 0.7, size = Nsource)

    return xs, ys, zs

def get_source_grid(xmin, xmax, ymin, ymax, nx, ny):
    
    Nsource = nx * ny
    xedge = np.linspace(xmin, xmax, nx+1)
    yedge = np.linspace(ymin, ymax, ny+1)
    xmid = (xedge[1:] + xedge[:-1])/2.
    ymid = (yedge[1:] + yedge[:-1])/2.
    x, y = np.meshgrid(xmid, ymid, indexing = 'xy')
    xs = x.ravel()
    ys = y.ravel()
    
    # Random source redshift
    #zs = np.random.uniform(0.3, 0.7, size = Nsource)
    # Fix source redshift(0.5) for plotting whisker map
    zs = np.ones(Nsource) * 0.5

    return xs, ys, zs
