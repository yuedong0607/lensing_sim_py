from math import pi
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from configuration import *
import jc_lib as jc

def random_three_vector(Nsub):
    phi = np.random.uniform(0, np.pi*2, Nsub)
    costheta = np.random.uniform(-1, 1, Nsub)

    theta = np.arccos( costheta)
    x = np.sin( theta) * np.cos( phi)
    y = np.sin( theta) * np.sin( phi)
    z = np.cos( theta )

    return (x, y, z)

def random_two_vector(Nsub):
    phi = np.random.uniform(0., np.pi*2, Nsub)
    costheta = np.random.uniform(-1, 1, Nsub)
    theta = np.arccos( costheta)
    x = np.cos(theta)
    y = np.sin(theta)

    return (x, y)

def get_zl():
    # TODO
    zl = 0.2
    return zl

def getM200():
    # TODO
    M200 = 1e14
    return M200

def getNsub(M200):
    # TODO (maybe implementing mass richness relation)
    maxNsub = 1
    Nsub = maxNsub
    #Nsub = np.random.random_integers(maxNsub)
    return Nsub

def getsubhalo(Nsub, R200):
    # TODO
    #msub = np.random.uniform(1e11, 1e12, Nsub) # Msol /h
    msub = np.random.uniform(1e12, 1e13, Nsub) # Msol /h

    # Calculate the location of the subhaloes
    # (here we forced them to be at least 0.1Mpc away from the center)
    #r = np.random.uniform(0.1, R200, Nsub)
    r = np.random.uniform(0.4, 0.6, Nsub)

    # A random subhalo in 3D
    #x, y, z = r * random_three_vector(Nsub)

    # A random subhlo in 2D
    x, y = r * random_two_vector(Nsub)

    # TODO implementing truncation radius

    return x, y, msub
    
    
    
