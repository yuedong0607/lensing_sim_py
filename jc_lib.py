from math import pi
import numpy as np
import scipy.integrate as integrate
import scipy.special as sp

from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates as mapc

# Planck cosmology
def cosmo():

    om_m = 0.3175
    om_L = 0.6825
    h = 0.6711
    n_s = 0.9624
    sig_8 = 0.8344
    om_b = 0.0490

    return[om_m, om_L, h, n_s, sig_8, om_b]

om_m, om_L, h, n_s, sig_8, om_b = cosmo()

################################################
def omega_matter(N):
    a = om_m * np.exp(-3.*N)
    return a / (a + om_L)

def linear_growth(N):
    d0 = [np.exp(N[0]), np.exp(N[0])]

    def func(d, N):
        omega = omega_matter(N)
        return [(3.*omega/2.-2.) * d[0] + 3.*omega/2.*d[1], d[0]]

    temp = integrate.odeint(func, d0, N)
    return temp[:,1]

# Critical energy density at z (M_sol h^2 / Mpc^3)
def critical_density(z):
    kg_per_Msol = 1.98892e30
    Msol_per_mw = 1e12
    km_per_mpc = 3.08568e19
    
    G = 6.67300e-11            # m^3 / (kg s^2)
    H_0 = 100.                 # km h / (s Mpc)

    convert_rho = 1e9 * km_per_mpc / (kg_per_Msol)
    rho_cr = 3./(8.*pi) * (H_0**2 / G) * convert_rho
    
    return rho_cr * (om_m * (1.+z)**3 + om_L)

################################
def comoving_dist(z):
    c = 3e5               # km/s
    H0 = 100.             # km h / (s Mpc)
    def dist_int(z):
        return 1./np.sqrt(om_m*(1.+z)**3 + 1.-om_m)
    res, err = integrate.quad(dist_int, 0., z)
    return [res * (c/H0) * 1000., err * (c/H0) * 1000.]      # kpc/h

# Work in physical distance along LOS
z_int = np.linspace(0.001, 2.001, 6001)
chi_int = []
for j in range(len(z_int)):
    res, err = comoving_dist(z_int[j])
    chi_int.append(res)
chi_int = np.array(chi_int)                         # kpc/h
    
chi_func = interp1d(z_int, chi_int, kind='linear')

def sigma_crit(z_l, z_s):
    kg_per_Msol = 1.98892e30
    km_per_mpc = 3.08568e19

    G = 6.67300e-11                                      # m^3 / (kg s^2)
    
    fac = (3e8)**2 / (4.*pi*G)                           # kg/m
    fac = fac *(1e3 * km_per_mpc/1e3) / kg_per_Msol      # Msol/kpc
    
    chi_s = chi_func(z_s)                       # kpc/h
    chi_l = chi_func(z_l)                       # kpc/h
    D_s = chi_s / (1.+z_s)
    D_l = chi_l / (1.+z_l)

    D_ls = (chi_s - chi_l) / (1.+z_s)
    
    sig_c = fac * D_s / (D_l*D_ls)                      # Msol h/kpc^2
    return sig_c / 1e6                                  # Msol h/pc^2
################################
