from math import pi
from matplotlib.pyplot import *
import pyfits as pf
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from configuration import *
import randomize_profile as rp
import jc_lib as jc

import os, sys
import os.path

import time

#############################
# Cosmological Parameters
#############################

om_m, om_L, h, n_s, sig_8, om_b = jc.cosmo()

# Constants
kg_per_Msol = 1.98892e30
km_per_mpc = 3.08568e19
G = 6.67300e-11     # m^3 / (kg s^2)
H_0 = 100.          # km h / (s Mpc)

def fill_subhalo(ixsub, iysub, sigma, sigma_sub):
    fx = [-1, 0, 1, 0]
    fy = [0, 1, 0, -1]
    xnode = np.zeros(nx*ny)
    ynode = np.zeros(nx*ny)
    xnode_sub = np.zeros(nx*ny)
    ynode_sub = np.zeros(nx*ny)
    xnode[0] = ixsub
    ynode[0] = iysub
    xnode_sub[0] = nx/2
    ynode_sub[0] = ny/2
    head = -1
    tail = 0
    check = np.zeros((nx, ny))
    check[nx/2, ny/2] = 1
    while(head < tail):
        head += 1
        x = xnode[head]
        y = ynode[head]
        xsub = xnode_sub[head]
        ysub = ynode_sub[head]
        if(sigma_sub[xsub, ysub] > sigma[x, y]):
            sigma[x, y] = sigma_sub[xsub, ysub]
            for k in range(4):
                xx = x + fx[k]
                yy = y + fy[k]
                xx_s = xsub + fx[k]
                yy_s = ysub + fy[k]
                if((xx >= 0) & (xx < nx) & (yy >= 0) & (yy < ny) & (xx_s >= 0) & (xx_s < nx) & (yy_s >= 0) & (yy_s < ny)):
                    if(check[xx_s, yy_s] == 0):
                        tail += 1   
                        xnode[tail] = xx
                        ynode[tail] = yy
                        xnode_sub[tail] = xx_s
                        ynode_sub[tail] = yy_s
                        check[xx_s, yy_s] = 1
    return sigma

#@profile
def get_lens(nx, ny, xmin, xmax, ymin, ymax, ilens):
    t1 = time.time()
    sigma = np.zeros((nx, ny))
    xedge = np.linspace(xmin, xmax, nx+1)
    yedge = np.linspace(ymin, ymax, ny+1)
    xmid = (xedge[1:] + xedge[:-1]) / 2.
    ymid = (yedge[1:] + yedge[:-1]) / 2.
    ds = (xedge[1] - xedge[0]) * (yedge[1] - yedge[0])
    redge = np.sqrt(xmax**2 + ymax**2)

    #############################
    # Point mass
    #############################
    if(test_case == 'point-mass'):
        Mpoint = 1e14                  # Msol /h
        sigma[xpt, ypt] = Mpoint /ds/(1e12)    # Msol h/pc^2
        zl = 0.2

    #############################
    # SIS
    #############################
    if(test_case == 'SIS'):
        print('Creating a SIS profile...')
        zl = 0.2
        v_c = 300 # km/s
        r_c = 0.001 # Mpc/h

        nr = 2000
        r_grid = np.linspace(0., np.sqrt(xmax**2 + ymax**2), nr)
        rlosmin = 0.00001
        rlosmax = 100.
        nrlos = 3000
        rlos = np.logspace(np.log10(rlosmin), np.log10(rlosmax), nrlos)
        sigma_temp = np.zeros(nr)
        
        for k in range(nr):
            x = np.sqrt(rlos**2 + r_grid[k]**2)
            rho = v_c**2/(2*pi*G*(x**2 + r_c**2))*(km_per_mpc)*(1e9)*(kg_per_Msol)**(-1.)
            sigma_temp[k] = 2.*integrate.simps(rho, rlos)

        f = interp1d(r_grid, sigma_temp)
        
        for i in range(len(xmid)):
            for j in range(len(ymid)):
                
                r = np.sqrt(xmid[i]**2 + ymid[j]**2)
                sigma[i,j] = f(r)

                if(r > redge):sigma[i,j] = 0.

        sigma = sigma / (1e12)     # Msol h/pc^2
        plot(xmid[(nx/2):],sigma[(nx/2):,(ny/2)], label = '\Sigma(R)')
        yscale('log')
        xscale('log')
        xlabel(r'$R [Mpc/h]$', size = 'x-large')
        ylabel(r'$\Sigma(R)$', size = 'x-large')
        plotname = 'sigma_testcase%s.pdf'%(test_case)
        savefig('plots/' + plotname)
    
    #############################
    # NFW
    #############################
    if(test_case == 'NFW'):
        print('Creating a NFW profile...')
        M200 = 1e14      # Msol /h
        zl = 0.2

        rho_cr = jc.critical_density(zl)         # Msol h^2 / Mpc^3
        rho_m = rho_cr * om_m * (1.+zl)**3
        R200 = (3.*M200 / (4.*200.*pi*rho_m))**(1./3.) # Mpc/h

        print('R200 = ', R200)

        # From Diemer and Kravtsov and references therein
        c_A = 11.39
        c_B = -0.107
        c_C = -1.16
        c = c_A * (M200 / 2e12)**c_B * (1.+zl)**c_C

        rs = R200 / c
        delc = (200./3.) * c**3 / (np.log(1.+c) - c/(1.+c))

        nr = 2000
        r_grid = np.linspace(0., np.sqrt(xmax**2 + ymax**2), nr)
        rlosmin = 0.00001
        rlosmax = 100.
        nrlos = 3000
        rlos = np.logspace(np.log10(rlosmin), np.log10(rlosmax), nrlos)
        sigma_temp = np.zeros(nr)

        for k in range(nr):
            x = np.sqrt(rlos**2 + r_grid[k]**2)
            rho = (delc*rho_m)/((x/rs)*(1.+x/rs)**2)
            sigma_temp[k] = 2.*integrate.simps(rho, rlos)

        f = interp1d(r_grid, sigma_temp)

        for i in range(len(xmid)):
            for j in range(len(ymid)):
                
                r = np.sqrt(xmid[i]**2 + ymid[j]**2)
                sigma[i,j] = f(r)
                if(r > redge):sigma[i,j] = 0.
        sigma = sigma / (1e12)     # Msol h/pc^2
        plot(xmid[(nx/2):],sigma[(nx/2):,(ny/2)], label = '\Sigma(R)')
        yscale('log')
        xscale('log')
        xlabel(r'$R [Mpc/h]$', size = 'x-large')
        ylabel(r'$\Sigma(R)$', size = 'x-large')
        plotname = 'sigma_testcase%s.pdf'%(test_case)
        savefig('plots/' + plotname)

    #############################
    # Subhalo
    #############################
    if(test_case == 'subhalo'):
        M200 = rp.getM200()
        zl = rp.get_zl()

        # Calculate the R200 for the host halo (assuming NFW)
        rho_cr = jc.critical_density(zl)         # Msol h^2 / Mpc^3
        rho_m = rho_cr * om_m * (1.+zl)**3
        R200 = (3.*M200 / (4.*200.*pi*rho_m))**(1./3.) # Mpc/h
        print('R200 = ', R200)
        
        Nsub = rp.getNsub(M200)
        xsub, ysub, msub = rp.getsubhalo(Nsub, R200)

        # Calculate the projected density profile for host halo
        #############################
        # From Diemer and Kravtsov and references therein
        c_A = 11.39
        c_B = -0.107
        c_C = -1.16
        c = c_A * (M200 / 2e12)**c_B * (1.+zl)**c_C
        rs = R200 / c
        delc = (200./3.) * c**3 / (np.log(1.+c) - c/(1.+c))
        nr = 500
        r_grid = np.linspace(0., np.sqrt(xmax**2 + ymax**2), nr)
        rlosmin = 0.00001
        rlosmax = 100.
        nrlos = 1000
        rlos = np.logspace(np.log10(rlosmin), np.log10(rlosmax), nrlos)
        sigma_r = np.zeros(nr)
        for k in range(nr):
            x = np.sqrt(rlos**2 + r_grid[k]**2)
            rho = (delc*rho_m)/((x/rs)*(1.+x/rs)**2)
            sigma_r[k] = 2.*integrate.simps(rho, rlos)

        f = interp1d(r_grid, sigma_r)

        xx, yy = np.meshgrid(xmid, ymid)
        r = np.sqrt(xx**2 + yy**2)
        sigma = f(r)
        sigma[r > redge] = 0.
        #for i in range(len(xmid)):
            #for j in range(len(ymid)):
                
                #r = np.sqrt(xmid[i]**2 + ymid[j]**2)
                #sigma[i,j] = f(r)
                #if(r > redge):sigma[i,j] = 0.
        sigma = sigma / (1e12)     # Msol h/pc^2
        #############################
        # Calculate the projected density profile for subhaloes
        #############################
        for isub in range(Nsub):
            print('generating subhalo #%d...'%(isub))
            sigma_sub = np.zeros((nx, ny))
            sigma_temp = np.zeros((nx, ny))
            # TODO different mass concentration relation for subhaloes?
            c = c_A * (msub[isub] / 2e12)**c_B * (1.+zl)**c_C
            r200 = (3.* msub[isub] / (4.*200.*pi*rho_m))**(1./3.) # Mpc/h
            rs = r200 / c
            delc = (200./3.) * c**3 / (np.log(1.+c) - c/(1.+c))
            nr = 500
            r_grid = np.linspace(0., np.sqrt(xmax**2 + ymax**2), nr)
            rlosmin = 0.00001
            rlosmax = 100.
            nrlos = 1000
            rlos = np.logspace(np.log10(rlosmin), np.log10(rlosmax), nrlos)
            sigma_sub_r = np.zeros(nr)
            for k in range(nr):
                x = np.sqrt(rlos**2 + r_grid[k]**2)
                rho = (delc*rho_m)/((x/rs)*(1.+x/rs)**2)
                sigma_sub_r[k] = 2.*integrate.simps(rho, rlos)
            f = interp1d(r_grid, sigma_sub_r)

            sigma_sub = f(r)
            sigma_sub[r > redge] = 0.
            #for i in range(len(xmid)):
                #for j in range(len(ymid)):
                    #r = np.sqrt(xmid[i]**2 + ymid[j]**2)
                    #sigma_sub[i,j] = f(r)
                    #if(r > redge):sigma_sub[i,j] = 0.
            sigma_sub = sigma_sub / (1e12) # Msol h/pc^2

            ix = np.searchsorted(xmid, xsub[isub])
            iy = np.searchsorted(ymid, ysub[isub])
            sigma = fill_subhalo(ix, iy, sigma, sigma_sub)
            #for i in range(len(xmid)):
            #    for j in range(len(ymid)):
            #        xx = (i - (nx/2) + ix)
            #        yy = (j - (ny/2) + iy)
            #        if((xx >= 0) & (xx < nx) & (yy >= 0) & (yy < ny)):
            #            sigma[xx,yy] = max(sigma[xx,yy], sigma_sub[i,j])
        subhalo_file = 'members-list_lens%dof%d.fits'%(ilens, Nlens)
        key = ['x', 'y']
        form = ['D', 'D']
        output_data = np.array([xsub, ysub])
        tmpcols = []
        for i in range(len(key)):
            tmpcols.append(pf.Column(name = key[i], format = form[i], array = output_data[i]))
        # Print tmpcols
        hdu_out = pf.PrimaryHDU()
        tbhdu = pf.new_table(tmpcols)
        thdulist = pf.HDUList([hdu_out, tbhdu])
         
        # Check if fit file exist
        if(os.path.exists(subhalo_output + subhalo_file)):
            os.remove(subhalo_output + subhalo_file)
             
        thdulist.writeto(subhalo_output + subhalo_file)
        thdulist.close()
        

    #############################
    # Plotting & return profile
    #############################
    xgrid, ygrid = np.meshgrid(xmid, ymid, indexing = 'xy')
    clevel = [15, 21, 31, 45, 65, 95, 135, 200, 500]
    cs = contour(xgrid, ygrid, sigma.transpose(), levels = clevel, lw = 5)
    plotname = 'sigma2D_testcase%s.pdf'%(test_case)
    savefig('plots/' + plotname)

    t2 = time.time()
    print('Time for get_lens = ',t2-t1)
    return [sigma, zl, ds, xmid, ymid]
