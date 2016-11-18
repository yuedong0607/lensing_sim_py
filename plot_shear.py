from matplotlib.pyplot import *
import pyfits as pf
import numpy as np

from configuration import *

import jc_lib as jc

R_edge = np.logspace(np.log10(Rmin), np.log10(Rmax), nR+1)
R_mid = np.sqrt(R_edge[1:]*R_edge[:-1])

print R_edge
print R_mid

gt_arr = np.zeros(nR)
gx_arr = np.zeros(nR)
ct_arr = np.zeros(nR)
dsig_arr = np.zeros(nR)

for ilens in range(Nlens):
    src_file = 'simulated_shear_case%s_lens%dof%d.fit' % (test_case, ilens+1, Nlens)
    hdu = pf.open(shear_output + src_file)
    columns = hdu[1].columns
    print columns

    xs = hdu[1].data.field('xs')
    ys = hdu[1].data.field('ys')
    zs = hdu[1].data.field('zs')
    g1 = hdu[1].data.field('g1')
    g2 = hdu[1].data.field('g2')

    hdu.close()
    # Read in lens redshift (TODO)
    zl = 0.2

    if(test_case == 'subhalo'):
        subhalo_file = 'members-list_lens%dof%d.fits'%(ilens, Nlens)
        hdu = pf.open(subhalo_output + subhalo_file)
        columns = hdu[1].columns
        print columns
        xsub = hdu[1].data.field('x')[0]
        ysub = hdu[1].data.field('y')[0]
        print 'xsub = ',xsub
        print 'ysub = ',ysub
        hdu.close()
        xs = xs - xsub
        ys = ys - ysub
        r_c = np.sqrt(xsub**2 + ysub**2)
        axvline(r_c)

    R_s= np.sqrt(xs**2 + ys**2)  # Mpc/h
    sigma_cr = jc.sigma_crit(zl, zs)

    # Calculate tangential shear
    psi = np.arctan2(ys, xs)

    gt = (-g1 * np.cos(2.*psi) - g2 * np.sin(2.*psi))
    gx = (g1 * np.sin(2.*psi) - g2 * np.cos(2.*psi))

    dsig_t = gt * sigma_cr

    for i in range(nR):
        ind = ((R_edge[i] <= R_s) & (R_s < R_edge[i+1]))
        gt_arr[i] = gt_arr[i] + np.sum(gt[ind])
        gx_arr[i] = gx_arr[i] + np.sum(gx[ind])
        ct_arr[i] = ct_arr[i] + np.sum(ind)
        dsig_arr[i] = dsig_arr[i] + np.sum(dsig_t[ind])

# For testing, without error estimation
gt_mean = np.zeros(nR)
gx_mean = np.zeros(nR)
dsig_mean = np.zeros(nR)

ind = (ct_arr != 0)
gt_mean[ind] = gt_arr[ind]/ct_arr[ind]
gx_mean[ind] = gx_arr[ind]/ct_arr[ind]
dsig_mean[ind] = dsig_arr[ind]/ct_arr[ind]

#plot(R_mid[ind], gt_mean[ind], 'ro', label = r'$\gamma_t$')
#plot(R_mid[ind], gx_mean[ind], 'rx', label = r'$\gamma_x$')
plot(R_mid[ind], dsig_mean[ind], 'ro', label = r'$\Delta\Sigma [M_\odot h / {\rm pc}^2]$')
legend(loc = 'upper right')
axhline(0, c = 'k', ls = ':')
xlabel(r'$R [Mpc/h]$', size = 'x-large')
xscale('log')
yscale('log')
xlim([0.01, 100.])
ylim([1., 1e3])
title(test_case, size = 'x-large')

# Save plot
plotname = 'observation_simulated_%s_case%s.pdf' % (plot_option, test_case)
savefig('plots/' + plotname)
