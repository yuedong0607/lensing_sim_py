from matplotlib.pyplot import *
from configuration import *
from scipy import optimize
from pylab import *
import pyfits as pf
import numpy as np

from lens import *
import jc_lib as jc

R_edge = np.logspace(np.log10(Rmin), np.log10(Rmax), nR+1)
R_mid = np.sqrt(R_edge[1:]*R_edge[:-1])

gt_arr = np.zeros(nR)
gx_arr = np.zeros(nR)
ct_arr = np.zeros(nR)

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

    if(test_case == 'point-mass'):
        sigma, zl, ds, xl, yl = get_lens(nx, ny, xmin, xmax, ymin, ymax)
        xs = xs - xl[xpt]
        ys = ys - yl[ypt]
    if(test_case == 'subhalo'):
        subhalo_file = 'members-list_lens%dof%d.fits'%(ilens, Nlens)
        hdu = pf.open(subhalo_output + subhalo_file)
        columns = hdu[1].columns
        print columns
        xsub = hdu[1].data.field('x')[0]
        ysub = hdu[1].data.field('y')[0]
        print xsub
        print ysub
        hdu.close()
        xs = xs - xsub
        ys = ys - ysub
        r_c = np.sqrt(xsub**2 + ysub**2)

    R_s= np.sqrt(xs**2 + ys**2)  # Mpc/h
    sigma_cr = jc.sigma_crit(zl, zs)

    # Calculate tangential shear
    psi = np.arctan2(ys, xs)

    gt = (-g1 * np.cos(2.*psi) - g2 * np.sin(2.*psi))
    gx = (g1 * np.sin(2.*psi) - g2 * np.cos(2.*psi))

    dsig_t = gt * sigma_cr

if(plot_option == 'gammat'):
    plot(R_s, gt, 'ro', label = r'$\gamma_t$')
    plot(R_s, gx, 'rx', label = r'$\gamma_x$')
elif(plot_option == 'dsig'):
    plot(R_s, dsig_t, 'ro', label = r'$\Delta\Sigma [M_\odot h / {\rm pc}^2]$')
    axvline(r_c)
    print r_c

# Test: fit a power law profile to delta_sigma measurement
ind = (dsig_t > 0)
logx = log10(R_s[ind])
logy = log10(dsig_t[ind])
fitfunc = lambda p, x: p[0] + (x*p[1])
errfunc = lambda p, x, y: (y - fitfunc(p, x))

pinit = [1.0, -1.0]
out = optimize.leastsq(errfunc, pinit, args=(logx, logy), full_output=1)

pfinal = out[0]

powerlaw = lambda x, amp, index: amp * (x**index)

#plot(R_s[ind], powerlaw(R_s[ind], 10.**pfinal[0], pfinal[1]), label = r'$y = %f * x^{%f}$'%(10.**pfinal[0], pfinal[1]))

# Ploting adjustment
legend(loc = 'upper right')
axhline(0, c = 'k', ls = ':')
xlabel(r'$R [Mpc/h]$', size = 'x-large')
xscale('log')
yscale('log')
xlim([0.05, 5.])
ylim([0.1,1000])
#ylim([0.1, max(dsig_t)])
title(test_case, size = 'x-large')

# Save plot
plotname = 'simulated_%s_case%s.pdf' % (plot_option, test_case)
savefig('plots/' + plotname)
