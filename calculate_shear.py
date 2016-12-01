import numpy as np
import pyfits as pf
from math import pi

from configuration import *
from lens import *
from source import *

import jc_lib as jc

import os, sys
import os.path

import cProfile, pstats
pr = cProfile.Profile()
pr.enable()

# Keys for shear fits file
key = ['xs', 'ys', 'g1', 'g2', 'zs']
form = ['D', 'D', 'D', 'D', 'D']

for ilens in range(Nlens):

    sigma, zl, ds, xl, yl = get_lens(nx, ny, xmin, xmax, ymin, ymax, ilens)
    #xs, ys, zs = get_source()
    # For testing subhalo
    if(test_case == 'subhalo'):
        subhalo_file = 'members-list_lens%dof%d.fits'%(ilens, Nlens)
        hdu = pf.open(subhalo_output + subhalo_file)
        # Test the first member
        xsub = hdu[1].data.field('x')[0]
        ysub = hdu[1].data.field('y')[0]
        print xsub, ' ', ysub
        hdu.close()
        # Test the sources around subhaloes
        #xs, ys, zs = get_source_around(xsub, ysub)

        # Get the sources grid for plotting a whisker map
        xs, ys, zs = get_source_grid(xmin_p, xmax_p, ymin_p, ymax_p, nx_p, ny_p)

    sigma_cr = jc.sigma_crit(zl, zs)
    #print 'sigma_cr = ',sigma_cr
    
    gamma = np.zeros(len(zs),dtype = complex)

    eps = (xmax - xmin)*np.sqrt(2)/nx
    print 'eps = ',eps
    
    xx, yy = np.meshgrid(xl, yl)
    for k in range(len(xs)):
        relx = xs[k] - xx
        rely = ys[k] - yy
        r = relx + 1j * rely

        tmp = 1./r
        tmp *= -ds
        tmp *= sigma
        tmp *= (1./sigma_cr[k])

        rsq = r.real**2
        rsq += r.imag**2
        ind = rsq > eps**2
        gamma[k] += np.sum(tmp[ind])
        #print 'tmp = ',tmp
        #tmp2 = np.sum((ds*sigma[ind]/sigma_cr[k]) * (-1) / (((xs[k] - xx[ind]) - 1.j*(ys[k] - yy[ind]))*((xs[k] - xx[ind]) - 1.j*(ys[k] - yy[ind]))))
        #print 'tmp2 = ',tmp2
        #print np.abs(np.sum(tmp-tmp2))
        #gamma[k] += np.sum((ds*sigma[ind]/sigma_cr[k]) * (-1) / (((xs[k] - xx[ind]) - 1.j*(ys[k] - yy[ind]))*((xs[k] - xx[ind]) - 1.j*(ys[k] - yy[ind]))))

    #for i in range(nx):
        #for j in range(ny):
            #r = np.sqrt((xs - xl[i]) ** 2 + (ys - yl[j])**2)
            #ind = (r > eps)
            #gamma = gamma + (ds*sigma[i,j]/sigma_cr) * (-1) / (((xs - xl[i]) - 1.j*(ys - yl[j]))*((xs - xl[i]) - 1.j*(ys - yl[j])))
            #gamma[ind] = gamma[ind] + (ds*sigma[i,j]/sigma_cr[ind]) * (-1) / (((xs[ind] - xl[i]) - 1.j*(ys[ind] - yl[j]))*((xs[ind] - xl[i]) - 1.j*(ys[ind] - yl[j])))
    print 'after loop over sources'

    gamma = gamma / pi
    g1 = gamma.real
    g2 = gamma.imag
    print 'g1,g2 = ',g1,g2

    # Write shear fits file
    src_fout = 'simulated_shear_case%s_lens%dof%d.fit' % (test_case, ilens+1, Nlens)
    output_data = np.array([xs, ys, g1, g2, zs])

    tmpcols = []
    for i in range(len(key)):
        tmpcols.append(pf.Column(name = key[i], format = form[i], array = output_data[i]))

    # Print tmpcols
    hdu_out = pf.PrimaryHDU()
    tbhdu = pf.new_table(tmpcols)
    thdulist = pf.HDUList([hdu_out, tbhdu])

    # Check if fit file exist
    if(os.path.exists(shear_output + src_fout)):
        os.remove(shear_output + src_fout)
    
    thdulist.writeto(shear_output + src_fout)
    thdulist.close()
    print 'wrote to ',src_fout

pr.disable()
ps = pstats.Stats(pr).sort_stats('time').reverse_order()
ps.print_stats()

