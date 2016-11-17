from matplotlib.pyplot import *
from matplotlib.pyplot import cm
from configuration import *
from scipy import optimize
from pylab import *
import pyfits as pf
import numpy as np

from lens import *
import jc_lib as jc

def get_major_angle(gam1, gam2):
    phi = 0.5 * (np.arctan2(gam2, gam1))
    return phi

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

    phi = get_major_angle(g1, g2)
    eabs = np.sqrt(g1**2 + g2**2)

    #U,V = eabs*np.cos(phi), eabs*np.sin(phi)
    xs = np.reshape(xs, (nx_p, ny_p))
    ys = np.reshape(ys, (nx_p, ny_p))
    U, V = np.cos(phi), np.sin(phi)
    U = np.reshape(U, (nx_p, ny_p))
    V = np.reshape(V, (nx_p, ny_p))
    eabs = np.reshape(eabs, (nx_p, ny_p))
    quiver(xs, ys, U, V, eabs, cmap = cm.seismic)
    colorbar()
    xlabel(r'$[Mpc/h]$', size = 'x-large')
    ylabel(r'$[Mpc/h]$', size = 'x-large')

    plotname = 'whisker_case%s.pdf'%(test_case)
    savefig('plots/' + plotname)
    

    
