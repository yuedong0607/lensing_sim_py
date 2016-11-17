#test_case = 'point-mass'
#test_case = 'SIS'
#test_case = 'NFW'
test_case = 'subhalo'

# point mass paprameters
Mpoint = 1e15
xpt = 100
ypt = 0
################################
Nlens = 1
nx = 1024
ny = 1024
xmin = -2. # Mpc/h
xmax = 2.  # Mpc/h
ymin = -2. # Mpc/h
ymax = 2.  # Mpc/h

shear_output = './shear_cat/'
subhalo_output = './members_cat/'

################################
# Ploting/measurement parameters
################################

#plot_option = 'gammat'
plot_option = 'dsig'
nR = 18
Rmin = 0.05
Rmax = 20.

# For whisker plot
xmin_p = -1. # Mpc/h
xmax_p = 1.  # Mpc/h
ymin_p = -1. # Mpc/h
ymax_p = 1.  # Mpc/h
nx_p = 25
ny_p = 25

