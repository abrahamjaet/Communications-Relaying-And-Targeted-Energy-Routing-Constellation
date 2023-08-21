"""
Script to demonstrate API for requesting periodic 3BP orbits from JPL database:
https://ssd.jpl.nasa.gov/tools/periodic_orbits.html#/periodic
"""


from csltk.jpl_query import JPLFamily
from csltk import cr3bp
from csltk.utilities import System

import numpy as np
import matplotlib.pyplot as plt

# define Earth-Moon system and orbit parameters
sys = System(mu=0.01215058560962404, lstar=389703.2648292776, tstar=382981.2891290545)
ax = sys.plot_system()

# supported values: lpo, dro, dragonfly, axial, halo, short, lyapunov, vertical, longp, dpo, resonant, butterfly
# pull JPL orbit initial conditions
fam = JPLFamily(sys='earth-moon', fam='axial', libr='2')
# fam = JPLFamily(sys='earth-moon', fam='resonant', libr='', branch = 'pq')
# fam.set_period_bounds(2, 2.3)
# fam.set_jacobi_bounds(2.75, 3.1)
# fam.set_stab_bounds(1, 30)
orbits = fam.request()
print('Number of orbits:', len(orbits))

# plot the orbits
num_to_plot = 1 # number of orbits in family to plot
for orbit in orbits[::int(len(orbits)/num_to_plot)]:

    # propagate orbit for plotting
    x0 = orbit.iState # initial conditions
    t_span = np.array([0, 6*orbit.T]) # time span
    sol = cr3bp.propagate(orbit.system, x0, t_span)

    # plot the orbit
    ax.plot(sol.y[0,:], sol.y[1,:], sol.y[2,:])
    ax.set_aspect('auto')

plt.show()
