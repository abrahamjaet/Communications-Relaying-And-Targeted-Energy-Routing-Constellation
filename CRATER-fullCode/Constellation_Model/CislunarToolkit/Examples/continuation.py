"""
Script to demonstrate continuation for findng orbit families. Neighboring orbits 
are found starting from an initial periodic orbit.
"""


from csltk import cr3bp
from csltk.utilities import Earth, Moon, System
from csltk.correctors import Continuation

import numpy as np
import matplotlib.pyplot as plt

# define Earth-Moon system and orbit parameters
sys = System(P1=Earth, P2=Moon)
ax = sys.plot_system(axlim=[-1.5, 1.5])

# initial conditions (x0, P)
x0c = np.array([4.8784941344943100E-1,	-6.7212754533055685E-1,	0.0000000000000000E+0,	2.4535031351852712E-1,	1.3805858489246356E-1,	0.0000000000000000E+0, 2.4071477560861904E+1])

# continuation
orbit_family = Continuation.psuedo_arclength(x0c, sys, 20, delta_s=1e-1)

# plot the system
for orbit in orbit_family:

    # propagate orbit for plotting
    x0,P = np.split(orbit, [-1])
    t_span = np.array([0, P]) # time span
    sol = cr3bp.propagate(sys, x0, t_span)

    # plot the orbit
    ax.plot(sol.y[0,:], sol.y[1,:], sol.y[2,:])

plt.show()
