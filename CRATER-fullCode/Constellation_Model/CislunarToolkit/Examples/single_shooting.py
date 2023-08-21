"""
Script to demonstrate single shooting corrector. Attempts to find a periodic orbit given 
an initial state and period guess.

SingleShooter.variable_time: period is free to vary
SingleShooter.fixed_time: period is fixed
"""

from csltk import cr3bp
from csltk.utilities import Earth, Moon, System
from csltk.correctors import SingleShooter

import numpy as np
import matplotlib.pyplot as plt

# define Earth-Moon system and orbit parameters
sys = System(P1=Earth, P2=Moon)
ax = sys.plot_system()
num_periods = 50 # number of periods for plotting

# initial state & period guess (x0, P)
ode_kw = dict(rel_tol=1E-13, abs_tol=1E-13) # integration parameters
x0_guess = np.array([1.07, 0.0, 0.2, 0.0, -0.19, 0.0, 2.29])
sol0 = cr3bp.propagate(sys, x0_guess[:6], np.array([0, num_periods*x0_guess[6]]), **ode_kw)

# correct to periodic orbit
x1,P1,_ = SingleShooter.fixed_time(x0_guess, sys, tol=1e-10) # period is fixed
# x1,P1,_ = SingleShooter.variable_time(x0_guess, sys, tol=1e-10) # period is free to vary
print('Corrected initial state: ', x1)
print('Corrected period: ', P1)
sol1 = cr3bp.propagate(sys, x1, np.array([0, num_periods*P1]), **ode_kw)

ax.plot(sol0.y[0,:], sol0.y[1,:], sol0.y[2,:], '--', label='Initial guess')
ax.plot(sol1.y[0,:], sol1.y[1,:], sol1.y[2,:], label='Corrected')
ax.legend()

plt.show()
