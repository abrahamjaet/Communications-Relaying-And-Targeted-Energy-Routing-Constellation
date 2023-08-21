import numpy as np
import pickle 
import orbitDict    
from csltk.jpl_query import JPLFamily
from csltk import cr3bp
from csltk.utilities import System
import matplotlib.pyplot as plt


sys = System(mu=0.01215058560962404, lstar=389703.2648292776, tstar=382981.2891290545)
ax = sys.plot_system()
for i in range(33):
    i = i+1 
    filename = 'trajectories/orbit' + str(i) + '.dat'
    orb = orbitDict.chosenOrbits.load(filename)
    ax.plot(orb.x,orb.y,orb.z)
    ax.set_aspect('auto')
    print(orb.family)

plt.show()
