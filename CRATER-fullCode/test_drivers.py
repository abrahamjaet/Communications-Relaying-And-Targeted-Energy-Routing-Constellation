import math
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from csltk.jpl_query import JPLFamily
from csltk import cr3bp
from csltk.utilities import System
import pickle
import itertools
from collections import Counter
from os import listdir
from os.path import isfile, join
import multiprocessing
import time
import sys

import orbitDict
chosenOrbits = orbitDict.chosenOrbits
from design import designConfiguration, orbitConfiguration, powerConfiguration, commsConfiguration
from ROI_Model.ROI_Driver import roi
roi = roi()
from Communications_Model.Comms_Driver import comms
comms = comms()
from Power_Model.Power_Driver import power
power = power()


orbitFiles = [f for f in listdir('trajectories') if isfile(join('trajectories', f))]
orbit = []
orbitNumber = []
numOrbits = 0
for i in range(len(orbitFiles)):
    filename = 'trajectories/' + orbitFiles[i]
    orbit.append(chosenOrbits.load(filename))
    orbitNumber.append(numOrbits)
    numOrbits += 1
num_to_plot = 2 # number of orbits in family to plot
orbit = orbit[::int(len(orbit)/num_to_plot)]
# orbitNumber = orbitNumber[::int(len(orbitNumber)/num_to_plot)]
orbitNumber = range(len(orbit))

print(range(len(orbit)))

# Intermediate Constants
Receiver_Reflectivity_Const = 0.1 # equal to 1 pane of glass
Receiver_Intensity_Cutoff = 1380*450 # [W/m^2]
P_per_kg = 1500 # [W/kg] (for Lithium Ion Battery)
E_per_kg = 200 # [Wh/kg] (for Lithium Ion Battery)

   
# Design Variables
diameterTxM = 1 # [m] antenna diameter of the receiver on the moon
batterySize_Range = np.linspace(200,5000,3)
solarPanelSize_Range = np.linspace(50, 1000, 3) # [m^2] solar panel area
laserPower_Range = [] # 1000  to  P_per_kg*LI_battery_mass [W]
apetureRad_range = np.linspace(0.01, 0.05, 3) # [m] radius of output lens on SC
powerReceiverRad_Range = np.linspace(1, 20, 3) # [m] radius of ground power receiver
diameterTxO_Range = np.linspace(1.3, 5, 3)
dataRate_Range = np.linspace(3e6, 40e6, 3) #list(np.arange(3e6, 40e6, 2.5e5)) # [bps] 
dataRate_ED_Range = dataRate_Range # [bps] desired data rate for earth downlink
numSat_max = 3
ID = 0


currOrbits = orbit
satDistr = [1,2]
batterySize = batterySize_Range[1]
solarPanelSize = solarPanelSize_Range[1]
laserPower = 1000
apetRad = apetureRad_range[1]
powerReceiverRad = powerReceiverRad_Range[1]
diameterTxO = diameterTxO_Range[1]
dataRate = dataRate_Range[1]
dataRate_ED = dataRate_ED_Range[1]


currDesign = designConfiguration(ID, currOrbits, satDistr, numSat_max, batterySize, solarPanelSize, laserPower, apetRad, powerReceiverRad, diameterTxM, diameterTxO, dataRate, dataRate_ED)

# %% Power
powerOutput = power.driver(currDesign)
print(powerOutput)

# %% Comms
commsOutput, commsConstraint = comms.driver(currDesign)
print("comms score done")
print(commsOutput)

# %% ROI
roiConstraint = roi.driver(currDesign, powerOutput, commsOutput)
print("roi score done")

print(currDesign.powerObj)
print(currDesign.commsObj)
print(currDesign.roiObj)