from typing import Tuple, Any

import matplotlib.pyplot as plt
import numpy as np
import time

from csltk.jpl_query import JPLFamily
from csltk import cr3bp
from csltk.utilities import System

import orbitDict

start = time.time()

LU = 389703  # 1 LU (distance from Earth to Moon in km)
moon_rad = 1740 / LU
step = 10  # 1/10 every step (100 seconds a step)
synodicPeriod = 29.523  # days
synodicS = synodicPeriod * 24 * 60 * 60  # s
loops = int(synodicS / (100 * step))

# filename1 = 'Trajectories/orbit1.dat'
# filename2 = 'Trajectories/orbit2.dat'

count = 0

distance = []
max_distance = []
for m in range (37):
    for l in range (37):

        count += 1

        print('------------ RUN', count, '------------ ')
        filename1 = 'Trajectories/orbit' + str(m) + '.dat'
        filename2 = 'Trajectories/orbit' + str(l) + '.dat'
        orb1 = orbitDict.chosenOrbits.load(filename1)
        orb2 = orbitDict.chosenOrbits.load(filename2)

        # - EVERY POSITION 2 SATELLITES - #
#        print(' # - 2 Satellites, every position compared to every position - #')
#        for i in range(int(np.floor(np.divide((len(orb1.x)), step)))):
#            for j in range(int(np.floor(np.divide((len(orb2.x)), step)))):
#                posX1 = orb1.x[step * i]
#                posY1 = orb1.y[step * i]
#                posZ1 = orb1.z[step * i]

#                posX2 = orb2.x[step * j]
#                posY2 = orb2.y[step * j]
#                posZ2 = orb2.z[step * j]

#                sat_to_satX = posX1 - posX2
#                sat_to_satY = posY1 - posY2
#                sat_to_satZ = posZ1 - posZ2

#                r = (sat_to_satX, sat_to_satY, sat_to_satZ)

#                distance.append(np.linalg.norm(r))

#        max_distance = max(distance)
#        print('Max distance POSSIBLE between satellites:', max_distance * LU)

        # - PHASING ATTEMPT - #

        print(' # - Phasing, 2 satellites - #')

        num_sat_orb1 = 1  # MAX OF 3
        num_sat_orb2 = 1  # MAX OF 3

        # 1 on first, 1 on second: max distance = 30189.32664616554
        # 2 on first, 1 on second: max distance = 30566.803838601976
        # 3 on first, 1 on second: max distance = 30566.803838601976
        # 1 on first, 2 on second: max distance = 30189.32664616554
        # 2 on first, 2 on second: max distance = 30566.803838601976
        # 3 on first, 2 on second: max distance = 30566.803838601976
        # 1 on first, 3 on second: max distance = 30189.32664616554
        # 2 on first, 3 on second: max distance = 30566.803838601976
        # 3 on first, 3 on second: max distance = 30566.803838601976
        # ^ avg of 30,441 km

        # keeping track of all combinations

        distance_this_run = []

        i = -1
        i1 = -1
        i2 = -1
        j = -1
        j1 = -1
        j2 = -1

        start_i = int(np.floor(np.divide((len(orb1.x)), num_sat_orb1)))
        i1 = start_i * (1 / num_sat_orb1)
        i2 = start_i * (2 / num_sat_orb1)
        start_j = int(np.floor(np.divide((len(orb2.x)), num_sat_orb2)))
        j1 = start_j * (1 / num_sat_orb2)
        j2 = start_j * (2 / num_sat_orb2)

        for k in range(loops):

            # for i in range(int(np.floor(np.divide((len(orb1.x)), step)))):
            i += 1  # remove this if using i loop
            i1 += 1
            i2 += 1

            if (step * i) > (len(orb1.x)):
                i = 0
            if (step * i1) > (len(orb1.x)):
                i1 = 0
            if (step * i2) > (len(orb1.x)):
                i2 = 0

            # for j in range(int(np.floor(np.divide((len(orb2.x)), step)))):
            j += 1
            j1 += 1
            j2 += 1

            if (step * j) > (len(orb2.x)):
                j = 0
            if (step * j1) > (len(orb2.x)):
                j1 = 0
            if (step * j2) > (len(orb2.x)):
                j2 = 0

            posX1 = orb1.x[(step * i)-1]
            posY1 = orb1.y[(step * i)-1]
            posZ1 = orb1.z[(step * i)-1]

            posX2 = orb2.x[(step * j)-1]
            posY2 = orb2.y[(step * j)-1]
            posZ2 = orb2.z[(step * j)-1]

            if num_sat_orb1 > 1:
                posX1_1 = orb1.x[step * i1]
                posY1_1 = orb1.y[step * i1]
                posZ1_1 = orb1.z[step * i1]

            if num_sat_orb2 > 1:
                posX2_1 = orb2.x[step * j1]
                posY2_1 = orb2.y[step * j1]
                posZ2_1 = orb2.z[step * j1]

            if num_sat_orb1 > 2:
                posX1_2 = orb1.x[step * i2]
                posY1_2 = orb1.y[step * i2]
                posZ1_2 = orb1.z[step * i2]

            if num_sat_orb2 > 2:
                posX2_2 = orb2.x[step * j2]
                posY2_2 = orb2.y[step * j2]
                posZ2_2 = orb2.z[step * j2]

            # r1
            sat_to_satX_1 = posX1 - posX2
            sat_to_satY_1 = posY1 - posY2
            sat_to_satZ_1 = posZ1 - posZ2
            r1 = (sat_to_satX_1, sat_to_satY_1, sat_to_satZ_1)
            distance.append(np.linalg.norm(r1))
            distance_this_run.append(np.linalg.norm(r1))

            # r2 & r3
            if num_sat_orb1 > 1:
                sat_to_satX_2 = posX1_1 - posX2
                sat_to_satY_2 = posY1_1 - posY2
                sat_to_satZ_2 = posZ1_1 - posZ2
                r2 = (sat_to_satX_2, sat_to_satY_2, sat_to_satZ_2)
                distance.append(np.linalg.norm(r2))
                distance_this_run.append(np.linalg.norm(r2))

                sat_to_satX_3 = posX1_1 - posX1
                sat_to_satY_3 = posY1_1 - posY1
                sat_to_satZ_3 = posZ1_1 - posZ1
                r3 = (sat_to_satX_3, sat_to_satY_3, sat_to_satZ_3)
                distance.append(np.linalg.norm(r3))
                distance_this_run.append(np.linalg.norm(r3))

            # r4 & r5
            if num_sat_orb2 > 1:
                sat_to_satX_4 = posX1 - posX2_1
                sat_to_satY_4 = posY1 - posY2_1
                sat_to_satZ_4 = posZ1 - posZ2_1
                r4 = (sat_to_satX_4, sat_to_satY_4, sat_to_satZ_4)
                distance.append(np.linalg.norm(r4))
                distance_this_run.append(np.linalg.norm(r4))

                sat_to_satX_5 = posX2 - posX2_1
                sat_to_satY_5 = posY2 - posY2_1
                sat_to_satZ_5 = posZ2 - posZ2_1
                r5 = (sat_to_satX_5, sat_to_satY_5, sat_to_satZ_5)
                distance.append(np.linalg.norm(r5))
                distance_this_run.append(np.linalg.norm(r5))

            # r6, r7 & r8
            if num_sat_orb1 > 2:
                sat_to_satX_6 = posX2 - posX1_2
                sat_to_satY_6 = posY2 - posY1_2
                sat_to_satZ_6 = posZ2 - posZ1_2
                r6 = (sat_to_satX_6, sat_to_satY_6, sat_to_satZ_6)
                distance.append(np.linalg.norm(r6))
                distance_this_run.append(np.linalg.norm(r6))

                sat_to_satX_7 = posX1 - posX1_2
                sat_to_satY_7 = posY1 - posY1_2
                sat_to_satZ_7 = posZ1 - posZ1_2
                r7 = (sat_to_satX_7, sat_to_satY_7, sat_to_satZ_7)
                distance.append(np.linalg.norm(r7))
                distance_this_run.append(np.linalg.norm(r7))

                sat_to_satX_8 = posX1_1 - posX1_2
                sat_to_satY_8 = posY1_1 - posY1_2
                sat_to_satZ_8 = posZ1_1 - posZ1_2
                r8 = (sat_to_satX_8, sat_to_satY_8, sat_to_satZ_8)
                distance.append(np.linalg.norm(r8))
                distance_this_run.append(np.linalg.norm(r8))

            # r9, r10 & r11
            if num_sat_orb2 > 2:
                sat_to_satX_9 = posX2 - posX2_2
                sat_to_satY_9 = posY2 - posY2_2
                sat_to_satZ_9 = posZ2 - posZ2_2
                r9 = (sat_to_satX_9, sat_to_satY_9, sat_to_satZ_9)
                distance.append(np.linalg.norm(r9))
                distance_this_run.append(np.linalg.norm(r9))

                sat_to_satX_10 = posX1 - posX2_2
                sat_to_satY_10 = posY1 - posY2_2
                sat_to_satZ_10 = posZ1 - posZ2_2
                r10 = (sat_to_satX_10, sat_to_satY_10, sat_to_satZ_10)
                distance.append(np.linalg.norm(r10))
                distance_this_run.append(np.linalg.norm(r10))

                sat_to_satX_11 = posX2_1 - posX2_2
                sat_to_satY_11 = posY2_1 - posY2_2
                sat_to_satZ_11 = posZ2_1 - posZ2_2
                r11 = (sat_to_satX_11, sat_to_satY_11, sat_to_satZ_11)
                distance.append(np.linalg.norm(r11))
                distance_this_run.append(np.linalg.norm(r11))

            # r12
            if num_sat_orb1 > 1 and num_sat_orb2 > 1:
                sat_to_satX_12 = posX2_1 - posX1_1
                sat_to_satY_12 = posY2_1 - posY1_1
                sat_to_satZ_12 = posZ2_1 - posZ1_1
                r12 = (sat_to_satX_12, sat_to_satY_12, sat_to_satZ_12)
                distance.append(np.linalg.norm(r12))
                distance_this_run.append(np.linalg.norm(r12))

            # r13
            if num_sat_orb1 > 1 and num_sat_orb2 > 2:
                sat_to_satX_13 = posX1_1 - posX2_2
                sat_to_satY_13 = posY1_1 - posY2_2
                sat_to_satZ_13 = posZ1_1 - posZ2_2
                r13 = (sat_to_satX_13, sat_to_satY_13, sat_to_satZ_13)
                distance.append(np.linalg.norm(r13))
                distance_this_run.append(np.linalg.norm(r13))

            # r14
            if num_sat_orb1 > 2 and num_sat_orb2 > 1:
                sat_to_satX_14 = posX1_2 - posX2_1
                sat_to_satY_14 = posY1_2 - posY2_1
                sat_to_satZ_14 = posZ1_2 - posZ2_1
                r14 = (sat_to_satX_14, sat_to_satY_14, sat_to_satZ_14)
                distance.append(np.linalg.norm(r14))
                distance_this_run.append(np.linalg.norm(r14))

            # r15
            if num_sat_orb1 > 2 and num_sat_orb2 > 2:
                sat_to_satX_15 = posX2_2 - posX2_2
                sat_to_satY_15 = posY2_2 - posY2_2
                sat_to_satZ_15 = posZ2_2 - posZ2_2
                r15 = (sat_to_satX_15, sat_to_satY_15, sat_to_satZ_15)
                distance.append(np.linalg.norm(r15))
                distance_this_run.append(np.linalg.norm(r15))

        if m != l:
            this_run_max_distance = max(distance_this_run)
            print('Run:', count, ' Max distance between satellites:', this_run_max_distance * LU, ' (for this run)')
            max_distance.append((max(distance)))
        else:
            print('m = l, comparing same orbit so distances were just 0')

max_max_distance = max(max_distance)
print('Max distance for all runs between satellites:', max_max_distance * LU)
avg_max_distance = np.mean(max_distance)
print('Average Max distance between satellites:', avg_max_distance * LU)

end = time.time()
print('run time:',(end-start)/60,'min')