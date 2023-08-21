from re import T
import numpy as np
from cmath import log10
from math import pi
import matplotlib.pyplot as plt
import time
from csltk.utilities import System


# Move path to main CRATER directory
import sys
import os
# getting the name of the directory where the this file is present.
current = os.path.dirname(os.path.realpath(__file__))
# Getting the parent directory name where the current directory is present.
parent = os.path.dirname(current) 
# adding the parent directory to the sys.path.
sys.path.append(parent)
 
from orbitDict import chosenOrbits
from design import designConfiguration, orbitConfiguration

class comms:
    def SpaceLoss(self, distance,wavelength):
        #Inputs: Distance in METERS & Wavelength in Hz
        #Outputs: Spaceloss in DB
        #Method: Using FSPL formula

        loss = -1*10*log10(((4*pi*distance)/wavelength)**2)
        return(loss.real) #return only the real part

    def AntennaParameters(self, D, wavelength):
        #Inputs: Diameter (meters) Wavelength (Hz)
        #Outputs: Gain in dB, area in m^2 (trivial calculation but sped up with function implementation)
        #Method: use a form to find gain that can be sampled via
        #diameter and wavelength

        k = 1.38e-23 #boltzman's constant
        area = pi*(D/2)**2
        beamwidth = 70 * (wavelength/D) 
        area_eff = area*.55 # = effective area simplified formula
        gain = 10*log10((6.9*area_eff)/(wavelength**2))
        return [gain.real, area, area_eff, beamwidth]

    def FOM(self, orbitDesign):

        def reorder(lst, newIdx):
            pos = newIdx #lst.index(first)
            return lst[pos:] + lst[:pos]
        
        orbits = orbitDesign.orbits
        numSats = orbitDesign.numSats
        # print('num orbits', len(orbits), 'num sats')
        # for sats in numSats: 
        #     print(sats)

        ## inputs are a design 
        ## takes in list of orbits with corresponding number of sats per orbit 
        
        # ax = plt.axes(projection = '3d')
        # sys = System(mu=0.01215058560962404, lstar=389703.2648292776, tstar=382981.2891290545)
        # ax = sys.plot_system()
        # ax = plt.axes(projection ="3d") 
        ###############################################
        ################# Constants ###################
        ###############################################

        # TU = 382981.2891290545 # seconds/TU
        # LU = 389703 # 1 LU (distance from Earth to Moon in km)
        # moon_rad = 1740/LU # radius of the Moon in [km]
        moonPos = np.array([0.98784942,0,0]) # LU 
        # lat_spacing = 10 # [deg]
        # max_points = 12 
        # latitude = np.linspace(90,-90,19)
        # gap = 360/max_points 
        synodicPeriod = 29.523 # days 
        synodicS = synodicPeriod*24*60*60 # s 
        # sToHr = (1/(60*60))

        ################################################
        ############## Grid point Calcluation ##########
        ################################################

        grid_point_coordinates = np.genfromtxt('Communications_Model/FOMLoad.csv', delimiter=',') # x, y, z, r,theta,phi,lat,long
        grid_point_coordinates[:,0] = grid_point_coordinates[:,0] + 0.98784942 
        grid_point_coordinates = grid_point_coordinates[:,0:3]
        totalPoints = len(grid_point_coordinates[:,0])

        # plot grid points 
        # ax = plt.axes(projection ="3d") 
        # ax.scatter3D(grid_point_coordinates[:,0],grid_point_coordinates[:,1],grid_point_coordinates[:,2])
        # plt.show()

        #################################################
        ############# Orbit Simulation #################
        #################################################


        # ax.set_xlabel('x')
        # ax.set_ylabel('y')
        # ax.set_zlabel('z')
       

        
        loops = round(synodicS/100)
        rows = totalPoints
        cols = loops
        
        coverageMain = np.zeros((rows,cols)) # each row corresponds to a grid point, each col corresponds to a time 
        counter = 0 
        # go through each orbit


        a = 1
        # sys = System(mu=0.01215058560962404, lstar=389703.2648292776, tstar=382981.2891290545)
        # ax = sys.plot_system()
        # for orb in orbits: 
        #     ax.plot(orb.x,orb.y,orb.z)
        #     ax.set_aspect('auto') 
        # # plt.show

        for orb in orbits:
           
            
            coverage = np.zeros((rows,cols)) # rows = points, cols = time 
            timeCounterCoverage = 0 # restart time counter goes with loops 
            counterLoops = 0 # restarting orbits counter 
        
            if not numSats[counter] == 0:
                for i in range(loops): 
                    
                    # allocate sat position 
                    if i >= (len(orb.x)*counterLoops + len(orb.x)):
                        counterLoops = counterLoops + 1
                    i = i - len(orb.x)*counterLoops
                    sat_xPos = orb.x[i]
                    sat_yPos = orb.y[i]
                    sat_zPos = orb.z[i]
            
                    ############# new stuff ##################

                    r_spacecraft = [sat_xPos, sat_yPos, sat_zPos] # vector from center of earth to satellite
                    r_spacecraftToPoint = np.zeros((len(grid_point_coordinates[:,0]), 3))
                    r_spacecraftToPoint[:,0] = grid_point_coordinates[:,0] - sat_xPos
                    r_spacecraftToPoint[:,1] = grid_point_coordinates[:,1] - sat_yPos
                    r_spacecraftToPoint[:,2] = grid_point_coordinates[:,2] - sat_zPos
                    r_moonToPoint = np.zeros((len(grid_point_coordinates[:,0]), 3))
                    r_moonToPoint[:,0] = grid_point_coordinates[:,0] - moonPos[0]
                    r_moonToPoint[:,1] = grid_point_coordinates[:,1] - moonPos[1]
                    r_moonToPoint[:,2] = grid_point_coordinates[:,2] - moonPos[2]
                    r_moonToSat = [r_spacecraft[0]-moonPos[0],r_spacecraft[1]-moonPos[1],r_spacecraft[2]-moonPos[2]] # vector from moon center to sat 
                    r_moonToSat = np.array(r_moonToSat)
                    
                    part1 = np.sum(np.array(r_spacecraftToPoint)*r_moonToPoint,axis=1)
                    part2 = np.sum(np.abs(r_spacecraftToPoint)**2,axis=-1)**(1./2)
                    part3 = np.sum(np.abs(r_moonToPoint)**2,axis=1)**(1./2)
                    part4 = part2*part3 
                    part5 = part1 / part4
                    angle1 = np.arccos(part5)
        
                    part1B = np.sum(r_moonToSat*r_moonToPoint,axis=1)
                    part2B = np.sum(np.abs(r_moonToSat)**2,axis=-1)**(1./2)
                    part3B = np.sum(np.abs(r_moonToPoint)**2,axis=1)**(1./2)
                    part4B = part2B*part3B 
                    part5B = part1B / part4B
                    angle2 = np.arccos(part5B)

                    idxAngle1 = np.where(angle1 > np.pi/2) # gives us the indices that are
                    idxAngle2 = np.where(angle2 < np.pi/2) # gives us the indices that are
                    idx = np.intersect1d(idxAngle1,idxAngle2)
                    
                    coverage[idx,timeCounterCoverage] = 1
            
                    timeCounterCoverage = timeCounterCoverage+1 

                ################ Phasing multiple satellites  ########################
        
                satellites = numSats[counter]
                satIDXstep = round(loops/satellites)
                for a in range(satellites):
                    for b in range(totalPoints): 
                        point = coverage[b,:]
                        point = reorder(list(point),a*satIDXstep)
                        for c in range(loops): 
                            if point[c] == 1: 
                                coverageMain[b,c] = 1
            
                if np.count_nonzero(coverageMain) == loops*totalPoints: 
                    print('working')
                    break          
                        
            counter = counter + 1
                    
        # print('Calculating FOM...')
        ###################################################
        ############## calculate FOM ######################     done once per design
        ###################################################

        coverageMainNP = np.array(coverageMain)
        if np.count_nonzero(coverageMainNP) == loops * totalPoints:
            percentCoverage = np.zeros(totalPoints) + 100  # % covered
            maxCoverageGap = np.zeros(totalPoints)  # longest coverage gap by each point
            meanCoverageGap = np.zeros(totalPoints)  # average coverage gap for each point
            timeAvgGap = np.zeros(totalPoints)  # time avg gap for each point
            meanResponseTime = np.zeros(totalPoints)
        else:
            percentCoverage = np.zeros(totalPoints)  # % covered
            maxCoverageGap = np.zeros(totalPoints)  # longest coverage gap by each point
            meanCoverageGap = np.zeros(totalPoints)  # average coverage gap for each point
            timeAvgGap = np.zeros(totalPoints)  # time avg gap for each point
            meanResponseTime = np.zeros(totalPoints)
            for i in range(totalPoints):
                arr = coverageMain[i, :]
                checkArr = arr
                checkArr = np.append(checkArr,1)  # adding a 1 to the start and end so the diff function works if theres no LOS to start
                checkArr = np.insert(checkArr, 0, 1)
                # find the indices where the array changes value from 0 to 1
                zero_indices = np.where(checkArr == 1)[0]  # Find the indices of all zeros
                zero_diffs = np.diff(zero_indices)  # Compute the differences between adjacent indices

                maxCoverageGap[i] = (max(np.concatenate(([0], zero_diffs)) - 1))*100  # Find the maximum number of repeated zeros
                zero_gaps = np.concatenate(([0], zero_diffs)) - 1

                ### Getting Mean Coverage Gap ########
                zero_list = np.delete(zero_gaps, [0])  # deleting the first appended element because it always comes to -1 due to line 28, i wont lose any data from this
                numGaps = np.count_nonzero(zero_list)
                if not numGaps == 0:
                    meanCoverageGap[i] = (np.sum(zero_list) / numGaps)*100
                else: 
                    meanCoverageGap[i] = 0
                # meanCoverageGap[i] = (np.sum(zero_list) / numGaps)*100

                ######### GETTING PERCENT COVERAGE ############
                percentCoverage[i] = (np.count_nonzero(arr) / (len(arr))) * 100

                ######## TIME AVG GAP #########
                timeAvgGap[i] = ((np.sum(zero_list ** 2)) / len(arr)) * 100

                ####### MEAN RESPONSE TIME #######
                # using formula for nth trinagular numbers
                # zero array * ((n + (n+1) / 2)
                mask = zero_list != 0
                zero_list_MRT = zero_list[mask]
                MRTvecA = zero_list_MRT + 1  # this is the (n + 1) step
                MRTvecB = MRTvecA * zero_list_MRT  # this is the (n * (n+1)) step
                MRTvecC = MRTvecB / 2  # final MRT vec # this the (n * (n+1)) /2 step

                meanResponseTime[i] = (sum(MRTvecC) / len(arr))*100
        avgPercentCoverage = np.mean(percentCoverage)
        avgMaxCoverageGap = np.mean(maxCoverageGap)
        avgMeanCoverageGap = np.mean(meanCoverageGap)
        avgTimeAvgGap = np.mean(timeAvgGap)
        avgMeanResponseTime = np.mean(meanResponseTime)
        ###############################################
        #################### Score / Return ####################
        ###############################################

        TimeLOS = synodicS*avgPercentCoverage/100
        print(' TimeLOS ', TimeLOS, ' avgPerCov ',avgPercentCoverage, ' avgMeanCovGap ',avgMeanCoverageGap,' avgTimeAvgGap ',avgTimeAvgGap, ' avgMaxCovGap ',avgMaxCoverageGap, ' avgMeanResponseTime ',avgMeanResponseTime)
        orbitDesign.FOM = [TimeLOS, avgPercentCoverage, avgMeanCoverageGap,avgTimeAvgGap, avgMaxCoverageGap, avgMeanResponseTime]
        return orbitDesign

    def driver(self, currDesign):
        #Inputs: Diameter_TxM = antenna diameter of the reciever on the moon in meters
        #DiameterTxO = antenna diameter of the reciever on the sattelite in meters
        #Freq = frequency in Hz
        #DataRate = desired data rate in bps
        #DataRateED = desired data rate for earth downlink (Keep this above 10e6)
        #Range = distance between Tx & Rx in m
        #Range_Sidelink = distance between sattelites in m
        #Time LOS = #of seconds over a target per period

        LU = 389703 # 1 LU (distance from Earth to Moon in km)
        distMoon = [0.98784942, 0, 0] # moon pos in LU
        rMoon = 1737.4/LU # radius of the moon in LU
        orbits = currDesign.orbits
        numSats = currDesign.numSats
        alt = 0
        for orbit in orbits:
            posX = orbit.x
            posY = orbit.y
            posZ = orbit.z
            currAlt = np.sqrt( (posX - distMoon[0])**2 + (posY - distMoon[1])**2 + (posZ - distMoon[2])**2) - rMoon # distance between satellite and moon center
            currAlt = max(currAlt)
            if currAlt > alt:
                alt = currAlt

        alt = alt*389703000 #convert from LU to m  

        Range = alt
        Range_Sidelink = 114039473.3696403 # m

        Diameter_TxM = currDesign.diameterTxM
        DiameterTxO = currDesign.diameterTxO
        DataRate = currDesign.dataRate
        DataRateED = currDesign.dataRate_ED
        FOM_input = currDesign.FOM

        # Call FOM

        TimeLOS = FOM_input[0]
        percentCoverage = FOM_input[1]
        meanCoverageGap = FOM_input[2]
        timeAvgGap = FOM_input[3]
        maxCoverageGap = FOM_input[4]
        meanResponseTime = FOM_input[5]
        coverage = [percentCoverage, meanCoverageGap, timeAvgGap, maxCoverageGap, meanResponseTime]

        #Outputs:
        #margin = 4 element of array of the link margin of all 4 links in order (Lunar Uplink, Lunar downlink, Cross Link, Earth downlink)
        #margin_check = 0 or 1 if the margin will be 6db or over, same order as margin
        #data amount = data amount of the link per day
        
        EbNo = 9.8 #dB change this if BER changes
        Efficiency = .4 #Ka deployable mesh antennas are not smooth
        DesignMargin = 6 #dB - hihg confidence with 6db Margin, complete confidence with 10db link margin
        boltzman = 1.38e-23

        #Preallocating arrays
        margin = []
        # margin_check = []
        DRS_list = []
        constraints = 1

        data_amount = DataRate * TimeLOS #this is the data amount that can be transmitted per period

        
        DiameterRx = DiameterTxO
        index = np.linspace(1,4, num=4)

        for Selection in index:
            #Chosing Selection Constants
            if Selection == 1: #For Lunar Uplink
                Text = 250 #k
                Tphys = 400 #k
                CableLoss = -.25 #dB
                AtmAtten = 0
                Freq = 3e10
                TransmitPower = 100
                Diameter_Tx = Diameter_TxM
                

            if Selection == 2: #For Lunar Downlink
                Text = 25
                Tphys = 400 #4
                CableLoss = -.5 #dB
                AtmAtten = 0
                Freq = 3e10
                TransmitPower = 100
                Diameter_Tx = DiameterTxO

            if Selection == 3: #For Crosslink
                Text = 25
                Tphys = 400 #k
                CableLoss = -.5 #dB
                AtmAtten = 0
                Freq = 3e10
                TransmitPower = 100
                Diameter_Tx = DiameterTxO
                Range = Range_Sidelink

            if Selection == 4: #For Earth Downlink
                Text = 300
                Tphys = 286.15 #k
                CableLoss = -.5 #dB
                AtmAtten = -7 #dB
                TransmitPower = 100
                Diameter_Tx = DiameterTxO
                Freq = 1.8e10
                Range = 3.84e8 - Range
                DataRate = DataRateED
                DiameterRx = 10 #Only link where reciever does not = transmitter

            wavelength = 3e8/Freq #convert freq to wavelength
            #Noise Constants
            Line_Loss = -.5 #dB
            #Antenna Constants
            Tref = 290 #k
            Tr = 289 #k - Noise Temp

            ################# DATA CALCULATIONS #################
            
            CNo_Min = EbNo + 10*log10(DataRate) + DesignMargin #Required CNoReq = Margin + EbnO + Datarate (in db)
            ################# NOISE CALCULATIONS #################
            Tant = Efficiency*Text + (1-Efficiency)*Tphys #Antenna temp in kelvin
            NF = 1 - (Tphys/(Tr)) #Noise Figure
            Ts = Tr + Tant + Text #Reciever System Noise Temp
            No = (10*log10(boltzman*Ts)).real #Reciever System Noise Power

            ################ NOISE PARAMETERS ################
            
            pointingLoss = -3 #dB
            spaceLoss = self.SpaceLoss(Range, wavelength)

            ################# RECEIVER PARAMETERS #################
            #Tx = Transmitter Rx = receiver
            gainTx, areaTx, areaEffRx, beamWidthRx = self.AntennaParameters(Diameter_Tx, wavelength)
            gainRx, areaRx, areaEffRx, beamWidthRx = self.AntennaParameters(DiameterRx, wavelength)

            powerTxDb = (10*log10(TransmitPower)).real
            
            EIRP = (powerTxDb + gainTx + Line_Loss).real
            
            ################## FINAL PARAMETERS ##################

            PropLosses = spaceLoss + AtmAtten 
            # print(PropLosses)
            # print(pointingLoss)
            # print(gainTx)

            ReceivedPower = EIRP + pointingLoss + PropLosses + gainRx
            CNo_received = ReceivedPower - No
            #print(CNo_Min.real)

            marg = (CNo_received - CNo_Min).real #calculating margin
            margin.append(marg)
            #print(marg)

            margin_checker = 0 #Default is 0, margin check of 1 means the link works
            
            if marg >= 3 and Range_Sidelink != 0: #assuming 3dB error correction is added
                margin_checker = 1 #do we have power and margin for the link to exist
            
            #Setting Weighted Average
            if Selection == 1:
                DRS = DataRate *  .35  * margin_checker #this is the data amount that can be transmitted per period
                #^^ sets amnt to 0 if not feasible
                
            if Selection == 2:
                DRS = DataRate *  .35 * margin_checker #this is the data amount that can be transmitted per period
                    
            if Selection == 3:
                DRS = DataRate * .20 * margin_checker #this is the data amount that can be transmitted per period
                
            if Selection == 4: #earth downlink is the least important rate
                DRS = DataRateED *  .10 * margin_checker #this is the data amount that can be transmitted per period
                
            DRS_list.append(DRS)

            if data_amount < 600e6 or margin_checker != 1: #if requirement is not satisfied seat feasibility to 0
                constraints = 0

        DataRateReturn = sum(DRS_list) #return weighted average of dataRates
        currDesign.add_commsObj(coverage, DataRateReturn)
        print("comms score done")
        return data_amount, constraints



