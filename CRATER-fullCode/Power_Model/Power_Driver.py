import numpy as np
import math
import mpmath
from os import listdir
from os.path import isfile, join



# Move path to main CRATER directory to import design and orbit classes
import sys
import os
# getting the name of the directory where the this file is present.
current = os.path.dirname(os.path.realpath(__file__))
# Getting the parent directory name where the current directory is present.
parent = os.path.dirname(current) 
# adding the parent directory to the sys.path.
sys.path.append(parent)
from orbitDict import chosenOrbits
from design import designConfiguration

class power:
    def position_eff_func(self, theta, pos_err, point_err, F_disp, h, r):
        # This function determines the efficiency associated with the error in pointing and position
        # error values given to this function are kept constant and positive, allowing those to be varied via calling this func
        # takes heritage valaues for attitude and position knowledge tolerances
        
        # define errors in each dimension
        xpos_err = pos_err[0];
        ypos_err = pos_err[1];
        hpos_err = pos_err[2];
        theta_err = point_err[0];
        phi_err = point_err[1];

        
        # Define/calculate other necessary values, taking errors in
            
        x = h*mpmath.tan(theta);
        d = math.sqrt((x+xpos_err)**2 + (h+hpos_err)**2)
        FOV = 2*(mpmath.atan((x+xpos_err+r)/(h+hpos_err))-mpmath.atan((x+xpos_err)/(h+hpos_err)))
        r_prime = d*mpmath.atan(FOV/2)
        
        # find total errors
        
        xtheta_err = h*mpmath.tan(theta+theta_err)-x;
        x_err = xpos_err + xtheta_err;
        
        yphi_err = d*mpmath.tan(phi_err);
        y_err = ypos_err + yphi_err;

        # define shell radii, & their respective flux %
        shell_num = len(F_disp[0])-1;
        
        r_shell = F_disp[0];
        Fperc_shell = F_disp[1];
        A_tot = np.pi*r_shell[-1]**2
        
        # Preallocation
        A_hit = np.zeros(shell_num);
        A_avail = np.zeros(shell_num);
        hit_eff = np.zeros(shell_num);
        shell_eff = np.zeros(shell_num);
        
        # Loop through shells, shells start at SECOND column for r_shell and F_disp
        
        for shell_index in range(0,shell_num):
            if shell_index == 0:
                r_outer = r_shell[shell_index+1];
                r_inner = 0;
            else:
                r_outer = r_shell[shell_index+1];
                r_inner = r_shell[shell_index];
        
            
            dx = r_outer/10
            dA = dx**2
                
            # loop through the total x of the shell
            for x in np.arange(-r_outer,r_outer,dx):
                
                # loop through the max and min y at x, not including within the shell
                y_max = math.sqrt(r_outer**2 - x**2);
                
                if abs(x) > abs(r_inner):
                    # before scanning through the inner shell
                
                    for y in np.arange(-y_max,y_max,dx):
                    
                        # check if point is on the receiver 1
                        hit_value = ((x-x_err)**2)/(r_prime**2) + ((y-y_err)**2)/(r**2) 
                    
                        if hit_value <= 1:
                            A_hit[shell_index] = A_hit[shell_index] + dA
                            A_avail[shell_index] = A_avail[shell_index] + dA
                        else:
                            A_avail[shell_index] = A_avail[shell_index] + dA
                else:
                    # scans top and bottom parts of the shell
                    
                    y_inner = math.sqrt(r_inner**2 - x**2)
                    
                    for y in np.arange(-y_max,-y_inner,dx):
                    
                        # check if point is on the receiver 2
                        hit_value = ((x-x_err)**2)/(r_prime**2) + ((y-y_err)**2)/(r**2) 
                    
                        if hit_value <= 1:
                            A_hit[shell_index] = A_hit[shell_index] + dA
                            A_avail[shell_index] = A_avail[shell_index] + dA
                        else:
                            A_avail[shell_index] = A_avail[shell_index] + dA
                            
                    for y in np.arange(y_inner,y_max,dx):
                                        
                        # check if point is on the receiver 3
                        hit_value = ((x-x_err)**2)/(r_prime**2) + ((y-y_err)**2)/(r**2) 
                    
                        if hit_value <= 1:
                            A_hit[shell_index] = A_hit[shell_index] + dA
                            A_avail[shell_index] = A_avail[shell_index] + dA
                        else:
                            A_avail[shell_index] = A_avail[shell_index] + dA
                    
            # now, per shell, find out efficiency
        
            if A_avail[shell_index] == 0: # This only happens when the area available to hit is so small the step size rounds it to zero
                # check if center of beam is on the ellips or not, and add to hit if it does
                if ((0-x_err)**2)/(r_prime**2) + ((0-y_err)**2)/(r**2) <= 1:
                    hit_eff[shell_index] = 1
                else:
                    hit_eff[shell_index] = 0
            else:
                hit_eff[shell_index] = A_hit[shell_index]/A_avail[shell_index] # if this is close to 1, then the shell is 100% on the receiver, even with error

            # now, figure out overall percentage of flux that 
            
            shell_eff[shell_index] = hit_eff[shell_index]*Fperc_shell[shell_index+1]

        # sum different flux percentages to get total percentage
        total_eff = sum(shell_eff)
        return total_eff

    def receiver_eff_func(self, theta, zero_loss_eff, b_0):
        # This function determines the efficiency of a receiver as a function of incident angle and flux   
        theta_cutoff = mpmath.acos(b_0/(1+b_0))
        if theta > theta_cutoff:
            n_rec = 0
        else:
            K_theta = 1-b_0*((1/mpmath.cos(theta))-1)
            n_rec = zero_loss_eff*K_theta
        
        return n_rec 

    def Current_Orbit_Values(self, h,t_end,t,V,r):
        # Constants
        r_m = 1737500; # radius of the moon, m
        
        alpha = abs(V*(t-t_end/2)/(r_m+h)) # angle sc and receiver make with center of moon, only care about abs to keep things acurate in the trig funcs below
        
        d = math.sqrt(r_m**2+(r_m+h)**2-2*r_m*(r_m+h)*mpmath.cos(alpha)) # linear distance to receiver
        
        theta_s = mpmath.acos( (d**2+(r_m+h)**2-r_m**2)/(2*d*(r_m+h)) ) # satellite view angle away from nadir (what laser pointing will need to do)
        
        theta_r = alpha + theta_s # incident angle of the receiver
        
        r_prime = r*mpmath.sin(np.pi/2-theta_r) # semi-minor axis of ellipse that receiver appears as to the sc
        r_prime = float(r_prime)
        
        FOV = 2*(mpmath.atan(r_prime/(d-math.sqrt(r**2-r_prime**2)))) # field of view of the receiver
        
        # adjust for negetives in second half of the orbit
        
        if t > t_end/2:
            theta_s = -theta_s
            
        return [d,theta_s,theta_r,FOV,r_prime]

    def gaussNcone_transmission_func(self, r_aperture, r_ave, d_ave, d, P0): 
        r = r_aperture/1.52 # 99% of the beam is in 1.52*r, where r is the radius of the beam IN THE GAUSSIAN FORMULATION
                            # dividing by 1.52 here allows for the calculation of a beam that contains 99% in r_aperature = r*1.52 
                   
        d_lens = 1 # hardcoded value for distance from the laser to the output lens

        # Get the flux distribution and radii for gaussian beam with radius of the output lens, at the output lens
        F_disp = self.gaussian_transmission_func(r, d_lens, d_lens, P0) # this does the gaussian formulation of the math to determine power in each shell
        
        # determine cone shape and radius ratio
        
        alpha = mpmath.atan((r_ave-r_aperture)/d_ave) # view angle, half of FOV basically
        focal_length = r_aperture/mpmath.tan(alpha)
                   
        for i in range(0,len(F_disp[0])): # Calculates the scaled radii of the larger 'spotlight'
            
            r_new = F_disp[0][i]*(d_ave/focal_length + 1) # radius of this ring at average distance
            
            F_disp[0][i] = r_new/d_ave*d # scale this new radius by current distance
        
        # CHECK: is the sum of the shell percentages close to 100%? throw error if not
        tot_shell_perc = sum(F_disp[1][:])-F_disp[1][0]    
        if tot_shell_perc < 0.98:
            # print('Total shell percentage = ',tot_shell_perc)
            raise ValueError('Total Shell Percentage doesnt sum to >98%. Check gaussian_transmission_func & calls.')
            
            
            
        return F_disp

    def gaussian_transmission_func(self, radius, av_d, curr_d, P0):
        # Chosen constants
        Lambda = 532 * 10**(-12) # meters (nanometers)
        
        # Waist of the beam is the size of the aperture
        w0 = radius

        # Calculate Rayleigh length needed for the average distance, wavelength, and target radius
        z_R = (w0**2 * np.pi / Lambda)
            
        # Calculate the actual radius of the beam at the changing distance (near average distance)
        w = w0 * np.sqrt(1 + (curr_d / z_R)**2)
        
        # Choose max radius to calculate intensity to, should be equal to radius of beam at surface & current time
        r_max = w*1.52

        # Initialize variables prior to for loop
        
        N = 10
        r_step = r_max/N
        P_within = []
        r = np.arange(0,r_max+r_step,r_step)
        r_vec = []
        
        # Loop through a range from 1 to selected end radius
        for i in range(0,len(r)):
            
            if i == 0:
                P_within.append(P0)
            elif i == 1:
                # Calculate Power within that radius, percent
                P_within.append((1-np.exp((-2 * r[i]**2) / (w)**2)))
            else:
                P_within.append((1-np.exp((-2 * r[i]**2) / (w)**2)) - (1-np.exp((-2 * r[i-1]**2) / (w)**2)))
                
            r_vec.append(r[i])

        # Convert lists to np arrays
        r = np.array(r)                                
                    
        return [r_vec, P_within]

    def driver(self, currDesign):
        
        # attitude & position errors - keep these at 0 unless analyzing error effects
        pos_err = [0,0,0]
        point_err = [0,0]

    #####################################################################################
    ################################   Constants    #####################################
    #####################################################################################

        # moon
        mu_m = 4.905E12; # Gravitational parameter of the moon, 
        LU = 389703000 # 1 LU to m
        TU = 806.80415 # 1 TU (in seconds)
        r_m = 1737500; # radius of the moon, m
        distMoon = [0.98784942, 0, 0] # moon pos in LU in J2000 frame

        # battery & pane constants
        satLife = float(10); # years
        degPYear = float(0.01); # 1% 
        thetaMax = float(0); # informs peak power production
        I_d = float(0.77); # inherent degradation (0.49-0.88)----------SMAD
        BOLEff = float(0.3); #Beginning of Life Efficiency, 30 % ----- https://www.nasa.gov/smallsat-institute/sst-soa/power 
        BOLPmp = float(400); # W/m^2 ----------------------------------https://www.nasa.gov/smallsat-institute/sst-soa/power 
        specPow = float(100); # W/kg ----------------------------------https://www.nasa.gov/smallsat-institute/sst-soa/power
        DoD = 0.4; # Depth pf Discharge
        LI_battery_upperBound = 0.15; # Battery can't allocate this capacity to anything else
        LI_battery_lowerBound = 0.15; # Battery can't allocate this capacity to anything else
        bounds = LI_battery_upperBound + LI_battery_lowerBound # total sum of the bounds
        SatSurvival = 0.05; # Battery dedicated to onboard computing
        LI_EOL_const = 0.85;  #0.85 is from EOL modeling
        Panel_Eff = 0.32 # solar panel efficiency within the assumed range of +/- 22.5 Degrees
        theta_panel = float(0.4); # Influences cosine loss 22.5 deg worst case -> 0.4 rad
        P_per_kg = 1500 # W/kg
        E_per_kg = 200 # Wh/kg

        # C = 7.5 -> relatively healthy discharge rate
        # battery power and energy density
        # Battery Type
            # specs : https://www.nasa.gov/smallsat-institute/sst-soa/power --> table 3-7
            # Lithium ion -> Power Density: 1500 W/kg
            #                Gravimetric Energy: 100-265 Wh/kg
            # BATTERY CYCLING SOURCE
            # https://iopscience.iop.org/article/10.1149/1945-7111/abf05f/pdf 

        # Comms
        Comm_Power = 100 # watts of constant power draw for comms system
        
        # laser
        laser_loss = 0.55 # percentage of loss of power in the laser itself
        
        # receiver
        rec_zleff = 0.30 # receiver's zero loss efficiency (normal / maximum efficiency of flux power to elec power) 
        rec_b_0 = 0.1 # reflectivity constant, 0.1 for 1 sheet of glass, 0.2 for 2
        rec_I_cutoff = 1380*450 # W/m^2 Max flux receiver can withstand, any higher flux and this is the power accepted. This caps the flux allowed.

    #####################################################################################
    #####################################################################################
    #####################################################################################antenna
        # design variable extraction
        ID = currDesign.ID
        orbits = currDesign.orbits # List of all orbits (family, trajectory, velocity, period, percent eclipsed) in current design
        numSats = currDesign.numSats # Number of satellites on each orbit as list
        totSats = currDesign.totSats # Total number of satellites in constellation 
        solarPanelSize = currDesign.solarPanelSize # Solar panel area [m^2]
        batterySize = currDesign.batterySize # Battery mass [kg]
        laserPower = currDesign.laserPower # Wattage required to power the laser [W]
        r_aperture = currDesign.apetRad # radius of output lens on SC [m]
        r = currDesign.receiverRad_power # radius of ground receiver [m]
        diameterTxM = currDesign.diameterTxM # antenna diameter of the receiver on the moon [m]
        diameterTxO = currDesign.diameterTxO # antenna diameter of the receiver on the satellite [m]
        dataRate = currDesign.dataRate # desired lunar data rate [bps]
        dataRate_ED = currDesign.dataRate_ED # desired data rate for earth downlink [bps]

    #####################################################################################
    #####################################################################################
    #####################################################################################
        feasibility = []
        E_R_tot_lst = []
        for orb in range(len(orbits)):
            orbit = orbits[orb]
            sats = numSats[orb]
            if sats == 0:
                continue
            # orbitDict value extraction
            x = orbit.x
            y = orbit.y
            z = orbit.z
            vx = orbit.vx
            vy = orbit.vy
            vz = orbit.vz
            eclipse_percent = orbit.eclipse
            Period = orbit.T

            southern_hem = [0]
            transmitTime = []
            transmitIdx_lst = []
            trajLen = range(len(z))
            for i in trajLen:
                idx1 = 0
                idx2 = 0
                if z[i] < 0:
                    southern_hem.append(1)
                else:
                    southern_hem.append(0)
                if i == len(z) - 1:
                    southern_hem[i] = 0
                if i > 0:
                    if southern_hem[i] > southern_hem[i-1]:
                        idx1 = i
                    if southern_hem[i] < southern_hem[i-1]:
                        idx2 = i
                        transmitTime.append(abs(idx2-idx1))
                        transmitIdx_lst.append(round(abs((idx2-idx1)/2)))
            if not transmitIdx_lst:
                continue
            transmitIdx = transmitIdx_lst[transmitTime.index(max(transmitTime))]   
            h = np.sqrt( (x[transmitIdx] - distMoon[0])**2 + (y[transmitIdx] - distMoon[1])**2 + (z[transmitIdx] - distMoon[2])**2)*LU - r_m
            V = np.sqrt( (vx[transmitIdx] - distMoon[0])**2 + (vy[transmitIdx] - distMoon[1])**2 + (vz[transmitIdx] - distMoon[2])**2)*(LU/TU)


            

        ###### Calculate power generated by the solar panels ######

            L_d = (1-degPYear)**satLife # % (How much the satellite degrades over a lifetime)

            P_eol = Panel_Eff* BOLPmp * L_d * math.cos(theta_panel) # Specific power at end of life
            P_0 = P_eol*solarPanelSize # power available at end of life, assume this is the power available during the whole lifetime

            P_0_noComm = P_0 - Comm_Power

        ###### Battery losses and allocation ######

            LI_usablebattery_mass = batterySize * 0.5; # Redundancy: Makes sure that there is a secondary battery system if the first fails for some reason

            LI_battery_capacity_total = LI_usablebattery_mass*E_per_kg* LI_EOL_const *(1-SatSurvival-bounds); # [Wh] # same assumption of end of life power output of panels; battery at end of life has LI_EOL_const amount of initial value,
            LI_battery_discharge = LI_usablebattery_mass*P_per_kg * LI_EOL_const;     # [W]  # so entire lifetime we assume we are operating at end of life conditions

            LI_battery_capacity_laser = DoD*LI_battery_capacity_total # energy capacity for the laser

            # feasibility checks
            Feasible = 1
            if (LI_battery_discharge < laserPower): 
                Feasible = 0; # not enough mass of batt to power high wattage laser
                continue

            #Satellite can charge panels in half an orbital period
            E2Batt = P_0_noComm * Period * (1-eclipse_percent); # [Wh] Assume battery charges for half of the orbit period
            if E2Batt < LI_battery_capacity_total:
                Feasible = 0; 
                continue   

            feasibility.append(Feasible)
        ###### calculate max transmission time for current battery specs & receiver specs ######


            # laser loss and maximum discharge time
            L_W = laserPower*(1-laser_loss) # Laser Wattage, this is the battery/capaciter AVERAGE watt output possible, minus the power loss of the laser
            t_max_battery = LI_battery_capacity_laser/laserPower*3600; # max discharge time, equal to maximum transmission time for this battery

            # receiver maximum trasnmission time given receiver efficiency, this takes time of total view into account
            theta_r_max = mpmath.acos(rec_b_0/(1+rec_b_0))
            theta_s_max_receiver = mpmath.asin( r_m* mpmath.sin(np.pi-theta_r_max)/(r_m+h) )
            alpha_max_receiver = theta_r_max-theta_s_max_receiver
            t_max_receiver = 2*(r_m+h)*alpha_max_receiver/V

            t_end = min([t_max_battery,t_max_receiver]) # choose the smallest maximum time possible with given orbit & battery 


        ###### do one pass simulation, calculate orbit averages and focal length to define beam conditions ######


            # time step accuracy
            t_step = t_end/1000; # size of time step, s
            N = t_end/t_step; # number of elements in every model array ALWAYS USE THIS FOR LOOPS
            N = int(N)

            # preallocate
            t = np.zeros(N);
            d = np.zeros(N);
            theta_s = np.zeros(N);
            theta_r = np.zeros(N);
            FOV = np.zeros(N);
            r_prime = np.zeros(N);
            dtheta_s = np.zeros(N);
            dtheta_s_approx = np.zeros(N);
            ddtheta_s_approx = np.zeros(N);
            # loop through time period to figure out average distance and average size of the receiver
            for i in range(0,N):
                t[i] = i*t_step; # calculate current time
                
                try:
                    current_sich = self.Current_Orbit_Values(h,t_end,t[i],V,r)
                except:
                    return 0
                
                
                # split up the output from Current_Orbit_Values into useful values to save
                d[i] = current_sich[0]
                theta_s[i] = current_sich[1]
                theta_r[i] = current_sich[2]
                FOV[i] = current_sich[3]
                r_prime[i] = current_sich[4]
                
            # find pointing speed, approximatly
            if i == 0 :
                dtheta_s_approx[i] = 0
                ddtheta_s_approx[i] = 0
            else :
                dtheta_s_approx[i] = abs(theta_s[i] - theta_s[i-1])
                ddtheta_s_approx[i] = dtheta_s_approx[i] - dtheta_s_approx[i-1]
                
            
            # calculate average distance and r_prime for part 3, and find focal length

            d_ave = np.mean(d)
            r_b = np.mean(r_prime)

            alpha_ave = mpmath.atan((r_b-r_aperture)/d_ave) # ideal view angle, this is the defining angle of the shape of the beam
            focal_length = -r_aperture/mpmath.tan(alpha_ave) # focal length of necessary diverging lens, neg cause diverging
            
        ###### using defined beam conditions, simulate the beam at every point on reciever during transmission ######

            # preallocations
            F_disp = [];
            P_T = [];
            UA_F_disp = [];
            I_ave = np.zeros([N-1,100]);
            I_max = np.zeros(N-1)
            
            # loop through the transmission period 
            # this section also accounts for the max intensity, but also measures the UA = UNADJUSTED flux dispursion to determin soley position error down the line
            for i in range(0,N-1):
            
                current_disp = self.gaussNcone_transmission_func(r_aperture,r_b,d_ave,d[i],L_W)
                current_disp = np.array(current_disp)
            
                # check if intensity of shell is above the maximum, adjust the percentage within to keep I_ave[shell] < rec_I_cutoff
                for j in range(1,len(current_disp[1,:])):
                    P_perc_old = current_disp[1,j]
                    A_shell = (np.pi*(current_disp[0,j]**2-current_disp[0,j-1]**2)) # area of the shell rn
                    I_ave[i,j] = (current_disp[1,0]*current_disp[1,j])/A_shell # (total power * percent of power within shell) / area of shell
            
                    if I_ave[i,j] >= rec_I_cutoff: # we need to reassign the second row of current_disp to rec_I_cutoff = P_within / A_shell
                        P_allowed = rec_I_cutoff*A_shell
                        P_perc_new = P_allowed/current_disp[1,0]
                        current_disp[1,j] = P_perc_new
            
                I_max[i] = max(I_ave[i,:])  
                F_disp.append(current_disp) 
                P_T.append(current_disp[1,0])

            F_disp = np.array(F_disp) # flux dispursion for each time step in a matrix
            UA_F_disp = np.array(UA_F_disp) # this is the incident flux, unadjusted (UA) for the receiver's max intensity
            
            
        ###### Using the flux dispursion at every time step, determine efficiency at every step, and in total ######
            
            # preallocations
            n_rec = np.zeros(N);
            n_pos = np.zeros(N);
            UA_n_pos = np.zeros(N);
            E_R = np.zeros(N);
            E_T = np.zeros(N);
            
            for i in range(0,N-1):

                # this function calculates the efficiency associated with incidence angle and receiver efficiency
                n_rec[i] = self.receiver_eff_func(theta_r[i], rec_zleff, rec_b_0);
                
            
                # this is the position error, without taking max intensity into account -> useful for checking position error effects
                UA_F_disp = self.gaussNcone_transmission_func(r_aperture,r_b,d_ave,d[i],L_W)
                UA_n_pos[i] = self.position_eff_func(theta_r[i], pos_err, point_err, UA_F_disp, h, r); 
                
                # this function will determine the efficiency associated with the pointing and position error of the satellite
                # This also incorperates lost energy from changing apperent reciever size
                n_pos[i] = self.position_eff_func(theta_r[i], pos_err, point_err, F_disp[i,:,:], h, r);

                E_R[i] = t_step*P_T[i]*n_rec[i]*n_pos[i] # total energy recieved, per time step
                E_T[i] = t_step*P_T[i] # Total energy transmit, per time step
            
                
            E_R_tot = sats*sum(E_R); # This is in Joules ->
            E_R_tot_lst.append(E_R_tot*2.7778*10**-7) # kWh

            E_T_tot = sats*sum(E_T); # This is in Joules ->
            # E_T_tot = E_T_tot*2.7778*10**-7 # kWh

            Total_eff = E_R_tot/E_T_tot*100

        # if (not E_R_tot_lst) or (1 not in feasibility):
        #     print("not feasible")
        #     return 0

        if (not E_R_tot_lst):
            return 0

        if (1 not in feasibility):
            return 0
        print("power score done")
        receivedE = max(E_R_tot_lst)
        bestOrbit = orbits[E_R_tot_lst.index(receivedE)] 
        currDesign.add_powerObj(receivedE)
        return [receivedE, bestOrbit.T]
        # total recieved energy


# p = power()   

# allDesigns=[]
# designScores = []
# # # Import orbitDict files
# # files = [131073, 131074, 131075, 131076, 131077, 131078, 131079, 131080, 131081, 131082, 131083, 131084, 131085, 131086, 131087, 131088, 131089, 131090, 131091, 131092, 131093, 131094, 131095, 131096, 131097, 131098, 131099, 131100, 131101, 131102, 131103, 131104, 131105, 131106, 131107, 131108, 211201, 221201, 611201, 621201, 781201]
# # lengthOfFiles = len(files)
# # print(lengthOfFiles, 'designs being tested')
# # for i in files:
# #     filename = 'Bullshit/design' + str(i) + '.dat'
# #     allDesigns.append(chosenOrbits.load(filename))

# # files = [f for f in listdir('Bullshit') if isfile(join('Bullshit', f))]
# # for i in range(len(files)):
# #     filename = 'Bullshit/' + files[i]
# #     if filename != 'Bullshit/.DS_Store':
# #         allDesigns.append(chosenOrbits.load(filename))

# allDesigns.append(chosenOrbits.load(filename))

# lengthOfFiles = len(files)-1

# ###
# ### CDR DESIGN ID: 42
# ###

# design_ID = 42

# #
# # ## Spring semester i dont remember anything
# # ############################################################
# # print('R E M E M B E R : : ID =', allDesigns[design_ID].ID)
# # print('R E M E M B E R : : orbits =', len(allDesigns[design_ID].orbits))
# # print('R E M E M B E R : : numSats =', allDesigns[design_ID].numSats)
# # print('R E M E M B E R : : totSats =', allDesigns[design_ID].totSats)
# # print('R E M E M B E R : : solarPanelSize =', allDesigns[design_ID].solarPanelSize)
# # print('R E M E M B E R : : batterySize =', allDesigns[design_ID].batterySize)
# # print('R E M E M B E R : : laserPower =', allDesigns[design_ID].laserPower)
# # print('R E M E M B E R : : apetRad =', allDesigns[design_ID].apetRad)
# # print('R E M E M B E R : : receiverRad_power =', allDesigns[design_ID].receiverRad_power)
# # print('R E M E M B E R : : diameterTxM =', allDesigns[design_ID].diameterTxM)
# # print('R E M E M B E R : : diameterTxO =', allDesigns[design_ID].diameterTxO)
# # print('R E M E M B E R : : dataRate =', allDesigns[design_ID].dataRate)
# # print('R E M E M B E R : : dataRate_ED =', allDesigns[design_ID].dataRate_ED)
# # print('R E M E M B E R : : commsObj =', allDesigns[design_ID].commsObj)
# # print('R E M E M B E R : : powerObj =', allDesigns[design_ID].powerObj)
# # print('R E M E M B E R : : roiObj =', allDesigns[design_ID].roiObj)
# # print('R E M E M B E R : : constraint =', allDesigns[design_ID].constraint)
# # ############################################################
  


 
# ################### now start the actual power driver ################


# E_tot = []
# N_start = 0
# N_end = lengthOfFiles-1
# #N_end = 1
# print(N_end-N_start, 'designs being tested')
# for i in range(N_start,N_end):
#     design_ID = i

#     design = allDesigns[design_ID]
#     print("______________________________________________________")
#     print('Design',design_ID,':\n')
#     E_tot.append(p.driver(design))

#     # print(p.driver(design))

# max_E = max(E_tot)
# max_E_ID = E_tot.index(max_E)

# print("______________________________________________________")
# print('\nFinal Results for All Designs:')
# print('Maximum Energy achieved by a design:',max_E,'kWh/24h')
# print('Maximum Energy Orbit ID:',max_E_ID)

