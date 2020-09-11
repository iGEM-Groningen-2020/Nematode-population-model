# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 21:41:50 2020

@author: Sietse
"""

import numpy as np
import matplotlib.pyplot as plt

#The population model function defines the differential equations and runs the model, storing the data of each time point in a nested list.

default_parameters = [12, 4000000, 0.5, 300, 0.00826, 0.00128, 0.001, 0.002, 0.0099, 0.00215, 0.0099, 0, 1, 0.5]

#The population model has the optional arguments of the parameters so you do not need to add all the parameters each time you want to make the model
#You only specify the parameter that is different from the default settings
def population_model(eggs_initial, dt, soil_type, K = default_parameters[1], b = default_parameters[2], eggs_per_cyst = default_parameters[3], Rae = default_parameters[4], Res = default_parameters[5], Rsc = default_parameters[6], Rcp = default_parameters[7], Ds = default_parameters[8], Rpa = default_parameters[9], Da = default_parameters[10], N = default_parameters[11], Rcr = default_parameters[12], Rrs = default_parameters[13], Dummy = default_parameters[0]):
    soil_properties = {'Sand':[1650, -0.0029, 0.1115, 0.5281],'Silt':[1500, -0.0017, 0.075, 0.6667],'Clay':[1350, -0.0017, 0.1083, -0.25],'Loam':[1500, -0.0033, 0.1833, -1],'Sandy_loam':[1550, -0.0029, 0.1115, 0.5281],'Silty_loam':[1500, -0.0033, 0.1833, -1],'Clay_loam':[1450, -0.0017, 0.075, 0.6667]} #We have included the 7 broad soil types with their densities in kg/m^3
    Ws = 400 * soil_properties[soil_type][0] #we have a plot of 40 x 100 x 100 cm = 0.4 m^3, going from kg to g means multiplying by 1000 so we multiply the density in kg/m^3 by 400 to get the weight of our plot in g
    E = eggs_initial
    S = 0
    C = 0
    P = 0
    A = 0
    X = 0
    R = 0
    MF = 1 #Multiplication factor = ratio of eggs to starting density of eggs; used to compare to existing literature
    #This list will be updated at each point and is what the model returns at the end.
    data_points = [[0],[E],[S],[C],[P],[R],[A],[X],[MF]]

    total_time = int(150/dt) #The total time is always 150 days, which is the time between planting and harvesting.
    for i in range(total_time):
        t = round((i + 1) * dt, 1) #i ranges from 0 to total_time - 1 so I add 1, round is done to reduce the number to one decimal place
        #Because temperature and moisture are not a single constant parameter it is not possible to do the sensitivity analysis in this state.
        #To do sensitivity analysis for temperature and moisture, replace the relevant equation with the Dummy parameter and use that as your variable.
        T = 4.75 + 14.09 * np.sin(0.01057 * t + 0.3018)
        moisture_percent = 30 + 10 * np.sin((2 / 243) * np.pi * t) #The first half year can be approximated by a sine curve that oscillates between 20% and 40% with a period of eight months = 243 days
        motility_index = soil_properties[soil_type][1] * moisture_percent**2 + soil_properties[soil_type][2] * moisture_percent + soil_properties[soil_type][3] #This is a parabolic curve that is dependent on the soil retention curve of the specific soil type
        #Define the differential equations 
        dX = 0.25 * Ws * (A + b * P)/K
        dE = eggs_per_cyst * Rae * T * A - Res * T * E * (1 - X) - 0.02 * E
        dS = Res * T * E * (1 - X) - Ds * S - Rsc * motility_index * T * S + Rrs * R
        dC = Rsc * motility_index * T * S - Rcp * T * C * (1 - X) - Ds * C - Rcr * N * C
        dR = Rcr * N * C - Rrs * R - Ds * R
        dP = Rcp * T * C * (1 - X) - Rpa * T * P * (1 - X)
        dA = 0.5 * Rpa * T * P * (1 - X) - Da * A - Rae * T * A
        #Update values of the state variables
        if X >= 1:
            X = 1 #If the entire plant is infected up to carrying capacity it can no longer support nematodes and it cannot be damaged any further --> X = 1
        else:
            X += dX * dt
        E += dE * dt
        S += dS * dt
        C += dC * dt
        P += dP * dt
        R += dR * dt
        A += dA * dt
        MF = E / eggs_initial

        new_data = [t, E, S, C, R, P, A, X, MF]
        #append the new data to the 2D list for graphing later
        for j in range(9):
            data_points[j].append(new_data[j]) #9 columns: t, E, S, C, P, R, A, X, MF
        
    #yield_loss = (1 - 0.7 * X) * 100
    #print(X)
    # print('Remaining yield is {}%'.format((1 - 0.7 * X) * 100))  
    return data_points
    
#print(population_model(1, 0.1)[0])

def plot_dataset(x_values, y_values, x_label, y_label):
    plt.plot(x_values, y_values)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()
    
# plot_dataset(population_model(5, 0.1, 'Sand')[0], population_model(5, 0.1, 'Sand')[7], 'Time', 'Eggs')

#This function compares initial egg densities to final egg densities to compare our model with Ward (1985)

def final_vs_initial(start_range, end_range, steps, soil_type):
    data_points = [[],[]]
    for i in np.arange(start_range, end_range, steps):
        new_data = [i, (population_model(i, 0.1, soil_type)[1][1500])/i]
        data_points[0].append(new_data[0])
        data_points[1].append(new_data[1])
            
    return data_points

# plot_dataset(final_vs_initial(0, 50, 1, 'Sand')[0], final_vs_initial(0, 50, 1, 'Sand')[1], 'Initial', 'Final')

#This function is used for sensitivity analysis. The model is repeatedly run while changing the value of one parameter and the final result of each run is stored.
#The parameter that you want to analyze must be input as 'parameter' = i
#To see the final egg density choose column 1 of the population model, to see total damage inflicted choose column 7.           

def test_parameter_influence(start_range, end_range, steps, soil_type):
    data_points = [[],[]]
    for i in np.arange(start_range, end_range, steps):
        new_data = [i, population_model(100, 0.1, soil_type, N = i)[1][1500]]
        data_points[0].append(new_data[0])
        data_points[1].append(new_data[1])
        
    return data_points
plot_dataset(test_parameter_influence(0, 5, 0.05, 'Loam')[0], test_parameter_influence(0, 5, 0.05, 'Loam')[1], 'NLP', 'damage')

#The NLP bifurcation function tests the required RootPatch efficacy to keep the damage inflicted upon the potato plant below a certain threshold for a given rate of NLP clearance by the nematodes.

def nlp_bifurcation(start_range, end_range, steps, steps_nlp, eggs_initial, dt, soil_type): 
    #This function will test the required rate of C-individuals becoming R-individuals (= N * Rcr) to let X < 0.5 for a given rate of R-individuals becoming S-individuals (= Rrs)
    #!!! IMPORTANT The starting value of nlp should be very close to, but always lower than, the value of N * Rcr for which an X < 0.5 can be achieved. Otherwise the runtime will be extremely long.
    #!!! IMPORTANT To find a good starting value of nlp, use the test_parameter_influence function to create a graph of X as a function of N and pick the value that closely approximates X = 0.5.
    #!!! IMPORTANT Rrs * dt should NEVER be larger than 1. This will result in more than 100% of repelled nematodes to become unrepelled per time step, leading to negative values of R. With these conditions it cannot be guaranteed that the equations can be solved, leading to runtimes of hours before the programme has exhausted all possible values of N.
    data_points = [[],[]]
    nlp = 0.05
    for i in np.arange(start_range, end_range, steps): #Different values of Rrs are used similar to the sensitivity analysis function.
        while nlp <= 10: 
            if population_model(eggs_initial, dt, soil_type, N = nlp, Rcr = 1, Rrs = i)[7][1500] >= 0.5: #The model is repeatedly run with increasing values for the RootPatch efficacy (N * Rcr) until a certain threshold is met, in this case a damage reduction of 50%
                nlp += steps_nlp #Because the required RootPatch efficacy can only increase as the clearance rate increases, we can ignore any values for N that are lower than the one where the threshold was met for the previous value of Rrs. This massively reduces the search space and runtime.
            else:
                data_points[0].append(i)
                data_points[1].append(nlp)
                #If the threshold is met the function writes down the coordinates and moves to the next value of Rrs.
                break
    return data_points

# plot_dataset(nlp_bifurcation(0, 1, 0.01, 0.0001, 10, 0.1, 'Sand')[0], nlp_bifurcation(0, 1, 0.01, 0.0001, 10, 0.1, 'Sand')[1], 'Recovery rate', 'NLP efficacy')
            
        
        
