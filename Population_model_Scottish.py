# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 21:41:50 2020

@author: Sietse
"""

import numpy as np
import matplotlib.pyplot as plt

default_parameters = [600000, 4000000, 0.5, 300, 0.00826, 0.00128, 0.001, 0.002, 0.0099, 0.00215, 0.0099]

#I give the population model the optional arguments of the parameters so you do not need to add all the parameters each time you want to make the model
#You only specify the parameter that is different from the default settings
def population_model(eggs_initial, dt, soil_type, K = default_parameters[1], b = default_parameters[2], eggs_per_cyst = default_parameters[3], Rae = default_parameters[4], Res = default_parameters[5], Rsc = default_parameters[6], Rcp = default_parameters[7], Ds = default_parameters[8], Rpa = default_parameters[9], Da = default_parameters[10]):
    soil_properties = {'Sand':[1650, -0.0029, 0.1115, 0.5281],'Silt':[1500, -0.0017, 0.075, 0.6667],'Clay':[1350, -0.0017, 0.1083, -0.25],'Loam':[1500, -0.0033, 0.1833, -1],'Sandy_loam':[1550, -0.0029, 0.1115, 0.5281],'Silty_loam':[1500, -0.0033, 0.1833, -1],'Clay_loam':[1450, -0.0017, 0.075, 0.6667], 'Scottish':[1650, -0.0029, 0.1115, 0.5281]} 
    #Scottish data comes from Aitkenhead & Coull (2019) from the James Hutton institute
    #The median density was 550 kg/m^3 while the first and third quartiles were 200 and 1000, which is quite a bit lower than what Ilya found as the bulk densities of general soil types
    #We have included the 7 broad soil types with their densities in kg/m^3
    Ws = 400 * soil_properties[soil_type][0] #we have a plot of 40 x 100 x 100 cm = 0.4 m^3, going from kg to g means multiplying by 1000 so we multiply the density in kg/m^3 by 400 to get the weight of our plot in g
    E = eggs_initial
    S = 0
    C = 0
    P = 0
    A = 0
    X = 0
    MF = 1
    data_points = [[0],[E],[S],[C],[P],[A],[X],[MF]]
    

    total_time = int(150/dt)
    for i in range(total_time):
        t = round((i + 1) * dt, 1) #i ranges from 0 to total_time - 1 so I add 1, round is done to reduce the number to one decimal place
        T = 10 + 5 * np.sin((2 / 365) * np.pi * t) #I got my data from metoffice.uk.gov with data from Eastern Scotland since the James Hutton institute is located in Eastern Scotland
        moisture_percent = 17.5 - 7.5 * np.sin((2/365) * np.pi * t) #I got the data from COSMOS-UK 2019, I got the data from Bunny Park since that had a similar moisture regime as Eastern Scotland
        motility_index = soil_properties[soil_type][1] * moisture_percent**2 + soil_properties[soil_type][2] * moisture_percent + soil_properties[soil_type][3]
        # print(motility_index)
        #Define the differential equations 
        dX = 0.15 * Ws * (A + b * P)/K
        dE = eggs_per_cyst * Rae * T * A - Res * T * E * (1 - X) - 0.02 * E
        dS = Res * T * E * (1 - X) - Ds * S - Rsc * motility_index * T * S 
        dC = Rsc * motility_index * T * S - Rcp * T * C * (1 - X) - Ds * C
        dP = Rcp * T * C * (1 - X) - Rpa * T * P * (1 - X)
        dA = 0.5 * Rpa * T * P * (1 - X) - Da * A - Rae * T * A
        #Update values of the state variables
        if X >= 1:
            X = 1 #If the entire plant is infected up to carrying capacity it is dead --> X = 1
        else:
            X += dX * dt
        E += dE * dt
        S += dS * dt
        C += dC * dt
        P += dP * dt
        A += dA * dt
        MF = E / eggs_initial

        new_data = [t, E, S, C, P, A, X, MF]
        #append the new data to the 2D list for graphing later
        for j in range(8):
            data_points[j].append(new_data[j]) #7 columns: t, E, S, P, A, X
        
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
    
# plot_dataset(population_model(10, 0.1, 'Scottish')[0], population_model(10, 0.1, 'Scottish')[7], 'Time', 'Eggs')

#This function compares initial egg densities to final egg densities to compare our model with Ward (1985)
def final_vs_initial(start_range, end_range, steps, soil_type):
    data_points = [[],[]]
    for i in np.arange(start_range, end_range, steps):
        new_data = [i, (population_model(i, 0.1, soil_type)[1][1500])/i]
        data_points[0].append(new_data[0])
        data_points[1].append(new_data[1])
            
    return data_points

plot_dataset(final_vs_initial(0, 50, 0.1, 'Scottish')[0], final_vs_initial(0, 50, 0.1, 'Scottish')[1], 'Initial', 'Final')
            
def test_parameter_influence(parameter, start_range, end_range, steps):
    data_points = [[],[]]
    for i in np.arange(start_range, end_range, steps):
        new_data = [i, population_model(1, 0.1, parameter = i)[1][1500]]
        data_points[0].append(new_data[0])
        data_points[1].append(new_data[1])
        
    return data_points
#plot_dataset(population_model(10,0.1)[0],population_model(10,0.1)[5],'x','y')
#plot_dataset(test_parameter_influence(b, 0, 1, 0.1)[0], test_parameter_influence(b, 0, 1, 0.1)[1], 'b', 'E')
    