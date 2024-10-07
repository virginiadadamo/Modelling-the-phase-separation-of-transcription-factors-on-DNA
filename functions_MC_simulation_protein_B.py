# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:16:22 2024

@author: virgi
"""
import numpy as np
import time 
import random




def step_MC (time_step, list_DNA, list_A, list_B, k, list_empty_DNA, E_ad, E_aa, E_ab, E_bb, residence_times, times_variables):
    
   
    random.seed(time.time())
    random_A = np.random.randint(0, len(list_A))
    random_B = np.random.randint(0, len(list_B))
    
    random_site_A = np.random.randint(0, len(list_empty_DNA))# Return random integers from low (inclusive) to high (exclusive). 
    random_site_B = np.random.randint(0, k)
    
    empty_random_site = list_empty_DNA [random_site_A] #select an empty site 
    
    #Randomly select between AddA or RemoveA   
    random_event = np.random.random()  #draw a random number between 0 and 1
    if random_event < 0.5:  
        #The adding event is selected 
        if list_A[random_A] == (-1): #A is free
        
            list_DNA [empty_random_site] = 1  #always add because it lowers the energy 
            list_A [random_A] = empty_random_site 
            residence_times[random_A] = time_step
            list_empty_DNA.remove (empty_random_site)
            times_variables[random_A]['Count binding events'] += 1
            
            
            
                
    else:
        
        #Remove event is selected 
        if list_A [random_A] != (-1):  #if the site is occupied 
            random_binding = np.random.random()  #draw a random number between 0 and 1
            energy = energy_function (list_A [random_A], list_DNA, E_ad, E_aa)
            
            if random_binding < 1/np.exp(energy):
                list_DNA [list_A [random_A]] = 0 #A is removed 
                list_empty_DNA.append(list_A [random_A])
                list_A [random_A] = -1 #A becomes unbound state
                time_binding = residence_times[random_A]
                residence_times[random_A] = time_step-time_binding
                times_variables[random_A]['Residence times'].append(residence_times[random_A] )
                residence_times[random_A]  = 0 
                
                
              
        
    time_step = time_step + 1
    return time_step, list_DNA, list_A, list_empty_DNA, times_variables 
        
def energy_function (index,list_DNA, E_ad, E_aa): 
    
    
    if index == 0: #first site
        energy = list_DNA [index + 1] * E_aa + list_DNA [index]*E_ad
        
    elif index == (len(list_DNA)-1) : #last site 
        energy = list_DNA [index-1]*E_aa + list_DNA [index]*E_ad
    
    else: 
        energy = list_DNA [index + 1] * E_aa + list_DNA [index-1]*E_aa + list_DNA [index]*E_ad
    return energy 
    
