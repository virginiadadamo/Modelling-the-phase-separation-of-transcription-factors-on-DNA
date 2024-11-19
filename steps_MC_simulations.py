# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:16:22 2024

@author: virgi
"""
import numpy as np
import time 
import random
import events_for_MC_steps

def step_MC_protein_A (time_step, list_DNA, list_A, list_B, list_empty_DNA, E_ad, E_aa, E_ab, residence_times, times_variables):
    
    random.seed(time.time())
    random_A = np.random.randint(0, list_A.shape[0])
    random_site = np.random.randint(0, len(list_empty_DNA))# Return random integers from low (inclusive) to high (exclusive).
    empty_random_site = list_empty_DNA [random_site] #select an empty site 
    
    
    #Randomly select between AddA or RemoveA   
    random_event = np.random.random()  #draw a random number between 0 and 1
    if random_event < 0.5:  
        #The adding event is selected 
        list_DNA,list_A, residence_times, list_empty_DNA, times_variables = events_for_MC_steps.add_A (list_DNA, empty_random_site, list_A , random_A, residence_times, list_empty_DNA, times_variables, time_step)
               
    else:
        
        #Remove event is selected 
        list_DNA, list_A, random_A,list_empty_DNA,residence_time  = events_for_MC_steps.remove_A (list_DNA, list_empty_DNA, list_A , random_A, list_B, residence_times, times_variables, E_ad, E_aa, E_ab, time_step)
       
    time_step = time_step + 1
    
    return time_step, list_DNA, list_A, list_empty_DNA, times_variables, residence_times
   
def step_MC_proteins_A_B (time_step, list_DNA, list_A, list_B, list_empty_DNA, L, E_ad, E_aa, E_ba, E_ab, residence_times, times_variables, does_B_bind):
    
    random.seed(time.time())
    random_A = np.random.randint(0, list_A.shape[0])
    random_site = np.random.randint(0, len(list_empty_DNA))# Return random integers from low (inclusive) to high (exclusive).   
   
    empty_random_site = list_empty_DNA [random_site] #select an empty site 
    
    random_event = np.random.random()  #draw a random number between 0 and 1
    
    
    if random_event < 0.5:  
         #The adding event is selected 
             list_DNA,list_A, residence_times, list_empty_DNA, times_variables = events_for_MC_steps.add_A (list_DNA, empty_random_site, list_A , random_A, residence_times, list_empty_DNA, times_variables, time_step)
    else: 
         #Removing event is selected 
        list_DNA, list_A, random_A,list_empty_DNA,residence_time  = events_for_MC_steps.remove_A (list_DNA, list_empty_DNA, list_A , random_A, list_B, residence_times, times_variables, E_ad, E_aa, E_ab, time_step)
    
    # ### B PART ###
    random_B = np.random.randint(0, list_B.shape[0])#choose random B means choosing one between the rows of the nB x K matrix representing Bs  
    
    
    B_bound_sites = np.where (list_B[random_B, :] != -1)[0] #selecting all binding sites which are not free
    
    if len(B_bound_sites) == 0:
        number_of_bound_sites = 0
    else:
        number_of_bound_sites = len(B_bound_sites)
    
    probability_removal = events_for_MC_steps.compute_probability_removal(number_of_bound_sites)
    #print ('probability_removal', probability_removal)
    
    if probability_removal <= 0.5:  #Adding B event is selected 
        
        if (list_B[random_B, :] == -1).any() :#if there is at least one empty binding site
              list_DNA, list_A, list_B , does_B_bind= events_for_MC_steps.add_B_event( list_DNA, list_A, list_B, random_B, L, does_B_bind)
    
    else: #Removing B event is selected 
        if (list_B[random_B, :] != -1).any(): #if there is at least one occupied site 
            #print ('Is there a free site ? ', list_B[random_B,:])
    
            list_A, list_B = events_for_MC_steps.remove_B_event(list_A, list_B,random_B, E_ba)
    
     
    time_step = time_step + 1
    return time_step, list_DNA, list_A, list_B, list_empty_DNA, times_variables, residence_times, does_B_bind



        
