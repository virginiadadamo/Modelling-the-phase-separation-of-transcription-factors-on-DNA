# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:16:22 2024

@author: virgi
"""
import numpy as np
import time 
import random


def add_A (list_DNA, empty_random_site, list_A , random_A, residence_times, list_empty_DNA, times_variables ):
    list_DNA [empty_random_site] = 1  #always add because it lowers the energy 
    list_A [random_A] = empty_random_site 
    residence_times[random_A] = time_step
    list_empty_DNA.remove (empty_random_site)
    times_variables[random_A]['Count binding events'] += 1
    
def remove_A (list_DNA, empty_random_site, list_A , random_A, residence_times, times_variables ):
    if list_A [random_A] != (-1):  #if A is not free 
        if list_DNA [list_A [random_A]] == 1 : #only A is bound to the DNA

            random_binding = np.random.random()  #draw a random number between 0 and 1
            energy = energy_function (list_A [random_A], list_DNA, E_ad, E_aa)
            
            if random_binding < 1/np.exp(energy):
                remove_A_not_bound_B()
                
        if list_DNA [list_A [random_A]] == 2 : # B protein is also bound 
            #Energy condition 
            remove_A_bound_B ()
        
def remove_A_bound_B (list_DNA, list_A, random_A):
     

def remove_A_not_bound_B (list_DNA, list_A, random_A,list_empty_DNA,residence_times,time_step): 
    list_DNA [list_A [random_A]] = 0 #A is removed 
    list_empty_DNA.append(list_A [random_A])
    list_A [random_A] = -1 #A becomes unbound state
    time_binding = residence_times[random_A]
    residence_times[random_A] = time_step-time_binding
    times_variables[random_A]['Residence times'].append(residence_times[random_A] )
    residence_times[random_A]  = 0 
    
            
def add_B ():

    
def remove_B ():


def step_MC_protein_B (time_step, list_empty_DNA ):
    
    random.seed(time.time())
    
    ### A PART ###
    random_A = np.random.randint(0, len(list_A)) #choose random A
    random_event = np.random.random()  # choose random event - draw a random number between 0 and 1
    random_site_A = np.random.randint(0, len(list_empty_DNA))# Return random integers from low (inclusive) to high (exclusive).
    empty_random_site = list_empty_DNA [random_site_A] #select an empty site 
    if random_event < 0.5:  
        #The adding event is selected 
        if list_A[random_A] == (-1): #A is free
            add_A (list_DNA, empty_random_site, list_A , random_A, residence_times, list_empty_DNA, list_empty_DNA )
    else: 
        #Removing event is selected 
        remove_A ()
    
    ### B PART ###
    random_B = np.random.randint(0, len(list_B))#choose random B 
    random_event = np.random.random()  # choose random event - draw a random number between 0 and 1
    
    
    if random_event < 0.5:  
        #Adding B event is selected 
         add_B()
    
    else:
        remove_B()
    
     
    time_step = time_step + 1
    return time_step,
        
