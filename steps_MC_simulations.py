# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:16:22 2024

@author: virgi
"""
import numpy as np
import time 
import random
import events_for_MC_steps

###STEP OF MC SIMULATION IN CASE THERE IS JUST THE A PROTEIN AND NOT THE B

def step_MC_protein_A (time_step, list_DNA, list_A, list_B, list_empty_DNA, E_ad, E_aa, E_ab, residence_times, times_variables):
    
    """
    Performs a single Monte Carlo (MC) step for the binding or unbinding of an A protein on the DNA.

    Parameters:
    ----------
    time_step : int
        The current simulation time step.
    
    list_DNA : numpy.ndarray (1xN)
        A 1D array representing the DNA strand, where each element indicates whether a site is occupied.
    
    list_A : numpy.ndarray (nAx2)
        A 2D array where each row represents an A protein. The first column contains the DNA site bound by A, 
        and the second column stores the index of the B protein it is bound to (-1 if unbound).
    
    list_B : numpy.ndarray (nBxk)
        A 2D array where each row represents a B protein, with columns storing its binding states 
        (the index of the A where it's bound, the corresponding row in list_A').
    
    list_empty_DNA : list
        A list of DNA site indices that are currently unoccupied by A proteins.
    
    E_ad : int
        The energy contribution from A binding directly to DNA.
    
    E_aa : int
        The energy contribution from A-A interactions in the system.
    
    E_ab : int
        The energy contribution from B-A interactions in the system. (energy effect of B on A)
    
    residence_times : numpy.ndarray
        An array tracking the residence times of A proteins on DNA.
    
    times_variables : numpy.ndarray
        An array tracking various time-dependent variables for the simulation.

    Returns:
    -------
    time_step : int
        The updated simulation time step.
    
    list_DNA : numpy.ndarray
        Updated DNA occupancy after the MC step.
    
    list_A : numpy.ndarray
        Updated list of A proteins after the MC step.
    
    list_empty_DNA : list
        Updated list of empty DNA sites.
    
    times_variables : numpy.ndarray
        Updated time-dependent variables.
    
    residence_times : numpy.ndarray
        Updated residence times for A proteins.

    Notes:
    ------
    - A random A protein is selected for the step.
    - A random DNA site is chosen from the available empty sites.
    - A random event is selected: 
        - With 50% probability, an A protein is added to an empty site.
        - Otherwise, an A protein is removed from the DNA.
    - Calls `events_for_MC_steps.add_A()` for adding events.
    - Calls `events_for_MC_steps.remove_A()` for removal events.
    - Increments `time_step` after executing the event.
    """

    
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
        list_DNA, list_A, random_A,list_empty_DNA,residence_time  = events_for_MC_steps.remove_A (list_DNA, list_empty_DNA, list_A , random_A, list_B, residence_times, times_variables, E_ad, E_aa, time_step)
       
    time_step = time_step + 1
    
    return time_step, list_DNA, list_A, list_empty_DNA, times_variables, residence_times


## STEP OF THE SIMULATION WHEN BOTH ARE PRESENT 
def step_MC_proteins_A_B (time_step, list_DNA, list_A, list_B, list_empty_DNA, L, E_ad, E_aa, E_ba, residence_times, times_variables, does_B_bind, k):
    
    """
    Performs a single Monte Carlo (MC) step for the binding or unbinding of A and B proteins on the DNA.
    
    Parameters:
    ----------
    time_step : int
        The current simulation time step.
    
    list_DNA : numpy.ndarray (1xN)
        A 1D array representing the DNA strand, where each element indicates whether a site is occupied.
    
    list_A : numpy.ndarray (nAx2)
        A 2D array where each row represents an A protein. The first column contains the DNA site bound by A, 
        and the second column stores the index of the B protein it is bound to (-1 if unbound).
    
    list_B : numpy.ndarray (nBxk)
        A 2D array where each row represents a B protein, with columns storing its binding states (the index of the A where it's bound, the corresponding row in list_A, -1 if free).    
    list_empty_DNA : list
        A list of DNA site indices that are currently unoccupied by A proteins.
    
    L : int
        The interaction distance used for determining B protein binding availability.
    
    E_ad : float
        The energy contribution from A binding directly to DNA.
    
    E_aa : float
        The energy contribution from A-A interactions in the system.
    
    E_ba : float
        The energy contribution from B-A interactions in the system. (energy effect of B on A)
    
    residence_times : numpy.ndarray
        An array tracking the residence times of A proteins on DNA.
    
    times_variables : numpy.ndarray
        An array tracking various time-dependent variables for the simulation.
    
    does_B_bind : numpy.ndarray
        A 2D array indicating whether B proteins successfully bind during each event.
    
    k : int
        A parameter related to the number of total binding sites for B.
    
    Returns:
    -------
    time_step : int
        The updated simulation time step.
    
    list_DNA : numpy.ndarray
        Updated DNA occupancy after the MC step.
    
    list_A : numpy.ndarray
        Updated list of A proteins after the MC step.
    
    list_B : numpy.ndarray
        Updated list of B proteins after the MC step.
    
    list_empty_DNA : list
        Updated list of empty DNA sites.
    
    times_variables : numpy.ndarray
        Updated time-dependent variables.
    
    residence_times : numpy.ndarray
        Updated residence times for A proteins.
    
    does_B_bind : numpy.ndarray
        Updated binding states of B proteins.
    
    Notes:
    ------
    - A random A protein is selected for the step, and a random DNA site is chosen from available empty sites.
    - A random event is selected:
        - With 50% probability, an A protein is added to an empty site.
        - Otherwise, an A protein is removed from the DNA.
    - Calls `events_for_MC_steps.add_A()` for adding A proteins.
    - Calls `events_for_MC_steps.remove_A()` for removing A proteins.
    - A random B protein is selected, and another event is chosen:
        - With 50% probability, an attempt is made to add a B protein if any binding sites are free.
        - Otherwise, a B protein is removed if any sites are occupied.
    - Calls `events_for_MC_steps.add_B_event()` for adding B proteins.
    - Calls `events_for_MC_steps.remove_B_event()` for removing B proteins.
    - The time step is incremented after executing the event.
    """
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
        list_DNA, list_A, random_A,list_empty_DNA,residence_time  = events_for_MC_steps.remove_A (list_DNA, list_empty_DNA, list_A , random_A, list_B, residence_times, times_variables, E_ad, E_aa, time_step)
    
    # ### B PART ###
    random_B = np.random.randint(0, list_B.shape[0])#choose random B means choosing one between the rows of the nB x K matrix representing Bs  
    
    random_event = np.random.random() 
    
    if random_event < 0.5:  #Adding B event is selected 
        
        if (list_B[random_B, :] == -1).any() :#if there is at least one empty binding site
              list_DNA, list_A, list_B , does_B_bind= events_for_MC_steps.add_B_event( list_DNA, list_A, list_B, random_B, L, does_B_bind)
    
    else: #Removing B event is selected 
    
    
        if (list_B[random_B, :] != -1).any(): #if there is at least one occupied site 
            #print ('Is there a free site ? ', list_B[random_B,:])
    
            list_A, list_B = events_for_MC_steps.remove_B_event(list_DNA, list_A, list_B,random_B, L, E_ba)
    
     
    time_step = time_step + 1
    return time_step, list_DNA, list_A, list_B, list_empty_DNA, times_variables, residence_times, does_B_bind



        
