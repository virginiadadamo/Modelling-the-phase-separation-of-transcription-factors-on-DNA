# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 15:56:31 2024

@author: virgi
"""

import numpy as np
import time 
import random

def step_MC (time_step, list_DNA, list_A, list_empty_DNA, E_ad, E_aa, residence_times, times_variables):
    
   
    random.seed(time.time())
    random_A = np.random.randint(0, len(list_A))
    random_site = np.random.randint(0, len(list_empty_DNA))# Return random integers from low (inclusive) to high (exclusive).   
    empty_random_site = list_empty_DNA [random_site] #select an empty site 
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
    

def count_consecutive_ones(DNA_list, return_only_nA= True):
    '''
    Find the size of clusters and the position where they start
    
    Parameters
    ----------
    DNA_list : list of integers (1 or 0)
        list of 1s and 0s. The state = 1 corresponds at the binding of the TF; the state = 0 when it is not bound 

    Returns
    -------
    group_sizes : list of integers
        list containing the size of each cluster
    max_count : integer
        number of clusters 
    positions_first_clusters : list of integers 
        correspond to position in DNA_list of first elements of clusters. Follows Python convention (starting at 0)

    '''
    group_sizes = []
    current_count = 0
    
    positions_first_clusters = []
    
    for i in range(len(DNA_list)):
        if DNA_list[i] == 1:
            if current_count == 0: 
                positions_first_clusters.append(i)
                
            current_count += 1
            
        else:
            if current_count > 0:
                group_sizes.append(current_count)
                current_count = 0

    if current_count > 0:
        group_sizes.append(current_count)

    max_count = max(group_sizes) if group_sizes else 0
    
    if return_only_nA:
        return sum (group_sizes)
    else:
        return group_sizes, max_count, positions_first_clusters, sum (group_sizes)

def count_A (list_A):
    return len([x for x in list_A if x != (-1)])


def take_sample (list_DNA, list_A, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time,max_cluster_sizes, rate_counter, all_group_sizes_histogram):
    
    return_only_nA= False #we want to take samples of all the variables
    
    group_sizes, max_count, positions_first_clusters, nA_bound = count_consecutive_ones(list_DNA, return_only_nA)
    if not group_sizes: #if the group_sizes is empty 
        group_sizes = [0]
    nA_bound_snapshots.append(nA_bound)
    group_sizes_snapshots.append(group_sizes)
    all_group_sizes_histogram = all_group_sizes_histogram + group_sizes
    
    mean_cluster_sizes_over_time.append(np.mean(group_sizes))
    max_cluster_sizes.append(max_count)
    
    rate_counter =1 #everytime take a sample put the counter back to one 
    
    # print(f"Group sizes: {group_sizes}") 
    # print(f"Max count: {max_count}") 
    # print(f"Position of the first member of each cluster: {positions_first_clusters}") 
    
    
    
    
   
    return  group_sizes, max_count, positions_first_clusters, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time, max_cluster_sizes, rate_counter, all_group_sizes_histogram

#def correlation_times:
    
 