# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 10:08:27 2024

@author: virgi
"""

import numpy as np
import matplotlib.pyplot as plt
import os 
import functions_MC_simulation_protein_A

#PARAMETERS 

nA = 5000 #number of As
N = 10000 #total number of binding sites in the DNA 
stop_time = 100

ignoring_steps = 1
m = 10

#AGGIUNGERE ERRORE

folder_name = 'Simulations_protein_A'

#nA_values = range(0, 100000, 10000)
E_ad_values = range(0, 101, 10) 

#E_ad_values = np.arange(0, 2, 0.5) # E_aa goes from 0 to 1000
E_aa = 0
#E_aa_values = range(0, 101, 10)  # E_ad goes from 0 to 1000

mean_cluster_sizes = [] #mean cluster sizes for each energy value
mean_max_cluster_sizes = [] #mean of the maximum value for cluster for each energy value 

# Iterate over each combination of E_aa and E_ad

for E_ad in E_ad_values:
    
    #for E_aa in E_aa_values:
    
    nA_bound = 0 
    list_DNA = [0] * N #at the initial state there is no bound A to the DNA 
    time_step = 0
    time_step_sampled = [] #To store the time step that are being sampled 
    nA_bound_snapshots = [] # To store the number of A bound for each time step
    group_sizes_snapshots = [] #To store at each time step the corresponding group sizes 
    mean_cluster_sizes_over_time = [] #To store the mean of the cluster sizes at each time step
    max_cluster_sizes = [] #To store at each time step the value of the highest value of cluster size
    all_group_sizes_histogram = [] #union of all the groups sizes list in one list, that will be used to compute the final histogram of distribution of cluster sizes 
    rate_counter = 0 #to count the steps after the first sampled one
    
    
    while time_step < stop_time:
        
        time_step, list_DNA = functions_MC_simulation_protein_A.step_MC(time_step, list_DNA, E_aa, E_ad, nA_bound, nA)
        
        nA_bound= functions_MC_simulation_protein_A.count_consecutive_ones(list_DNA)
        
        if time_step > ignoring_steps:
            
            rate_counter = rate_counter + 1
            
            if rate_counter == m: 
            
                group_sizes, max_count, positions_first_clusters, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time, max_cluster_sizes, rate_counter, all_group_sizes_histogram = functions_MC_simulation_protein_A.take_sample(list_DNA, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time,max_cluster_sizes, rate_counter, all_group_sizes_histogram)
                time_step_sampled.append(time_step)
                
             
        
        elif time_step == ignoring_steps :#or time_step == 0: 
        
            group_sizes, max_count, positions_first_clusters, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time, max_cluster_sizes, rate_counter, all_group_sizes_histogram = functions_MC_simulation_protein_A.take_sample(list_DNA, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time, max_cluster_sizes, rate_counter, all_group_sizes_histogram)
            time_step_sampled.append(time_step)
       
    
    mean_cluster_sizes.append(np.mean(mean_cluster_sizes_over_time)) #taking the mean of custer sizes for each energy value
    mean_max_cluster_sizes.append(np.mean(max_cluster_sizes)) #taking the mean of maximum sized cluster for each energy value 

    
    bin_width = 1
    max_value = max(all_group_sizes_histogram)
    
    bin_edges = np.arange(0, max_value + bin_width, bin_width)  # Bins of width 10
    counts, bins = np.histogram(all_group_sizes_histogram, bins=bin_edges)
    counts = counts / len(time_step_sampled)
    plt.stairs(counts, bins, fill = True)
    plt.title (f'cluster histogram Eaa {E_aa} Ead {E_ad}')
    plt.show()
    
    
   
    # # Save the histogram plot
    # plot_filename = os.path.join(folder_name, f'cluster_histogram_Eaa_{E_aa}_Ead_{E_ad}.png')
        # plt.savefig(plot_filename)
        
       
    # Plotting nA that are bound at each time step 
    time_steps = range(len(nA_bound_snapshots))  
    
    plt.figure(figsize=(8, 6))
    plt.plot(time_step_sampled, nA_bound_snapshots, label='nA_bound_snapshots', color='b', marker='o')
    plt.xlabel('Time Steps')
    plt.ylabel('Number of A bound to DNA')
    plt.title(f'Number of A bound to DNA vs. Time Steps (E_aa={E_aa}, E_ad={E_ad})')
    plt.legend()
    plt.show()
       
    # plt.figure(figsize=(8, 6))
    # plt.plot(time_step_sampled, nA_bound_snapshots, label='nA_bound_snapshots', color='b', marker='o')
    # plt.xlim(2000, max(time_step_sampled))  # Set the x-axis limit from 2000 to the last time step
    # plt.ylim(900, 1010)
    # plt.xlabel('Time Steps')
    # plt.ylabel('Number of A bound to DNA')
    # plt.title(f'Zoomed-in: Number of A bound to DNA vs. Time Steps (E_aa={E_aa}, E_ad={E_ad})')
    # plt.legend()
    # plt.show()    

    print(f"Saved plots for E_aa={E_aa} and E_ad={E_ad}")
    

plt.figure(figsize=(8, 6))
plt.plot(E_ad_values, mean_cluster_sizes, color='r')
plt.xlabel('E_ad')
plt.ylabel('Mean Cluster Size')
plt.title('Mean Cluster Size vs. E_ad')
plt.ylim(0, max(mean_cluster_sizes))  # Set y-axis limits from 0 to max value
plt.grid(True)
plt.show()

# Plot for Mean Max Cluster Size vs. E_ad
plt.figure(figsize=(8, 6))
plt.plot(E_ad_values, mean_max_cluster_sizes, marker='o', color='r')
plt.xlabel('E_ad')
plt.ylabel('Mean Max Cluster Size')
plt.title('Mean Max Cluster Size vs. E_ad')
plt.ylim(0, max(mean_max_cluster_sizes))  # Set y-axis limits from 0 to max value
plt.grid(True)
plt.show()
  
        

        
 
        
     
      