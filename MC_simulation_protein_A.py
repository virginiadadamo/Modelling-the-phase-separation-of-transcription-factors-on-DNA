# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 21:09:15 2024

@author: virgi
"""

import numpy as np
import matplotlib.pyplot as plt
import os 
import functions_MC_simulation_protein_A

#PARAMETERS 

nA = 5000 #number of As
N = 10000 #total number of binding sites in the DNA 
stop_time = 60000

ignoring_steps = 1000
m = 100

if (stop_time - ignoring_steps) % m != 0 :
    raise ValueError(f"Error: m ({m}) is not a divisor of {stop_time - ignoring_steps}")

folder_name = 'Simulations_protein_A'

#nA_values = range(0, 100000, 10000)
E_ad_values = range(0, 5, 1) 

E_aa = 8
no_change_time = 500 #time from which I start computing mean and standard deviation #ricorda che in verità lo calcoli SAMPLING TIME 
nA_bound_for_different_energies = []
mean_cluster_sizes = [] #mean cluster sizes for each energy value
mean_max_cluster_sizes = [] #mean of the maximum value for cluster for each energy value 
mean_no_changes = []
stdv_no_changes = []

for E_ad in E_ad_values:
    
    
    nA_bound_DNA = 0 # nA that are bound to the DNA, should be equal to the number of A that are bound in the list of A
    nA_bound_list_A = 0 #As that are bound 
    list_DNA = [0] * N #at the initial state there is no bound A to the DNA 
    list_A = [-1]*nA #-1 for unbound 1 for bound 
    list_empty_DNA  = list(range (0,len(list_DNA),1)) #list containing all the indexes of empty sites in the DNA
    time_step = 0 #initial time step
    time_step_sampled = [] #To store the time step that are being sampled 
    nA_bound_snapshots = [] # To store the number of A bound for each time step
    group_sizes_snapshots = [] #To store at each time step the corresponding group sizes 
    mean_cluster_sizes_over_time = [] #To store the mean of the cluster sizes at each time step
    max_cluster_sizes = [] #To store at each time step the value of the highest value of cluster size
    all_group_sizes_histogram = [] #union of all the groups sizes list in one list, that will be used to compute the final histogram of distribution of cluster sizes 
    rate_counter = 0 #to count the steps after the first sampled one
    
    
    while time_step < stop_time:
        
                        
        time_step, list_DNA, list_A, list_empty_DNA = functions_MC_simulation_protein_A.step_MC(time_step, list_DNA, list_A, list_empty_DNA, E_ad, E_aa)
        
        nA_bound_list_A= functions_MC_simulation_protein_A.count_A(list_A)
        nA_bound_DNA = functions_MC_simulation_protein_A.count_consecutive_ones(list_DNA)
        
        if nA_bound_DNA != nA_bound_list_A:
            raise ValueError(f"Error: nA_bound_DNA ({nA_bound_DNA}) and nA_bound ({nA_bound_list_A}) are not equal")
        
        if time_step > ignoring_steps:
            
            rate_counter = rate_counter + 1
            
            if rate_counter == m: 
            
                group_sizes, max_count, positions_first_clusters, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time, max_cluster_sizes, rate_counter, all_group_sizes_histogram = functions_MC_simulation_protein_A.take_sample(list_DNA, list_A, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time,max_cluster_sizes, rate_counter, all_group_sizes_histogram)
                time_step_sampled.append(time_step)
                
             
        
        elif time_step == ignoring_steps :
        
            group_sizes, max_count, positions_first_clusters, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time, max_cluster_sizes, rate_counter, all_group_sizes_histogram = functions_MC_simulation_protein_A.take_sample(list_DNA, list_A, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time, max_cluster_sizes, rate_counter, all_group_sizes_histogram)
            time_step_sampled.append(time_step)
       
    
    nA_bound_for_different_energies.append(nA_bound_snapshots)
    mean_cluster_sizes.append(np.mean(mean_cluster_sizes_over_time)) #taking the mean of custer sizes for each energy value
    mean_max_cluster_sizes.append(np.mean(max_cluster_sizes)) #taking the mean of maximum sized cluster for each energy value 
    
    
    stdv = np.std(nA_bound_snapshots[no_change_time:])
    mean = np.mean(nA_bound_snapshots[no_change_time:])
    print ((f'Standard deviation, when there are no more visible changes for Eaa {E_aa} and Ead {E_ad} is ', stdv))
    print ((f'Mean, when there are no more visible changes for Eaa {E_aa} and Ead {E_ad} is ', mean))
    
    mean_no_changes.append(mean)
    stdv_no_changes.append(stdv)

    
    bin_width = 1
    max_value = max(all_group_sizes_histogram)
    
    bin_edges = np.arange(1, max_value + bin_width, bin_width)  # Bins of width 10
    counts, bins = np.histogram(all_group_sizes_histogram, bins=bin_edges)
    counts = counts / len(time_step_sampled)
    plt.stairs(counts, bins, fill = True)
    plt.title (f'cluster histogram Eaa {E_aa} Ead {E_ad}')
    plot_filename = os.path.join(folder_name, f'nA_{nA}_cluster_histogram_Eaa_{E_aa}_Ead_{E_ad}.png')
    plt.savefig(plot_filename)
    plt.show()    

cmap = plt.get_cmap('viridis')
plt.figure(figsize=(8, 6))

for i, (nA_bound_snapshots, E_ad) in enumerate(zip(nA_bound_for_different_energies, E_ad_values)):

    color = cmap(i / len(E_ad_values))
    plt.plot(time_step_sampled, nA_bound_snapshots, label=f'E_ad={E_ad}', color=color, marker='o')


plt.xlabel('Time Steps')
plt.ylabel('Number of A bound to DNA')
plt.title(f'Number of A bound to DNA vs. Time Steps (E_aa={E_aa})')
plt.legend()
plot_filename = os.path.join(folder_name, f'nA_{nA}_bound_in_time_Eaa_{E_aa}_Ead_{E_ad}.png')
plt.savefig(plot_filename)
plt.show()

    

plt.figure(figsize=(8, 6))
plt.plot(E_ad_values, mean_cluster_sizes, color='r')
plt.xlabel('E_ad')
plt.ylabel('Mean Cluster Size')
plt.title('Mean Cluster Size vs. E_ad')
plt.ylim(0, max(mean_cluster_sizes))  # Set y-axis limits from 0 to max value
plt.grid(True)
plt.show()


plt.figure(figsize=(8, 6))
plt.plot(E_ad_values, mean_max_cluster_sizes, marker='o', color='r')
plt.xlabel('E_ad')
plt.ylabel('Mean Max Cluster Size')
plt.title('Mean Max Cluster Size vs. E_ad')
plt.ylim(0, max(mean_max_cluster_sizes))  # Set y-axis limits from 0 to max value
plt.grid(True)
plt.show()
  

plt.figure(figsize=(8, 6))
plt.plot(E_ad_values, stdv_no_changes, color='r')
plt.xlabel('E_ad')
plt.ylabel('Stdv')
plt.title(f'Stdv of A bound after time {no_change_time} ')
plt.ylim(0, max(stdv_no_changes))  # Set y-axis limits from 0 to max value
plt.grid(True)
plt.show()


plt.figure(figsize=(8, 6))
plt.plot(E_ad_values, mean_no_changes, color='r')
plt.xlabel('E_ad')
plt.ylabel('Mean')
plt.title(f'Mean of A bound after time {no_change_time} ')
plt.ylim(0, max(mean_no_changes))  # Set y-axis limits from 0 to max value
plt.grid(True)
plt.show()        

        
 