# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:16:30 2024

@author: virgi
"""
import os 
import numpy as np
import matplotlib.pyplot as plt

def create_folders (folder_name, alfa): 


    subfolder_name = f'alfa_{alfa}'
    subfolder_path = os.path.join(folder_name, subfolder_name)
    os.makedirs(subfolder_path, exist_ok=True)
    return subfolder_path

def create_txt_parameters (subfolder_path, alfa, stop_time, ignoring_steps):

    txt_filename = os.path.join(subfolder_path, 'simulation_parameters.txt')

    with open(txt_filename, 'a') as f:
        f.write("# Simulation Parameters\n")
        f.write(f"alfa = {alfa}  # Ratio between nA and N\n")
        f.write(f"stop_time = {stop_time}  # Simulation stop time\n")
        f.write(f"ignoring_steps = {ignoring_steps}  # Steps to ignore\n")

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
    

def plot_histogram (list_to_plot, title, legend, subfolder_path, name_to_save, time_step_sampled, mean = False, bin_width = 1): 
    
    max_value = max(list_to_plot)
    min_value = min(list_to_plot)
    
    bin_edges = np.arange(min_value, max_value + bin_width, bin_width)  
     
    counts, bins = np.histogram(list_to_plot, bins=bin_edges)
    if mean:
            counts = counts / len(time_step_sampled)
    plt.stairs(counts, bins, fill = True)
    
    plt.title (title)
    plt.text(0.95, 0.05, legend, 
              horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)
    
    plot_filename = os.path.join(subfolder_path, name_to_save)
    plt.savefig(plot_filename)
    plt.show()

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