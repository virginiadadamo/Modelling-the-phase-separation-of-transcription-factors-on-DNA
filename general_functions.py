# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:16:30 2024

@author: virgi
"""
import os 
import numpy as np
import matplotlib.pyplot as plt
import cluster_class

###FUNCTIONS CREATE FOLDERS AND TXT FILES ###

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
        

### SAMPLING FUNCTIONS ###

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
    
    #positions_first_clusters = []
    clusters = []
    
    for i in range(DNA_list.shape[1]):
        if DNA_list[0,i] == 1:
            if current_count == 0: 
                #positions_first_clusters.append(i)
                cluster = cluster_class.Cluster (i)
                
            current_count += 1
            
        else:
            if current_count > 0:
                group_sizes.append(current_count)
                cluster.set_size (current_count)
                clusters.append(cluster)
                current_count = 0

    if current_count > 0:
        group_sizes.append(current_count)
        cluster.set_size (current_count)
        clusters.append(cluster)

    max_count = max(group_sizes) if group_sizes else 0
    
    if return_only_nA:
        return sum (group_sizes)
    else:
        return group_sizes, max_count, clusters, sum (group_sizes)
    
def count_A (list_A):
    return len([x for x in list_A[0,:] if x != (-1)])


def take_sample (time_step,list_DNA, list_A, nA_bound_snapshots, average_cluster_sizes,max_cluster_sizes, rate_counter, all_group_sizes_histogram, number_previously_sampled, time_step_sampled):
    
    return_only_nA= False #we want to take samples of all the variables
    
    group_sizes, max_count, clusters, nA_bound = count_consecutive_ones(list_DNA, return_only_nA)
    if not group_sizes: #if the group_sizes is empty 
        group_sizes = [0]
    nA_bound_snapshots[0,number_previously_sampled] = nA_bound
    #group_sizes_snapshots.append(group_sizes)
    all_group_sizes_histogram = all_group_sizes_histogram + group_sizes
    
    average_cluster_sizes[0,number_previously_sampled] = np.mean(group_sizes)
    max_cluster_sizes[0,number_previously_sampled] =max_count
    #clusters_each_time_sampled.append(clusters)
    rate_counter =0 #everytime take a sample put the counter back to zero
    
    # print(f"Group sizes: {group_sizes}") 
    # print(f"Max count: {max_count}") 
    # print(f"Position of the first member of each cluster: {positions_first_clusters}") 
    time_step_sampled[0,number_previously_sampled] = time_step
    number_previously_sampled = number_previously_sampled +1 
    return  group_sizes, max_count, nA_bound_snapshots, average_cluster_sizes, max_cluster_sizes, rate_counter, all_group_sizes_histogram , number_previously_sampled, time_step_sampled    
    
    
    
### PLOTS FUNCTIONS ###
    

def plot_histogram (list_to_plot, title, legend, subfolder_path, x_label, y_label, name_to_save, time_step_sampled, mean = False, bin_width = 1): 
    
    max_value = max(list_to_plot)
    min_value = min(list_to_plot)
    
    bin_edges = np.arange(min_value, max_value + bin_width, bin_width)  
     
    counts, bins = np.histogram(list_to_plot, bins=bin_edges)
    if mean:
        
        if len(time_step_sampled) > 0:
            counts = counts / len(time_step_sampled)
        else:
            print("Warning: time_step_sampled is empty. No mean normalization applied.")
    plt.stairs(counts, bins, fill = True)
    
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title (title)
    plt.text(1, 1, legend, 
         horizontalalignment='right', verticalalignment='top', 
         transform=plt.gca().transAxes)
    
    plot_filename = os.path.join(subfolder_path, name_to_save)
    plt.savefig(plot_filename)
    plt.clf()
    plt.close() 


def scatter_plot (x,y, legend, subfolder_path, x_label, y_label, title, saving_name):
    plt.scatter(x, y, color='r', s=5) 
    plt.xlim(min(x), max(x))
    
    #mean_y = np.mean(y)
    
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    
    
    # plt.text(0.05, 0.95, f'Mean of y axis: {mean_y:.2f}', 
    #          horizontalalignment='left', verticalalignment='top', 
    #          transform=plt.gca().transAxes, fontsize=10, color='blue')
    
    plt.text(1, 1, legend, 
         horizontalalignment='right', verticalalignment='top', 
         transform=plt.gca().transAxes)
    
    plot_filename = os.path.join(subfolder_path, saving_name )
    plt.savefig(plot_filename)
    
    plt.clf()
    plt.close()

def plot_different_Ead_in_time (element_for_different_energies, E_ad_values, time_step_sampled, xlabel, ylabel,title,legend,subfolder_path,saving_name, inverse = False):
        
    cmap = plt.get_cmap('viridis')
    plt.figure(figsize=(8, 6))
    #print(len (element_for_different_energies))
    #print(len (element_for_different_energies[0]))
    time_step_sampled = time_step_sampled.flatten()
    for i, (element, E_ad) in enumerate(zip(element_for_different_energies, E_ad_values)):
        
        color = cmap(i / len(E_ad_values))
        element = element.flatten()
        
        
        if inverse :
            inverse_element = 1 / element
            inverse_time = 1 / time_step_sampled 
            plt.plot(inverse_time, inverse_element, label=f'E_ad={E_ad}', color=color, marker='o')
            
        else:
            plt.plot(time_step_sampled, element, label=f'E_ad={E_ad}', color=color, marker='o')
        
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.text(0.95, 0.05, legend, 
              horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)
    plt.legend()
    plot_filename = os.path.join(subfolder_path,saving_name)
    plt.savefig(plot_filename)
    plt.clf()
    plt.close()
    
def plot_histogram_distribution_different_E_ad (list_to_plot,E_ad_values, E_aa, xlabel, ylabel,title,legend,subfolder_path,saving_name):

    plt.figure(figsize=(8, 6))
    cmap = plt.get_cmap('viridis')

    
    for i, (list_, E_ad) in enumerate(zip(list_to_plot, E_ad_values)):

        bin_width = 1
        bin_edges = np.arange(min(list_), max(list_) + bin_width, bin_width)
        counts, bins = np.histogram(list_, bins=bin_edges)
        
        color = cmap(i/len(E_ad_values))  
    
        plt.stairs(counts, bins, fill=True, color=color, label = f'E_aa={E_aa}, E_ad={E_ad}')
    
    plt.title(title)
    plt.text(0.95, 0.05, legend, 
              horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    
    plot_filename = os.path.join(subfolder_path, saving_name)
    plt.savefig(plot_filename)

    plt.clf()
    plt.close()

def comparison_plots_different_E_aa(list_to_plot, E_aa_values, E_ad_values, xlabel,ylabel,title,legend,subfolder_path,saving_name):
    cmap = plt.get_cmap('viridis')
    plt.figure(figsize=(12,8))
    
    for i, (element, E_aa) in enumerate(zip(list_to_plot, E_aa_values)):
    
        color = cmap(i / len(E_aa_values))
        
        
        plt.plot(E_ad_values, element, 
                  label=f'E_aa={E_aa}', color=color, marker='o')
        
        for x, y in zip(E_ad_values, element):
            plt.text(x, y, f'{y:.3f}', fontsize=8, ha='left', va='bottom', color=color)  #f'{y:.3f}' converts the value of y into a string with exactly two digits after the decimal point.
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    
    
    plt.text(0.5, 0.95, legend, 
              horizontalalignment='center', verticalalignment='top', transform=plt.gca().transAxes)
    
    plt.legend()
    plot_filename = os.path.join(subfolder_path, saving_name)
    plt.savefig(plot_filename)
    
    plt.clf()
    plt.close()

def plot_figure (x,y,xlabel,ylabel,title,subfolder_path,saving_name) :
    
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, color='r')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.ylim(0, max(y))
    plot_filename = os.path.join(subfolder_path, saving_name)
    plt.savefig(plot_filename)
    plt.grid(True)
    plt.clf()
    plt.close()

    
