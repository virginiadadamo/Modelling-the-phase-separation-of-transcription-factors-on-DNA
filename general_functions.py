# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:16:30 2024

@author: virgi
"""
import os 
import numpy as np
import matplotlib.pyplot as plt
import cluster_class

from scipy.stats import norm, expon, gamma  # Import common distributions


###FUNCTIONS CREATE FOLDERS AND TXT FILES ###


def create_folders(folder_name, alfa, nB, k, L):
    
    """
    Creates a hierarchical folder structure for simulations and returns the path to the innermost folder.
    
    The function first creates a subfolder named `Final_simulations_alfa_{alfa}` inside the specified 
    `folder_name`. Then, within this subfolder, it creates another folder named `nB_{nB}_K_{k}_L_{L}`. 
    If the folders already exist, they are not recreated (to avoid overwriting).
    
    Args:
        folder_name (str): The root folder where the folder structure will be created.
        alfa (float or int): A parameter used to name the first-level subfolder.
        nB (int): A parameter used to name the second-level subfolder.
        k (int): A parameter used to name the second-level subfolder.
        L (int): A parameter used to name the second-level subfolder.
    
    Returns:
        str: The path to the innermost folder `nB_{nB}_K_{k}_L_{L}`.
    """

    # Create the subfolder 
    subfolder_name = f'Final_simulations_alfa_{alfa}'
    subfolder_path = os.path.join(folder_name, subfolder_name)
    
    # Create the subfolder if it doesn't already exist
    os.makedirs(subfolder_path, exist_ok=True)
    
    # Now create a subfolder inside 'alfa_{alfa}'
    L_folder_name = f'nB_{nB}_K_{k}_L_{L}'  
    L_folder_path = os.path.join(subfolder_path, L_folder_name)

    os.makedirs(L_folder_path, exist_ok=True)
    
    # Return the path to the L folder
    return L_folder_path

    
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
    

def gaps_and_consecutives(idx_B_on_DNA):
    
    """
    Analyzes a list representing the DNA sites with B bound to them.
    
    This function takes a list of indices (`idx_B_on_DNA`) and:
        
    1. Identifies groups of consecutives DNA sites.
    2. Measures the lengths of these groups (ignoring groups of length 1).
    3. Calculates the gaps (the number of positions) between groups.
    
    Args:
        idx_B_on_DNA (list[int]): A list of integer indices representing positions on DNA.
    
    Returns:
        tuple:
            - consecutive_B (list[int]): A list of lengths of consecutive groups of indices 
              (groups of length 1 are ignored).
            - gaps (list[int]): A list of gap sizes between consecutive groups of indices.
    
    Example:
        >>> idx_B_on_DNA = [1, 2, 3, 7, 8, 12]
        >>> gaps_and_consecutives(idx_B_on_DNA)
        ([3, 2], [3, 3])
        
        Explanation:
        - There are two consecutive groups: [1, 2, 3] (length 3) and [7, 8] (length 2).
        - There are two gaps: 7 - 3 - 1 = 3 and 12 - 8 - 1 = 3.
    
    Notes:
        - A group is defined as a sequence of consecutive indices.
        - A gap is the difference between two non-consecutive indices, minus 1.
    """
    
    consecutive_B = []
    gaps = []
    start = 0
    
    # Loop through the list to find consecutive groups
    for i in range(1, len(idx_B_on_DNA)):
        if idx_B_on_DNA[i] != idx_B_on_DNA[i-1] + 1:
            # If not consecutive, capture the length of the sequence
            group_length = i - start
            if group_length > 1:  # Skip groups of length 1
                consecutive_B.append(group_length)
            # Calculate the gap difference
            gaps.append((idx_B_on_DNA[i] - idx_B_on_DNA[i-1]-1))
            #The start of the next group
            start = i
    # Add the last consecutive group length if > 1
    group_length = len(idx_B_on_DNA) - start
    if group_length > 1:
        consecutive_B.append(group_length)
    
    return consecutive_B, gaps



def count_fraction_occupied_sites_B(B_list):
    
    """
    Calculate the fraction of occupied sites (non -1 entries) for B in the input array (corresponding to the B matrix).
    
    This function takes a 2D array (B_list), where each row represents a set of sites, and 
    counts the fraction of sites that are not marked as -1. It returns a 1D array containing 
    the fraction of non -1 entries for each row.
    
    Parameters:
    B_list (numpy.ndarray): A 2D NumPy array where each row contains site data, 
                             with -1 indicating an unoccupied site and other values indicating occupied sites.
    
    Returns:
    numpy.ndarray: A 1D array where each element corresponds to the fraction of occupied sites 
                    (i.e., entries that are not -1) in the respective row of B_list.
    """
   
    count = np.count_nonzero(B_list != -1, axis=1)
    count = count / B_list.shape[1]
    
    return np.reshape (count, ((1, B_list.shape[0] )))


def take_sample (time_step,list_DNA, nA_bound_snapshots, average_cluster_sizes, stdv_cluster_sizes, max_cluster_sizes, rate_counter, all_group_sizes_histogram, number_previously_sampled, time_step_sampled):
    
    """
    Takes a sample of the system's state at a given time step, calculating various properties
    of the clusters formed by the DNA sequences, and updates several statistics.
    
    This function computes the group sizes, the maximum group size, and the number of bound 'nA' sites, 
    and updates the provided arrays that track these values over time. It also updates the overall group sizes 
    histogram, as well as the statistics (average, standard deviation, maximum) for cluster sizes.
    
    Parameters:
    time_step (int): The current time step of the simulation.
    list_DNA (list): A list of DNA sequences representing the state of the system at the current time step
    nA_bound_snapshots (numpy.ndarray): A 2D array to store the number of bound 'nA' sites at each sampled time step.
    average_cluster_sizes (numpy.ndarray): A 2D array to store the average cluster size at each sampled time step.
    stdv_cluster_sizes (numpy.ndarray): A 2D array to store the standard deviation of cluster sizes at each sampled time step.
    max_cluster_sizes (numpy.ndarray): A 2D array to store the maximum cluster size at each sampled time step.
    rate_counter (int): A counter tracking the sampling, reset each time a sample is taken.
    all_group_sizes_histogram (list): accumulates the sizes of all groups encountered.
    number_previously_sampled (int): used to index the snapshots.
    time_step_sampled (numpy.ndarray): A 2D array to store the time steps when samples were taken.
    
    Returns:
    tuple: A tuple containing updated values for:
        - group_sizes (list): The sizes of the groups (clusters) found at the current time step.
        - max_count (int): The size of the largest group.
        - nA_bound_snapshots (numpy.ndarray): The updated array of 'nA' site counts over time.
        - average_cluster_sizes (numpy.ndarray): The updated array of average cluster sizes over time.
        - stdv_cluster_sizes (numpy.ndarray): The updated array of standard deviation of cluster sizes over time.
        - max_cluster_sizes (numpy.ndarray): The updated array of maximum cluster sizes over time.
        - rate_counter (int): Reset rate counter (set to 0).
        - all_group_sizes_histogram (list): The updated histogram of group sizes.
        - number_previously_sampled (int): The updated number of previously sampled time steps.
        - time_step_sampled (numpy.ndarray): The updated array of sampled time steps.
    """
    
    return_only_nA= False #we want to take samples of all the variables
    
    group_sizes, max_count, clusters, nA_bound = count_consecutive_ones(list_DNA, return_only_nA)
    if not group_sizes: #if the group_sizes is empty 
        group_sizes = [0]
    
    nA_bound_snapshots[0,number_previously_sampled] = nA_bound
    #group_sizes_snapshots.append(group_sizes)
    all_group_sizes_histogram = all_group_sizes_histogram + group_sizes
    
    average_cluster_sizes[0,number_previously_sampled] = np.mean(group_sizes)
    stdv_cluster_sizes [0,number_previously_sampled] = np.std(group_sizes)
    max_cluster_sizes[0,number_previously_sampled] =max_count
    #clusters_each_time_sampled.append(clusters)
    rate_counter =0 #everytime take a sample put the counter back to zero
    
    # print(f"Max count: {max_count}") 
    # print(f"Position of the first member of each cluster: {positions_first_clusters}") 
    time_step_sampled[0,number_previously_sampled] = time_step
    number_previously_sampled = number_previously_sampled +1 
    return  group_sizes, max_count, nA_bound_snapshots, average_cluster_sizes, stdv_cluster_sizes, max_cluster_sizes, rate_counter, all_group_sizes_histogram , number_previously_sampled, time_step_sampled    
    
    
    
### PLOTS FUNCTIONS ###

def plot_histogram(
    list_to_plot, title, legend, subfolder_path, x_label, y_label, 
    name_to_save, time_step_sampled, mean=False, bin_width=1, 
    fit_distribution=None  
):
    max_value = max(list_to_plot)+1
    min_value = min(list_to_plot)
    
    bin_edges = np.arange(min_value, max_value + bin_width, bin_width)
    
    counts, bins = np.histogram(list_to_plot, bins=bin_edges)
    
    if mean:
        if len(time_step_sampled) > 0:
            
            counts = counts / len(time_step_sampled)
            
    
    # Plot histogram
    plt.stairs(counts, bins, fill=True)

    x_data = bin_edges
    y_data = np.append(counts, 0)  

    # Fit and overlay distribution if specified
    if fit_distribution:
        x_vals = np.linspace(min_value, max_value, 1000)  # Fine grid for PDF
        
        # Fit the specified distribution to data
        if fit_distribution == 'normal':
            mu, std = norm.fit(list_to_plot)
            pdf_vals = norm.pdf(x_vals, mu, std)
            label = f"Normal Fit: μ={mu:.2f}, σ={std:.2f}"
        
        elif fit_distribution == 'exponential':
            loc, scale = expon.fit(list_to_plot)
            pdf_vals = expon.pdf(x_vals, loc, scale)
            label = f"Exponential Fit: λ={1/scale:.2f}"
        
        elif fit_distribution == 'gamma':
            shape, loc, scale = gamma.fit(list_to_plot)
            pdf_vals = gamma.pdf(x_vals, shape, loc, scale)
            label = f"Gamma Fit: shape={shape:.2f}, scale={scale:.2f}"
        
        else:
            raise ValueError(f"Unsupported distribution: {fit_distribution}")
        
        # Normalize PDF to match histogram scale
        pdf_vals *= (bin_width * len(list_to_plot) if not mean else bin_width)
        
        # Plot PDF
        plt.plot(x_vals, pdf_vals, label=label, color="red", linewidth=2)

        # Save PDF values along with histogram if fit_distribution is specified
        x_data = np.concatenate((x_data, x_vals))  # Combine histogram and fit x values
        y_data = np.concatenate((y_data, pdf_vals))  # Combine histogram and fit y values

    # Add labels and title
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend(loc='best')  
    
    # Add legend text (annotation) at the bottom of the figure
    plt.text(
    0.5, -0.11, legend,  # Position text below the plot
    horizontalalignment='center',
    verticalalignment='center',
    transform=plt.gcf().transFigure,
    fontsize=10
)

    
    # Save plot
    plot_filename = os.path.join(subfolder_path, name_to_save + '.png')
    plt.savefig(plot_filename)
    plt.show()
    plt.clf()
    plt.close()
    
    # Save x and y data to txt file
    data_filename = os.path.join(subfolder_path, name_to_save + '.txt')
    with open(data_filename, 'w') as f:
        f.write("x\ty\n")  # Header
        for x_val, y_val in zip(x_data, y_data):
            f.write(f"{x_val}\t{y_val}\n")



def scatter_plot (x,y, legend, subfolder_path, x_label, y_label, title, saving_name):
    plt.scatter(x, y, color='r', s=5) 
   
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
    
    plot_filename = os.path.join(subfolder_path, saving_name  + '.png' )
    plt.savefig(plot_filename)
    
    plt.show()
    plt.clf()
    plt.close()
    
    # Save x and y data to txt file
    data_filename = os.path.join(subfolder_path, saving_name + '.txt')
    with open(data_filename, 'w') as f:
        f.write("x\ty\n")  # Header
        for x_val, y_val in zip(x, y):
            f.write(f"{x_val}\t{y_val}\n")
            
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
    plt.legend(loc='best')
    plot_filename = os.path.join(subfolder_path,saving_name)
    plt.savefig(plot_filename)
    plt.show()
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
    plt.legend(loc='best')
    
    plot_filename = os.path.join(subfolder_path, saving_name)
    plt.savefig(plot_filename)
    plt.show()
    plt.clf()
    plt.close()

def comparison_plots_different_E_aa(list_to_plot, E_aa_values, E_ad_values, xlabel,ylabel,title,legend,subfolder_path,saving_name):
    cmap = plt.get_cmap('viridis')
    plt.figure(figsize=(8,6))
    
    for i, (element, E_aa) in enumerate(zip(list_to_plot, E_aa_values)):
    
        color = cmap(i / len(E_aa_values))
        
        
        plt.plot(E_ad_values, element, 
                  label=f'E_aa={E_aa}', color=color, marker='o')
        
        for x, y in zip(E_ad_values, element):
            plt.text(x, y, f'{y:.3f}', fontsize=8, ha='left', va='bottom', color=color)  #f'{y:.3f}' converts the value of y into a string with exactly two digits after the decimal point.
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    
    
    plt.text(
    0.5, -0.06, legend,  # Position text below the plot
    horizontalalignment='center',
    verticalalignment='center',
    transform=plt.gcf().transFigure,
    fontsize=10
)
    
    plt.legend(loc='best')
    plot_filename = os.path.join(subfolder_path, saving_name)
    plt.savefig(plot_filename)
    plt.show()
    plt.clf()
    
    plt.close()

def plot_figure (x,y,xlabel,ylabel,title,subfolder_path,saving_name, legend, stdv = None, yerr = False) :
    
    plt.figure(figsize=(8, 6))
    if yerr :
       
        plt.errorbar(x[::1000], y[::1000], yerr=stdv[::1000], fmt='o', capsize=5, linestyle='-', color='b', label='Data with error bars')
       
    else: 
        plt.plot(x, y, color='r')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    # Add legend text (annotation) at the bottom of the figure
    plt.text(
    0.5, -0.02, legend,  # Position text below the plot
    horizontalalignment='center',
    verticalalignment='center',
    transform=plt.gcf().transFigure,
    fontsize=10
)
    plot_filename = os.path.join(subfolder_path, saving_name)
    plt.savefig(plot_filename)
    plt.grid(True)
    plt.show()
    plt.clf()
    plt.close()

    
