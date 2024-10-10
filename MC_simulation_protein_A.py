# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 21:09:15 2024

@author: virgi
"""

import numpy as np
import matplotlib.pyplot as plt
import os 
import functions_MC_simulation_protein_A
import functions_MC_simulation_both 

#PARAMETERS 


alfa = 0.15 #ratio between nA/N 
N = 6000 #total number of binding sites in the DNA
nA = int (N*alfa) #number of As


nB = 0 #number of Bs
k = 0 #number of B interacting sites with A  

#TIME PARAMETERS
stop_time = 2000000
ignoring_steps = 10000
m = 50
no_change_time = 500 #time from which I start computing mean and standard deviation #ricorda che in verità lo calcoli SAMPLING TIME 2000

if (stop_time - ignoring_steps) % m != 0 :
    raise ValueError(f"Error: m ({m}) is not a divisor of {stop_time - ignoring_steps}")

#ENERGY PARAMETERS 

E_ad_values = np.arange(0, 4, 1)
E_aa_values = [0, 2.5]

#CREATING FOLDERS 
folder_name = 'Simulations_protein_A'
subfolder_path = functions_MC_simulation_both.create_folders(folder_name, alfa)
functions_MC_simulation_both.create_txt_parameters(subfolder_path, alfa, stop_time, ignoring_steps)

legend = f'stop_time={stop_time}\nignoring_steps={ignoring_steps}\nm={m}\nnA={nA}\n{N}'


residence_times_for_different_E_aa = []
std_devs_for_different_E_aa = []


for E_aa in E_aa_values: 
    
    nA_bound_for_different_energies = []
    mean_cluster_sizes = [] #mean cluster sizes for each energy value
    mean_max_cluster_sizes = [] #mean of the maximum value for cluster for each energy value 
    mean_no_changes = []
    stdv_no_changes = []
    mean_residence_times = []
    max_residence_times = []
    binding_events_lists = []
    
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
        residence_times = [0]*nA
        times_variables = [{'Index of Transcription Factor': i, 
                        'Residence times': [],
                        'Mean residence time' : 0 ,
                        'Standard deviation':0,
                        'Max residence time' : 0,
                        'Count binding events': 0} for i in range (nA)]
        
        average_cluster_sizes = []
        index_tfs = []
        stdv_each_tf = []
        mean_each_tf = []
        
        while time_step < stop_time:
            
                            
            time_step, list_DNA, list_A, list_empty_DNA, times_variables = functions_MC_simulation_protein_A.step_MC(time_step, list_DNA, list_A, list_empty_DNA, E_ad, E_aa, residence_times, times_variables)
            
            nA_bound_list_A= functions_MC_simulation_both.count_A(list_A)
            nA_bound_DNA = functions_MC_simulation_both.count_consecutive_ones(list_DNA)
            
            if nA_bound_DNA != nA_bound_list_A:
                raise ValueError(f"Error: nA_bound_DNA ({nA_bound_DNA}) and nA_bound ({nA_bound_list_A}) are not equal")
            
            if time_step > ignoring_steps:
                
                rate_counter = rate_counter + 1
                
                if rate_counter == m: 
                
                    group_sizes, max_count, positions_first_clusters, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time, max_cluster_sizes, rate_counter, all_group_sizes_histogram = functions_MC_simulation_both.take_sample(list_DNA, list_A, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time,max_cluster_sizes, rate_counter, all_group_sizes_histogram)
                    time_step_sampled.append(time_step)
                    average_cluster_sizes.append(np.mean (group_sizes))
                    
                 
            
            elif time_step == ignoring_steps :
            
                group_sizes, max_count, positions_first_clusters, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time, max_cluster_sizes, rate_counter, all_group_sizes_histogram = functions_MC_simulation_both.take_sample(list_DNA, list_A, nA_bound_snapshots, group_sizes_snapshots, mean_cluster_sizes_over_time, max_cluster_sizes, rate_counter, all_group_sizes_histogram)
                time_step_sampled.append(time_step)
                average_cluster_sizes.append(np.mean (group_sizes))
        
        #print ((f'Residence times for Eaa {E_aa} and Ead {E_ad} is ', times_variables))
        
        count_never_unbind = np.count_nonzero(residence_times)
        times_never_unbind =  [x for x in residence_times if x != 0]
        
        # with open(txt_filename, 'a') as f:
        #     f.write(f" For Eaa = {E_aa} and Ead = {E_ad} The number of Tfs that never unbinds is = {count_never_unbind}  \n") 
        #     f.write(f" For Eaa = {E_aa} and Ead = {E_ad} Time steps at which it happens = {times_never_unbind}  \n") 
        
        # print ('The number of Tfs that never unbinds is', count_never_unbind )
        # print ('Time steps at which it happens', times_never_unbind)
        
        mean_residence_time = 0
        count_residence_time= 0 
        mean_residence_list = [] 
        count_binding_events_list = []
        
        
        for i in range (nA):
            count_binding_events_list.append (times_variables[i]['Count binding events'])
            if len (times_variables[i]['Residence times']) != 0: 
                times_variables[i]['Mean residence time'] = np.mean (times_variables[i]['Residence times'])
                mean_each_tf.append(times_variables[i]['Mean residence time'])
                times_variables[i]['Standard deviation'] = np.std (times_variables[i]['Residence times'])
                stdv_each_tf.append (times_variables[i]['Standard deviation'])
                times_variables[i]['Max residence time'] = np.max (times_variables[i]['Residence times'])
                
                index_tfs.append (times_variables[i]['Index of Transcription Factor'])
            if times_variables[i]['Mean residence time'] != 0 :
                mean_residence_time += times_variables[i]['Mean residence time'] 
                count_residence_time = count_residence_time +1
                mean_residence_list.append (times_variables[i]['Mean residence time'] )
        
        
        binding_events_lists.append (count_binding_events_list)
        std_mean_residence_times = np.std(mean_residence_list)
        std_devs_for_different_E_aa.append(std_mean_residence_times)
        #print ('Max residence times', [times_variables[i]['Max residence time'] for i in range(nA)])
        max_residence_time = np.max([times_variables[i]['Max residence time'] for i in range(nA)])
        #print('Max of max', max_residence_time)
        
        mean_residence_times.append(mean_residence_time/ count_residence_time)
        max_residence_times.append (max_residence_time)
        nA_bound_for_different_energies.append(nA_bound_snapshots)
        mean_cluster_sizes.append(np.mean(mean_cluster_sizes_over_time)) #taking the mean of custer sizes for each energy value
        mean_max_cluster_sizes.append(np.mean(max_cluster_sizes)) #taking the mean of maximum sized cluster for each energy value 
        
        unique_idx = np.unique(index_tfs)
        if len(unique_idx) != len(index_tfs):
            raise ValueError("Error: the index of TFS are not unique")
        
        plt.scatter(index_tfs, stdv_each_tf, color='r', label=f'TF Index {i}')
        
        plt.xlabel('Index of TF')
        plt.ylabel('Standard Deviation of Residence Times')
        plt.title(f'Standard deviation of Residence times for different TFs (E_aa={E_aa}, E_ad={E_ad})')
        
        plt.text(0.95, 0.05, legend, 
                 horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)
        
        plot_filename = os.path.join(subfolder_path, f'nA_{nA}_n_{N}_stdv_rt_tf__Eaa_{E_aa}_Ead_{E_ad}.png')
        plt.savefig(plot_filename)
        
        plt.show()
        
        plt.clf()
        
        plt.scatter(index_tfs, mean_each_tf, color='r', label=f'TF Index {i}')
        
        plt.xlabel('Index of TF')
        plt.ylabel('Mean of Residence Times')
        plt.title(f'Mean of Residence times for different TFs (E_aa={E_aa}, E_ad={E_ad})')
        
        plt.text(0.95, 0.05, legend, 
                 horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)
        
        
        
        plot_filename = os.path.join(subfolder_path, f'nA_{nA}_n_{N}_mean_rt_tf__Eaa_{E_aa}_Ead_{E_ad}.png')
        plt.savefig(plot_filename)
        
        plt.clf()
        plt.show()
        
        histogram_title = f'Histogram of distribution of binding events Eaa {E_aa} Ead {E_ad}'
        saving_histogram_name = f'nA_{nA}_n_{N}_histogram_binding_events_Eaa_{E_aa}_Ead_{E_ad}.png'
        functions_MC_simulation_both.plot_histogram(count_binding_events_list, histogram_title, legend, subfolder_path, saving_histogram_name, time_step_sampled )
        
                
        # stdv = np.std(nA_bound_snapshots[no_change_time:])
        # mean = np.mean(nA_bound_snapshots[no_change_time:])
        # print ((f'Standard deviation, when there are no more visible changes for Eaa {E_aa} and Ead {E_ad} is ', stdv))
        # print ((f'Mean, when there are no more visible changes for Eaa {E_aa} and Ead {E_ad} is ', mean))
        
        # mean_no_changes.append(mean)
        # stdv_no_changes.append(stdv)
    
        histogram_title = f'Cluster histogram Eaa {E_aa} Ead {E_ad}'
        saving_histogram_name = f'nA_{nA}_n_{N}_cluster_histogram_Eaa_{E_aa}_Ead_{E_ad}.png'
        functions_MC_simulation_both.plot_histogram(count_binding_events_list, histogram_title, legend, subfolder_path, saving_histogram_name, time_step_sampled, True )
        
      
           
    
    residence_times_for_different_E_aa.append (mean_residence_times)
    
    cmap = plt.get_cmap('viridis')
    plt.figure(figsize=(8, 6))
    
    for i, (nA_bound_snapshots, E_ad) in enumerate(zip(nA_bound_for_different_energies, E_ad_values)):
    
        color = cmap(i / len(E_ad_values))
        plt.plot(time_step_sampled, nA_bound_snapshots, label=f'E_ad={E_ad}', color=color, marker='o')
    
    
    
    
    plt.xlabel('Time Steps')
    plt.ylabel('Number of A bound to DNA')
    plt.title(f'Number of A bound to DNA vs. Time Steps (E_aa={E_aa})')
    plt.text(0.95, 0.05, legend, 
              horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)
    plt.legend()
    plot_filename = os.path.join(subfolder_path, f'nA_{nA}_n_{N}_bound_in_time_Eaa_{E_aa}_Ead_{E_ad}.png')
    plt.savefig(plot_filename)
    plt.show()
    
    
    plt.figure(figsize=(8, 6))
    cmap = plt.get_cmap('viridis')

    
    for i, (count_binding_events_list, E_ad) in enumerate(zip(binding_events_lists, E_ad_values)):

        bin_width = 1
        bin_edges = np.arange(min(count_binding_events_list), max(count_binding_events_list) + bin_width, bin_width)
        counts, bins = np.histogram(count_binding_events_list, bins=bin_edges)
        
        color = cmap(i/len(E_ad_values))  
    
        plt.stairs(counts, bins, fill=True, color=color, label=f'E_aa={E_aa}, E_ad={E_ad}')
    
    plt.title('Histogram of Distribution of Binding Events for Different E_aa and E_ad')
    plt.text(0.95, 0.05, legend, 
             horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)
    
    plt.xlabel('Binding Events')
    plt.ylabel('Frequency')
    plt.legend()
    
    plot_filename = os.path.join(subfolder_path, f'nA_{nA}_n_{N}_combined_histogram_binding_events.png')
    plt.savefig(plot_filename)

    plt.show()
    
    # cmap = plt.get_cmap('viridis')
    # plt.figure(figsize=(8, 6))
    
    # for i, (nA_bound_snapshots, E_ad) in enumerate(zip(nA_bound_for_different_energies, E_ad_values)):
        
    #     inverse_nA =  [1 / nA for nA in nA_bound_snapshots]
    #     inverse_time = [1 / time for time in time_step_sampled ]
    
    #     color = cmap(i / len(E_ad_values))
    #     plt.plot(inverse_time, inverse_nA, label=f'E_ad={E_ad}', color=color, marker='o') 
        
    # plt.xlabel('Inverse Time Steps')
    # plt.ylabel('Inverse Number of A bound to DNA')
    # plt.title(f'Lineweaver–Burk plot (E_aa={E_aa})')
    # plt.text(0.95, 0.05, f'stop_time={stop_time}\nignoring_steps={ignoring_steps}\nm={m}\nnA={nA}\n{N}', 
    #           horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)
    # plot_filename = os.path.join(subfolder_path, f'nA_{nA}_n_{N}_bound_in_time_Eaa_{E_aa}_Ead_{E_ad}.png')
    # plt.legend()
    # plt.savefig(plot_filename)
    # plt.show()
    
    plt.figure(figsize=(8, 6))
    plt.plot(max_residence_times, mean_cluster_sizes, color='r')
    plt.xlabel('Max residence times')
    plt.ylabel('Mean Cluster Size')
    plt.title('Mean Cluster Size vs Max Cluster size')
    plt.ylim(min (mean_cluster_sizes), max(mean_cluster_sizes))  
    plt.grid(True)
    plt.show()
    
    
    # plt.figure(figsize=(8, 6))
    # plt.plot(E_ad_values, mean_cluster_sizes, color='r')
    # plt.xlabel('E_ad')
    # plt.ylabel('Mean Cluster Size')
    # plt.title('Mean Cluster Size vs. E_ad for E_aa = {E_aa}')
    # plt.ylim(0, max(mean_cluster_sizes))  # Set y-axis limits from 0 to max value
    # plt.grid(True)
    # plt.show()
    
    
    # plt.figure(figsize=(8, 6))
    # plt.plot(E_ad_values, mean_max_cluster_sizes, marker='o', color='r')
    # plt.xlabel('E_ad')
    # plt.ylabel('Mean Max Cluster Size')
    # plt.title('Mean Max Cluster Size vs. E_ad for E_aa = {E_aa}')
    # plt.ylim(0, max(mean_max_cluster_sizes))  # Set y-axis limits from 0 to max value
    # plt.grid(True)
    # plt.show()
      
    
    # plt.figure(figsize=(8, 6))
    # plt.plot(E_ad_values, stdv_no_changes, color='r')
    # plt.xlabel('E_ad')
    # plt.ylabel('Stdv')
    # plt.title(f'Stdv of A bound after time {no_change_time} for E_aa = {E_aa} ')
    # plt.ylim(0, max(stdv_no_changes))  # Set y-axis limits from 0 to max value
    # plt.grid(True)
    # plt.show()
    
    
    # plt.figure(figsize=(8, 6))
    # plt.plot(E_ad_values, mean_no_changes, color='r')
    # plt.xlabel('E_ad')
    # plt.ylabel('Mean')
    # plt.title(f'Mean of A bound after time {no_change_time} for E_aa = {E_aa} ')
    # plt.ylim(0, max(mean_no_changes))  # Set y-axis limits from 0 to max value
    # plt.grid(True)
    # plt.show() 
    
           
    

cmap = plt.get_cmap('viridis')
plt.figure(figsize=(8, 6))

for i, (mean_residence_times, E_aa) in enumerate(zip(residence_times_for_different_E_aa, E_aa_values)):

    color = cmap(i / len(E_aa_values))
    print ('Log value of mean residence times: ', np.log(mean_residence_times) )
    plt.plot(E_ad_values, np.log(mean_residence_times), 
             label=f'E_aa={E_aa}', color=color, marker='o')

plt.xlabel('E ad values')
plt.ylabel('Log (Mean Residence Time)')
plt.title('Log Mean Residence Time vs. E ad values')

plt.text(0.5, 0.95, f'stop_time={stop_time}\nignoring_steps={ignoring_steps}\nm={m}\nnA={nA}\n{N}', 
         horizontalalignment='center', verticalalignment='top', transform=plt.gca().transAxes)

plt.legend()

plot_filename = os.path.join(subfolder_path, f'nA_{nA}_n_{N}_log_mean_residence_time_Eaa_{E_aa_values[0]}_Ead_{E_ad_values[0]}.png')
plt.savefig(plot_filename)

plt.show()