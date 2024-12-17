# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 21:09:15 2024

@author: virgi
"""

import numpy as np
import os
import steps_MC_simulations
import general_functions


'''
ASSUMPTIONS:
    -For the residence times: Difference between unbinding and binding event of the transcription factor. If a TF binds but never
    unbinds the residence time will not be computed. These trasncription factors will be ignored when computing the mean residence times over all transcription factors 
    -If an A is free and there is an empty site of DNA, protein A will always find that site 
    -E_ab is now set to infinity, id est if a B is bound to an A it will never leave 
'''

# Fare E_ab_values 
#Mettere i plots delle varie Average occupied faction in una sola per diverse E_ba_values 

###PARAMETERS###

alfa = 0.15  #ratio between nA/N 
N = 6000 #total number of binding sites in the DNA
nA = int (N*alfa) #number of As


nB = 200 #number of Bs 
k = 2# 2 #, 5, 10]  #number of B interacting sites with A  
#beta fraction B over As
#Adding protein B in the simulation (True if you want to add, False otherwise)
protein_B= True     


#Time parameters
stop_time = 3000000
ignoring_steps = 1000000
m = 50

number_of_time_steps_sampled = int ((stop_time - ignoring_steps) /m)

if (stop_time - ignoring_steps) % m != 0 :
    raise ValueError(f"Error: m ({m}) is not a divisor of {stop_time - ignoring_steps}")
    
number_of_time_steps_sampled = int ((stop_time - ignoring_steps) /m) 

#Energy parameters 

#E_ad_values = np.arange(0, 4, 1)
E_ad_values = [1]
E_aa_values = [2]


#B parameters 
if protein_B:
    #E_ab = 7
    E_ba_values = [1,4]
    L_values = [10] 
    
else: #put these same parameters to 0 
    E_ab = [0]
    E_ba_values = [0]
    L_values = [0]
    nB=0
    k=0


###PLOTS' TAGS ### - select the plots by putting the corresponding value to true 

#For each combination of Eaa and Ead
if protein_B== False  :
    average_B_fraction = False #plot to identify the average (over time steps sampled) of the fraction of bound sites for each B protein
    plot_B_bound_final_time_step = False 
else: 
    average_B_fraction = True
    plot_B_bound_final_time_step =  True 
    
plot_never_unbind = False #plot to identify the time at which the TFs that never leave the DNA bind
histogram_never_unbind = False #histogram to identify the distribution of Binding Times of the TFs that never leave the DNA
plot_cluster_sizes_over_time = False #To plot the mean cluster size at each time step, with error bars with the corresponding standard deviation 
histogram_mean_residence_time = True # To plot the corresponding distribution of mean residence times of the transcription factors
scatter_plot_std = False #To plot the standard deviation of the residence times of each transcription factors 
scatter_plot_mean = False #To plot the mean of the residence times of each transcription factors
histogram_binding_events = True #To plot the distribution of binding events of the transcription factors 
histogram_cluster_size = True #To plot the distribution of cluster sizes 

#For each E_aa
plot_nA_bound = True #Plot the number of transcription factors bound in time
plot_inverse = False #Plot the inverse of the number of transcription factors vs the inverse of the corresponding time steps 
histogram_binding_events_E_aa = False #Plot the distribution of binding events of the transcription factors for each E_aa
plot_mean_cluster_size_max_cluster_size = False #Plot Mean Cluster Size vs Max Cluster size
plot_mean_cluster_size_E_ad = False #Mean Cluster Size vs. E_ad 
plot_mean_cluster_size_max_cluster_size = False #Mean Cluster Size vs Max Cluster size

#Comparison between different E_aa
plot_log_mean_residence_time = True #Plot for different E_aa the ln of Mean Residence Time vs E_ad values 
plot_stdv_residence_time = False #Plot for different E_aa the Standard Deviation vs E_ad values 
plot_ratio_mean_var = True  #Plot for different E_aa the ratio of Mean Residence Time and Standard Deviation vs E_ad values 


###CREATING FOLDERS###
if protein_B :
    folder_name = 'Simulations_proteins_A_B'
else:
    folder_name = 'Simulations_protein_A'
    



###MONTE CARLO SIMULATION###

results_dict_B = {}

for L  in L_values:

    for E_ba in E_ba_values:
        
                mean_residence_times_for_different_E_aa = [] #To store a list for each E_aa containing the Mean Residence Times for each E_ad values (computed by taking the mean over all the transcription factors over the mean residence time of each transcription factor) - see 
                std_devs_for_different_E_aa = [] #To store a list for each E_aa containing the Mean Standard deviation of Residence Times for each E_ad values (computed by taking the mean over the transcription factors over the standard deviation of residence times of each transcription factor) -see stdv_list
    
                for E_aa in E_aa_values: 
                    
                    nA_bound_for_different_energies = []#To store a list for each E_ad containing the number of transcription factors bound at each time step sampled 
                    mean_cluster_sizes = [] #To store a list for each E_ad containing the mean cluster size (computed by taking the mean over time steps sampled of the mean cluster size at each time sampled) 
                    mean_max_cluster_sizes = [] #To store a list for each E_ad containing the mean of the max cluster sizes (computed by taking the mean over time steps sampled of the max cluster size at each time sampled) 
                    max_residence_times = [] #To store for each E_ad the maximum over all transcription factor between the maximum Residence Times of each transcription factor 
                    binding_events_lists = [] #To store for each E_Ad the number of binding events of each transcription factor 
                    
                    
                    mean_residence_times = [] #To store for each E_ad the mean Residence time (computed by taking the mean over the transcription factors of the mean residence time of each transcription factor)
                    std_devs = [] #To store for each E_ad the Mean Standard deviation of Residence Times (computed by taking the mean over the transcription factors of the standard deviation of residence times of each transcription factor)
                    
                    
                    for E_ad in E_ad_values:
                        
                                
                                subfolder_path = general_functions.create_folders(folder_name, alfa, nB, k, L)
                                legend = (
                                    f"n = {N}\n"
                                    f"E_aa = {E_aa}, E_ad = {E_ad}\n"
                                    f"E_ba = {E_ba}\n"
                                    f"L = {L}\n"
                                    f"k = {k}\n"
                                    f"stop_time = {stop_time}\n"
                                    f"ignoring_steps = {ignoring_steps}"
                                    )
                                
                                ### PARAMETERS FOR THE INITIAL STATE###
                                #Time parameters 
                                time_step = 0 #initial time step
                                number_previously_sampled = 0
                                time_step_sampled = np.zeros((1,number_of_time_steps_sampled))
                                rate_counter = 0 #To count the steps after the first sampled one
                                residence_times = np.zeros((1,nA)) # To store for each transcription the residence times and binding events. It is initialised at 0, if a binding event occurs the time of binding will be store. If an unbinding event then happens the time stored will be used to compute the Residence time and the value will be put to 0 again
                                does_B_bind = np.ones((1,nA)) #To study if a transcription factor never sees a B. Start all at 1, put a 0 if a B binds to that A 
                                times_variables = [{'Index of Transcription Factor': i, 
                                                'Residence times': [],
                                                'Mean residence time' : 0 ,
                                                'Standard deviation':0,
                                                'Max residence time' : 0,
                                                'Count binding events': 0} for i in range (nA)] #Dictionary to store for each Transcription factors the time variables indicated 
                                
                                
                                index_tfs = [] # To store the indexes of the Transcription factor with a least one Residence Time
                                stdv_each_tf = [] # np.zeros((1,nA)) #[] # To store for each transcription factor the Standard Deviation of Residence Times 
                                mean_each_tf = [] # np.zeros((1,nA)) # [] # To store for each transcription factor the Mean Residence time 
                                count_binding_events_list = [] # np.zeros((1,nA))#[] #To store for each transcription factor the number of Binding Events
                                mean_residence_time = 0 #To compute the mean of the residence times over all transcription factors
                                
                                #DNA parameters 
                                list_DNA = np.zeros((1,N)) #Array representing the DNA sites: 0 corresponds to an empty site, 1 to a site with an transcription factor bound to it- at the initial state there is no bound A to the DNA 
                                
                        
                                list_empty_DNA  = list(range (0,list_DNA.shape[1],1)) #list containing all the indexes of empty sites in the DNA
                                nA_bound_DNA = 0 # number of Transcription factors that are bound to the DNA
                                
                                
                                average_cluster_sizes = np.zeros((1,number_of_time_steps_sampled))#[] #To store the mean of the cluster sizes at each time step
                                stdv_cluster_sizes = np.zeros((1,number_of_time_steps_sampled)) #To store the standard deviation of cluster sizes at each time step 
                                max_cluster_sizes = np.zeros((1,number_of_time_steps_sampled))#[] #To store at each time step the value of the highest value of cluster size
                                all_group_sizes_histogram = [] #union of all the cluster sizes list in one list, that will be used to compute the final histogram of distribution of cluster sizes 
                                
                                #Transcription factors parameters 
                                list_A = np.full((nA, 2), -1) #Array representing the transcription factors: -1 for unbound, will store the position of the site on the DNA when bound. The second column will be -1 if the A is not bound to a B, otherwise will store the index of the B to which it is bound
                                
                                nA_bound_snapshots = np.zeros((1,number_of_time_steps_sampled))
                                list_B = np.full((nB, k), -1) #-1 if it's unbound, otherwhise store the A to which it is bound 
                                fraction_occupied_sites = np.zeros ((number_of_time_steps_sampled,nB))
                                
                                #File to store the DNA over time sampled
                                DNA_over_time_file = os.path.join(subfolder_path, f"nA_{nA}_n_{N}_DNA_list_Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}.txt")
                                    
                                while time_step < stop_time:
                                    
                                    if protein_B:
                                        time_step, list_DNA, list_A, list_B, list_empty_DNA, times_variables, residence_times, does_B_bind = steps_MC_simulations.step_MC_proteins_A_B(time_step, list_DNA, list_A, list_B, list_empty_DNA, L, E_ad, E_aa, E_ba, residence_times, times_variables, does_B_bind, k)
                                        
                                        if time_step == stop_time: #last time step
                                            print('Final time step', time_step) 
                                            # print ('List_A', list_A)
                                            # print ('List_B', list_B)
                                            # print ('List_DNA', list_DNA)
                                            indices = np.where(list_A[:, 1] != -1)[0]
                                            idx_DNA_B_bound_final_time_step = []
                                            #print ('Indices', indices)
                                            for idx in indices:
                                                idx_DNA_B_bound_final_time_step.append(list_A[idx, 0])
                                
                                    else:
                                        time_step, list_DNA, list_A, list_empty_DNA, times_variables,residence_times = steps_MC_simulations.step_MC_protein_A(time_step, list_DNA, list_A, list_B, list_empty_DNA, E_ad, E_aa, E_ab, residence_times, times_variables)
                                
                                    
                                    nA_bound_DNA = general_functions.count_consecutive_ones(list_DNA)
                                    
                                    if time_step >= ignoring_steps: 
                                        
                                        rate_counter = rate_counter + 1
                                        
                                        if rate_counter == m: #sampling occurs
                                           
                                            
                                            # with open(DNA_over_time_file, 'a') as file:
                                            
                                                
                                            #     indices = np.where(list_A[:, 1] != -1)[0]  # Returns indices of As bound to B 
                                                
                                            #     for idx in indices:
                                            #         list_DNA[0,list_A[idx, 0]] = 2  #Put the corresponding sites as 2, marked by B protein
                                            #     file.write(f"Time step {time_step}: DNA = ")
                                            
                                            #     # Loop through each element in the DNA and write each element on the same line
                                            #     for row in list_DNA:
                                            #         for element in row:
                                            #             file.write(f"{int(element)} ")  # Write each element, converted to an integer, followed by a space
                                                
                                            #     file.write("\n")  # Newline after writing all elements for the current time step
                                                
                                            #     for idx in indices: #Put everything back as it was before sampling 
                                            #         list_DNA[0,list_A[idx, 0]] = 1
                                            
                                            group_sizes, max_count, nA_bound_snapshots, average_cluster_sizes, stdv_cluster_sizes, max_cluster_sizes, rate_counter, all_group_sizes_histogram, number_previously_sampled, time_step_sampled = general_functions.take_sample(time_step, list_DNA, list_A, nA_bound_snapshots, average_cluster_sizes, stdv_cluster_sizes, max_cluster_sizes, rate_counter, all_group_sizes_histogram, number_previously_sampled, time_step_sampled)
                                            
                                    
                                            if protein_B:
                                                fraction_occupied_sites[number_previously_sampled-1, :]  = general_functions.count_fraction_occupied_sites_B(list_B)
                                                
                                                 
                                            
                                                    #print ('Index of B bound', idx_DNA_B_bound_final_time_step )
                                times_never_unbind =  [x for x in residence_times[0] if x != 0]
                                
                                index_never_unbind = [i for i, x in enumerate(residence_times[0]) if x != 0]
                                
                                index_A_never_bind_B = np.where(does_B_bind[0,:] == 1)[0]
                            
                                file_path = os.path.join(subfolder_path, f'mean_residence_time_A_never_bind_to_B_E_ba={E_ba}_L{L}_E_aa={E_aa}_E_ad={E_ad}.txt')
                        
                                # Open the file in append mode to store mean residence times
                                with open(file_path, 'a') as file:
                        
                                
                                 
                                    for i in range (nA):
                                        
                                        count_binding_events_list.append (times_variables[i]['Count binding events'])
                                        
                                        if len (times_variables[i]['Residence times']) != 0: #if there is at leat one value of Residence time
                                            
                                            index_tfs.append (times_variables[i]['Index of Transcription Factor'])
                                            
                                            times_variables[i]['Mean residence time'] = np.mean (times_variables[i]['Residence times'])
                                            
                                            
                                            mean_each_tf.append(times_variables[i]['Mean residence time'])
                                            
                                            
                                            times_variables[i]['Standard deviation'] = np.std (times_variables[i]['Residence times'])
                                            stdv_each_tf.append (times_variables[i]['Standard deviation'])
                                            
                                            times_variables[i]['Max residence time'] = np.max (times_variables[i]['Residence times'])
                                            
                                            
                                            # Check if the index `i` is in the list of A that never binds to B
                                            if i in index_A_never_bind_B:
                                                # Write the mean residence time to the file
                                                file.write(f"Index {i}: Mean Residence Time = {times_variables[i]['Mean residence time']}\n")
                                        
                                        
                                
                                mean_residence_times.append(np.mean(mean_each_tf)) 
                                std_devs.append(np.mean(stdv_each_tf))
                                binding_events_lists.append (count_binding_events_list) 
                                
                                max_residence_time = np.max([times_variables[i]['Max residence time'] for i in range(nA)])
                                max_residence_times.append (max_residence_time)
                                
                                
                                nA_bound_for_different_energies.append(nA_bound_snapshots) 
                                mean_cluster_sizes.append(np.mean(average_cluster_sizes)) #taking the mean of cluster sizes for each energy value
                                mean_max_cluster_sizes.append(np.mean(max_cluster_sizes)) #taking the mean of maximum sized cluster for each energy value 
                                
                                if protein_B: 
                                    average_occupied_B_fraction_over_time = np.mean(fraction_occupied_sites, axis=0)
                                    
                                    ### PLOTS FOR EACH COMBINATION OF PARAMETERS ###
                                    
                                    consecutive_B, gaps = general_functions.gaps_and_consecutives(sorted(idx_DNA_B_bound_final_time_step))
                                    
                                    key = f"L = {L}, E_ba = {E_ba}, E_aa = {E_aa}, E_ad = {E_ad}"
                    
                                    # Store results in the dictionary
                                    results_dict_B[key] = {
                                        "Mean number of consecutives B": np.mean(consecutive_B),
                                        "Mean width of the gaps between the Bs": np.mean (gaps)
                                    }
                                    
                                    # Optionally log or print progress
                                    print(f"Processed {key}: ncons = {np.mean(consecutive_B)}, ngaps = {np.mean (gaps)}")

                                
                                if plot_B_bound_final_time_step: 
                                    bin_width_values = [1,5,10]
                                    for bin_width in bin_width_values: 
                                        xlabel_plot_B_bound_final_time_step = 'DNA site'
                                        title_plot_B_bound_final_time_step = f'B bound at the last time step (E_aa={E_aa}, E_ad={E_ad})'
                                        saving_name_plot_B_bound_final_time_step = f'nA_{nA}_n_{N}_B_bound_last_time_step_Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}_L_{L}_k_{k}_bin_width_{bin_width}.png' 
                                        general_functions.plot_histogram(idx_DNA_B_bound_final_time_step, title_plot_B_bound_final_time_step, legend, subfolder_path, xlabel_plot_B_bound_final_time_step, 'Frequency', saving_name_plot_B_bound_final_time_step, time_step_sampled, False, bin_width=bin_width) 
                                        
                                if average_B_fraction :
                                    
                                    plot_average_B_fraction_title = f'Avarage over time of occupied sites for each B(E_aa={E_aa}, E_ad={E_ad})'
                                    saving_average_B_fraction_name = f'nA_{nA}_n_{N}_av_frac_bound_sites_Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}_L_{L}_k_{k}.png'
                                    x_label_average_B_fraction = 'Index of B'
                                    y_label_average_B_fraction = 'Avearge number of occupied sites'
                                    
                                    general_functions.plot_figure (np.arange(nB), average_occupied_B_fraction_over_time, x_label_average_B_fraction ,y_label_average_B_fraction,plot_average_B_fraction_title,subfolder_path,saving_average_B_fraction_name)
                                
                               
                                
                                
                                if plot_never_unbind :
                                    plot_never_unbind_title = f'Transcription factor that bind but never unbind (E_aa={E_aa}, E_ad={E_ad})'
                                    saving_never_unbind_name = f'nA_{nA}_n_{N}_never_unbind_Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}_L_{L}_k_{k}.png'
                                    x_label_plot_never_unbind = 'Index of TF'
                                    y_label_plot_never_unbind = 'Binding Time'
                                    
                                    general_functions.plot_figure (index_never_unbind,times_never_unbind, x_label_plot_never_unbind ,y_label_plot_never_unbind,plot_never_unbind_title,subfolder_path,saving_never_unbind_name)
                                
                               
                                if histogram_never_unbind :
                                    histogram_never_unbind_title = f'Histogram for binding times for TFs that bind but never unbind (E_aa={E_aa}, E_ad={E_ad})'
                                    saving_histogram_never_unbind_name = f'nA_{nA}_n_{N}_histo_never_unbind_Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}_L_{L}_k_{k}.png'
                                    x_label_histogram_never_unbind = 'Binding Times'
                                    y_label_histogram_never_unbind = 'Frequency'
                                    general_functions.plot_histogram(times_never_unbind, histogram_never_unbind_title, legend, subfolder_path, x_label_histogram_never_unbind, 'Frequency', saving_histogram_never_unbind_name, time_step_sampled, False, 100) 
                                
                                if plot_cluster_sizes_over_time :
                                    plot_cluster_size_title = f'Mean Cluster size over time (E_aa={E_aa}, E_ad={E_ad})'
                                    saving_cluster_size_name = f'nA_{nA}_n_{N}_mean_cluster_size_Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}_L_{L}_k_{k}.png'
                                    x_label_plot_cluster_size = 'Time_step'
                                    y_label_plot_cluster_size = 'Mean Cluster size'
                                    
                                    general_functions.plot_figure (time_step_sampled[0],average_cluster_sizes[0],x_label_plot_cluster_size,y_label_plot_cluster_size,plot_cluster_size_title,subfolder_path,saving_cluster_size_name, stdv_cluster_sizes[0], True)
                                    
                                    
                                
                                
                                if histogram_mean_residence_time:
                                    histogram_title_mean = f'Histogram of distribution of mean residence times (E_aa={E_aa}, E_ad={E_ad})'
                                    saving_histogram_name = f'nA_{nA}_n_{N}_histo_mean_resident_times_Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}_L_{L}_k_{k}'
                                    x_label_histogram_mean = 'Mean Residence Time'
                                    
                                    general_functions.plot_histogram(mean_each_tf, histogram_title_mean, legend, subfolder_path, x_label_histogram_mean, 'Frequency', saving_histogram_name, time_step_sampled, False, 100, 'normal')
                                
                                unique_idx = np.unique(index_tfs)
                                
                                if len(unique_idx) != len(index_tfs):
                                    raise ValueError("Error: the index of TFS are not unique")
                                    
                                if scatter_plot_std:
                                    xlabel_scatter = 'Index of TF'
                                    ylabel_scatter_stdv = 'Standard Deviation of Residence Times'
                                    title_scatter_stdv = f'Standard deviation of Residence times for different TFs (E_aa={E_aa}, E_ad={E_ad})'
                                    saving_name_scatter_stdv = f'nA_{nA}_n_{N}_stdv_res_times_tf__Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}_L_{L}_k_{k}.png' 
                                    
                                    general_functions.scatter_plot (index_tfs,stdv_each_tf, legend, subfolder_path, xlabel_scatter, ylabel_scatter_stdv, title_scatter_stdv, saving_name_scatter_stdv)
                                
                                if scatter_plot_mean:
                                    xlabel_scatter = 'Index of TF'
                                    ylabel_scatter_mean = 'Mean of Residence Times'
                                    title_scatter_mean = f'Mean of Residence times for different TFs (E_aa={E_aa}, E_ad={E_ad})'
                                    saving_name_scatter_mean = f'nA_{nA}_n_{N}_mean_res_times_tf__Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}_L_{L}_k_{k}.png'
                                    
                                    general_functions.scatter_plot (index_tfs,mean_each_tf, legend, subfolder_path, xlabel_scatter, ylabel_scatter_mean, title_scatter_mean, saving_name_scatter_mean)
                        
                                if histogram_binding_events:
                                    histogram_title_be = f'Histogram of distribution of binding events (Eaa {E_aa}, Ead {E_ad})'
                                    saving_histogram_name_be = f'nA_{nA}_n_{N}_histo_binding_events_Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}_L_{L}_k_{k}'
                                    x_label_histogram_be = 'Binding Events'
                                    
                                    general_functions.plot_histogram(count_binding_events_list, histogram_title_be, legend, subfolder_path,x_label_histogram_be,'Frequency', saving_histogram_name_be, time_step_sampled[0], False, 1, 'normal' )
                                    del count_binding_events_list
                                    
                                
                                if histogram_cluster_size:
                                    histogram_title_cluster = f'Cluster histogram (Eaa {E_aa}, Ead {E_ad})'
                                    saving_histogram_name_cluster = f'nA_{nA}_n_{N}_cluster_histo_Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}_L_{L}_k_{k}'
                                    x_label_histogram_cluster = 'Cluster size'
                                    
                                    general_functions.plot_histogram(all_group_sizes_histogram, histogram_title_cluster, legend, subfolder_path,x_label_histogram_cluster, 'Frequency', saving_histogram_name_cluster, time_step_sampled[0], True )
                                    del all_group_sizes_histogram       
                            
                        
                            
                    std_devs_for_different_E_aa.append(std_devs)       
                    mean_residence_times_for_different_E_aa.append (mean_residence_times)
                    del index_tfs,std_devs, mean_residence_times
                        
                    ### PLOTS FOR EACH E_AA 
                    if plot_nA_bound:
                        xlabel_nA_bound_snapshots = 'Time Steps'
                        ylabel_nA_bound_snapshots = 'Number of A bound to DNA'
                        title_nA_bound_snapshots = f'Number of A bound to DNA vs. Time Steps (E_aa={E_aa})'
                        saving_name_nA_bound_snapshots = f'nA_{nA}_n_{N}_bound_in_time_Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}.png'
                        
                        general_functions.plot_different_Ead_in_time (nA_bound_for_different_energies, E_ad_values, time_step_sampled[0], xlabel_nA_bound_snapshots, ylabel_nA_bound_snapshots,title_nA_bound_snapshots,legend,subfolder_path,saving_name_nA_bound_snapshots)
                        del nA_bound_for_different_energies
                        
                    if plot_inverse:
                        xlabel_lineweaver_burk = 'Inverse Time Steps'
                        ylabel_lineweaver_burk = 'Inverse Number of A bound to DNA'
                        title_lineweaver_burk = f'Lineweaver–Burk plot (E_aa={E_aa})'
                        saving_name_lineweaver_burk = f'nA_{nA}_n_{N}_Lineweaver_Burk_Eaa_{E_aa}_Ead_{E_ad}_E_ba_{E_ba}.png'
                        
                        general_functions.plot_different_Ead_in_time (nA_bound_for_different_energies, E_ad_values, time_step_sampled[0], xlabel_lineweaver_burk, ylabel_lineweaver_burk,title_lineweaver_burk,legend,subfolder_path,saving_name_lineweaver_burk, True)
                        del nA_bound_for_different_energies
                    
                    if histogram_binding_events_E_aa:
                        xlabel_binding_events = 'Binding Events'
                        ylabel_binding_events = 'Frequency'
                        title_binding_events = f'Histogram of Distribution of Binding Events for different E_ad (E_aa={E_aa})'
                        saving_name_binding_events = f'nA_{nA}_n_{N}_combined_histo_binding_events_E_ba_{E_ba}.png'
                      
                        general_functions.plot_histogram_distribution_different_E_ad (binding_events_lists,E_ad_values, E_aa, xlabel_binding_events, ylabel_binding_events,title_binding_events,legend,subfolder_path,saving_name_binding_events)
                        del binding_events_lists
                    
                    if plot_mean_cluster_size_max_cluster_size:
                        general_functions.plot_figure (max_residence_times,mean_cluster_sizes,'Max residence times','Mean Cluster Size','Max Residence Time vs Max Cluster size (E_aa={E_aa})',subfolder_path,f'nA_{nA}_n_{N}_max_rt_mean_cluster_size.png') 
                    if plot_mean_cluster_size_E_ad:
                        general_functions.plot_figure (E_ad_values,mean_cluster_sizes,'E_ad','Mean Cluster Size',f'Mean Cluster Size vs. E_ad (E_aa={E_aa})',subfolder_path,f'nA_{nA}_n_{N}_E_ad_mean_cluster_size.png') 
                   
                ###COMPARISON BETWEEN DIFFERENT E_AA ###
                x_label_E_aa_comparison = 'E ad values'
                
                if plot_log_mean_residence_time:
                    log_list = [np.log(element) for element in mean_residence_times_for_different_E_aa]
                    ylabel_log = 'Log (Mean Residence Time)'
                    title_log = 'Log Mean Residence Time vs. E ad values'
                    saving_name_log = f'nA_{nA}_n_{N}_log_mean_residence_time_E_ba_{E_ba}.png'
                    
                    general_functions.comparison_plots_different_E_aa(log_list, E_aa_values, E_ad_values, x_label_E_aa_comparison, ylabel_log,title_log,legend,subfolder_path,saving_name_log)
                
                if plot_stdv_residence_time:
                    ylabel_stdv = 'Stdv residence times'
                    title_stdv = 'Mean standard deviation of residence times vs. E ad values'
                    saving_name_stdv = f'nA_{nA}_n_{N}_mean_stdv_residence_times_E_ba_{E_ba}.png'
                    
                    general_functions.comparison_plots_different_E_aa(std_devs_for_different_E_aa, E_aa_values, E_ad_values, x_label_E_aa_comparison, ylabel_stdv,title_stdv,legend,subfolder_path,saving_name_stdv)
                 
                 
                if plot_ratio_mean_var:
                    ratio_mean_stdev = [[mean / stdv for mean, stdv in zip(res_times, std_devs)]
                                        for res_times, std_devs in zip(mean_residence_times_for_different_E_aa, std_devs_for_different_E_aa)]
                    
                    ylabel_ratio = 'Mean/stdv devation'
                    title_ratio = 'Ratio Mean and Stdev of residence times vs. E ad values'
                    saving_name_ratio = f'nA_{nA}_n_{N}_ratio_mean_stdev_residence_times_E_ba_{E_ba}.png'
                    
                    general_functions.comparison_plots_different_E_aa(ratio_mean_stdev, E_aa_values, E_ad_values, x_label_E_aa_comparison, ylabel_ratio,title_ratio,legend,subfolder_path,saving_name_ratio)
            
# import matplotlib.pyplot as plt
# from collections import defaultdict
# # Group results by parameters (E_ba, E_aa, E_ad) while keeping track of L

# main_folder = 'Simulations_proteins_A_B'
# subfolder =  f'alfa_{alfa}'
# full_folder_path = os.path.join(folder_name, subfolder)



# grouped_results = defaultdict(list)
# for key, value in results_dict_B.items():
#     # Extract parameters
#     params = key.split(", ")
#     L = int(params[0].split(" = ")[1])  # Extract L value
#     E_ba = params[1]
#     E_aa = params[2]
#     E_ad = params[3]
    
#     # Group data
#     parameter_combination = f"{E_ba}, {E_aa}, {E_ad}"
#     grouped_results[parameter_combination].append((L, value["Mean number of consecutives B"], value["Mean width of the gaps between the Bs"]))


# # Plot for each combination
# for params, data in grouped_results.items():
#     # Sort data by L for consistent plotting
#     data.sort(key=lambda x: x[0])  # Sort by L
#     L_values = [item[0] for item in data]
#     ncons_values = [item[1] for item in data]
#     ngaps_values = [item[2] for item in data]
    
#     # Create a plot for ncons vs L
#     plt.figure(figsize=(8, 6))
#     plt.plot(L_values, ncons_values, label="Mean Consecutive", marker="o")
#     plt.plot(L_values, ngaps_values, label="Mean Gaps", marker="s", linestyle="--")
#     plt.title(f"Mean Consecutive and Gaps vs L\nParameters: {params}")
#     plt.xlabel("L")
#     plt.ylabel("Mean Values")
#     plt.legend()
#     plt.grid(True)
#     plt.tight_layout()
    
#     # Save the plot in the specified folder
#     safe_params = params.replace(", ", "_").replace(" = ", "_")  # Make filename safe
#     plot_path = os.path.join(full_folder_path, f"width_gaps_vsL_k_{k}_{safe_params}.png")
#     plt.savefig(plot_path)
#     plt.show()
#     plt.close()  # Close the plot to avoid memory issues
