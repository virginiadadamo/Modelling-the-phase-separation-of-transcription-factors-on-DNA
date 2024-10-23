# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 21:09:15 2024

@author: virgi
"""

import numpy as np
import functions_MC_simulation_protein_A
import functions_MC_simulation_both 


'''
ASSUMPTIONS:
    -For the residence times: Difference between unbinding and binding event of the transcription factor. If a TF binds but never
    unbinds the residence time will not be computed. These trasncription factors will be ignored when computing the mean residence times over all transcription factors 
    -If an A is free and there is an empty site of DNA, protein A will always find that site 

'''

###PARAMETERS###

alfa = 0.15 #ratio between nA/N 
N = 10000#10000 #total number of binding sites in the DNA
nA = int (N*alfa) #number of As


nB = 0 #number of Bs
k = 0 #number of B interacting sites with A  

#Time parameters
stop_time = 2000000
ignoring_steps = 10000
m = 50

if (stop_time - ignoring_steps) % m != 0 :
    raise ValueError(f"Error: m ({m}) is not a divisor of {stop_time - ignoring_steps}")

#Energy parameters 

E_ad_values = np.arange(0, 4, 1)
E_aa_values = [0, 2.5]

###PLOTS' TAGS ### - select the plots by putting the corresponding value to true 

#For each combination of Eaa and Ead
histogram_mean_residence_time = True # To plot the corresponding distribution of mean residence times of the transcription factors
scatter_plot_std = True #To plot the standard deviation of the residence times of each transcription factors 
scatter_plot_mean = True #To plot the mean of the residence times of each transcription factors
histogram_binding_events = True #To plot the distribution of binding events of the transcription factors 
histogram_cluster_size = True #To plot the distribution of cluster sizes 

#For each E_aa
plot_nA_bound = True #Plot the number of transcription factors bound in time
plot_inverse = True #Plot the inverse of the number of transcription factors vs the inverse of the corresponding time steps 
histogram_binding_events_E_aa = True #Plot the distribution of binding events of the transcription factors for each E_aa
plot_mean_cluster_size_max_cluster_size = False #Plot Mean Cluster Size vs Max Cluster size
plot_mean_cluster_size_E_ad = False #Mean Cluster Size vs. E_ad 
plot_mean_cluster_size_max_cluster_size = False #Mean Cluster Size vs Max Cluster size

#Comparison between different E_aa
plot_log_mean_residence_time = True #Plot for different E_aa the ln of Mean Residence Time vs E_ad values 
plot_stdv_residence_time = True #Plot for different E_aa the Standard Deviation vs E_ad values 
plot_ratio_mean_stdv = True  #Plot for different E_aa the ratio of Mean Residence Time and Standard Deviation vs E_ad values 


###CREATING FOLDERS###
folder_name = 'Simulations_protein_A'
subfolder_path = functions_MC_simulation_both.create_folders(folder_name, alfa)
functions_MC_simulation_both.create_txt_parameters(subfolder_path, alfa, stop_time, ignoring_steps)

legend = f'stop_time={stop_time}\nignoring_steps={ignoring_steps}\nm={m}\nnA={nA}\nn={N}'


###MONTE CARLO SIMULATION###
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
        
        
        ### PARAMETERS FOR THE INITIAL STATE###
        #Time parameters 
        time_step = 0 #initial time step
        time_step_sampled = [] #To store the time steps that are being sampled
        rate_counter = 0 #To count the steps after the first sampled one
        residence_times = np.zeros((1,nA)) # To store for each transcription the residence times and binding events. It is initialised at 0, if a binding event occurs the time of binding will be store. If an unbinding event then happens the time stored will be used to compute the Residence time and the value will be put to 0 again
        times_variables = [{'Index of Transcription Factor': i, 
                        'Residence times': [],
                        'Mean residence time' : 0 ,
                        'Standard deviation':0,
                        'Max residence time' : 0,
                        'Count binding events': 0} for i in range (nA)] #Dictionary to store for each Transcription factors the time variables indicated 
        
        index_tfs = [] # To store the indexes of the Transcription factor with a least one Residence Time
        stdv_each_tf = [] # To store for each transcription factor the Standard Deviation of Residence Times 
        mean_each_tf = [] # To store for each transcription factor the Mean Residence time 
        count_binding_events_list = [] #To store for each transcription factor the number of Binding Events
        mean_residence_time = 0 #To compute the mean of the residence times over all transcription factors
        
        #DNA parameters 
        list_DNA = np.zeros((1,N)) #Array representing the DNA sites: 0 corresponds to an empty site, 1 to a site with an transcription factor bound to it, 2 if a B is bound - at the initial state there is no bound A nor B to the DNA 
        list_empty_DNA  = list(range (0,list_DNA.shape[1],1)) #list containing all the indexes of empty sites in the DNA
        nA_bound_DNA = 0 # number of Transcription factors that are bound to the DNA
        
        group_sizes_snapshots = [] #To store at each time step the corresponding cluster sizes 
        average_cluster_sizes = [] #To store the mean of the cluster sizes at each time step
        max_cluster_sizes = [] #To store at each time step the value of the highest value of cluster size
        all_group_sizes_histogram = [] #union of all the cluster sizes list in one list, that will be used to compute the final histogram of distribution of cluster sizes 
        
        clusters_each_time_sampled = [] #List containing the cluster objects at each time sampled 
        
        #Transcription factors parameters 
        list_A = np.full(nA, -1).reshape(1, nA) #Array representing the transcription factors: -1 for unbound, will store the position of the site on the DNA when bound 
        nA_bound_list_A = 0 #As that are bound in the list of A- will be used to check that the number of transcription factors bound in the list of transcription factors is the same as the number of bound transcription factors in the DNA
        nA_bound_snapshots = [] # To store the number of A bound for each time step
        
        if nB > 0 :
            list_B = np.array ((nB,k))
            #do SIMULATION FOR PROTEIN B 
            
        
        
        
        while time_step < stop_time:
                        
            time_step, list_DNA, list_A, list_empty_DNA, times_variables = functions_MC_simulation_protein_A.step_MC(time_step, list_DNA, list_A, list_empty_DNA, E_ad, E_aa, residence_times, times_variables)
            
            nA_bound_list_A= functions_MC_simulation_both.count_A(list_A)
            nA_bound_DNA = functions_MC_simulation_both.count_consecutive_ones(list_DNA)
            
            if nA_bound_DNA != nA_bound_list_A:
                raise ValueError(f"Error: nA_bound_DNA ({nA_bound_DNA}) and nA_bound ({nA_bound_list_A}) are not equal")
            
            if time_step >= ignoring_steps: 
                
                rate_counter = rate_counter + 1
                
                if rate_counter == m: #sampling occurs
                
                    group_sizes, max_count, nA_bound_snapshots, group_sizes_snapshots, average_cluster_sizes, max_cluster_sizes, rate_counter, all_group_sizes_histogram, clusters_each_time_sampled = functions_MC_simulation_both.take_sample(list_DNA, list_A, nA_bound_snapshots, group_sizes_snapshots, average_cluster_sizes,max_cluster_sizes, rate_counter, all_group_sizes_histogram, clusters_each_time_sampled)
                    time_step_sampled.append(time_step)
                    
                    
        # count_never_unbind = np.count_nonzero(residence_times)
        # times_never_unbind =  [x for x in residence_times if x != 0]
        
         
        for i in range (nA):
            
            count_binding_events_list.append (times_variables[i]['Count binding events'])
            
            if len (times_variables[i]['Residence times']) != 0: #if there is at leat one value of Residence time
                
                index_tfs.append (times_variables[i]['Index of Transcription Factor'])
                
                times_variables[i]['Mean residence time'] = np.mean (times_variables[i]['Residence times'])
                mean_each_tf.append(times_variables[i]['Mean residence time'])
                
                
                times_variables[i]['Standard deviation'] = np.std (times_variables[i]['Residence times'])
                stdv_each_tf.append (times_variables[i]['Standard deviation'])
                
                times_variables[i]['Max residence time'] = np.max (times_variables[i]['Residence times'])
                
                
        
    
        mean_residence_times.append(np.mean(mean_each_tf)) 
        std_devs.append(np.mean(stdv_each_tf))
        binding_events_lists.append (count_binding_events_list) 
        
        max_residence_time = np.max([times_variables[i]['Max residence time'] for i in range(nA)])
        max_residence_times.append (max_residence_time)
        
        nA_bound_for_different_energies.append(nA_bound_snapshots) 
        mean_cluster_sizes.append(np.mean(average_cluster_sizes)) #taking the mean of cluster sizes for each energy value
        mean_max_cluster_sizes.append(np.mean(max_cluster_sizes)) #taking the mean of maximum sized cluster for each energy value 
        
        ### PLOTS FOR EACH COMBINATION OF E_AA AND E_AD ###
        
        
        if histogram_mean_residence_time:
            histogram_title_mean = f'Histogram of distribution of mean residence times (E_aa={E_aa}, E_ad={E_ad})'
            saving_histogram_name = f'nA_{nA}_n_{N}_histo_mean_resident_times_Eaa_{E_aa}_Ead_{E_ad}.png'
            x_label_histogram_mean = 'Mean Residence Time'
            
            functions_MC_simulation_both.plot_histogram(mean_each_tf, histogram_title_mean, legend, subfolder_path, x_label_histogram_mean, 'Frequency', saving_histogram_name, time_step_sampled, False, 100)
        
        unique_idx = np.unique(index_tfs)
        
        if len(unique_idx) != len(index_tfs):
            raise ValueError("Error: the index of TFS are not unique")
            
        if scatter_plot_std:
            xlabel_scatter = 'Index of TF'
            ylabel_scatter_stdv = 'Standard Deviation of Residence Times'
            title_scatter_stdv = f'Standard deviation of Residence times for different TFs (E_aa={E_aa}, E_ad={E_ad})'
            saving_name_scatter_stdv = f'nA_{nA}_n_{N}_stdv_res_times_tf__Eaa_{E_aa}_Ead_{E_ad}.png' 
            
            functions_MC_simulation_both.scatter_plot (index_tfs,stdv_each_tf, legend, subfolder_path, xlabel_scatter, ylabel_scatter_stdv, title_scatter_stdv, saving_name_scatter_stdv)
        
        if scatter_plot_mean:
            xlabel_scatter = 'Index of TF'
            ylabel_scatter_mean = 'Mean of Residence Times'
            title_scatter_mean = f'Mean of Residence times for different TFs (E_aa={E_aa}, E_ad={E_ad})'
            saving_name_scatter_mean = f'nA_{nA}_n_{N}_mean_res_times_tf__Eaa_{E_aa}_Ead_{E_ad}.png'
            
            functions_MC_simulation_both.scatter_plot (index_tfs,mean_each_tf, legend, subfolder_path, xlabel_scatter, ylabel_scatter_mean, title_scatter_mean, saving_name_scatter_mean)

        if histogram_binding_events:
            histogram_title_be = f'Histogram of distribution of binding events (Eaa {E_aa}, Ead {E_ad})'
            saving_histogram_name_be = f'nA_{nA}_n_{N}_histo_binding_events_Eaa_{E_aa}_Ead_{E_ad}.png'
            x_label_histogram_be = 'Binding Events'
            
            functions_MC_simulation_both.plot_histogram(count_binding_events_list, histogram_title_be, legend, subfolder_path,x_label_histogram_be,'Frequency', saving_histogram_name_be, time_step_sampled )
        
        if histogram_cluster_size:
            histogram_title_cluster = f'Cluster histogram (Eaa {E_aa}, Ead {E_ad})'
            saving_histogram_name_cluster = f'nA_{nA}_n_{N}_cluster_histo_Eaa_{E_aa}_Ead_{E_ad}.png'
            x_label_histogram_cluster = 'Cluster size'
            
            functions_MC_simulation_both.plot_histogram(all_group_sizes_histogram, histogram_title_cluster, legend, subfolder_path,x_label_histogram_cluster, 'Frequency', saving_histogram_name_cluster, time_step_sampled, True )
                    
      
    std_devs_for_different_E_aa.append(std_devs)       
    mean_residence_times_for_different_E_aa.append (mean_residence_times)
    
    ### PLOTS FOR EACH E_AA 
    if plot_nA_bound:
        xlabel_nA_bound_snapshots = 'Time Steps'
        ylabel_nA_bound_snapshots = 'Number of A bound to DNA'
        title_nA_bound_snapshots = f'Number of A bound to DNA vs. Time Steps (E_aa={E_aa})'
        saving_name_nA_bound_snapshots = f'nA_{nA}_n_{N}_bound_in_time_Eaa_{E_aa}_Ead_{E_ad}.png'
        
        functions_MC_simulation_both.plot_different_Ead_in_time (nA_bound_for_different_energies, E_ad_values, time_step_sampled, xlabel_nA_bound_snapshots, ylabel_nA_bound_snapshots,title_nA_bound_snapshots,legend,subfolder_path,saving_name_nA_bound_snapshots)
        
    if plot_inverse:
        xlabel_lineweaver_burk = 'Inverse Time Steps'
        ylabel_lineweaver_burk = 'Inverse Number of A bound to DNA'
        title_lineweaver_burk = f'Lineweaver–Burk plot (E_aa={E_aa})'
        saving_name_lineweaver_burk = f'nA_{nA}_n_{N}_Lineweaver_Burk_Eaa_{E_aa}_Ead_{E_ad}.png'
        
        functions_MC_simulation_both.plot_different_Ead_in_time (nA_bound_for_different_energies, E_ad_values, time_step_sampled, xlabel_lineweaver_burk, ylabel_lineweaver_burk,title_lineweaver_burk,legend,subfolder_path,saving_name_lineweaver_burk, True)
        
    
    if histogram_binding_events_E_aa:
        xlabel_binding_events = 'Binding Events'
        ylabel_binding_events = 'Frequency'
        title_binding_events = f'Histogram of Distribution of Binding Events for different E_ad (E_aa={E_aa})'
        saving_name_binding_events = f'nA_{nA}_n_{N}_combined_histo_binding_events.png'
      
        functions_MC_simulation_both.plot_histogram_distribution_different_E_ad (binding_events_lists,E_ad_values, E_aa, xlabel_binding_events, ylabel_binding_events,title_binding_events,legend,subfolder_path,saving_name_binding_events)
        
    
    if plot_mean_cluster_size_max_cluster_size:
        functions_MC_simulation_both.plot_figure (max_residence_times,mean_cluster_sizes,'Max residence times','Mean Cluster Size','Max Residence Time vs Max Cluster size (E_aa={E_aa})',subfolder_path,f'nA_{nA}_n_{N}_max_rt_mean_cluster_size.png') 
    if plot_mean_cluster_size_E_ad:
        functions_MC_simulation_both.plot_figure (E_ad_values,mean_cluster_sizes,'E_ad','Mean Cluster Size',f'Mean Cluster Size vs. E_ad (E_aa={E_aa})',subfolder_path,f'nA_{nA}_n_{N}_E_ad_mean_cluster_size.png') 
   
###COMPARISON BETWEEN DIFFERENT E_AA ###
x_label_E_aa_comparison = 'E ad values'
    
if plot_log_mean_residence_time:
    log_list = [np.log(element) for element in mean_residence_times_for_different_E_aa]
    ylabel_log = 'Log (Mean Residence Time)'
    title_log = 'Log Mean Residence Time vs. E ad values'
    saving_name_log = f'nA_{nA}_n_{N}_log_mean_residence_time.png'
    
    functions_MC_simulation_both.comparison_plots_different_E_aa(log_list, E_aa_values, E_ad_values, x_label_E_aa_comparison, ylabel_log,title_log,legend,subfolder_path,saving_name_log)

if plot_stdv_residence_time:
    ylabel_stdv = 'Stdv residence times'
    title_stdv = 'Mean standard deviation of residence times vs. E ad values'
    saving_name_stdv = f'nA_{nA}_n_{N}_mean_stdv_residence_times.png'
    
    functions_MC_simulation_both.comparison_plots_different_E_aa(std_devs_for_different_E_aa, E_aa_values, E_ad_values, x_label_E_aa_comparison, ylabel_stdv,title_stdv,legend,subfolder_path,saving_name_stdv)
    
if plot_ratio_mean_stdv:
    ratio_mean_stdv = [[mean / stdv for mean, stdv in zip(res_times, std_devs)]
                        for res_times, std_devs in zip(mean_residence_times_for_different_E_aa, std_devs_for_different_E_aa)]
    
    ylabel_ratio = 'Mean/standard devation'
    title_ratio = 'Ratio Mean and standard devation of residence times vs. E ad values'
    saving_name_ratio = f'nA_{nA}_n_{N}_ratio_mean_stdv_residence_times.png'
    
    functions_MC_simulation_both.comparison_plots_different_E_aa(ratio_mean_stdv, E_aa_values, E_ad_values, x_label_E_aa_comparison, ylabel_ratio,title_ratio,legend,subfolder_path,saving_name_ratio)


