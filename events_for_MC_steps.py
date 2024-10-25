# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 13:45:42 2024

@author: virgi
"""
import numpy as np



def add_A (list_DNA, empty_random_site, list_A , random_A, residence_times, list_empty_DNA, times_variables, time_step):
    if list_A[0,random_A] == (-1): #A is free
        list_DNA [0,empty_random_site] = 1  #always add because it lowers the energy 
        list_A [0,random_A] = empty_random_site 
        residence_times[0,random_A] = time_step
        list_empty_DNA.remove (empty_random_site)
        times_variables[random_A]['Count binding events'] += 1
    return list_DNA,list_A, residence_times, list_empty_DNA, times_variables
    
def remove_A (list_DNA, list_empty_DNA, empty_random_site, list_A , random_A, residence_times, times_variables, E_ad, E_aa, time_step):
    if list_A [0,random_A] != (-1):  #if A is not free 
        if list_DNA [0,list_A [0,random_A]] == 1 : #only A is bound to the DNA

            random_binding = np.random.random()  #draw a random number between 0 and 1
            energy = energy_function (list_A [random_A], list_DNA, E_ad, E_aa)
            
            if random_binding < 1/np.exp(energy):
                list_DNA, list_A, random_A,list_empty_DNA,residence_times, times_variables = remove_A_from_DNA_site(list_DNA, list_A, random_A,list_empty_DNA,residence_times,time_step, times_variables)
                
        if list_DNA [0, list_A [0,random_A]] == 2 : # B protein is also bound 
            #Energy condition 
            list_DNA, list_A, random_A,list_empty_DNA,residence_times, times_variables = remove_A_from_DNA_site (list_DNA, list_A, random_A,list_empty_DNA,residence_times,time_step, times_variables)
            #IF A REMOVE AN A FROM A SITE DO I ALSO REMOCE THE B ? 
    return list_DNA, list_A, random_A,list_empty_DNA,residence_times 


def remove_A_from_DNA_site (list_DNA, list_A, random_A,list_empty_DNA,residence_times,time_step, times_variables): 
    list_DNA [list_A [random_A]] = 0 #A is removed 
    list_empty_DNA.append(list_A [random_A])
    list_A [random_A] = -1 #A becomes unbound state
    time_binding = residence_times[random_A]
    residence_times[random_A] = time_step-time_binding
    times_variables[random_A]['Residence times'].append(residence_times[random_A] )
    residence_times[random_A]  = 0 
    return list_DNA, list_A, random_A,list_empty_DNA,residence_times, times_variables
    
def add_B_to_DNA_site (list_DNA, list_B, sites_possible_for_binding, random_binding_site, random_B ):
    
    ########THIS IS ASSUMING B BINDS INDEPENTLY FROM ENERGY 
    sites_with_bound_A = np.where (sites_possible_for_binding == 1)
    random_site_with_A_to_bind = sites_with_bound_A[np.random.randint(0, len (sites_with_bound_A))]
    list_DNA[0,random_site_with_A_to_bind] = 2
    list_B[random_B, random_binding_site] = random_site_with_A_to_bind
    #changes da fare con A ???
    return list_DNA, list_B
    
            
def add_B_event (list_B,random_B, list_DNA, L):
    if list_B[random_B,:].all() == -1: #check if B is free, i.e. all the binding sites are empty (representing by -1 in the original matrix)
        random_binding_site = np.random.randint(0, list_B.shape[1]) #choose randomly one of the sites
        list_DNA, list_B = add_B_to_DNA_site (list_DNA, list_B, list_DNA, random_binding_site, random_B )
        
    else: #there is at least one Binding site in B protein that is already bound 
        bound_B_sites = np.where (list_B[random_B, :] != -1) #list of all the DNA sites where a binding site of the B protein is already bound 
        random_B_bound_site = bound_B_sites[np.random.randint(0, len (bound_B_sites))] #pick randomly one of the bound sites 
        range_sites_within_L = list_DNA [0, random_B_bound_site-L: random_B_bound_site+L+1] #look within a maximum distance L for an A bound 
        list_DNA, list_B = add_B_to_DNA_site (list_DNA, list_B, range_sites_within_L, random_B_bound_site, random_B )
        ##Energy condition
                                                              
    return list_DNA, list_B
        
def remove_B_event (list_B,random_B, list_DNA):
 bound_B_sites =np.where (list_B[random_B, :] != -1) #get the binding sites that are bound to the DNA
 index_of_B_binding_site = np.random.randint(0, len (bound_B_sites))
 random_B_bound_site = bound_B_sites[index_of_B_binding_site] #pick randomly one of the bound sites 
 #energy condition
 #if succeed then remove 
 list_DNA, list_B=remove_B_from_DNA_site(random_B_bound_site, list_DNA, list_B, index_of_B_binding_site, random_B)
 return list_DNA, list_B
    
def remove_B_from_DNA_site (random_B_bound_site, list_DNA, list_B, index_of_B_binding_site, random_B):
    list_DNA[0,random_B_bound_site] = 1 #removing B from DNA site 
    list_B[random_B,index_of_B_binding_site] = -1 #setting B as free in the list of Bs
    return list_DNA, list_B
    
def energy_function (index,list_DNA, E_ad, E_aa): 
    
    
    if index == 0: #first site
        energy = list_DNA[0, index + 1] * E_aa + list_DNA [0, index]*E_ad
        
    elif index == (list_DNA.shape[1]-1) : #last site 
        energy = list_DNA [0, index-1]*E_aa + list_DNA [0, index]*E_ad
    
    else: 
        energy = list_DNA [0, index + 1] * E_aa + list_DNA [0, index-1]*E_aa + list_DNA [0, index]*E_ad
    return energy 