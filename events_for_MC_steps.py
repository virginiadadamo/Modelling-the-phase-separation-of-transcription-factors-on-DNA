# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 13:45:42 2024

@author: virgi
"""
import numpy as np
import random



def add_A (list_DNA, empty_random_site, list_A , random_A, residence_times, list_empty_DNA, times_variables, time_step):
    
    if list_A[random_A,0] == (-1): #A is free
        #print (f'Step: adding of a free A : (time step={time_step}, A selected is ={random_A})')
        list_DNA [0,empty_random_site] = 1  #always add because it lowers the energy 
        list_A [random_A,0] = empty_random_site #change the value in list A, instead of (-1) to the site where A is bound 
        #print (f'The site to which the A was added is {empty_random_site}, the DNA is  ={list_DNA})')
        residence_times[0,random_A] = time_step
        list_empty_DNA.remove (empty_random_site)
        times_variables[random_A]['Count binding events'] += 1
        
    return list_DNA,list_A, residence_times, list_empty_DNA, times_variables
    
def remove_A (list_DNA, list_empty_DNA, list_A , random_A, list_B, residence_times, times_variables, E_ad, E_aa, E_ab, time_step):
    
    random_binding = np.random.random()  #draw a random number between 0 and 1
    A = list_A [random_A,:] 
    if A[0] != (-1):  #if A is not free 
    
        energy = energy_function_unbinding_A (list_DNA, A, list_B, E_ad, E_aa, E_ab)#compute unbinding energy
        
        if random_binding < 1/np.exp(energy): #if energy condition then it will remove
        
            list_DNA, list_A,list_empty_DNA,residence_times, times_variables = remove_A_from_DNA_site(list_DNA, list_A, random_A,list_empty_DNA,residence_times,time_step, times_variables)
            #print (f'Removal of A succceed: list_DNA {list_DNA}, list_A  ={list_A}, list_empty_DNA={list_empty_DNA} ')
            
            if A[1] != (-1): #Also a B protein is bound, it must be freed
                #print ('Step: remove A but also B is bound')
                B_bound = A[1]
                B_bound_site = np.where(list_B[B_bound,:] == random_A)[0] #find the site 
                #print('List a and B before:', list_A, list_B)
                list_B[B_bound,B_bound_site] = -1 #setting B as free in the list of Bs
                list_A [random_A,1] = -1 #mark as free the corresponding site on A 
                #print('List a and B after:', list_A, list_B)
    return list_DNA, list_A, random_A,list_empty_DNA,residence_times 


def remove_A_from_DNA_site (list_DNA, list_A, random_A,list_empty_DNA,residence_times,time_step, times_variables):  #CHANG ETHE NAME (RESET)
    
    #print (f'Before starting removing the list of A is {list_A}, random A is {random_A}')
    list_DNA [0,list_A [random_A, 0]] = 0 #A is removed 
    list_empty_DNA.append(list_A [random_A,0]) #the site is added to the list of empty sites 
    list_A [random_A,0] = -1 #A becomes unbound state
    
    #Now the A unbinds we compute the residence time
    time_binding = residence_times[0,random_A]
    residence_times[0,random_A] = time_step-time_binding
    times_variables[random_A]['Residence times'].append(residence_times[0,random_A] )
    residence_times[0,random_A]  = 0 #time put back to zero, to compute the next residence time when A will bind again
    
    return list_DNA, list_A,list_empty_DNA,residence_times, times_variables
    
            
def add_B_event ( list_DNA, list_A, list_B, random_B, L):
    
    #print ('Adding B event selected', random_B)
    
    if (list_B[random_B, :] == -1).all():  #check if B is free, i.e. all the binding sites are empty (representing by -1 in the original matrix)
        #print ('All random sites !')
        random_binding_site = np.random.randint(0, list_B.shape[1]) #choose randomly one of the sites
        #print ('Random binding site', random_binding_site)
        #print ('List A and List B before having added B, ', list_A, list_B)
        starting_site = 0
        ending_site = len(list_DNA[0]) 
        list_A, list_B = add_B_to_DNA_site (list_A, list_B, list_DNA, starting_site, ending_site, random_binding_site, random_B ) #add B to one random site within all DNA
        #print ('List A and List B after having added B, ', list_A, list_B)
        
    else: #there is at least one Binding site in B protein that is already bound 
       
        bound_B_sites = np.where (list_B[random_B, :] != -1) #list of all the DNA sites where a binding site of the B protein is already bound 
        free_B_sites = np.where (list_B[random_B, :] == -1) #liste of all free B sites
        random_B_bound_site = random.choice(bound_B_sites[0]) #pick randomly one of the bound sites 
        random_B_free_site = random.choice(free_B_sites[0])  #pick randomly one of the free sites 
        DNA_index_where_B_bound = list_A[list_B[random_B,random_B_bound_site], 0]
      
        if DNA_index_where_B_bound-L > 0:
            starting_site = int (DNA_index_where_B_bound-L)
        else:
            starting_site = 0
        if DNA_index_where_B_bound+L > len(list_DNA[0]):
            ending_site = len(list_DNA[0])
        else: 
            ending_site = int(DNA_index_where_B_bound+L)
            
        #print ('Dna index where the B is bound', DNA_index_where_B_bound)
        
        #print ('List A and List B before having added B, ', list_A, list_B)
        list_A, list_B = add_B_to_DNA_site (list_A, list_B, list_DNA, starting_site, ending_site, random_B_free_site, random_B )
        #print ('List A and List B after having added B, ', list_A, list_B)
        
                                                              
    return list_DNA, list_A, list_B

def add_B_to_DNA_site (list_A, list_B, list_DNA, starting_site, ending_site, random_binding_site, random_B ):
    
    #print ('Starting site', starting_site)
    #print ('Ending site', ending_site)
    
    #The problem if the random A selcted is already bound to B 
    range_sites_within_L = list_DNA [0, starting_site:ending_site+1] #look within a maximum distance L for an A bound 
    #print('Range sites within a distance L', range_sites_within_L)
    sites_with_bound_A = list (np.where (range_sites_within_L == 1)[0]) #find the sites in the specified distance
    #print('sites_with_bound_A', sites_with_bound_A)
    if starting_site != 0:
    # Convert relative indices to original indices in list_DNA
        sites_with_bound_A = [i + starting_site for i in sites_with_bound_A]
    #print('sites_with_bound_A sfter conversion', sites_with_bound_A)    
    
    for site in sites_with_bound_A[:]:  # Use a slice to create a copy of the list
     corresponding_A_to_site = np.where(list_A[:, 0] == site)[0]
     #print('Site', site, 'corresponding_A_to_site', corresponding_A_to_site)
     if list_A[corresponding_A_to_site, 1] != -1:  # If that A is already bound to a B
        #print('Before removing the site', sites_with_bound_A)
        sites_with_bound_A.remove(site)  # Remove the site
        #print('After removing the site', sites_with_bound_A)

            
    if len(sites_with_bound_A) > 0 : #if there is at leas ìt one bound A in the given lenght 
        #print ('Sites with bound A', sites_with_bound_A )
        random_site_with_A_to_bind = random.choice (sites_with_bound_A)#[np.random.randint(0, len (sites_with_bound_A))]#choose randomly one of these sites. We assume B always bind independently of the energy
        #print ('Random site with bound A', random_site_with_A_to_bind)
        del sites_with_bound_A #clearing memory after usage
        #Get the A to which that occupied site correspond to 
        A_to_bind = np.where(list_A[:, 0] == random_site_with_A_to_bind)[0]
        
        #print('A to bind:', A_to_bind)
        
        list_B[random_B, random_binding_site] = A_to_bind #mark B as occupied by replacing -1 with the A to which it is bound
        list_A[A_to_bind, 1] = random_B #mark also that the A is bound to B 
    
    return list_A, list_B
    
        
def remove_B_event (list_A, list_B,random_B, Eba):

 #print ('Removing B event')
 bound_B_sites =np.where (list_B[random_B, :] != -1)[0] #get the binding sites that are bound to the DNA
 #print ('Bound B sites', bound_B_sites)
 # randomly choose one index within the bound sites 
 print (bound_B_sites)
 random_B_bound_site = random.choice(bound_B_sites) #pick randomly one of the bound sites 
 #print ('Random B bound site', random_B_bound_site )
 random_binding = np.random.random()  #draw a random number between 0 and 1
 energy = energy_unbind_function_B (list_B[random_B, :], Eba)#compute unbinding energy
 #print (energy)
 
 if random_binding < 1/np.exp(energy): #if energy condition succeed then it will remove
     #print ('Removal succeeded')
     #print ('List A and B before removing', list_A, list_B)
     list_A, list_B=remove_B_from_DNA_site(random_B_bound_site, list_A, list_B, random_B)
     #print ('List A and B after removing', list_A, list_B)
     
 return list_A, list_B
    
def remove_B_from_DNA_site (random_B_bound_site, list_A, list_B, random_B):

    index_bound_A = list_B[random_B,random_B_bound_site]
    #print (index_bound_A)
    list_A[index_bound_A, 1] = -1 #marking that the B is no longer bound to the A
    list_B[random_B,random_B_bound_site] = -1 #setting B as free in the list of Bs
    
    return list_A, list_B

def energy_unbind_function_B (B, Eba): 
    
    #computing the number of sites in B protein that are bound
    energy = compute_energy_B_binding (B, Eba)
    return energy

def energy_function_unbinding_A (list_DNA, A, list_B, E_ad, E_aa, E_ab): 
    DNA_site_index = A[0] #find the index to which A is bound 
    if A[1] != (-1 ): #There is a B bound to that A
            B = list_B[A[1],:] #all binding sites for that B 
            B_energy = compute_energy_B_binding (B, E_ab) #compute the energy of B protein bound to A, depending on the number of B sites 
    else: 
            B_energy = 0 #if there are no B bound to A the corresponding energy will be zero 
            
    A_energy = compute_energy_A_binding (DNA_site_index, list_DNA, E_ad, E_aa)
    return A_energy + B_energy


def compute_energy_A_binding (index, list_DNA, E_ad, E_aa):
    '''
    Functions to compute the energy of a As protein on the DNA
    The energy is equal to E_ad for the site to which the TF binds, E_aa to account for the energy contribution of the neighbours 
    '''
    if index == 0: #first site
        energy = list_DNA[0, index + 1] * E_aa + list_DNA [0, index]*E_ad
        
    elif index == (list_DNA.shape[1]-1) : #last site 
        energy = list_DNA [0, index-1]*E_aa + list_DNA [0, index]*E_ad
    
    else: 
        energy = list_DNA [0, index + 1] * E_aa + list_DNA [0, index-1]*E_aa + list_DNA [0, index]*E_ad
    return energy 
    

def compute_energy_B_binding (B, energy):
    '''
    Functions to compute the energy of a B protein
    Assumption: this energy is directly proportional on the number of B proteins that are bound
    '''
    bound_B_sites =np.where (B != -1)[0]
    B_energy = energy * len(bound_B_sites)
    return B_energy 
    