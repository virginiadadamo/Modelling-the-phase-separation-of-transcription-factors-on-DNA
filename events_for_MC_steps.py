# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 13:45:42 2024

@author: virgi
"""
import numpy as np
import random



def add_A (list_DNA, empty_random_site, list_A , random_A, residence_times, list_empty_DNA, times_variables, time_step):
    
    """
    Simulates the binding of A to an empty site along the DNA.

    This function handles the process of binding an entity from `list_A` to an available site 
    in `list_DNA`. If the selected entity (A) is free, it binds to the specified empty site, 
    updates the relevant data structures, and tracks the binding event.

    Parameters:
    ----------
    list_DNA : numpy.ndarray (1xN)
        Array that simulate the DNA where a value of `1` indicates a bound site and `0` 
        indicates a free site.
    empty_random_site : int
        The index of the empty site in `list_DNA` to which A will bind.
    list_A : numpy.ndarray (nA x 2)
        An array representing As. A value of `-1` in the first column indicates a free entity, 
        while a positive integer indicates the index of the site to which the entity is bound.
    random_A : int
        The index of the randomly selected A from `list_A` attempting to bind.
    residence_times : numpy.ndarray
        A 1D array tracking the time steps at which A was last bound.
    list_empty_DNA : list
        A list of indices representing currently empty sites in `list_DNA`.
    times_variables : dict
        A dictionary tracking various event counts and times for each A, such as 
        the number of binding events.
    time_step : int
        The current simulation time step.

    Updates:
    -------
    - `list_DNA`: Marks the specified site as bound (`1`).
    - `list_A`: Updates the binding state of the selected entity to the site index.
    - `residence_times`: Records the current time step for the binding event.
    - `list_empty_DNA`: Removes the bound site from the list of empty sites.
    - `times_variables`: Increments the binding event counter for the selected entity.

    Returns:
    -------
    tuple:
        - list_DNA : numpy.ndarray
            Updated representation of DNA-like sites.
        - list_A : numpy.ndarray
            Updated state of entities (A).
        - residence_times : numpy.ndarray
            Updated residence times array.
        - list_empty_DNA : list
            Updated list of empty sites.
        - times_variables : dict
            Updated event tracking dictionary.

    Notes:
    ------
    - Binding is always successful if the selected entity (A) is free (`-1`).
    - The function assumes `empty_random_site` and `random_A` are valid indices.
    """
    
    if list_A[random_A,0] == (-1): #A is free
        list_DNA [0,empty_random_site] = 1  #always add because it lowers the energy 
        list_A [random_A,0] = empty_random_site #change the value in list A, instead of (-1) to the site where A is bound 
        
        residence_times[0,random_A] = time_step
        list_empty_DNA.remove (empty_random_site)
        times_variables[random_A]['Count binding events'] += 1
        
        
    return list_DNA,list_A, residence_times, list_empty_DNA, times_variables
    
def remove_A (list_DNA, list_empty_DNA, list_A , random_A, list_B, residence_times, times_variables, E_ad, E_aa, time_step):
    
    """
    Simulates the removal of A from a site on the DNA.
    
    Parla delle energie

    Parameters:
    ----------
    list_DNA : numpy.ndarray (1xN)
        Array that simulate the DNA where a value of `1` indicates a bound site and `0` 
        indicates a free site. 
    list_empty_DNA : list
        A list of indices representing currently empty sites in `list_DNA`.
    
    list_A : numpy.ndarray (nA x 2)
        An array representing As. A value of `-1` in the first column indicates a free entity, 
        while a positive integer indicates the index of the site to which the entity is bound.
    random_A : int
        The index of the randomly selected A from `list_A` for removal.
    residence_times : numpy.ndarray
        A 1D array tracking the time steps at which A was last bound..
    times_variables : dict
        A dictionary tracking various event counts and times for each A, such as 
        the number of binding events.
    E_ad : int 
        Represent the binding energy of an A to the DNA 
    E_aa: int  
        Energy that reflects the influence of neighboring A molecules on the selected A being removed.
    k: int (default set to 0)
    time_step : int
        The current simulation time step.

    Updates:
    -------
    - `list_DNA`: Marks the specified site as bound (`1`).
    - `list_A`: Updates the binding state of the selected entity to the site index.
    - `residence_times`: Records the current time step for the binding event.
    - `list_empty_DNA`: Removes the bound site from the list of empty sites.
    - `times_variables`: Increments the binding event counter for the selected entity.

    Returns:
    -------
    tuple:
        - list_DNA : numpy.ndarray
            Updated representation of DNA-like sites.
        - list_A : numpy.ndarray
            Updated state of entities (A).
        - residence_times : numpy.ndarray
            Updated residence times array.
        - list_empty_DNA : list
            Updated list of empty sites.

    Notes:
    ------
    - 
    """

    random_binding = np.random.random()  #draw a random number between 0 and 1
    A = list_A [random_A,:] 
    if A[0] != (-1):  #if A is not free 
    
        energy = energy_function_unbinding_A (list_DNA,list_A, A, list_B, E_ad, E_aa)#compute unbinding energy
        
        if random_binding < 1/np.exp(energy): #if energy condition then it will remove
            
            
            list_DNA, list_A,list_empty_DNA,residence_times, times_variables = remove_A_from_DNA_site(list_DNA, list_A, random_A,list_empty_DNA,residence_times,time_step, times_variables)
            
            if A[1] != (-1): #Also a B protein is bound, it must be freed
                B_bound = A[1]
                B_bound_site = np.where(list_B[B_bound,:] == random_A)[0] #find the site 
                list_B[B_bound,B_bound_site] = -1 #setting B as free in the list of Bs
                list_A [random_A,1] = -1 #mark as free the corresponding site on A 
                
                
    return list_DNA, list_A, random_A,list_empty_DNA,residence_times 


def remove_A_from_DNA_site (list_DNA, list_A, random_A,list_empty_DNA,residence_times,time_step, times_variables):  

    """
    Handles the removal of an A from the DNA and updates associated data structures.

    This function removes a specified A from its bound site in `list_DNA`, marks the site 
    as empty, and computes the residence time of the entity. The residence time is recorded, and 
    the state of the entity is reset for future binding events.

    Parameters:
    ----------
    list_DNA : numpy.ndarray (1xN)
        Array that simulate the DNA where a value of `1` indicates a bound site and `0` 
        indicates a free site.
    list_A : numpy.ndarray (nA x 2)
        An array representing As. A value of `-1` in the first column indicates a free entity, 
        while a positive integer indicates the index of the site to which the entity is bound.
    random_A : int
        The index of the A being removed from its bound site.
    list_empty_DNA : list
        A list of indices representing currently empty sites in `list_DNA`.
    residence_times : numpy.ndarray
        A 1D array tracking the time steps at which each A was last bound.
    time_step : int
        The current simulation time step.
    times_variables : dict
        A dictionary tracking various event counts and times for each A, including residence times.

    Updates:
    -------
    - `list_DNA`: Marks the previously occupied site as empty (`0`).
    - `list_A`: Sets the state of the specified entity to unbound (`-1`).
    - `list_empty_DNA`: Adds the unbound site to the list of empty sites.
    - `residence_times`: Computes and records the residence time for the entity.
    - `times_variables`: Appends the residence time to the entity's history and resets its tracking.

    Returns:
    -------
    tuple:
        - list_DNA : numpy.ndarray
            Updated representation of DNA-like sites.
        - list_A : numpy.ndarray
            Updated state of As.
        - list_empty_DNA : list
            Updated list of empty sites.
        - residence_times : numpy.ndarray
            Updated residence times array.
        - times_variables : dict
            Updated event tracking dictionary.

    Notes:
    ------
    - Residence time is calculated as the difference between the current time step and the time 
      the entity was bound.
    - The function assumes `random_A` corresponds to a valid bound entity in `list_A`.
    """
    
    
    list_DNA [0,list_A [random_A, 0]] = 0 #A is removed 
    list_empty_DNA.append(list_A [random_A,0]) #the site is added to the list of empty sites 
    list_A [random_A,0] = -1 #A becomes unbound state
    
    #Now the A unbinds we compute the residence time
    time_binding = residence_times[0,random_A]
    residence_times[0,random_A] = time_step-time_binding
    
    times_variables[random_A]['Residence times'].append(residence_times[0,random_A] )
    #print ('Times variables',  times_variables )
    residence_times[0,random_A]  = 0 #time put back to zero, to compute the next residence time when A will bind again
    
    return list_DNA, list_A,list_empty_DNA,residence_times, times_variables
    
            
def add_B_event ( list_DNA, list_A, list_B, random_B, L, does_B_bind):

    """
    Handles the event of adding a B protein to the DNA, updating the relevant data structures.

    This function determines whether a B protein can bind to the DNA. If the selected B protein is free,
    it randomly chooses a binding site. If the B protein is partially bound, it selects a free site 
    within a defined region around an already bound site for binding.

    Parameters:
    ----------
    list_DNA : numpy.ndarray
        A representation of the DNA where binding events can occur. A value of `1` indicates a bound site, 
        and `0` indicates a free site.
    list_A : numpy.ndarray (nA x 2)
        A matrix tracking the state of A proteins. Contains the indices of the DNA sites where each A is bound.
    list_B : numpy.ndarray
        A matrix tracking the state of B proteins. Rows correspond to individual B proteins, and columns represent
        their multiple binding sites. A value of `-1` indicates a free binding site.
    random_B : int
        The index of the B protein selected for the binding event.
    L : int
        The range around a bound DNA site to consider for additional B bindings.
    does_B_bind : 
        

    Updates:
    -------
    - `list_A`: Updates the DNA binding states of A proteins if B proteins interact with sites already bound by A.
    - `list_B`: Updates the binding states of the selected B protein, marking the new site as bound.
    - `list_DNA`: Tracks the DNA sites occupied as a result of the binding event.
    - `does_B_bind`: Indicates whether the B protein successfully bound during the current event.

    Returns:
    -------
    tuple:
        - list_DNA : numpy.ndarray
            Updated DNA representation.
        - list_A : numpy.ndarray
            Updated state of A proteins.
        - list_B : numpy.ndarray
            Updated state of B proteins.
        - does_B_bind : bool
            Indicates if the binding event was successful.

    Logic:
    ------
    1. If the selected B protein is entirely free, a random site is chosen for binding.
    2. If the B protein has some binding sites already occupied:
        - A random bound site is selected, and its position determines the region (`L`) for a new binding.
        - A random free site is selected within the defined region for the binding event.
    3. The `add_B_to_DNA_site` function is called to execute the actual binding logic and update structures.

    Notes:
    ------
    - The function assumes that the dimensions of `list_DNA`, `list_A`, and `list_B` are consistent with
      the intended modeling.
    - Binding events consider the spatial constraints defined by `L`, limiting new bindings to a neighborhood
      of already bound sites.
    """
    
    #print ('Adding B event selected', random_B)
    
    if (list_B[random_B, :] == -1).all():  #check if B is free, i.e. all the binding sites are empty (representing by -1 in the original matrix)
        #print ('All random sites !')
        random_binding_site = np.random.randint(0, list_B.shape[1]) #choose randomly one of the sites
        #print ('Random binding site', random_binding_site)
        #print ('List A and List B before having added B, ', list_A, list_B)
        starting_site = 0
        ending_site = len(list_DNA[0]) 
        list_A, list_B, does_B_bind = add_B_to_DNA_site (list_A, list_B, list_DNA, starting_site, ending_site, random_binding_site, random_B, does_B_bind ) #add B to one random site within all DNA
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
        list_A, list_B, does_B_bind = add_B_to_DNA_site (list_A, list_B, list_DNA, starting_site, ending_site, random_B_free_site, random_B, does_B_bind )
        #print ('List A and List B after having added B, ', list_A, list_B)
        
                                                              
    return list_DNA, list_A, list_B, does_B_bind

def add_B_to_DNA_site (list_A, list_B, list_DNA, starting_site, ending_site, random_binding_site, random_B, does_B_bind ):
    
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
        does_B_bind[0,A_to_bind] = 0 #putting the corresponding equal to 0 to mark that a B has bound 
    return list_A, list_B, does_B_bind
    
        
def remove_B_event (list_DNA, list_A, list_B,random_B, L, Eba):
#def remove_B_event (list_A, list_B,random_B, Eba):
 #print ('Removing B event')
 bound_B_sites =np.where (list_B[random_B, :] != -1)[0] #get the binding sites that are bound to the DNA
 #print ('Bound B sites', bound_B_sites)
 # randomly choose one index within the bound sites 
 random_B_bound_site = random.choice(bound_B_sites) #pick randomly one of the bound sites 
 #print ('Random B bound site', random_B_bound_site )
 random_binding = np.random.random()  #draw a random number between 0 and 1
 #energy = energy_unbind_function_B (list_B[random_B, :], Eba)#compute unbinding energy
 DNA_site_to_which_A_bound = list_A[list_B[random_B, random_B_bound_site], 0]
 energy = energy_unbind_function_B_adjacent (list_DNA, list_A, DNA_site_to_which_A_bound, L, Eba)
 
 
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
    #print ('energy unbind B', B, energy  )
    return energy


def energy_unbind_function_B_adjacent (list_DNA, list_A, A_corresponding_B_to_unbind, L, Eba): 
    
    #computing the number of sites in B protein that are bound
    energy = compute_energy_B_binding_adjacent (list_DNA, list_A, A_corresponding_B_to_unbind, L, Eba) 
    #print ('energy unbind B', B, energy  )
    return energy


def energy_function_unbinding_A (list_DNA, list_A, A, list_B, E_ad, E_aa): 
    DNA_site_index = A[0] #find the index to which A is bound 
    if A[1] != (-1 ): #There is a B bound to that A
            #B = list_B[A[1],:] #all binding sites for that B 
            #print (B)
            B_energy = np.nan #set this energy to infinity so it will never leave
            #print (energy_function_unbinding_A, B_energy)
            #B_energy = compute_energy_B_binding_adjacent (list_DNA, list_A,DNA_site_index, k, E_ab) 
            #B_energy = compute_energy_B_binding (B, E_ab) #compute the energy of B protein bound to A, depending on the number of B sites 
            #print ('energy unbind a', B, B_energy  )
    else: 
            B_energy = 0 #if there are no B bound to A the corresponding energy will be zero 
            #print ('energy unbind a', B_energy  )  
    if B_energy != (np.nan):
        A_energy = compute_energy_A_binding (DNA_site_index, list_DNA, E_ad, E_aa)
        #print ('DNA', DNA_site_index, list_DNA)
        total_energy = A_energy + B_energy
    #print ('energy unbind a', A_energy ) 
    else:
        total_energy = B_energy
    #print ('Total energy', total_energy)
    return total_energy


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


def compute_energy_B_binding_adjacent (list_DNA, list_A, DNA_site, L, energy):
    
    if DNA_site-L < 0:
        starting_site = 0
    else:
        starting_site = DNA_site-L
    if DNA_site+ L > list_DNA.shape[1]:
        ending_site = list_DNA.shape[1]
    else:
        ending_site = DNA_site+ L
        
        
    range_sites_within_L = list_DNA [0, starting_site:ending_site+1] 
    
    sites_with_bound_A = list (np.where (range_sites_within_L == 1)[0]) #find the sites in the specified distance
  
    
    if starting_site != 0 : 
        # Convert relative indices to original indices in list_DNA
        sites_with_bound_A = [i + starting_site for i in sites_with_bound_A]
   
    sites_with_bound_B = [] #to count the adjacent sites that are bound to a B
    
    for site in sites_with_bound_A[:]:  # Use a slice to create a copy of the list
        corresponding_A_to_site = np.where(list_A[:, 0] == site)[0]
        
     
        if list_A[corresponding_A_to_site, 1] != -1:  # If that A is already bound to a B
            
            sites_with_bound_B.append(site)  # Add the site
            

    if len (sites_with_bound_B) > 0:      
        number_adjacent_bound_Bs = len (sites_with_bound_B)  
    else: 
        number_adjacent_bound_Bs = 0
    
    
    p = number_adjacent_bound_Bs / (2*L)
    B_energy = (p)*energy
    return B_energy


def compute_probability_removal (k_bound):
#function that takes the number of B bound sites and computes the probability of removal event for B 
#We want to be 0 if all the sites are free and 1 if all the sites are occupied 
    p = k_bound/(1+k_bound)
    return p     