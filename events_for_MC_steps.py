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
    If the A is free it is added otherwhise the move is ignored 
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
    If the A is free the move is ignored, otherwise there is a probalistic condition depending on an energy condition (computed by the function energy_function_unbinding_A)
    If the energy condition is met then remove_A_from_DNA_site is called.
    If a B is bound to an A, and an A is removed also the corresponding B is removed.

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
        The second column is set to -1 if the A is not bound to a B, or the index of the corresponding B (the row in the matrix list_B) otherwise 
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
    time_step : int
        The current simulation time step.

    Updates:
    -------
    The variables are updtated in the function 'remove_A_from_DNA_site'

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
    It should be noted that the energy of unbinding if the A is bound to a B is now set to infinity, making it impossible
    for an A to unbind if it is bound to a B. If in future this assumption is to be removed, you can only change the energy of unbinding
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
        A matrix tracking the state of A proteins. Contains the indices of the DNA sites where each A is bound in the first column and if the A is bound to a B in the second column.
    list_B : numpy.ndarray
        A matrix tracking the state of B proteins. Rows correspond to individual B proteins, and columns represent
        their multiple binding sites. A value of `-1` indicates a free binding site, othetwhise the index (row) of the corresponding A is tracked
    random_B : int
        The index of the B protein selected for the binding event.
    L : int
        The range around a bound DNA site to consider for additional B bindings.
    does_B_bind : 
        

    Updates:
    -------
    The updates are done by the function add_B_to_DNA_site
    
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
    """
    Attempts to bind a B protein to a DNA site that is already occupied by an A protein, updating the binding states accordingly.
    
    Parameters:
    ----------
    list_A : numpy.ndarray
        A 2D array where each row represents an A protein. The first column contains the DNA site bound by A, 
        and the second column stores the B protein it is bound to (-1 if unbound).
    list_B : numpy.ndarray
        A 2D array where each row represents a B protein, and the columns store the binding state of each B.
    list_DNA : numpy.ndarray
        A 1D array representing the DNA strand, where 1 indicates an A protein is bound and 0 means unoccupied.
    starting_site : int
        The start index of the DNA region where B attempts to bind.
    ending_site : int
        The end index of the DNA region where B attempts to bind.
    random_binding_site : int
        A randomly chosen binding site within the possible binding region.
    random_B : int
        The index of the B protein attempting to bind.
    does_B_bind : numpy.ndarray
        A 2D array that tracks whether a B protein has successfully bound to an A protein.
    
    Returns:
    -------
    list_A : numpy.ndarray
        Updated array reflecting A proteins that have been bound by B.
    list_B : numpy.ndarray
        Updated array reflecting B proteins that have successfully bound to A.
    does_B_bind : numpy.ndarray
        Updated binding state indicating which A proteins have been bound by B.
    
    Notes:
    ------
    - The function first checks for A proteins already bound within the specified DNA region.
    - If an A protein in this region is already bound to another B, it is ignored.
    - If one or more unbound A proteins are found, a random one is selected for binding.
    - The binding updates the corresponding entries in `list_A`, `list_B`, and `does_B_bind`.
    """
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
    
     """
    Attempts to remove a B protein from its bound DNA site based on an energy-dependent probability.
    
    Parameters:
    ----------
    list_DNA : numpy.ndarray
        A 1D array representing the DNA strand, indicating which sites are occupied.
    list_A : numpy.ndarray
        A 2D array where each row represents an A protein. The first column contains the DNA site bound by A, 
        and the second column stores the B protein it is bound to (-1 if unbound).
    list_B : numpy.ndarray
        A 2D array where each row represents a B protein, and the columns store the binding state of each B.
    random_B : int
        The index of the B protein attempting to unbind.
    L : int
        The interaction range used to calculate energy adjustments.
    Eba : float
        The binding energy between A and B, used to determine the unbinding probability.
    
    Returns:
    -------
    list_A : numpy.ndarray
        Updated array reflecting the removal of B proteins from A.
    list_B : numpy.ndarray
        Updated array reflecting the unbound state of B proteins.
    
    Notes:
    ------
    - The function first identifies the sites where the selected B protein is bound.
    - A random site among the bound locations is chosen for the unbinding attempt.
    - An unbinding probability is calculated using an energy-based function.
    - If the unbinding condition is met, the B protein is removed from the A protein it was bound to.
    - The removal updates the corresponding entries in `list_A` and `list_B`.
      """
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
    
    """
    Removes a B protein from its bound DNA site, updating the binding states of both A and B proteins.
    
    Parameters:
    ----------
    random_B_bound_site : int
        The specific site where the B protein is currently bound.
    list_A : numpy.ndarray
        A 2D array where each row represents an A protein. The first column contains the DNA site bound by A, 
        and the second column stores the B protein it is bound to (-1 if unbound).
    list_B : numpy.ndarray
        A 2D array where each row represents a B protein, and the columns store the binding state of each B.
    random_B : int
        The index of the B protein being removed.
    
    Returns:
    -------
    list_A : numpy.ndarray
        Updated array reflecting that the A protein is no longer bound to B.
    list_B : numpy.ndarray
        Updated array indicating that the B protein is now unbound.
    
    Notes:
    ------
    - The function identifies the A protein that was bound to the selected B.
    - It updates `list_A` to mark the A protein as unbound.
    - It updates `list_B` to indicate that the B protein is now free.
    """

    index_bound_A = list_B[random_B,random_B_bound_site]
    #print (index_bound_A)
    list_A[index_bound_A, 1] = -1 #marking that the B is no longer bound to the A
    list_B[random_B,random_B_bound_site] = -1 #setting B as free in the list of Bs
    
    return list_A, list_B

###ENERGY FUNCTIONS 
'''
Below there are the functions used to compute the energies. For the unbinding of B two energies functions have been investigated.
The energy effect of an A to the unbiding of a B is denoted Eba 
1. A first approach was to compute the energy unbinding of the B proportional to the number of site bound of that same B
2. Then, to investigate more the 'clustering effect of the B' we computed the energy proportional to the number of bound B within a distance |L|

'''

def energy_unbind_function_B (B, Eba): 
    
    '''
    Function to compute the energy binding for the unbinding of a B, proportional to the number of site bound of that same B
    The function compute_energy_B_binding is called for computing the energy value
    
    Parameters:
    ----------
    B : list
        The B (all the binding sites) for which the energy must be computed
    Eba : int
        energy effect of an A on the unbinding of a B
    
    Returns:
    -------
    energy: int 
        total energy value used for compute if the removing event of a B succeeds 
    
    '''
    
    #computing the number of sites in B protein that are bound
    energy = compute_energy_B_binding (B, Eba)
    #print ('energy unbind B', B, energy  )
    return energy


def energy_unbind_function_B_adjacent (list_DNA, list_A, A_corresponding_B_to_unbind, L, Eba): 
    
    
    '''
    Function to compute the energy associated with the unbinding of a B protein, proportional to the number of B proteins bound 
    within a distance |L|. It calls the `compute_energy_B_binding_adjacent` function to calculate the energy value.

    Parameters:
    ----------
    list_DNA : numpy.ndarray
        A 1D array representing the DNA strand, where each element indicates whether a site is occupied or not.
        
    list_A : numpy.ndarray
        A 2D array where each row represents an A protein. The first column contains the DNA site bound by A, 
        and the second column stores the index of the B protein it is bound to (-1 if unbound).
        
    A_corresponding_B_to_unbind : int
        The index of the A protein corresponding to the B protein that is attempting to unbind.
        
    L : int
        The maximum distance within which the influence of other B proteins is considered in the energy calculation.
        
    Eba : float
        The binding energy between an A protein and a B protein used in the energy calculation.

    Returns:
    -------
    energy : float
        The total energy value that determines whether the unbinding event of the B protein succeeds.
    
    Notes:
    ------
    - The function computes the energy for unbinding based on the number of B proteins bound within the range |L|.
    - The function depends on the `compute_energy_B_binding_adjacent` function to compute the energy value.
    '''
    
    #computing the number of sites in B protein that are bound
    energy = compute_energy_B_binding_adjacent (list_DNA, list_A, A_corresponding_B_to_unbind, L, Eba) 
    #print ('energy unbind B', B, energy  )
    return energy


def energy_function_unbinding_A (list_DNA, list_A, A, list_B, E_ad, E_aa): 
    
    """
    Computes the energy required for unbinding an A protein from a DNA site, considering its interactions with B proteins and DNA.
    
    Parameters:
    ----------
    list_DNA : numpy.ndarray
        A 1D array representing the DNA strand, where each element indicates whether a site is occupied.
    
    list_A : numpy.ndarray
        A 2D array where each row represents an A protein. The first column contains the DNA site bound by A, 
        and the second column stores the index of the B protein it is bound to (-1 if unbound).
    
    A : numpy.ndarray
        A 1D array representing a single A protein, where:
        - `A[0]` is the DNA site index to which A is bound.
        - `A[1]` is the index of the B protein it is bound to (-1 if unbound).
    
    list_B : numpy.ndarray
        A 2D array where each row represents a B protein, with columns storing its binding states.
    
    E_ad : float
        Energy contribution from A binding directly to DNA.
    
    E_aa : float
        Energy contribution from A-A interactions in the system.
    
    Returns:
    -------
    total_energy : float
        The total energy required for A to unbind from the DNA site.
    
    Notes:
    ------
    - If A is bound to a B protein, the energy is set to NaN (effectively preventing unbinding).
    - Otherwise, the unbinding energy is computed as the sum of:
      - `compute_energy_A_binding()`, which accounts for A-DNA and A-A interactions.
      - A B-related energy term (set to 0 if A is not bound to any B protein).
    """
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
    """
    Computes the binding energy of an A protein at a given DNA site, considering both direct DNA binding 
    and contributions from neighboring sites.

    Parameters:
    ----------
    index : int
        The DNA site index where the A protein is bound.
    
    list_DNA : numpy.ndarray
        A 2D array representing the DNA strand, where the first row (`list_DNA[0, :]`) indicates site occupancy.
    
    E_ad : float
        The energy contribution from the direct binding of A to its DNA site.
    
    E_aa : float
        The energy contribution from neighboring A proteins at adjacent sites.

    Returns:
    -------
    energy : float
        The total energy associated with the A protein at the specified DNA site.

    Notes:
    ------
    - If the A protein is at the first site (`index == 0`), only the right neighbor contributes.
    - If the A protein is at the last site, only the left neighbor contributes.
    - Otherwise, both left and right neighbors contribute.
    - The direct binding energy (`E_ad`) is always included in the calculation.
    """
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
    '''
    Functions to compute the energy of a B protein
    Assumption: this energy is directly proportional proportional to the number of bound B within a distance |L|
    '''
    
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