# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 15:55:20 2024

@author: virgi
"""

#CLUSTERING ALGORITHM 

def count_consecutive_ones(DNA_list):
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

    return group_sizes, max_count, positions_first_clusters


# Generate a large test list
#import random
#DNA_test = [random.choice([0, 1]) for x in range(10**6)]

lst = [1, 1, 0, 1, 1, 1, 0, 1, 1]
group_sizes, max_count, positions_first_clusters = count_consecutive_ones(lst)

print(f"Group sizes: {group_sizes}") 
print(f"Max count: {max_count}") 
print(f"Position of the first member of each cluster: {positions_first_clusters}") 

