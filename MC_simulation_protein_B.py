# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:13:54 2024

@author: virgi
"""


###PARAMETERS###

alfa = 0.15 #ratio between nA/N 
N = 10000 #total number of binding sites in the DNA
nA = int (N*alfa) #number of As


nB = 0 #number of Bs
k = 0 #number of B interacting sites with A 






#DNA parameters 
list_DNA = [0] * N #List representing the DNA sites: 0 corresponds to an empty site, 1 to a site with an transcription factor bound to it, 2 when also a B protein is bound- at the initial state there is no bound A nor B to the DNA 
list_empty_DNA  = list(range (0,len(list_DNA),1)) #list containing all the indexes of empty sites in the DNA
list_B = [[-1]*k for _ in range(nB)]# -1 if B is free, will store the position on the DNA when bound 


#Transcription factors parameters 
list_A = [-1]*nA #List representing the transcription factors: -1 for unbound, will store the position of the site on the DNA when bound 
