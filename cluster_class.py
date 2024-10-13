# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 09:11:45 2024

@author: virgi
"""

class Cluster:
    def __init__(self, position, size= 0 ):
        self.position_first_tf = position
        self.size = size #initial size will be 0 and then changed when the counter stops counting 
    def set_size (self, size):
        self.size = size 
    def get_size (self):
        return self.size
        
    