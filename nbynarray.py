# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 12:27:39 2022

@author: sdixo
"""

import scipy as cp
import numpy as np
import random
import matplotlib.pyplot as plt
plt.style.use("ggplot")


a = np.random
def matrix_testing(prevalence: float, n: int):
    
    test_matrix = np.zeros(n)
    number_sick = int(prevalence * len(test_matrix)) # upper bound for sick people in population
    
    for i in range(number_sick): # populate test matrix with sick people
        random_index = random.randint(0, len(test_matrix)-1)
        test_matrix[random_index] = 1
        
    a = test_matrix.reshape(int(n**(0.5)), int(n**(0.5)))
    b = a.transpose()
    

    infected_indices = [[], []] # [horizontal indices positive][vertical indices positive]
    for indx2, row_vector in enumerate(a):
        pooled_test_h = sum(row_vector)
    
        if pooled_test_h > 0: # someone in pool is positive
            infected_indices[0].append(indx2)
    
    for indx1, vertical_vector in enumerate(b):
        pooled_test_v = sum(vertical_vector)
      
        if pooled_test_v > 0: # someone in pool is positive
            infected_indices[1].append(indx1)
            
                
    return number_sick, max(len(infected_indices[0]), len(infected_indices[1]))


total = np.array(())
caught = np.array(())
for j in range(10):
    c = matrix_testing(0.1, 100)
    total = np.append(total, c[0])
    caught = np.append(caught, c[1])
    
plt.bar([i for i in range(1,len(total)+1)], (total-caught))
plt.xlabel("Test Number")
plt.ylabel("Cases Missed with COMP Scheme")
plt.title("10 Tests using COMP 10x10 array with groups of 10 pooling 5% prevalence")


        
    