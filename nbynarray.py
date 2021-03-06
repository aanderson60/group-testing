# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 12:27:39 2022

@author: sdixo
"""

from scipy.optimize import linprog
from statistics import mean
import scipy as cp
import numpy as np
import random
import matplotlib.pyplot as plt
plt.style.use("ggplot")


# a = np.random
# def matrix_testing(prevalence: float, n: int):

#     test_matrix = np.zeros(n)
#     number_sick = int(prevalence * len(test_matrix)) # upper bound for sick people in population

#     for i in range(number_sick): # populate test matrix with sick people
#         random_index = random.randint(0, len(test_matrix)-1)
#         test_matrix[random_index] = 1

#     a = test_matrix.reshape(int(n**(0.5)), int(n**(0.5)))
#     b = a.transpose()


#     infected_indices = [[], []] # [horizontal indices positive][vertical indices positive]
#     for indx2, row_vector in enumerate(a):
#         pooled_test_h = sum(row_vector)

#         if pooled_test_h > 0: # someone in pool is positive
#             infected_indices[0].append(indx2)

#     for indx1, vertical_vector in enumerate(b):
#         pooled_test_v = sum(vertical_vector)

#         if pooled_test_v > 0: # someone in pool is positive
#             infected_indices[1].append(indx1)


#     return number_sick, max(len(infected_indices[0]), len(infected_indices[1]))


# total = np.array(())
# caught = np.array(())
# for j in range(10):
#     c = matrix_testing(0.1, 100)
#     total = np.append(total, c[0])
#     caught = np.append(caught, c[1])

# plt.bar([i for i in range(1,len(total)+1)], (total-caught))
# plt.xlabel("Test Number")
# plt.ylabel("Cases Missed with COMP Scheme")
# plt.title("10 Tests using COMP 10x10 array with groups of 10 pooling 5% prevalence")

# my first thougts above for coding a model for this problem


# Grabbing Alex's Code and working from here


# Goals:
#
# Generate an array of size nxn and fill each index with a 1 or 0 depending on the prevalance
# "Pool" along rows and columns and determine whether or not the group is positive
# Attempt to determine which individuals are positive and check results
#
# Repeat this process many times and return final results
# TODO: Allow support for other types of testing (currently only COMP)

# The desired prevalance
P = 0.04

# The number of individuals
I = 100

# The matrix size
n = 10


# Function that will build the test matrix and run the testing scheme

def matrixSimulation():
    # The test matrix M
    M = []

    # Array with positive individual p
    truePositives = []

    # Generate test matrix according to defined prevalance
    '''
	for i in range(n):
		row = []
		for j in range(n):
			# Generate a random number with prevalance P
			ran = random()
			if (ran > P):
				row.append(0)
			else:
				row.append(1)
				truePositives.append(i*n+j)
		M.append(row)
	'''
    # Generate positive indices with a fixed prevalance
    for j in range(int(P*I)):
        posIndex = random.randint(0, I-1)
        while posIndex in truePositives:
            # Generate a new number (prevent duplicates)
            posIndex = random.randint(0, I-1)
        truePositives.append(posIndex)
    # Generate test matrix
    for i in range(n):
        row = []
        for j in range(n):
            if i*n+j in truePositives:
                row.append(1)
            else:
                row.append(0)
        M.append(row)
    printM(M)
    rowTests = []
    colTests = []
    diagTests = []

    # Perform 'tests' on rows and columns (and diags)
    # If any individual is positive in the row/col, the entire test will return a 'positive' result

    for row in range(n):
        rowTests.append(0)
        for col in range(n):
            if (M[row][col] == 1):
                rowTests[row] = 1

    for col in range(n):
        colTests.append(0)
        for row in range(n):
            if (M[row][col] == 1):
                colTests[col] = 1

    for diag in range(n):
        diagTests.append(0)
        combo = np.concatenate(
            (np.diagonal(M, diag), np.diagonal(M, -(n-diag))))
        for idx in combo:
            if idx == 1:
                diagTests[diag] = 1

    # printM(M,rowTests,colTests)

    # Perform recovery algorithm to determine DND/PDs
    testPositives = COMP(rowTests, colTests, diagTests)
    # LP Stuff
    '''
	floatResult = linear_prog(rowTests,colTests)
	testPositives = []
	for i,val in enumerate(floatResult):
		# Assume positive result if val is >10^-2
		if val > 0.01:
			testPositives.append(i)
	'''
    # Compare results
    # Define c = proportion of individuals incorrect to those correctly identified
    numIncorrect = 0
    totalCount = len(testPositives)
    # Check for false positives
    for item in testPositives:
        if not item in truePositives:
            numIncorrect += 1
    # Check for false negatives
    for item in truePositives:
        if not item in testPositives:
            numIncorrect += 1
    # if numIncorrect > 0:
        # Going to write this to file in the future
        # printM(M,[],[])
    if totalCount == 0:  # Case with no positive individuals
        c = 0.0
    else:
        c = numIncorrect/I
    return numIncorrect

# Debugging function for printing the test matrix and row/col tests


def printM(M, rowTests=[], colTests=[]):
    print("Test Matrix:")
    for row in range(n):
        print(M[row])
    print()
    if len(rowTests) > 0 and len(colTests) > 0:
        print("Row Tests:", rowTests)
        print("Column Tests:", colTests)

# Recovery algorithm using the COMP algorithm


def COMP(rowTests, colTests, diagTests):
    # Use COMP algorithm to find set of defectives
    # The possible defectives
    PD = []
    # The definite non-defectives
    DND = []
    # First determine the definite non-defectives
    for i in range(len(rowTests)):
        if rowTests[i] == 0:
            # Add the entire row to the DND list
            for x in range(len(rowTests)):
                DND.append(i*n+x)
    for j in range(len(colTests)):
        if colTests[j] == 0:
            # Add the entire col to the DND list (if not already present)
            for y in range(len(colTests)):
                if not((y*n+j) in DND):
                    DND.append(y*n+j)
    for k in range(len(diagTests)):
        n_array = []
        for x in range(n):
            row = []
            for y in range(n):
                row.append(x*n+y)
            n_array.append(row)
        if diagTests[k] == 0:
            # Add entire diag to DND list
            diag = np.concatenate(
                (np.diagonal(n_array, k), np.diagonal(n_array, -(n-k))))
            for item in diag:
                if not item in DND:
                    DND.append(item)

    # Now determine the PDs where the PDs are whatever does not exist in the DNDs
    for m in range(I):
        if not (m in DND):
            PD.append(m)

    # Assume the PDs are the Definite Defectives (DD)
    return(PD)

# Linear programming algorithm using scipy


def linear_prog(rowTests, colTests):
    # Minimize z (Vector length n*n)
    # If test is positive, each x element in that test Xti subject to, Xti*Zi >= 1
    # If test is negative, each x element in that test Xti subject to, Xti*Zi = 0
    # And each Zi is 0 or 1
    #
    # Scipy parameters:
    # c = output matrix (coefficients of z)
    # A_ub = 2D array, inequality constraint - each row is one constraint set in the form A_ub*x <= b_ub
    # b_ub = 1D array, upper bound on value of A_ub @ x
    # A_eq = 2D array, equality constraint
    # b_eq = 1D array, equality constraint vector
    # Returns: scipy.optimize.OptimizeResult object

    # Min z1,z2,...,zn where n is the number of individuals
    c = np.ones(I)

    # Generate upper bound vector (is actually a lower bound in our specific case)
    # >= 1 is the same as <= -1 in this case (multiply both sides by -1)
    y_t = np.append(rowTests, colTests)
    b_ub_t = []
    b_eq_t = []
    for item in y_t:
        if item:
            b_ub_t.append(-1)
            b_eq_t.append(1)
        else:
            b_ub_t.append(0)
            b_eq_t.append(0)
    # Array with test number on rows and individual on the cols
    # If an index is 1, it means the test corresponding to row contained the individual corresp to that col
    A_ub = np.zeros((2*n, I))
    # Fill row tests with indivs
    for t in range(0, n):
        for i in range(n*t, n*t+n):
            A_ub[t][i] = -1

    # Fill in col tests
    for t in range(n, 2*n):
        for i in range(t-n, I, n):
            A_ub[t][i] = -1

    # Generate equality constraint vector (zero when test is 0)

    A_eq = np.negative(A_ub)

    # Generate bounds of (0,1) for each z
    z_bound = (0, 1)
    bounds = []
    # Fill for every zi
    for i in range(I):
        bounds.append(z_bound)
    result = linprog(c, A_ub, b_ub_t, A_eq, b_eq_t, bounds)
    return result["x"]


# Output results
def outputResults(results, mc):
    print("Results for simulation of size ", mc, ", a prevalance of ", P,
          " and a matrix size of ", n, "x", n, " : ", round(mean(results), 4), sep="")

# Monte-Carlo simulation


def monteCarlo(mc):
    results = []
    for i in range(mc):
        results.append(matrixSimulation())
    # outputResults(results,mc)
    return results

# ------- MAIN --------


plt.style.use('ggplot')
a = matrixSimulation()
