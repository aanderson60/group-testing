# Goals:
#
# Generate an array of size nxn and fill each index with a 1 or 0 depending on the prevalance
# "Pool" along rows and columns and determine whether or not the group is positive
# Attempt to determine which individuals are positive and check results
#
# Repeat this process many times and return final results
# TODO: Allow support for other types of testing (currently only COMP)

# The desired prevalance
P = 0.01

# The number of individuals
I = 100

# The matrix size
n = 10

from random import random
from statistics import mean
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linprog

# Function that will build the test matrix and run the testing scheme
def matrixSimulation():
	# The test matrix M
	M = []

	# Array with positive individual p
	truePositives = []

	# Generate test matrix according to defined prevalance
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
		combo = np.concatenate((np.diagonal(M,diag),np.diagonal(M,-(n-diag))))
		for idx in combo:
			if idx == 1:
				diagTests[diag] = 1
		

	#printM(M,rowTests,colTests)

	# Perform recovery algorithm to determine DND/PDs
	#testPositives = COMP(rowTests,colTests,diagTests)
	# LP Stuff
	
	floatResult = linear_prog(rowTests,colTests)
	testPositives = []
	for i,val in enumerate(floatResult):
		# Assume positive result if val is >10^-2
		if val > 0.01:
			testPositives.append(i)
	
	# Compare results
	# Define c = proportion of individuals incorrect to those correctly identified
	numIncorrect = 0
	totalCount = len(testPositives)
	# Check for false positives
	for item in testPositives:
		if not item in truePositives:
			numIncorrect+=1
	# Check for false negatives
	for item in truePositives:
		if not item in testPositives:
			numIncorrect+=1

	if totalCount == 0: # Case with no positive individuals
		c = 0.0
	else:
		c = numIncorrect/I
	return c

# Debugging function for printing the test matrix and row/col tests
def printM(M,rowTests,colTests):
	print("Test Matrix:")
	for row in range(n):
		print(M[row])
	print()
	print("Row Tests:",rowTests)
	print("Column Tests:",colTests)

# Recovery algorithm using the COMP algorithm
def COMP(rowTests,colTests,diagTests):
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
			diag = np.concatenate((np.diagonal(n_array,k),np.diagonal(n_array,-(n-k))))
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
def linear_prog(rowTests,colTests):
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
	y_t = np.append(rowTests,colTests)
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
	A_ub = np.zeros((2*n,I))
	# Fill row tests with indivs
	for t in range(0,n):
		for i in range(n*t,n*t+n):
			A_ub[t][i] = -1

	# Fill in col tests
	for t in range(n,2*n):
		for i in range(t-n,I,n):
			A_ub[t][i] = -1

	# Generate equality constraint vector (zero when test is 0)

	A_eq = np.negative(A_ub)
	

	# Generate bounds of (0,1) for each z
	z_bound = (0,1)
	bounds = []
	# Fill for every zi
	for i in range(I):
		bounds.append(z_bound)
	result = linprog(c,A_ub,b_ub_t,A_eq,b_eq_t,bounds)
	return result["x"]
	

# Output results
def outputResults(results, mc):
	print("Results for simulation of size ",mc,", a prevalance of ", P, " and a matrix size of ",n,"x",n," : ",round(mean(results),4),sep="")

# Monte-Carlo simulation
def monteCarlo(mc):
	results = []
	for i in range(mc):
		results.append(matrixSimulation())
	outputResults(results,mc)
	return mean(results)

# ------- MAIN --------
'''
print("Group Testing Simulation (COMP)")
print("")
print("Please input the matrix size n (where M is nxn): ")
n = int(input())
I=n*n
print("Please input the prevalance as a decimal value (0 to 1): ")
P = float(input())
print("Please input the number of simulations to be run: ")
mc = int(input())
monteCarlo(mc)
'''
plt.style.use('ggplot')

sizes = []
prevs = []
ns = []
results = []
prev = 0
n=10
while prev <= 0.05:
	prevs.append(prev)
	ns.append(n)
	P = prev
	I = n*n
	val = monteCarlo(100) 
	results.append(val)
	prev += 0.0025

plt.plot(prevs,results)
plt.title("False Positive Rates for 10x10 LP Testing Scheme")
plt.xlabel("Prevalences")
plt.ylabel("Proportion of False Positives")
plt.show()
