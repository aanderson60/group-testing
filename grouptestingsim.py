# Goals:
#
# Generate an array of size nxn and fill each index with a 1 or 0 depending on the prevalance
# "Pool" along rows and columns and determine whether or not the group is positive
# Attempt to determine which individuals are positive and check results
#
# Repeat this process many times and return final results
# TODO: Allow support for other types of testing (currently only COMP)

# The desired prevalance
P = 0.0055555

# The number of individuals
I = 900

# The matrix size
n = 30

import random
from statistics import mean
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linprog
import math
import seaborn as sns
import matplotlib.pylab as plb
from itertools import combinations
from matplotlib import animation

# Function that will build the test matrix and run the testing scheme
def matrixSimulation(desiredPositives=None):
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
	if desiredPositives == None:
		for j in range(round(P*I)):
			posIndex = random.randint(0,I-1)
			while posIndex in truePositives:
				# Generate a new number (prevent duplicates)
				posIndex = random.randint(0,I-1)
			truePositives.append(posIndex)
		# Generate test matrix
		for i in range(n):
			row=[]
			for j in range(n):
				if i*n+j in truePositives:
					row.append(1)
				else:
					row.append(0)
			M.append(row)
	else:
		for item in desiredPositives:
			truePositives.append(item)
		for i in range(n):
			row=[]
			for j in range(n):
				if i*n+j in truePositives:
					row.append(1)
				else:
					row.append(0)
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
		

	#printM(M,rowTests,colTests,diagTests)

	# Perform recovery algorithm to determine DND/PDs
	testPositives = COMP(rowTests,colTests,diagTests)
	# LP Stuff
	'''
	floatResult = linear_prog(rowTests,colTests)
	testPositives = []
	for i,val in enumerate(floatResult):
		# Assume positive result if val is >10^-2
		if val > 0.01:
			testPositives.append(i)
	'''
	falsePositives = []
	# Compare results
	# Define c = proportion of individuals incorrect to those correctly identified
	numIncorrect = 0
	totalCount = len(testPositives)
	# Check for false positives
	for item in testPositives:
		if not item in truePositives:
			falsePositives.append(item)
			numIncorrect+=1
	# Check for false negatives
	for item in truePositives:
		if not item in testPositives:
			numIncorrect+=1
	#if numIncorrect > 0:
		# Going to write this to file in the future
		#printM(M,[],[])
	if len(falsePositives) != numIncorrect:
		raise Exception("Number of false positives not equal to number incorrect (false negative detected; something went wrong!)")
	if totalCount == 0: # Case with no positive individuals
		c = 0.0
	else:
		c = numIncorrect/I
	# Return the specificity ratio and which items were 'false positive'
	return ((I-numIncorrect)/I,len(testPositives))

# Debugging function for printing the test matrix and row/col tests
def printM(M,rowTests=[],colTests=[],diagTests=[]):
	print("Test Matrix:")
	for row in range(n):
		print(M[row])
	print()
	if len(rowTests) > 0 and len(colTests) > 0 and len(diagTests) > 0:
		print("Row Tests:",rowTests)
		print("Column Tests:",colTests)
		print("Diag Tests: ",diagTests)

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
	n_array = np.arange(0, I).reshape(n, n)
	for k in range(len(diagTests)):
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
	
# WIP Function
def DD(rowTests, colTests, diagTests):
	# Use DD algorithm to find defectives
	# 1) Any item in a negative test is DND
	# 2) If an item is the only one in a positive test it is DD
	# The possible defectives
	PD = []
	# The definite non-defectives
	DND = []
	posTestItems = []
	# First determine the definite non-defectives
	for i in range(len(rowTests)):
		if rowTests[i] == 0:
			# Add the entire row to the DND list
			for x in range(len(rowTests)):
				DND.append(i*n+x)
		if rowTests[i] == 1:
			for x in range(len(rowTests)):
				posTestItems.append(i*n+x)
	for j in range(len(colTests)):
		if colTests[j] == 0:
			# Add the entire col to the DND list (if not already present)
			for y in range(len(colTests)):
				if not((y*n+j) in DND):
					DND.append(y*n+j)
		if colTests[j] == 1:
			for y in range(len(colTests)):
				if not ((y*n+j) in posTestItems):
					posTestItems.append(y*n+j)
	n_array = np.arange(0, I).reshape(n, n)
	for k in range(len(diagTests)):
		if diagTests[k] == 0:
			# Add entire diag to DND list
			diag = np.concatenate((np.diagonal(n_array,k),np.diagonal(n_array,-(n-k))))
			for item in diag:
				if not item in DND:
					DND.append(item)
		if diagTests[k] == 1:
			diag = np.concatenate((np.diagonal(n_array,k),np.diagonal(n_array,-(n-k))))
			for item in diag:
				if not item in posTestItems:
					posTestItems.append(item)
	# Now determine the PDs where the PDs are whatever does not exist in the DNDs
	for m in range(I):
		if not (m in DND):
			PD.append(m)
		
	

# Output results
def outputResults(results, mc):
	print("Results for simulation of size ",mc,", a prevalance of ", P, " and a matrix size of ",n,"x",n," : ",round(mean(results),4),sep="")

# Monte-Carlo simulation
def monteCarlo(mc):
	failureRates = []
	falsePositives = []
	for i in range(mc):
		results = matrixSimulation()
		falsePositives.append(results[1])
		failureRates.append(results[0])
	outputResults(failureRates,mc)
	return (mean(failureRates),falsePositives)

# Function to exhaust every possibility for a 10x10 matrix with 3 positive individuals
def desiredPositivesSweep(start,stop):
	failureRates = []
	falsePositives = []
	arr = np.linspace(start,stop,10,dtype=int)
	comb = combinations(arr,3)
	for j,i in enumerate(list(comb)):
		result = matrixSimulation(i)
		for item in result[1]:
			falsePositives.append(item)
		failureRates.append(result[0])
	falsePositives = np.array(result[1])
	outputM = np.zeros((n,n))
	# Must convert false positive individual numbers into nxn array
	for item in falsePositives:
		outputM[int(np.floor((item-1)/n))][(item-1)%n] += 1
	return outputM


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
plt.style.use('seaborn-paper')

sizes = []
prevs = []
ns = []
tests = []
results100 = []
results1k = []
results10k = []
prev = 0.02
n=10

while n <= 50:
	prevs.append(prev)
	ns.append(n)
	tests.append(3*n)
	P = prev
	I = n*n
	#results100.append(monteCarlo(100))
	#results1k.append(monteCarlo(1000))
	results10k.append(monteCarlo(1000)[0])
	n += 5

#result = monteCarlo(10000)[1]
#plt.plot(prevs,results100,label='100')
#plt.plot(prevs,results1k,label='1k')

plt.plot(tests,results10k)
plt.title("Specificity of 10x10 COMP3 Scheme for Varying Number of Samples")
plt.xlabel("Number of Tests")
plt.ylabel("Specificity")
plt.grid(linestyle='--',linewidth=0.5)
plt.show()
plt.plot(ns,results10k)
plt.grid(linestyle='--',linewidth=0.5)
plt.show()
'''
result = np.array(result)
bins = np.arange(result.min(), result.max() + 1.5) - 0.5
fig,ax = plt.subplots()
_= ax.hist(result,bins,rwidth=0.75)
ax.set_xticks(bins+0.5)
plt.show()
'''
# Heat map of false positives
'''
outputM = desiredPositivesSweep()
ax = sns.heatmap(outputM)
plt.show()
'''