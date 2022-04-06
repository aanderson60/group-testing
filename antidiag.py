# local copy to make changes and get plots



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

import random
from statistics import mean
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linprog


# Function that will build the test matrix and run the testing scheme
def matrixSimulation(a=False):
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


    rowTests = []
    colTests = []
    diagTests = []
    adiagTests = []

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
    
   

    if a: # do antidiag testing
        flipM = np.fliplr(M)
        for adiag in range(n):
            adiagTests.append(0)
            combo = np.concatenate((np.diagonal(flipM,adiag),np.diagonal(flipM,-(n-adiag))))
            for idx in combo:
                if idx == 1:
                    adiagTests[adiag] = 1
        
        testPositives = COMP(rowTests,colTests,diagTests,adiagTests)
    else:
        testPositives = COMP(rowTests,colTests,diagTests)
        
    # Perform recovery algorithm to determine DND/PDs
    
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
            numIncorrect+=1
    # Check for false negatives
    for item in truePositives:
        if not item in testPositives:
            numIncorrect+=1
    #if numIncorrect > 0:
        # Going to write this to file in the future
        #printM(M,[],[])
    if totalCount == 0: # Case with no positive individuals
        c = 0.0
    else:
        c = numIncorrect/I
    return numIncorrect

# Debugging function for printing the test matrix and row/col tests
def printM(M,rowTests=[],colTests=[],diagTests=[], adiagTests=[]):
    print("Test Matrix:")
    for row in range(n):
        print(M[row])
    print()
    if len(rowTests) > 0 and len(colTests) > 0:
        print("Row Tests:",rowTests)
        print("Column Tests:",colTests)
    if len(diagTests) > 0:
        print("Diagonal Tests:",diagTests)
    if len(adiagTests) > 0:
        print("Anti Diagonal Tests",adiagTests)

# Recovery algorithm using the COMP algorithm
def COMP(rowTests,colTests,diagTests,adiagTests=[]):
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
                    
    # n_array is representing each sample as a unique integer
    # from 0-99 ( 100 total )
    
    
    n_array = np.arange(0, 100).reshape(10, 10)    
            
    for k in range(len(diagTests)):
        if diagTests[k] == 0:
            # Add entire diag to DND list
            diag = np.concatenate((np.diagonal(n_array,k),np.diagonal(n_array,-(n-k))))
            for item in diag:
                if not item in DND:
                    DND.append(item)

    for k in range(len(adiagTests)):
        if adiagTests[k] == 0:
            # Add entire diag to DND list
            diag = np.concatenate((np.diagonal(np.fliplr(n_array),k),np.diagonal(np.fliplr(n_array),-(n-k))))
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
def monteCarlo(mc, a=False):
    results = []
    for i in range(mc):
        results.append(matrixSimulation(a))
    #outputResults(results,mc)
    return results



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

"""
******
Plots below for Histograms with varying prevalence
******
"""


'''
plt.plot(prevs,results10k,label='10k')
plt.title("False Positive Rates for 10x10 COMP Testing Scheme")
plt.xlabel("Prevalences")
plt.ylabel("Proportion of False Positives")
plt.legend()
plt.show()
'''

'''P = 0.02
result = monteCarlo(10000)
labels, counts = np.unique(result, return_counts=True)
plt.bar(labels, counts, align='center')
plt.title("False Positives for 10000 Simulations with 2% Prevalence")
plt.gca().set_xticks(labels)
plt.savefig("Histogram10k2padiag")
plt.show()


P = 0.03
result = monteCarlo(10000)
plt.hist(result)
plt.title("False Positives for 10000 Simulations with 3% Prevalence")
labels, counts = np.unique(result, return_counts=True)
plt.bar(labels, counts, align='center')
plt.gca().set_xticks(labels)
plt.savefig("Histogram10k3padiag")
plt.show()


P = 0.04
result = monteCarlo(10000)
labels, counts = np.unique(result, return_counts=True)
plt.bar(labels, counts, align='center')
plt.title("False Positives for 10000 Simulations with 4% Prevalence")
plt.gca().set_xticks(labels)
plt.savefig("Histogram10k4padiag")
plt.show()

P = 0.05
result = monteCarlo(10000)
labels, counts = np.unique(result, return_counts=True)
plt.bar(labels, counts, align='center')
plt.title("False Positives for 10000 Simulations with 5% Prevalence")
plt.gca().set_xticks(labels)
plt.savefig("Histogram10k5padiag")
plt.show()'''



"""
Plots below for average 
false positives for varying prevalances for
4 total pools (COMP4):
    anti-diag
    diag
    rows
    cols
3 total pools (COMP3):
    diag
    rows
    cols
"""

plt.style.use('ggplot')

P = 0.02
COMP4 = []
COMP3 = []
prevs = []

while P <= 0.05:
    prevs.append(P)
    COMP4.append(mean(monteCarlo(1000, True)))
    P += 0.01
    
P = 0.02

while P <= 0.05:
    COMP3.append(mean(monteCarlo(10000)))
    P += 0.01
    
plt.plot(prevs, COMP3, label="COMP5")
plt.plot(prevs, COMP4, label="COMP4")
plt.xlabel("Prevalences")
plt.ylabel("Average Percentage of False Positives")
plt.legend()
plt.title("10k tests for varying prevalences with 2 different COMP pooling schemes")
plt.show()