import numpy as np
# Recovery algorithm for group testing for a nxn test matrix using the COMP2/3 algorithm [Aldridge et. al.]
# Optional parameter for diagonal tests
# Args:
#   rowTests : array of length n, the results of the row pooling top to bottom
#   colTests : array of length n, the results of the column pooling left to right
#   n : integer representing dimensions of nxn test matrix where the number of individuals I=n*n
#   diagTests: array of length n, optional, the results of the diagonal pooling left to right
# Returns:
#   Tuple of size 3, np.array objects
#   DND : the definite non-defective individuals
#   PD : the probable defective individuals
#   DD : the definite defective individuals
# Note: all individuals are numbered starting at 0, left to right, top to bottom (i.e. for n=10, individuals are numbered 0-99)
# (see n_array for specific format)
#
def COMP(rowTests,colTests,n,diagTests=None):
    # Use COMP algorithm to find set of defectives
    I = n*n
    # The possible defectives
    PD = []
    # The definite non-defectives
    DND = []
    # The definite defectives
    DD = []
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
    # Generate indexing array for diagonals
    n_array = np.arange(0, 100).reshape(10, 10)
    if diagTests != None:
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

    # Determine Definite Defectives (DDs) using the following fact:
    # COMP3 cannot generate false positives for prevelance <= 2%, thus if we receive a number of positives
    #   corresponding to a 1% of 2% prevelance, we know every element in the PDs is a DD.
    #   Additionally, we know that a number of positives corresponding to a 3% prevelance would suggest 
    #   a case where no false positives were generated as well - and that all PDs are once again DDs.
    #   This method could be problematic for non exact prevelances (i.e. not a multiple of 0.01) - use tolerances for now

    supPrev = len(PD)/I
    tolerance = 1e-3

    if ((supPrev - 0.01) < tolerance or (supPrev - 0.02) < tolerance or (supPrev - 0.03) < tolerance):
        for item in PD:
            DD.append(item)
        PD = []
    
    DND = np.array(DND)
    PD = np.array(PD)
    DD = np.array(DD)

    np.sort(DND)
    np.sort(PD)
    np.sort(DD)

    return((DND,PD,DD))

# ---------- MAIN -----------
rowTests = []
colTests = []
diagTests = []
diags = False
# Getting I and n from stdin
print("Enter the total number of individuals: ",end="")
I = float(input())
n = np.sqrt(I)
# Check that matrix is square
if not n == int(n):
    print("ERROR: Number of individuals must be in a square matrix")
    exit()
n = int(n)
# See if diagonal tests should be included
print("Does this testing scheme use diagonal tests? (Y/N): ",end="")
ans = input()
diags = True if ans=="Y" else False
# Make sure answer is in correct format
if ans != "Y" and ans != "N":
    print("ERROR: Please enter Y/N only")
    exit()
print("COMP Testing Recovery for ",n,"x",n," matrix:",sep='')
print("")
# Input loop for row tests
print("Please enter the results of the row tests:")
for i in range(n):
    print("Please enter the test result (0/1) for test R",i+1,": ",sep='',end="")
    rowTests.append(int(input()))
# Input loop for col tests
print("Please enter the results of the column tests:")
for i in range(n):
    print("Please enter the test result (0/1) for test C",i+1,": ",sep='',end="")
    colTests.append(int(input()))
# Input loop for diag tests if present
if diags:
    print("Please enter the results of the diagonal tests:")
    for i in range(n):
        print("Please enter the test result (0/1) for test D",i+1,": ",sep='',end="")
        diagTests.append(int(input()))
if diags:
    result = COMP(rowTests,colTests,n,diagTests)
else:
    result = COMP(rowTests,colTests,n)
# Add one to every index to get correct indiv number
for arr in result:
    arr += 1
print()
print("RESULTS:")
print("Definite Non-Defectives: ",result[0],sep='\n')
print("Definite Defectives: ",result[2])
print("Indeterminate Individuals: ",result[1])
print()
outputM = np.zeros((n,n))
# Fill in output matrix from the positive individual list using some math
print("Definite Defective Matrix: ")
for item in result[2]:
    outputM[int(np.floor((item-1)/n))][(item-1)%n] = 1

print("Output Testing Matrix:")
file = open("COMPRecoveryOutput.txt",'w')
for row in outputM:
    file.write(str(row)+'\n')
    print(row)
