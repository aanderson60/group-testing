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
    n_array = np.arange(0, I).reshape(n, n)
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

    # Unexplained tests
    # These unexplained items may help indicate where false negatives are (from lab/human error)
    # However, many of these will be unnecessary in working trials
    # These are only applicable if the number of tests are incorrect (i.e. false negative lab error) - how to check for this??

    posIndivs = []
    unexplainedIndivs = []

    # These loops basically create a list of individuals in positive tests, and if they appear in more than
    #   one positive test, flags them as unexplained - this is kind of a variation on the SCOMP algorithm [Aldridge et al.]
    for i,item in enumerate(rowTests):
        if item == 1:
            for j in range(len(rowTests)):
                if i*n+j not in posIndivs:
                    posIndivs.append(i*n+j)
                else:
                    if i*n+j not in DD and i*n+j not in PD:
                        unexplainedIndivs.append(i*n+j)
    for i,item in enumerate(colTests):
        if item == 1:
            for j in range(len(colTests)):
                if j*n+i not in posIndivs:
                    posIndivs.append(j*n+i)
                else:
                    if j*n+i not in DD and j*n+i not in PD:
                        unexplainedIndivs.append(j*n+i)
    if diagTests != None:
        for i,item in enumerate(diagTests):
            if item == 1:
                diag = np.concatenate((np.diagonal(n_array,i),np.diagonal(n_array,-(n-i))))
                for item in diag:
                    if item not in posIndivs:
                        posIndivs.append(item)  
                    else:
                        if item not in DD and item not in PD:
                            unexplainedIndivs.append(item)

    # Convert to Numpy arrays and sort so that we can manipulate these later
    DND = np.array(DND)
    PD = np.array(PD)
    DD = np.array(DD)
    unexplainedIndivs = np.array(unexplainedIndivs)

    np.sort(DND)
    np.sort(PD)
    np.sort(DD)
    np.sort(unexplainedIndivs)

    return((DND,PD,DD,unexplainedIndivs))

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

print(rowTests,colTests,diagTests)
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
print("Unexplained Individuals (PD): ",result[3])
print()
outputM = np.chararray((n,n),unicode=True)
outputM[:][:] = '0'
# Fill in output matrix from the positive individual list using some math
# Also fill in unexplained individuals with a '?' (note that this may or may not indicate a positive sample (see above))
print("Definite Defective Matrix: ")
for item in result[2]:
    outputM[int(np.floor((item-1)/n))][(item-1)%n] = 1
for item in result[3]:
    outputM[int(np.floor((item-1)/n))][(item-1)%n] = '?'

print("Output Testing Matrix:")
file = open("COMPRecoveryOutput.txt",'w')
for row in outputM:
    file.write(str(row)+'\n')
    print("[ ",end='')
    for item in row:
        print(item+' ',end='')
    print("]",end='\n')
