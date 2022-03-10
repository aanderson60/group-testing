import numpy as np
# Recovery algorithm using the COMP algorithm
# Optional parameter for diagonal tests
def COMP(rowTests,colTests,n,diagTests=None):
    # Use COMP algorithm to find set of defectives
    I = n*n
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
    if not diagTests == None:
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

# ---------- MAIN -----------
rowTests = []
colTests = []
diagTests = []
diags = False
print("Enter the total number of individuals: ",end="")
I = float(input())
n = np.sqrt(I)
if not n == int(n):
    print("ERROR: Number of individuals must be in a square matrix")
    exit()
n = int(n)
print("Does this testing scheme use diagonal tests? (Y/N): ",end="")
ans = input()
diags = True if ans=="Y" else False
if ans != "Y" and ans != "N":
    print("ERROR: Please enter Y/N only")
    exit()
print("COMP Testing Recovery for ",n,"x",n," matrix:",sep='')
print("")
print("Please enter the results of the row tests:")
for i in range(n):
    print("Please enter the test result (0/1) for test R",i+1,": ",sep='',end="")
    rowTests.append(int(input()))
print("Please enter the results of the column tests:")
for i in range(n):
    print("Please enter the test result (0/1) for test C",i+1,": ",sep='',end="")
    colTests.append(int(input()))
if diags:
    print("Please enter the results of the diagonal tests:")
    for i in range(n):
        print("Please enter the test result (0/1) for test D",i+1,": ",sep='',end="")
        diagTests.append(int(input()))
if diags:
    result = COMP(rowTests,colTests,n,diagTests)
else:
    result = COMP(rowTests,colTests,n)
result = np.array(result)
# Add one to every index to get correct indiv number
result += 1
print()
print("RESULTS:")
print("Positive individuals: ",result,sep='')
outputM = np.zeros((n,n))
for item in result:
    outputM[int(np.floor((item-1)/n))][(item-1)%n] = 1

print("Output Testing Matrix:")
for row in outputM:
    print(row)
