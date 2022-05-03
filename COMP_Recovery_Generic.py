import numpy as np

# Function to decode sensing matrix information from file
# Args:
#   filename: the string filename of the input file (including filetype)
#       Format of file: first line contains a single integer, n, the number of individuals
#                       remaining n lines contain a zero or one seperated by spaces to denote which tests each individual will participate in
# Returns:
#   The sensing matrix, a numpy array, and the integer n, number of individuals
#
def sensingMatrixInput(filename):
    sensingMatrix = []
    with open(filename,"r") as file:
        n = int(file.readline())
        for line in file:
            row = [int(x) for x in line.split()]
            sensingMatrix.append(row)
    sensingMatrix = np.array(sensingMatrix)
    return (sensingMatrix, n)


# Function to decode test result input from a file
# Args:
#   filename: the string filename of the input file (including filetype)
#       Format of the file: a single line containing the binary test results (1 or 0) seperated by a space
#       The number of test results should equal the number of columns in the sensing matrix
# Returns:
#   The array of test results, length equal to the number of tests
#
def testInput(filename):
    with open(filename,"r") as file:
        tests = [int(x) for x in file.readline().split()]
    tests = np.array(tests)
    return tests

# Function to apply COMP algorithm recovery to sensing matrix and tests
# Args:
#   sensingMatrix: the sensing matrix of the scheme
#   tests: the pooled test results, corresponding to the sensing matrix order
#   n: the number of items
# Returns:
#   The probable defectives, the definite non-defectives
#
def COMP(sensingMatrix, n, tests):
    # The probable defectives
    PDs = []
    # The definite non defectives
    DNDs = []

    # First determine DNDs
    for i,result in enumerate(tests):
        # Test is negative
        if result == 0:
            # Check sensing matrix for which items participated in this test
            for j,row in enumerate(sensingMatrix):
                if row[i] == 1 and j not in DNDs:
                    DNDs.append(j)
    
    # Now determine the PDs
    for x in range(n):
        if not (x in DNDs):
            PDs.append(x)

    # Cast arrays to numpy arrays
    PDs = np.array(PDs)
    DNDs = np.array(DNDs)

    # Add 1 to every item number to normalize them to start at 1 instead of 0
    PDs += 1
    DNDs += 1

    return (PDs,DNDs)


# --------- MAIN ----------
filename1 = input("Enter the filename for the sensing matrix: ")
filename2 = input("Enter the filename for the pooled test results: ")
tests = testInput(filename2)
sm, n = sensingMatrixInput(filename1)
PD, DND = COMP(sm,n,tests)
print("Definite Non Defectives:")
print(DND)
print("Probable Defectives:")
print(PD)