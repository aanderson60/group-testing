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
def positiveTestInput(filename):
    with open(filename,"r") as file:
        tests = [int(x) for x in file.readline().split()]
    tests = np.array(tests)
    print(tests)
    return tests


def COMP(sensingMatrix, tests, n):
    

positiveTestInput("testresults.txt")