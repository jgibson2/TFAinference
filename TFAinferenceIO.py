import random
import gurobipy
from gurobipy import *

"""
Input:
  three nested lists of lists of equal length (matrices)
  a single variance explained value
  four strings for filenames
Output:
  none

saves the input matrices and value into csv files according to the filenames
currently not in use
"""


def saveFactorResults(C, A, fVals, var, outFileName1, outFileName2, outFileName3, outFileName4):
    writeMatrixToFile(C, outFileName1)
    writeMatrixToFile(A, outFileName2)
    writeMatrixToFile(fVals, outFileName3)
    outfile = open(outFileName4, "a")
    outfile.write(str(var))
    outfile.write("\n")


"""
Input:
  two nested lists of lists of equal length (matrices)
  a single variance explained value
  three strings for filenames
Output:
  none

saves the input matrices and value into csv files according to the filenames
"""


def saveResults(C, A, var, outFileName1, outFileName2, outFileName3):
    writeMatrixToFile(C, outFileName1)
    writeMatrixToFile(A, outFileName2)
    outfile = open(outFileName3, "a")
    outfile.write(str(var))
    outfile.write("\n")


"""
Input:
  a nested list of lists of equal length (matrix)
  a string for file name
Output:
  none

saves the matrix to the end of the file in csv format
"""


def writeMatrixToFile(M, outFileName):
    outfile = open(outFileName, "a")
    for row in M:
        toSave = str(row[0])
        for valIndex in range(1, len(row)):
            toSave += "," + str(row[valIndex])
        toSave += "\n"
        outfile.write(toSave)
    outfile.close()


"""
Input:
  a nested lists of lists of equal length (matrix)
Output:
  none

prints the matrix to the screen
assumes that all values are gurobipy.Var or values with a default string conversion
"""


def matrixPrintOther(M):
    for row in M:
        print('{ ', end=' ')
        for val in row:
            if type(val) == gurobipy.Var:
                print("," + val.getAttr('VarName'), end=' ')
            else:
                print("     ," + str(val), end=' ')
        print(" }")


"""
Input:
  a nested lists of lists of equal length (matrix)
Output:
  none

prints the matrix to the screen
assumes that all values are numerical
"""


def matrixPrint(M):
    for row in M:
        toPrint = "{ "
        for val in row:
            toPrint += "{:<10}".format(str(round(val, 2))) + ","
        toPrint += " }"
        print(toPrint)


"""
Input:
  a string for filename
Output:
  a nested list of lists of numerical values

assumes that all values are floats
assumes that the file is csv
"""


def readMatrixFromFile(fileName):
    M = []
    fileMatrix = open(fileName, 'r')
    for line in fileMatrix:
        row = []
        values = line.strip().split(',')
        for val in values:
            row.append(float(val.strip(',')))
        M.append(row)
    fileMatrix.close()
    return M


"""
Input:
  two strings for filenames
  a boolean for cs sign constraint
Output:
  none

makes a random start point for cs matrix and saves as a csv file
"""


def makeRandStart(binaryFilename, outFilename, signFlag):
    C = readMatrixFromFile(binaryFilename)
    startingCS = []
    for i in range(len(C)):
        startingCS.append([])
        for j in range(len(C[i])):
            randVal = random.uniform(-10, 10)
            if signFlag:
                startingCS[i].append(abs(randVal) * C[i][j])
            else:
                startingCS[i].append(randVal * C[i][j])
        startingCS[i].append(random.uniform(0.01, 20))
    writeMatrixToFile(startingCS, outFilename)
