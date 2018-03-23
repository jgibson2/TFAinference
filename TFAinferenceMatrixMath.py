from gurobipy import *
import numpy as np

"""
matrix multiplication
Someday I will import numpy and this won't be used anymore ...

UNUSED 2/27/18 JG
"""


def matrixMultiply(C, A):
    if len(C[0]) == len(A):
        result = []
        for i in range(len(C)):
            result.append([])
            for j in range(len(A[0])):
                val = 0
                for k in range(len(A)):
                    val += C[i][k] * A[k][j]
                result[i].append(val)
        return result
    else:
        print("matrices are the wrong dimensions for multiplication")
        print(len(C[0]), len(A))
        return []


"""
Input:
  activity matrix, and a binary matrix of factor labels (factor x sample matrix)
  a list of lists of factor values
Output:
  activity matrix with factor values applied

This is for datasets such that certain samples can be grouped as having some factor
such that the factor is expected to have a consistent effect on TF activity
Currently not used anywhere
"""


def factorMultiplyA(A, F, fVals):
    newA = [x[:] for x in A]
    for i in range(len(fVals[0])):  # TF
        for j in range(len(A[i])):  # sample
            for fIndex in range(len(F)):  # factor
                if F[fIndex][j] == 1:
                    newA[i][j] *= fVals[fIndex][i]
    return newA


"""
Input:
  two matrices, where the first matrices #col = #row of second
Output:
  two rescaled matrices, such that the non-zero values of each row of the second is 1

one possible way to normalize the learned CS and activity matrices
currently not used anywhere
"""


def rescale(C, A):
    xMatrix = []
    xMatrixInverse = []
    for i in range(len(A)):
        xMatrix.append([0] * len(A))
        xMatrixInverse.append([0] * len(A))
        xMatrix[i][i] = 1.0 / (sum(abs(x) for x in A[i]) / len([x for x in A[i] if x != 0.0]))
        xMatrixInverse[i][i] = (sum(abs(x) for x in A[i]) / len([x for x in A[i] if x != 0.0]))

    Crescaled = np.multiply(C, xMatrixInverse)
    Arescaled = np.multiply(xMatrix, A)
    return [Crescaled, Arescaled]


"""
Input:
  two data matrices
  the first is assumed to be true data
  the second is assumed to be predicted/learned/fitted data
Output:
  list of var explained and SSE

calculates total variance of data, SSE of dataLearned-data
and variance explained of dataLearned
"""


def calcError(data, dataLearned, printFlag):
    numerator = 0
    denominator = 0
    error = 0
    numValues = 0
    for i in range(len(dataLearned)):
        # geneMean = sum(data[i]) / len(data[i])
        # for j in range(len(dataLearned[i])):
        #     numerator += (data[i][j] - dataLearned[i][j]) * (data[i][j] - dataLearned[i][j])
        #     denominator += (data[i][j] - geneMean) * (data[i][j] - geneMean)
        #     numValues += 1
        geneMean = np.mean(data[i]) # temp value to calculate error
        numerator += sum(np.square(np.subtract(data[i], dataLearned[i])))
        denominator += sum(np.square(np.apply_along_axis(lambda x: x-geneMean, 0, data[i])))
        numValues += len(data[i])
    if printFlag:
        print("variance: ", denominator / numValues)
        print("variance explained: ", 1 - (numerator / denominator))
        print("Sum of Squared Errors: ", numerator)
    return 1 - (numerator / denominator), numerator


"""
Input:
    matrix A, matrix B
    l, optional pseudocount value to add to each parameter when calculating fold change
Output:
    maximum absolute fold change across all parameters in matricies A and B
"""
def foldChange(A, B, l=0):
    # calculate the fold-change between A and B
    ab = np.multiply(np.array(A) + l, 1.0 / (np.array(B) + l)) - 1.0
    # calculate the fold-change between B and A
    ba = np.multiply(np.array(B) + l, 1.0 / (np.array(A) + l)) - 1.0
    # """
    # Here, we calculate the fold-change of each element between the matricies. This is currently unnecessary
    # but might be helpful in the future
    # d = np.fmax(np.absolute(ab), np.absolute(ba))

    # opt = np.get_printoptions()
    # np.set_printoptions(threshold=np.inf)
    # print(d)
    # np.set_printoptions(**opt)

    # """
    # calculate maximum absolute fold change
    # FIX nan RESULTS
    e = max(map(lambda x: max(abs(np.nanmax(x)), abs(np.nanmin(x))), (ab, ba)))

    # print(list(map(lambda x: max(abs(np.nanmax(x)), abs(np.nanmin(x))), (ab, ba))))

    return e


