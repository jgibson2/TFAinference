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
    # 100x speedup... WOW
    denominator = np.sum(np.square(data - np.repeat(np.mean(data, axis=1).reshape(-1, 1), data.shape[1], axis=1)))
    numValues = data.shape[0] * data.shape[1]
    numerator = np.sum(np.square(data - dataLearned))
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
def foldChange(A, B, l=0.0, method='max'):
    # calculate the fold-change between A and B
    ab = np.multiply(np.array(A) + l, 1.0 / (np.array(B) + l)) # - 1.0
    # calculate the fold-change between B and A
    ba = np.multiply(np.array(B) + l, 1.0 / (np.array(A) + l)) # - 1.0
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
    e = np.inf
    if method == 'max':
        e = max(map(lambda x: max(abs(np.nanmax(x)), abs(np.nanmin(x))), (ab, ba)))
    elif method == 'mean':
        e = np.nanmean(np.fmax(np.absolute(ab), np.absolute(ba)))
    elif method == 'median':
        e = np.nanmedian(np.fmax(np.absolute(ab), np.absolute(ba)))
    else:
        raise ValueError('Unrecognized method!')
    return e


"""
Input:
    matrix A, matrix B
    l, optional pseudocount value to add to each parameter when calculating fold change
Output:
    maximum log2 fold change across all parameters in matricies A and B
"""
def logFoldChange(A, B, l=0.0, method='max'):
    # calculate the fold-change between A and B
    # a = np.log2(np.array(A) + l)
    # calculate the fold-change between B and A
    # b = np.log2(np.array(B) + l)
    log_matrix = np.fmax(np.log(A + l) - np.log(B + l), np.log(B + l) - np.log(A + l))
    e = np.inf
    if method == 'max':
        e = np.nanmax(log_matrix)
    elif method == 'mean':
        e = np.nanmean(log_matrix)
    elif method == 'median':
        e = np.nanmedian(log_matrix)
    else:
        raise ValueError('Unrecognized method!')
    return e


"""
Calculates finite differences for each parameter for each gene to get significance of each parameter for weighing parameters

"""
def finiteDifferences(cs, tfa, data, epsilon=0.001):
    cs = np.array(cs, dtype=np.float64)
    tfa = np.array(tfa, dtype=np.float64)
    error = calcError(data, np.dot(cs, tfa), False)[0]
    csDerivatives = np.zeros(shape=cs.shape)
    tfaDerivatives = np.zeros(shape=tfa.shape)
    (i, j) = cs.shape
    (j, k) = tfa.shape
    for a in range(0, i):
        for b in range(0, j):
            mask = np.zeros(shape=cs.shape, dtype=np.float64)
            mask[a][b] += epsilon
            ######### BEGIN QUANTIFICATION METHOD HERE ###############
            # first derivative
            csDerivatives[a][b] = (calcError(data, np.dot(cs + mask, tfa), False)[0] -
                                    calcError(data, np.dot(cs - mask, tfa), False)[0]) / (2 * epsilon)
            # second derivative
            # csDerivatives[a][b] = ((calcError(data, np.dot(cs + mask, tfa), False)[0] - error) -
            #                         (calcError(data, np.dot(cs - mask, tfa), False)[0]) - error) / (2 * epsilon)
            ######### END QUANTIFICATION METHOD ######################
    for a in range(0, j):
        for b in range(0, k):
            mask = np.zeros(shape=tfa.shape, dtype=np.float64)
            mask[a][b] += epsilon
            derivs = (np.dot(cs, tfa + mask) - np.dot(cs, tfa - mask)) / (2*epsilon)
            ######### BEGIN QUANTIFICATION METHOD HERE ###############
            # first derivative
            tfaDerivatives[a][b] = (calcError(data, np.dot(cs, tfa + mask), False)[0] -
                                    calcError(data, np.dot(cs, tfa - mask), False)[0]) / (2*epsilon)
            # second derivative
            # tfaDerivatives[a][b] = ((calcError(data, np.dot(cs, tfa + mask), False)[0] - error) -
            #                         (calcError(data, np.dot(cs, tfa - mask), False)[0]) - error) / (2 * epsilon)
            ######### END QUANTIFICATION METHOD ######################
    return csDerivatives, tfaDerivatives


"""
Computes pseudocount for CS matrix (10% of average value)
"""

def computeCSPseudocount(cs):
    l = 0.1 * np.nanmean(cs, axis=0)
    return np.repeat(l.reshape(1, -1), cs.shape[0], axis=0)


'''
Returns boolean which states whether two CS matricies are converging
'''
def csConverges(A, B, absolute_limit=0.01, log_limit=1):
    converges = np.ones(A.shape,dtype=np.bool)
    # compute the pseudocount for each TF in the CS matrix
    l = computeCSPseudocount(A)
    # get the sign matrix of the first matrix
    sign_mat = np.sign(A)
    # check that we have the same signs
    same_signs = np.equal(sign_mat, np.sign(B))
    it = np.nditer(same_signs, flags=['multi_index'])
    while not it.finished:
        # print('Here!' + str(it.multi_index), flush=True)
        if it[0]:
            l[it.multi_index] *= sign_mat[it.multi_index]
        else: # False
            if abs(A[it.multi_index] - B[it.multi_index]) <= absolute_limit:
                pass
            else:
                converges[it.multi_index] = False
        it.iternext()
    # # calculate the fold-change between A and B
    # ab = np.multiply(np.array(A) + l, 1.0 / (np.array(B) + l))
    # # calculate the fold-change between B and A
    # ba = np.multiply(np.array(B) + l, 1.0 / (np.array(A) + l))
    # # calculate the log of the fold change
    # d = np.log(np.fmax(np.absolute(ab), np.absolute(ba)))
    # this is equivalent lol
    d = np.fmax(np.log(A + l) - np.log(B + l), np.log(B + l) - np.log(A + l))
    it = np.nditer(d, flags=['multi_index'])
    while not it.finished:
        # print('Here2!' + str(it.multi_index), flush=True)
        if it[0] >= log_limit:
            converges[it.multi_index] = False
        it.iternext()
    # basically check all the criteria for each i,j and generate the convergence conditional matrix,
    # then check if they all converge -- inefficient? yes, we could just return False at the first one -- this could be
    # an optimization. However, just in case we want to check which parameters aren't converging, this is how to do it
    # so I've left it here (for now) [I have a feeling we'll use the finite differences to check which are significant params]
    return converges.all()

'''
Similar to csConverges, check if a TFA matrix converges
However, we already know that the activities mean to 1, so we take 10% (0.1) as our pseudocount
In addition, all activities have a sign constraint -- they are positive
'''
def tfaConverges(A, B, l=0.1, absolute_limit=0.01, log_limit=1):
    converges = np.ones(A.shape, dtype=np.bool)
    d = np.fmax(np.log(A + l) - np.log(B + l), np.log(B + l) - np.log(A + l))

    it = np.nditer(A, flags=['multi_index'])
    while not it.finished:
        # print('Here!' + str(it.multi_index), flush=True)
        if abs(A[it.multi_index] - B[it.multi_index]) <= absolute_limit:
            pass
        else:  # False
            if d[it.multi_index] <= log_limit:
                pass
            else:
                converges[it.multi_index] = False
        it.iternext()
    return converges.all()
