import argparse
from TFAinferenceIO import *
from TFAinferenceMatrixMath import *
import numpy as np


def read_all_matricies(filename, iterations):
    data = readMatrixFromFile(filename)
    if not (len(data) / iterations) % 1 == 0.0:
        raise ValueError('Data length does not match iterations!')
    mat_length = int(len(data) / iterations)
    return np.array([data[mat_length*i:mat_length*(i+1)] for i in range(iterations)])


def check_cs_convergence(cs_mats, in_a_row=3, log_limit=0.1):
    iar = 0
    for i in range(0, len(cs_mats) - 1):
        if csConverges(cs_mats[i], cs_mats[i+1], log_limit=log_limit):
            iar += 1
            if iar >= in_a_row:
                return True, i+1
    return False, -1


def check_geneexpr_convergence(mats, in_a_row=3, log_limit=0.1):
    iar = 0
    for i in range(0, len(mats) - 1):
        if csConverges(mats[i], mats[i+1], log_limit=log_limit): # just use CS convergence criteria
            iar += 1
            if iar >= in_a_row:
                return True, i+1
    return False, -1


def check_tfa_convergence(tfa_mats, in_a_row=3, log_limit=0.1):
    iar = 0
    for i in range(0, len(tfa_mats) - 1):
        if tfaConverges(tfa_mats[i], tfa_mats[i+1], log_limit=log_limit):
            iar += 1
            if iar >= in_a_row:
                return True, i+1
    return False, -1


def main():
    cs_matricies = list()
    cs_names = dict()
    tfa_matricies = list()
    tfa_names = dict()
    varsExplained = list()
    cs_convergence_points = list()
    tfa_convergence_points = list()

    parser = argparse.ArgumentParser(description='Process matricies to test convergence')
    parser.add_argument('--cs', '-c', type=str, nargs='+', help='List of control strength matrix files', required=True)
    parser.add_argument('--tfa', '-t', type=str, nargs='+', help='List of transcription factor activity matricies', required=True)
    parser.add_argument('--variances', '-v', type=str, nargs='+', help='List of variance explained files', required=True)
    parser.add_argument('--iterations', '-i', type=int, action='store', default=100, help='Number of iterations in log files')
    parser.add_argument('--breaks', '-b', type=int, action='store', default=20, help='Number of iterations to consider at a time')
    parser.add_argument('--std_dev', action='store', type=str, help='Calculate standard deviations matrix for final iterations and write to file')
    parser.add_argument('--log_limit', action='store', type=int, default=0.5, help='Log limit')
    args = parser.parse_args()
    for index, file in enumerate(args.cs):
        cs_names[index] = file
        cs_matricies.append(read_all_matricies(file, args.iterations))
        cs_convergence_points.append(-1)
    for index, file in enumerate(args.tfa):
        tfa_names[index] = file
        tfa_matricies.append(read_all_matricies(file, args.iterations))
        tfa_convergence_points.append(-1)
    for file in args.variances:
        vars_tmp = []
        with open(file) as f:
            for line in f:
                line = line.strip()
                vars_tmp.append(float(line))
        varsExplained.append(vars_tmp)
    geneexpr_matricies = [[np.dot(C,A) for C,A in zip(cs, tfa)] for cs,tfa in zip(cs_matricies, tfa_matricies)]
    geneexpr_convergence_points = [-1 for i in range(len(geneexpr_matricies))]

    for i in np.arange(0, args.iterations, args.breaks):
        cs = [x[i:i+args.breaks] for x in cs_matricies]
        tfa = [x[i:i+args.breaks] for x in tfa_matricies]
        var = [x[i:i+args.breaks] for x in varsExplained]
        for index, (c, t) in enumerate(zip(cs, tfa)):
            if cs_convergence_points[index] == -1:
                cs_converges, cs_conv_pt = check_cs_convergence(c, log_limit=args.log_limit)
                if cs_converges:
                    cs_convergence_points[index] = cs_conv_pt + i
            if tfa_convergence_points[index] == -1:
                tfa_converges, tfa_conv_pt = check_tfa_convergence(t, log_limit=args.log_limit)
                if tfa_converges:
                    tfa_convergence_points[index] = tfa_conv_pt + i
            if geneexpr_convergence_points[index] == -1:
                geneexpr_converges, geneexpr_conv_pt = check_geneexpr_convergence(t, log_limit=args.log_limit)
                if geneexpr_converges:
                    geneexpr_convergence_points[index] = geneexpr_conv_pt + i

    cs_convergences = [[] for i in range(len(cs_matricies))]
    tfa_convergences = [[] for i in range(len(tfa_matricies))]
    geneexpr_convergences = [[] for i in range(len(geneexpr_matricies))]

    for index,conv_pt in enumerate(cs_convergence_points):
        if conv_pt != -1:
            for index2, conv_pt2 in enumerate(cs_convergence_points):
                if index2 != index:
                    if csConverges(cs_matricies[index][conv_pt], cs_matricies[index2][conv_pt2], log_limit=args.log_limit):
                        cs_convergences[index].append(index2)
                        
    for index,conv_pt in enumerate(tfa_convergence_points):
        if conv_pt != -1:
            for index2, conv_pt2 in enumerate(tfa_convergence_points):
                if index2 != index:
                    if tfaConverges(tfa_matricies[index][conv_pt], tfa_matricies[index2][conv_pt2], log_limit=args.log_limit):
                        tfa_convergences[index].append(index2)
                        
    for index,conv_pt in enumerate(geneexpr_convergence_points):
        if conv_pt != -1:
            for index2, conv_pt2 in enumerate(geneexpr_convergence_points):
                if index2 != index:
                    if csConverges(geneexpr_matricies[index][conv_pt], geneexpr_matricies[index2][conv_pt2], log_limit=args.log_limit): #use CS criteria
                        geneexpr_convergences[index].append(index2)

    for index, other_matricies in enumerate(cs_convergences):
        print('{0}-th matrix from {1} converged to same solution as: '.format(cs_convergence_points[index], 
                                                                              cs_names[index]), end='')
        print(' '.join(
            ['{0}-th matrix from {1}'.format(cs_convergence_points[x], cs_names[x]) for x in other_matricies]))

    for index, other_matricies in enumerate(tfa_convergences):
        print('{0}-th matrix from {1} converged to same solution as: '.format(tfa_convergence_points[index],
                                                                              tfa_names[index]), end='')
        print(' '.join(
            ['{0}-th matrix from {1}'.format(tfa_convergence_points[x], tfa_names[x]) for x in other_matricies]))
    
    for index, other_matricies in enumerate(geneexpr_convergences):
        print('{0}-th matrix from {1} converged to same solution as: '.format(geneexpr_convergence_points[index],
                                                                              '{0}*{1}'.format(cs_names[index], tfa_names[index])), end='')
        print(' '.join(
            ['{0}-th matrix from {1}'.format(geneexpr_convergence_points[x], '{0}*{1}'.format(cs_names[index], tfa_names[index])) for x in other_matricies]))

    if args.std_dev:
        final_iter_matricies = [g[geneexpr_convergence_points[idx]] for idx,g in enumerate(geneexpr_matricies)]
        M = np.nanstd(np.array(final_iter_matricies), axis=0) # clever hack to form a matrix of std. devs of each parameter
        outfile = open(args.std_dev, "w")
        for row in M:
            toSave = str(row[0])
            for valIndex in range(1, len(row)):
                toSave += "," + str(row[valIndex])
            toSave += "\n"
            outfile.write(toSave)
        outfile.close()


if __name__ == '__main__':
    main()
