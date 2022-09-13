'''
These are some very basic tests in Ising data. 
I verify that the 1- and 2-point interactions are approximately correct, using the MCMC graph, and the true underlying graph of dependencies.
run `sh runTests.sh` before running `pytest integration_tests.py`. 
'''

import numpy as np
import pandas as pd
import pytest

def test_onePts_MCMC():
    '''
    Makes sure that most estimable 1-points are significant (they should be in the {0, 1} basis.)
    '''
    c1_CI_F = np.load('interactions_order1_MCMCgraph_CL01_100000Cells_0009Genes_expectations_CI_F.npy', allow_pickle=True).astype(np.float32)
    c1_inf = np.load('interactions_order1_MCMCgraph_CL01_100000Cells_0009Genes_expectations_inf.npy', allow_pickle=True).astype(np.float32)
    c1_undef = np.load('interactions_order1_MCMCgraph_CL01_100000Cells_0009Genes_expectations_undef.npy', allow_pickle=True).astype(np.float32)

    mask = (~np.isnan(c1_CI_F)) & (c1_inf==0) & (c1_undef==0)
    print(np.mean(c1_CI_F[mask]<0.1))

    assert(np.mean(c1_CI_F[mask]<0.1)>0.5)


def test_twoPts_MCMC():
    '''
    Makes sure that the F1 score of recognising the nearest neighbour structure is better than random.
    '''
    c2_CI_F = np.load('interactions_order2_MCMCgraph_CL01_100000Cells_0009Genes_expectations_CI_F.npy', allow_pickle=True).astype(np.float32)
    c2_inf = np.load('interactions_order2_MCMCgraph_CL01_100000Cells_0009Genes_expectations_inf.npy', allow_pickle=True).astype(np.float32)
    c2_undef = np.load('interactions_order2_MCMCgraph_CL01_100000Cells_0009Genes_expectations_undef.npy', allow_pickle=True).astype(np.float32)

    mask = (~np.isnan(c2_CI_F)) & (c2_inf==0) & (c2_undef==0)

    groundTruth = pd.read_csv('isingTestData_trueAdjMat.csv', index_col=0).values

    c2_CI_F_shuf = c2_CI_F.copy()
    F1 = []
    for i in range(1000):
        np.random.shuffle(c2_CI_F_shuf.flat)
        TP = np.sum((c2_CI_F_shuf<=0.1) & groundTruth & mask)
        FP = np.sum((c2_CI_F_shuf<=0.1) & ~groundTruth & mask)
        FN = np.sum((c2_CI_F_shuf>0.1) & groundTruth & mask)
        TN = np.sum((c2_CI_F_shuf>0.1) & ~groundTruth & mask)
        F1.append(2*TP/(2*TP + FP + FN))

    TP = np.sum((c2_CI_F<=0.1) & groundTruth & mask)
    FP = np.sum((c2_CI_F<=0.1) & ~groundTruth & mask)
    FN = np.sum((c2_CI_F>0.1) & groundTruth & mask)
    TN = np.sum((c2_CI_F>0.1) & ~groundTruth & mask)

    print(2*TP/(2*TP + FP + FN), np.array(F1).mean())
    assert(2*TP/(2*TP + FP + FN) > np.array(F1).mean())



# ******************* Using true adjacency graph *******************



def test_onePts_true():
    '''
    Makes sure that most estimable 1-points are significant (they should be in the {0, 1} basis.)
    '''
    c1_CI_F = np.load('interactions_order1_isingTestData_trueAdjMat_expectations_CI_F.npy', allow_pickle=True).astype(np.float32)
    c1_inf = np.load('interactions_order1_isingTestData_trueAdjMat_expectations_inf.npy', allow_pickle=True).astype(np.float32)
    c1_undef = np.load('interactions_order1_isingTestData_trueAdjMat_expectations_undef.npy', allow_pickle=True).astype(np.float32)

    mask = (~np.isnan(c1_CI_F)) & (c1_inf==0) & (c1_undef==0)
    print(np.mean(c1_CI_F[mask]<0.1))

    assert(np.mean(c1_CI_F[mask]<0.1)>0.5)


def test_twoPts_true():
    '''
    Makes sure that the F1 score of recognising the nearest neighbour structure is better than random.
    '''
    c2_CI_F = np.load('interactions_order2_isingTestData_trueAdjMat_expectations_CI_F.npy', allow_pickle=True).astype(np.float32)
    c2_inf = np.load('interactions_order2_isingTestData_trueAdjMat_expectations_inf.npy', allow_pickle=True).astype(np.float32)
    c2_undef = np.load('interactions_order2_isingTestData_trueAdjMat_expectations_undef.npy', allow_pickle=True).astype(np.float32)

    mask = (~np.isnan(c2_CI_F)) & (c2_inf==0) & (c2_undef==0)

    groundTruth = pd.read_csv('isingTestData_trueAdjMat.csv', index_col=0).values

    c2_CI_F_shuf = c2_CI_F.copy()
    F1 = []
    for i in range(1000):
        np.random.shuffle(c2_CI_F_shuf.flat)
        TP = np.sum((c2_CI_F_shuf<=0.1) & groundTruth & mask)
        FP = np.sum((c2_CI_F_shuf<=0.1) & ~groundTruth & mask)
        FN = np.sum((c2_CI_F_shuf>0.1) & groundTruth & mask)
        TN = np.sum((c2_CI_F_shuf>0.1) & ~groundTruth & mask)
        F1.append(2*TP/(2*TP + FP + FN))

    TP = np.sum((c2_CI_F<=0.1) & groundTruth & mask)
    FP = np.sum((c2_CI_F<=0.1) & ~groundTruth & mask)
    FN = np.sum((c2_CI_F>0.1) & groundTruth & mask)
    TN = np.sum((c2_CI_F>0.1) & ~groundTruth & mask)

    print(2*TP/(2*TP + FP + FN), np.array(F1).mean())
    assert(2*TP/(2*TP + FP + FN) > np.array(F1).mean())



