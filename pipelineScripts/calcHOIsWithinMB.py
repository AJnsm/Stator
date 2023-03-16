PrintBool=False

# Only print if in the __main__ call of the script
if(__name__=="__main__"):
    PrintBool=True

if PrintBool: print('Importing modules...')
import concurrent.futures
import numpy as np
import pandas as pd
import igraph as ig
import time
import sys
import argparse
import numba
from numba import njit
import scipy

# The utilities module holds the estimation methods etc. 
from utilities import *

if PrintBool: print('Modules imported \n')

# Parse all command line arguments
parser = argparse.ArgumentParser(description='Args for coupling estimation')
parser.add_argument("--dataPath", type=str, help="Path to training data")
parser.add_argument("--graphPath", type=str, help="Path to graph file")
parser.add_argument("--nResamps", type=int, help="Number of BS resamples")
parser.add_argument("--nCores", type=int, help="Number of cores")
parser.add_argument("--nRandoms", type=int, help="Number of random interactions to calculate")
parser.add_argument("--genesToOne", type=str, help="Path to list of genes that should be set to 1")
parser.add_argument("--dataDups", type=int, help="Number of data duplications. 0 is no duplication, and another value is the min binsize allowed (recommended to be 15). ")
parser.add_argument("--boundBool", type=int, help="Boolean that decided whether bounds should also be considered.")
parser.add_argument("--asympBool", type=int, help="Boolean to decide whether to use Bootstrap resampling (0) or asymptotic uncertainty estimation (1).")


args = parser.parse_args()
    
dataPath = args.dataPath
graphPath = args.graphPath
nResamps = args.nResamps
nCores = args.nCores
nRands = args.nRandoms
genesToOnePath = args.genesToOne
dataDups = args.dataDups
boundBool = args.boundBool
asympBool = args.asympBool

trainDat = pd.read_csv(dataPath)

# DSname copies the naming scheme from the graphs.
DSname = graphPath.split('.')[0]
adjMat = pd.read_csv(graphPath, index_col=0)
graph = ig.Graph.Adjacency(adjMat.values.tolist()) 


try:
    genesToOneIndices = pd.read_csv(genesToOnePath)
except:
    if PrintBool: print('NOTE: all genes conditioned on 0s.')
    genesToOneIndices = []


if PrintBool: print('data import and graph construction done')

def calcInteractionsAndWriteNPYs(ID, graph, trainDat, maxWorkers, estimator, nResamps=1000):
    
    if PrintBool: print(f'Starting with {ID}...')
    genes = trainDat.columns
    n = len(genes)

    # First, generate random 3-, 4-, and 5- tuples
    print('Generating random pairs...')
    randPairs = np.array([np.random.choice(np.arange(n), 2, replace=False) for i in range(nRands)]).astype(int)
    args_randPairs = [(pair, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for pair in randPairs]

    print('Generating random triplets...')
    randTrips = np.array([np.random.choice(np.arange(n), 3, replace=False) for i in range(nRands)]).astype(int)
    args_randTrips = [(triplet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for triplet in randTrips]

    print('Generating random quads...')
    randQuads = np.array([np.random.choice(np.arange(n), 4, replace=False) for i in range(nRands)]).astype(int)
    args_randQuads = [(quad, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for quad in randQuads]

    print('Generating random pents...')
    randPents = np.array([np.random.choice(np.arange(n), 5, replace=False) for i in range(nRands)]).astype(int)
    args_randPents = [(pent, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for pent in randPents]


    ints_2pt = []
    ints_3pt = []
    ints_4pt = []
    ints_5pt = []
    
    # Then iterate over Markov blankets and take intersections to have fully Markov-connected tuples
    print('Generating connected tuples...')
    for g1 in range(n):
        MB1 = findMarkovBlanket(g1, graph)

        ints_2pt.append([(g1, x) for x in MB1])

        for g2 in MB1:
            MB2 = findMarkovBlanket(g2, graph)
            MB1_MB2 = set(MB1).intersection(set(MB2))

            ints_3pt.append([(g1, g2, x) for x in MB1_MB2])

            for g3 in MB1_MB2:
                MB3 = findMarkovBlanket(g3, graph)
                MB1_MB2_MB3 = set(MB3).intersection(MB1_MB2)

                ints_4pt.append([(g1, g2, g3, x) for x in MB1_MB2_MB3])

                for g4 in MB1_MB2_MB3:
                    MB4 = findMarkovBlanket(g4, graph)
                    MB1_MB2_MB3_MB4 = set(MB4).intersection(MB1_MB2_MB3)

                    ints_5pt.append([(g1, g2, g3, g4, x) for x in MB1_MB2_MB3_MB4])

    print('Generated all connected 3-, 4-, 5-tuples')

    # To aid estimation, order the variables by the size of their Markov blanket so that the estiamtion uses the smallest one. 
    def onlySmallestMB(ar):
        ar = [tuple(sorted(genes)) for intList in ar for genes in intList]
        ar = np.unique(ar, axis=0)
        ar = np.array([sorted(genes, key=lambda x: len(findMarkovBlanket(x, graph))) for genes in ar])
        return ar

    ints_2pt = onlySmallestMB(ints_2pt)
    ints_3pt = onlySmallestMB(ints_3pt)
    ints_4pt = onlySmallestMB(ints_4pt)
    ints_5pt = onlySmallestMB(ints_5pt)

    ints_2pt = [(intSet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for intSet in ints_2pt]
    ints_3pt = [(intSet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for intSet in ints_3pt]
    ints_4pt = [(intSet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for intSet in ints_4pt]
    ints_5pt = [(intSet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for intSet in ints_5pt]

    if PrintBool:
        print(f'Markov-connected pairs: {len(ints_2pt)}, trips: {len(ints_3pt)}, Quads: {len(ints_4pt)}, Pents: {len(ints_5pt)}')
    

    for order, args in [['random_2pts', args_randPairs], ['random_3pts', args_randTrips], ['random_4pts', args_randQuads], ['random_5pts', args_randPents], ['withinMB_2pts', ints_2pt], ['withinMB_3pts', ints_3pt], ['withinMB_4pts', ints_4pt], ['withinMB_5pts', ints_5pt]]:

        
        start = time.perf_counter()
        with concurrent.futures.ProcessPoolExecutor(max_workers=maxWorkers) as executor:
            results = executor.map(calcInteraction_withCI_parallel  , args)  
        finish = time.perf_counter()
        if PrintBool: print(f'Time elapsed: {round(finish-start, 2)} secs')
        if PrintBool: print('calculation done, storing results...')

        resultArr = np.array(list(results), dtype=object)
        if PrintBool: print('writing files...')
        
        np.save(f'interactions_{order}_{ID}', resultArr, allow_pickle=True)

        if PrintBool: print(f'********** DONE with {order} {ID} **********\n')



def main():
    np.random.seed(0)
    notes = ''
    
    estimationMethod = 'expectations'
    estimator = calcInteraction_expectations_numba
    
    
    print('Starting calculation on ' )
    print('Using estimation method:  ', estimationMethod)
    print(f'With {nResamps} bootstrap resamples')
    print(f'Parallelised over {nCores} cores. ')
    print(f'Asymptotic variance estimation: {bool(asympBool)}')

    calcInteractionsAndWriteNPYs(DSname + '_' + estimationMethod+notes, graph, trainDat, maxWorkers=nCores, estimator = estimator, nResamps=nResamps)
    
    
    print('***********DONE***********')

if __name__=="__main__":
    main()
    


