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
parser.add_argument("--estimationMode", type=str, help="Can be set to MFI to condition on Markov blankets, or to LOR to ignore the Markov blanket and calculate log-odds ratios.")


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
estimationMode = args.estimationMode

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

def calcInteractionsAndWriteNPYs(ID, graph, trainDat, maxWorkers, estimator, nResamps=1000, mode='MFI'):
    
    if PrintBool: print(f'Starting with {ID}...')

    if mode=='MFI':
        if PrintBool: print(f'Calculating MFIs, so conditioning on Markov blankets')
        estimationGraph = ig.Graph.Adjacency(adjMat.values.tolist())
    elif mode=='LOR':
        if PrintBool: print(f'Calculating LORs, so NOT conditioning on Markov blankets')
        estimationGraph = ig.Graph.Adjacency(np.zeros_like(adjMat).tolist())
    else:
        if PrintBool: print('Invalid mode, so switching to MFI mode')
        if PrintBool: print(f'Calculating MFIs, so conditioning on Markov blankets')
        estimationGraph = ig.Graph.Adjacency(adjMat.values.tolist())

        

    genes = trainDat.columns
    n = len(genes)

    # First, generate random tuples
    print('Generating random pairs...')
    randPairs = np.array([np.random.choice(np.arange(n), 2, replace=False) for i in range(nRands)]).astype(int)
    args_randPairs = [(pair, estimationGraph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for pair in randPairs]

    print('Generating random triplets...')
    randTrips = np.array([np.random.choice(np.arange(n), 3, replace=False) for i in range(nRands)]).astype(int)
    args_randTrips = [(triplet, estimationGraph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for triplet in randTrips]

    print('Generating random quads...')
    randQuads = np.array([np.random.choice(np.arange(n), 4, replace=False) for i in range(nRands)]).astype(int)
    args_randQuads = [(quad, estimationGraph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for quad in randQuads]

    print('Generating random pents...')
    randPents = np.array([np.random.choice(np.arange(n), 5, replace=False) for i in range(nRands)]).astype(int)
    args_randPents = [(pent, estimationGraph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for pent in randPents]

    # Then, generate the Markov-connected tuples
    connected_2pts = []
    connected_3pts = []
    connected_4pts = []
    connected_5pts = []
    
    # Iterate over Markov blankets and take intersections to have fully Markov-connected tuples
    print('Generating connected tuples...')
    for g1 in range(n):
        MB1 = findMarkovBlanket(g1, graph)

        connected_2pts.append([(g1, x) for x in MB1])

        for g2 in MB1:
            MB2 = findMarkovBlanket(g2, graph)
            MB1_MB2 = set(MB1).intersection(set(MB2))

            connected_3pts.append([(g1, g2, x) for x in MB1_MB2])

            for g3 in MB1_MB2:
                MB3 = findMarkovBlanket(g3, graph)
                MB1_MB2_MB3 = set(MB3).intersection(MB1_MB2)

                connected_4pts.append([(g1, g2, g3, x) for x in MB1_MB2_MB3])

                for g4 in MB1_MB2_MB3:
                    MB4 = findMarkovBlanket(g4, graph)
                    MB1_MB2_MB3_MB4 = set(MB4).intersection(MB1_MB2_MB3)

                    connected_5pts.append([(g1, g2, g3, g4, x) for x in MB1_MB2_MB3_MB4])

    print('Generated all connected 3-, 4-, 5-tuples')

    # To aid estimation, order the variables by the size of their Markov blanket so that the estimation uses the smallest one. 
    def onlySmallestMB(ar_):
        # ar_ is of depth 3, so flatten first two levels, and order so that uniques can be kept:
        ar = [tuple(sorted(genes)) for intList in ar_ for genes in intList]
        ar = np.unique(ar, axis=0)
        ar = np.array([sorted(genes, key=lambda x: len(findMarkovBlanket(x, graph))) for genes in ar])
        return ar

    connected_2pts = onlySmallestMB(connected_2pts)
    connected_3pts = onlySmallestMB(connected_3pts)
    connected_4pts = onlySmallestMB(connected_4pts)
    connected_5pts = onlySmallestMB(connected_5pts)

    args_connected_2pts = [(intSet, estimationGraph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for intSet in connected_2pts]
    args_connected_3pts = [(intSet, estimationGraph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for intSet in connected_3pts]
    args_connected_4pts = [(intSet, estimationGraph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for intSet in connected_4pts]
    args_connected_5pts = [(intSet, estimationGraph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for intSet in connected_5pts]

    if PrintBool:
        print(f'Markov-connected pairs: {len(connected_2pts)}, trips: {len(connected_3pts)}, Quads: {len(connected_4pts)}, Pents: {len(connected_5pts)}')
    

    for order, args in [['random_2pts', args_randPairs], ['random_3pts', args_randTrips], ['random_4pts', args_randQuads], ['random_5pts', args_randPents], ['withinMB_2pts', args_connected_2pts], ['withinMB_3pts', args_connected_3pts], ['withinMB_4pts', args_connected_4pts], ['withinMB_5pts', args_connected_5pts]]:

        
        start = time.perf_counter()
        with concurrent.futures.ProcessPoolExecutor(max_workers=maxWorkers) as executor:
            results = executor.map(calcInteraction_withCI_parallel, args)  
        finish = time.perf_counter()
        if PrintBool: print(f'Time elapsed: {round(finish-start, 2)} secs')
        if PrintBool: print('calculation done, storing results...')

        resultArr = np.array(list(results), dtype=object)
        if PrintBool: print('writing files...')
        
        np.save(f'interactions_{order}_{estimationMode}_{ID}', resultArr, allow_pickle=True)

        if PrintBool: print(f'********** DONE with {order} {ID} **********\n')



def main():
    np.random.seed(0)
    notes = ''
    
    estimationMethod = 'expectations'
    estimator = calcInteraction_expectations_numba
    
    
    print('Starting calculation on ' )
    print(f'Calculating {estimationMode}s')
    print('Using estimation method:', estimationMethod)
    if asympBool:
        print(f'Using asymptotic variance estimation')
    else:
        print(f'With {nResamps} bootstrap resamples')
    print(f'Parallelised over {nCores} cores. ')


    calcInteractionsAndWriteNPYs(DSname + notes, graph, trainDat, maxWorkers=nCores, estimator = estimator, nResamps=nResamps, mode=estimationMode)
    
    
    print('***********DONE***********')

if __name__=="__main__":
    main()
    


