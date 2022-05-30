PrintBool=False
print('test')
if(__name__=="__main__"):
    PrintBool=True # Only print if in the __main__ call of the script

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
from utilities import *

if PrintBool: print('Modules imported \n')

parser = argparse.ArgumentParser(description='Args for coupling estimation')

parser.add_argument("--dataPath", type=str, nargs=1, help="Path to training data")
parser.add_argument("--graphPath", type=str, nargs=1, help="Path to graph file")
parser.add_argument("--nResamps", type=int, nargs=1, help="Number of BS resamples")
parser.add_argument("--nCores", type=int, nargs=1, help="Number of cores")
parser.add_argument("--nRandoms", type=int, nargs=1, help="Number of random interactions to calculate")
parser.add_argument("--genesToOne", type=str, nargs='?', help="Path to list of genes that should be set to 1")
parser.add_argument("--dataDups", type=int, nargs='?', help="Number of data duplications. 0 is no duplication, and another value is the min binsize allowed (recommended to be 15). ")
parser.add_argument("--boundBool", type=int, nargs='?', help="Boolean that decided whether bounds should also be considered.")

args = parser.parse_args()
    
dataPath = args.dataPath[0]
graphPath = args.graphPath[0]
nResamps = args.nResamps[0]
nCores = args.nCores[0]
nRands = args.nRandoms[0]
genesToOnePath = args.genesToOne
dataDups = args.dataDups
boundBool = args.boundBool

trainDat = pd.read_csv(dataPath)
adjMat = pd.read_csv(graphPath, index_col=0)
graph = ig.Graph.Adjacency(adjMat.values.tolist()) 

try:
    genesToOneIndices = pd.read_csv(genesToOnePath)
except:
    genesToOneIndices = []


if PrintBool: print('data import and graph construction done')

def calcInteractionsAndWriteNPYs(ID, graph, trainDat, maxWorkers, estimator, nResamps=1000):
    
    if PrintBool: print(f'Starting with {ID}...')
    genes = trainDat.columns
    n = len(genes)


    print('Generating random triplets...')
    randTrips = np.array([np.random.choice(np.arange(n), 3, replace=False) for i in range(nRands)]).astype(int)
    args_randTrips = [(triplet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool) for triplet in randTrips]

    print('Generating random quads...')
    randQuads = np.array([np.random.choice(np.arange(n), 4, replace=False) for i in range(nRands)]).astype(int)
    args_randQuads = [(quad, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool) for quad in randQuads]

    print('Generating random pents...')
    randPents = np.array([np.random.choice(np.arange(n), 5, replace=False) for i in range(nRands)]).astype(int)
    args_randPents = [(pent, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool) for pent in randPents]


    ints_2pt = []
    ints_3pt = []
    ints_4pt = []
    ints_5pt = []
    
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
    def onlySmallestMB(ar):
        ar = [tuple(sorted(genes)) for intList in ar for genes in intList]
        ar = np.unique(ar, axis=0)
        ar = np.array([sorted(genes, key=lambda x: len(findMarkovBlanket(x, graph))) for genes in ar])
        return ar

    ints_2pt = onlySmallestMB(ints_2pt)
    ints_3pt = onlySmallestMB(ints_3pt)
    ints_4pt = onlySmallestMB(ints_4pt)
    ints_5pt = onlySmallestMB(ints_5pt)

    ints_2pt = [(intSet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool) for intSet in ints_2pt]
    ints_3pt = [(intSet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool) for intSet in ints_3pt]
    ints_4pt = [(intSet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool) for intSet in ints_4pt]
    ints_5pt = [(intSet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool) for intSet in ints_5pt]

    if PrintBool:
        print(f'Connected trips: {len(ints_3pt)}, Quads: {len(ints_4pt)}, Pents: {len(ints_5pt)}')
        
    for order, args in [['random_3pts', args_randTrips], ['random_4pts', args_randQuads], ['random_5pts', args_randPents],
                     ['MB_3pts', ints_3pt], ['MB_4pts', ints_4pt], ['MB_5pts', ints_5pt]]:
        
                
        
        start = time.perf_counter()
        with concurrent.futures.ProcessPoolExecutor(max_workers=maxWorkers) as executor:
            results = executor.map(calcInteraction_withCI_parallel, args)  
        finish = time.perf_counter()
        if PrintBool: print(f'Time elapsed: {round(finish-start, 2)} secs')
        if PrintBool: print('calculation done, storing results...')

        resultArr = np.array(list(results), dtype=object)
        if PrintBool: print('writing files...')
        
        np.save(f'interactions_withinMB_{order}_{ID}', resultArr, allow_pickle=True)

        if PrintBool: print(f'********** DONE with {ID} **********\n')



def main():
    np.random.seed(0)
    notes = ''
    
    estimationMethod = 'expectations'
    estimator = calcInteraction_expectations_numba
    
    
    print('Starting calculation on ' )
    print('Using estimation method:  ', estimationMethod)

    print(f'With {nResamps} bootstrap resamples')
    print(f'Parallelised over {nCores} cores. ')



   
    calcInteractionsAndWriteNPYs(estimationMethod+notes, graph, trainDat, maxWorkers=nCores, estimator = estimator, nResamps=nResamps)
    


        
    
      
    
    print('***********DONE***********')

if __name__=="__main__":
    main()
    


