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
import scipy
import argparse

# The utilities module holds the estimation methods etc. 
from utilities import *

import numba
from numba import njit

if PrintBool: print('Modules imported \n')

# Parse all command line arguments
parser = argparse.ArgumentParser(description='Args for coupling estimation')
parser.add_argument("--dataPath", type=str, help="Path to training data")
parser.add_argument("--graphPath", type=str, help="Path to graph file")
parser.add_argument("--nCores", type=int, help="Number of cores")
parser.add_argument("--nResamps", type=int, help="Number of BS resamples")
parser.add_argument("--nRandoms", type=int, help="Number of random interactions to calculate")
parser.add_argument("--genesToOne", type=str, help="Path to list of genes that should be set to 1")
parser.add_argument("--dataDups", type=int, help="Number of data duplications. 0 is no duplication, and another value is the min binsize allowed (recommended to be 15). ")
parser.add_argument("--pathTo5pts", type=str, help="Path to calculated 5-point interactions")
parser.add_argument("--alpha5pts", type=float, help="Significance threshold on 5-pts to search for 6- and 7-pts.")
parser.add_argument("--boundBool", type=int, help="Boolean that decided whether bounds should also be considered.")
parser.add_argument("--asympBool", type=int, help="Boolean to decide whether to use Bootstrap resampling (0) or asymptotic uncertainty estimation (1).")


args = parser.parse_args()
    
dataPath = args.dataPath
graphPath = args.graphPath
nResamps = args.nResamps
nCores = args.nCores
nRands = args.nRandoms
pathTo5pts = args.pathTo5pts
genesToOnePath = args.genesToOne
boundBool = args.boundBool
dataDups = args.dataDups
asympBool = args.asympBool

trainDat = pd.read_csv(dataPath)

# DSname copies the naming scheme from the graphs.
DSname = graphPath.split('.')[0]
adjMat = pd.read_csv(graphPath, index_col=0)
graph = ig.Graph.Adjacency(adjMat.values.tolist()) 
graph.vs['label'] = adjMat.columns.values

try:
    genesToOneIndices = pd.read_csv(genesToOnePath)
except:
    if PrintBool: print('NOTE: all genes conditioned on 0s.')
    genesToOneIndices = []


ints = np.load(pathTo5pts, allow_pickle=True)

# Significance threshold that determines when a 5-point is interesting enough to base a 6- or 7-point estimation on. 
alpha=0.05

# Keep only significant 5-points that are not bounds, and have no undefined or divergent resamples. 
perfectSigEsts = list(map(lambda x: (((x[[4, 5, 6]]==0).all()) & (x[3]<=alpha)), ints))

# Sometimes there are no such 5-points, or not enought variables: terminate. 
try:
    HHOIs = ints[perfectSigEsts][:, [0, -1]]
except:
    print('Could not find 5-points to base estimation on -- terminating...')
    HHOIs = []
    sys.exit()

if len(trainDat.columns)<7:
    print('Too few genes to calculate 6- & 7-points -- terminating...')
    sys.exit()
        
        
def calcInteractionsAndWriteNPYs(ID, maxWorkers, estimator, nResamps=1000):
    
    if PrintBool: print(f'Starting with {ID}...')
        
    genes = trainDat.columns

    # Generate random 6- and 7-tuples. 
    print('Generating random sextuples...')
    randSexs = np.array([np.random.choice(np.arange(len(genes)), 6, replace=False) for i in range(nRands)]).astype(int)
    args_randSexs = [(Sextet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for Sextet in randSexs]

    print('Generating random Septuples...')
    randSepts = np.array([np.random.choice(np.arange(len(genes)), 7, replace=False) for i in range(nRands)]).astype(int)
    args_randSepts = [(Septet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for Septet in randSepts]

    ints_6pt = []
    ints_7pt = []


    # Generate 6- and 7-tuples from the intersection of the Markov blanket of the 5-point interactors.
    print('Generating connected tuples...')
    for sigInt in HHOIs[:, 1]:

        MBs = [set(findMarkovBlanket(g, graph)) for g in sigInt]
        MB_intersection = set.intersection(*MBs)

        interactors = [x for x in MB_intersection if not x in sigInt]

        ints_6pt.append([np.hstack([sigInt, x]) for x in interactors])
        ints_7pt.append([np.hstack([sigInt, x, y]) for x in interactors for y in interactors if x!=y])
    print('Generated the 6- and 7-tuples')
    
    # To aid estimation, order the variables by the size of their Markov blanket so that the estiamtion uses the smallest one. 
    def onlySmallestMB(ar):
        ar = [tuple(sorted(genes)) for intList in ar for genes in intList]
        ar = np.unique(ar, axis=0)
        ar = np.array([sorted(genes, key=lambda x: len(findMarkovBlanket(x, graph))) for genes in ar])
        return ar

    ints_6pt = onlySmallestMB(ints_6pt)
    ints_7pt = onlySmallestMB(ints_7pt)

    ints_6pt = [(intSet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for intSet in ints_6pt]
    ints_7pt = [(intSet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool, asympBool) for intSet in ints_7pt]

    if PrintBool:
        print(f'Connected sextuplets: {len(ints_6pt)}, Septuplets: {len(ints_7pt)}')
        
    for order, args in [['random_6pts', args_randSexs], ['random_7pts', args_randSepts],
                         ['withinMBAndInteracting_6pts', ints_6pt], ['withinMBAndInteracting_7pts', ints_7pt]]:
        
                
        
        start = time.perf_counter()
        with concurrent.futures.ProcessPoolExecutor(max_workers=maxWorkers) as executor:
            results = executor.map(calcInteraction_withCI_parallel, args)  
        finish = time.perf_counter()
        if PrintBool: print(f'Time elapsed: {round(finish-start, 2)} secs')
        if PrintBool: print('calculation done, storing results...')

        resultArr = np.array(list(results), dtype=object)
        if PrintBool: print('writing files...')
        
        np.save(f'interactions_{order}_{ID}', resultArr, allow_pickle=True)
        
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
    print(f'Asymptotic variance estimation: {bool(asympBool)}')


    calcInteractionsAndWriteNPYs(ID = DSname + '_' + notes, maxWorkers=nCores, estimator = estimator, nResamps=nResamps)

    
    print('***********DONE***********')

if __name__=="__main__":
    main()
    


