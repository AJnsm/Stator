
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
import scipy
import argparse
from utilities import *

import numba
from numba import njit

if PrintBool: print('Modules imported \n')

parser = argparse.ArgumentParser(description='Args for coupling estimation')

parser.add_argument("--dataPath", type=str, nargs=1, help="Path to training data")
parser.add_argument("--graphPath", type=str, nargs=1, help="Path to graph file")
parser.add_argument("--nCores", type=int, nargs=1, help="Number of cores")
parser.add_argument("--nResamps", type=int, nargs=1, help="Number of BS resamples")
parser.add_argument("--nRandoms", type=int, nargs=1, help="Number of random interactions to calculate")
parser.add_argument("--genesToOne", type=str, nargs='?', help="Path to list of genes that should be set to 1")
parser.add_argument("--dataDups", type=int, nargs='?', help="Number of data duplications. 0 is no duplication, and another value is the min binsize allowed (recommended to be 15). ")
parser.add_argument("--pathTo5pts", type=str, nargs='?', help="Path to calculated 5-point interactions")
parser.add_argument("--boundBool", type=int, nargs='?', help="Boolean that decided whether bounds should also be considered.")

args = parser.parse_args()
    
dataPath = args.dataPath[0]
graphPath = args.graphPath[0]
nResamps = args.nResamps[0]
nCores = args.nCores[0]
nRands = args.nRandoms[0]
pathTo5pts = args.pathTo5pts
genesToOnePath = args.genesToOne
boundBool = args.boundBool

trainDat = pd.read_csv(dataPath)
DSname = graphPath.split('.')[0]
adjMat = pd.read_csv(graphPath, index_col=0)
graph = ig.Graph.Adjacency(adjMat.values.tolist()) 

try:
    genesToOneIndices = pd.read_csv(genesToOnePath)
except:
    genesToOneIndices = []


HHOIs = {}
alpha=0.05

ints = np.load(pathTo5pts, allow_pickle=True)
if PrintBool: print(ints)
perfectSigEsts = list(map(lambda x: (((x[[4, 5, 6]]==0).all()) & (x[3]<=alpha)), ints))
if PrintBool: print(perfectSigEsts)

try:
    HHOIs = ints[perfectSigEsts][:, [0, -1]]
except:
    print('Could not find 5-points to base estimation on, terminating')
    sys.exit()

        
        
def calcInteractionsAndWriteNPYs(ID, maxWorkers, estimator, nResamps=1000):
    
    if PrintBool: print(f'Starting with {ID}...')
        
    genes = trainDat.columns

    print('Generating random sextuples...')
    randSexs = np.array([np.random.choice(np.arange(len(genes)), 6, replace=False) for i in range(nRands)]).astype(int)
    args_randSexs = [(Sextet, graph, trainDat, estimator, nResamps) for Sextet in randSexs]

    print('Generating random Septuples...')
    randSepts = np.array([np.random.choice(np.arange(len(genes)), 7, replace=False) for i in range(nRands)]).astype(int)
    args_randSepts = [(Septet, graph, trainDat, estimator, nResamps) for Septet in randSepts]

    ints_6pt = []
    ints_7pt = []

    print('Generating connected tuples...')

    for sigInt in HHOIs[:, 1]:

        MBs = [set(findMarkovBlanket(g, graph)) for g in sigInt]
        MB_intersection = set.intersection(*MBs)

        interactors = [x for x in MB_intersection if not x in sigInt]

        ints_6pt.append([np.hstack([sigInt, x]) for x in interactors])
        ints_7pt.append([np.hstack([sigInt, x, y]) for x in interactors for y in interactors if x!=y])
    print('Generated the 6- and 7-tuples')
    
    
    def onlySmallestMB(ar):
        ar = [tuple(sorted(genes)) for intList in ar for genes in intList]
        ar = np.unique(ar, axis=0)
        ar = np.array([sorted(genes, key=lambda x: len(findMarkovBlanket(x, graph))) for genes in ar])
        return ar

    ints_6pt = onlySmallestMB(ints_6pt)
    ints_7pt = onlySmallestMB(ints_7pt)

    ints_6pt = [(intSet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool) for intSet in ints_6pt]
    ints_7pt = [(intSet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool) for intSet in ints_7pt]

    if PrintBool:
        print(f'Connected sextuplets: {len(ints_6pt)}, Septuplets: {len(ints_7pt)}')
        
    for order, args in [['random_6pts', args_randSexs], ['random_7pts', args_randSepts],
                         ['MB_6pts', ints_6pt], ['MB_7pts', ints_7pt]]:
        
                
        
        start = time.perf_counter()
        with concurrent.futures.ProcessPoolExecutor(max_workers=maxWorkers) as executor:
            results = executor.map(calcInteraction_withCI_parallel, args)  
        finish = time.perf_counter()
        if PrintBool: print(f'Time elapsed: {round(finish-start, 2)} secs')
        if PrintBool: print('calculation done, storing results...')

        resultArr = np.array(list(results), dtype=object)
        if PrintBool: print('writing files...')
        
        np.save(f'interactions_withinMB_{order}_ID', resultArr, allow_pickle=True)
        
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


    calcInteractionsAndWriteNPYs(ID = DSname + '_' + estimationMethod+notes, maxWorkers=nCores, estimator = estimator, nResamps=nResamps)

    
    print('***********DONE***********')

if __name__=="__main__":
    main()
    


