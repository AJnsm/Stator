PrintBool=False

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
parser.add_argument("--dataPath", type=str, nargs='?', help="Path to training data")
parser.add_argument("--graphPath", type=str, nargs=1, help="Path to graph file")
parser.add_argument("--intOrder", type=int, nargs='?', help="order of interaction")
parser.add_argument("--nResamps", type=int, nargs=1, help="Number of BS resamples")
parser.add_argument("--nCores", type=int, nargs='?', help="Number of cores")
parser.add_argument("--estimationMethod", type=str, nargs=1, help="Estimation method to use")
parser.add_argument("--edgeListAlpha", type=float, nargs='?', help="Significance threshold for edge list inclusion")
parser.add_argument("--genesToOne", type=str, nargs='?', help="Path to list of genes that should be set to 1")
parser.add_argument("--dataDups", type=int, nargs='?', help="Number of data duplications. 0 is no duplication, and another value is the min binsize allowed (recommended to be 15). ")
parser.add_argument("--boundBool", type=int, nargs='?', help="Boolean that decided whether bounds should also be considered.")

args = parser.parse_args()

dataPath = args.dataPath
graphPath = args.graphPath[0]
intOrder = args.intOrder
nResamps = args.nResamps[0]
nCores = args.nCores
estimationMethod = args.estimationMethod[0]
edgeListAlpha = args.edgeListAlpha
genesToOnePath = args.genesToOne
dataDups = args.dataDups
boundBool = args.boundBool

trainDat = pd.read_csv(dataPath)
DSname = graphPath.split('.')[0]
adjMat = pd.read_csv(graphPath, index_col=0)
graph = ig.Graph.Adjacency(adjMat.values.tolist()) 

try:
    if PrintBool: print('Loading genes to condition on 1')
    genesToOne = pd.read_csv(genesToOnePath).columns.values
    genesToOneIndices = np.where([gene in genesToOne for gene in trainDat.columns.values])[0]
    if PrintBool: print(f'{len(genesToOneIndices)} genes will be conditioned on a 1')
except Exception as e:
    print(e)
    if PrintBool: print('NOTE: all genes conditioned on 0s.')
    genesToOneIndices = []

if PrintBool: 
    if boundBool==1:
        print('including bounds')

    if dataDups>0:
        print(f'Duplicating data up to {dataDups} times')
    else:
        print('no data duplication')


# Creating empty control graph
graph_ctrl = graph.copy()

for e in graph_ctrl.es():
    graph_ctrl.es.delete(e) 

    
    
if PrintBool: print('data import and graph construction done')

        
def calcInteractionsAndWriteNPYs(ID, graph, trainDat, maxWorkers, order, estimator, nResamps=1000):
    
    # if PrintBool: print(f'Starting with {ID}...')
    genes = trainDat.columns
    n = len(genes)


    #We're now doing every interaction multiple times, but that's ok since they come with different markov blankets
    #(As long as mode is not set to 'Min')
    
    if (order==1):
        args = [([x], graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool) for x in range(n)]

    if (order==2):
        args = [([x, y], graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool) for x in range(n) for y in range(n)]
    
    if (order==3):
        trips = []
        print('Generating all connected triplets...')
        for a in range(n):
            for b in range(n):
                if(b!=a):
                    for c in range(n):
                        if((c!=b) & (c!=a)):
                            if (int(a in set(graph.neighbors(b))) + int(b in set(graph.neighbors(c))) + int(c in set(graph.neighbors(a)))>1):
                                trips.append([a, b, c])
        print(f'{len(trips)} triplets generated')
        
#         # Generate random triplets:
#         from itertools import product 
#         from random import sample
        

#         tmp = np.array(sample(list(product(np.arange(100), repeat=3)), k=20000))

#         tmp2 = tmp[(tmp[:, 0]!=tmp[:, 1]) & (tmp[:, 0]!=tmp[:, 2]) & (tmp[:, 1]!= tmp[:, 2])]
#         trips = list(set(tuple([tuple(x) for x in tmp2])))[:10000]
#         trips = [list(trip) for trip in trips]
#         print(f'{len(trips)} triplets generated')
        
        args = [(triplet, graph, trainDat, estimator, nResamps, genesToOneIndices, dataDups, boundBool) for triplet in trips]
    
    start = time.perf_counter()
    with concurrent.futures.ProcessPoolExecutor(max_workers=maxWorkers) as executor:
        results = executor.map(calcInteraction_withCI_parallel, args)  
    finish = time.perf_counter()
    if PrintBool: print(f'Time elapsed: {round(finish-start, 2)} secs')
    if PrintBool: print('calculation done, storing results...')

    resultArr = np.array(list(results), dtype=object)
    if PrintBool: print('writing files...')
    
    if (order==1):
        TLcoups = resultArr[:, 0]
        TLcoups_LB = resultArr[:, 1]
        TLcoups_UB = resultArr[:, 2]
        TLcoups_nonZero = resultArr[:, 3]
        TLcoups_undef = resultArr[:, 4]
        TLcoups_inf = resultArr[:, 5]
        boundArr = resultArr[:, 6]

    if (order==2):
        TLcoups = resultArr[:, 0].reshape([n for i in range(order)])
        TLcoups_LB = resultArr[:, 1].reshape([n for i in range(order)])
        TLcoups_UB = resultArr[:, 2].reshape([n for i in range(order)])
        TLcoups_nonZero = resultArr[:, 3].reshape([n for i in range(order)])
        TLcoups_undef = resultArr[:, 4].reshape([n for i in range(order)])
        TLcoups_inf = resultArr[:, 5].reshape([n for i in range(order)])
        boundArr = resultArr[:, 6].reshape([n for i in range(order)])

    elif (order==3):
        TLcoups, TLcoups_LB, TLcoups_UB, TLcoups_nonZero, TLcoups_undef, TLcoups_inf, boundArr = np.empty((n, n, n)), np.empty((n, n, n)), np.empty((n, n, n)), np.empty((n, n, n)), np.empty((n, n, n)), np.empty((n, n, n)), np.empty((n, n, n))
        TLcoups[:], TLcoups_LB[:], TLcoups_UB[:], TLcoups_nonZero[:], TLcoups_undef[:], TLcoups_inf[:], boundArr[:] = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        for r in resultArr:
            TLcoups[tuple(r[-1])] = r[0]
            TLcoups_LB[tuple(r[-1])] = r[1]
            TLcoups_UB[tuple(r[-1])] = r[2]
            TLcoups_nonZero[tuple(r[-1])] = r[3]
            TLcoups_undef[tuple(r[-1])] = r[4]
            TLcoups_inf[tuple(r[-1])] = r[5]
            boundArr[tuple(r[-1])] = r[6]
            
    np.save(f'interactions_order{order}_{ID}', TLcoups)
    np.save(f'interactions_order{order}_{ID}_CI_LB', TLcoups_LB)
    np.save(f'interactions_order{order}_{ID}_CI_UB', TLcoups_UB)
    np.save(f'interactions_order{order}_{ID}_CI_F', TLcoups_nonZero)
    np.save(f'interactions_order{order}_{ID}_undef', TLcoups_undef)
    np.save(f'interactions_order{order}_{ID}_inf', TLcoups_inf)

    if boundBool:
        np.save(f'interactions_order{order}_{ID}_boundVal', boundArr)



    # ********** writing Cytoscape files ************
    
    def compTups(t1, t2):
        for i in range(len(t1)):
            if t1[i]!=t2[i]:
                return False
        else:
            return True

    def arr2SIF(coups, Fs, alpha = 0.05):    
        nanMask = (~np.isnan(np.array(coups).astype(float)))
        fMask = (Fs<=alpha)
        
        sigCoups = np.array(np.where(nanMask & fMask)).T
        
        
        return pd.DataFrame.from_dict({'genes' : [sigCoup for sigCoup in sigCoups],
                   'coup' : [coups[tuple(sigCoup)] for sigCoup in sigCoups],
                   'F' : [Fs[tuple(sigCoup)] for sigCoup in sigCoups]})

    def onlyUniques_mostSig(sigArr):
        us_mostSig = []
        trips = [tuple(np.sort(gs)) for gs in sigArr['genes'].values]
        us, inds, cs = np.unique(trips, axis=0, return_index=True, return_counts=True)
        for i, u in enumerate(us):
            dups = sigArr[[compTups(x, u) for x in trips]]
            mostSig = np.argmin(dups['F'])
            us_mostSig.append(dups.iloc[mostSig])
        df = pd.DataFrame(data = us_mostSig)
        return df


    if (order==2):
        with open(f"edgeList_interactions_order{order}_{ID}.csv", 'w', encoding = 'utf-8') as f:
            f.write('G1,G2,coup,1-F\n')
            for i, row in onlyUniques_mostSig(arr2SIF(TLcoups, TLcoups_nonZero, alpha = edgeListAlpha)).iterrows():
                s = f"{genes[row['genes'][0]]},{genes[row['genes'][1]]},"
                f.write(s)
                f.write(str(round(row['coup'], 5)) + ',')
                f.write(str(round(1-row['F'], 5)))
                f.write('\n')

    if (order==3):
        with open(f"edgeList_interactions_order{order}_{ID}.csv", 'w', encoding = 'utf-8') as f:
            f.write('G1,G2,coup,1-F\n')
            for i, row in onlyUniques_mostSig(arr2SIF(TLcoups, TLcoups_nonZero, alpha = edgeListAlpha)).iterrows():
                s = f"{genes[row['genes'][0]]},{genes[row['genes'][1]]},"
                f.write(s)
                f.write(str(round(row['coup'], 5)) + ',')
                f.write(str(round(1-row['F'], 5)))
                f.write('\n')
                
                s = f"{genes[row['genes'][1]]},{genes[row['genes'][2]]},"
                f.write(s)
                f.write(str(round(row['coup'], 5)) + ',')
                f.write(str(round(1-row['F'], 5)))
                f.write('\n')
                
                s = f"{genes[row['genes'][0]]},{genes[row['genes'][2]]},"
                f.write(s)
                f.write(str(round(row['coup'], 5)) + ',')
                f.write(str(round(1-row['F'], 5)))
                f.write('\n')


        with open(f"edgeList_interactions_order{order}_collapsed_{ID}.csv", 'w', encoding = 'utf-8') as f:
            f.write('S1,C1,S2,C2\n')
            arr = onlyUniques_mostSig(arr2SIF(TLcoups, TLcoups_nonZero, alpha = edgeListAlpha))
            if len(arr)>0:
                geneSets = [set(x) for x in arr['genes']]
                
                for i in range(len(geneSets)):
                    for j in range(i+1, len(geneSets)):
                        g1 = list(geneSets[i])
                        g2 = list(geneSets[j])
                        for k in range(len(geneSets[i].intersection(geneSets[j]))):
                            f.write(f'{genes[g1[0]]};{genes[g1[1]]};{genes[g1[2]]}')
                            f.write(',')
                            f.write(str(round(arr.iloc[i]['coup'], 5)))
                            f.write(',')
                            
                            f.write(f'{genes[g2[0]]};{genes[g2[1]]};{genes[g2[2]]}')
                            f.write(',')
                            f.write(str(round(arr.iloc[j]['coup'], 5)))
                            f.write('\n')



    if PrintBool: print(f'DONE with {ID}...\n')
    
def main():
    np.random.seed(0)
    # Arguments:
    # 0: used by sys
    # 1: training data
    # 2: graph
    # 3: order of interaction
    # 4: number of bootstrap resamples
    # 5: number of cores
    # 6: pVal table for Hartigan Dip test
    # 7: string to determine estimation method
    
    print('Starting calculation on ' + DSname)
    print('Using estimation method:  ', estimationMethod)

    print(f'Calculating interactions at order {intOrder}')
    print(f'With {nResamps} bootstrap resamples')
    print(f'Parallelised over {nCores} cores. ')



    notes = ''
    

    if estimationMethod == 'both':
        estimator = calcInteraction_binTrick
        calcInteractionsAndWriteNPYs(DSname+'_'+'probabilities'+notes, graph, trainDat, maxWorkers=nCores, order = intOrder, estimator = estimator, nResamps=nResamps)
        estimator = calcInteraction_expectations_numba
        calcInteractionsAndWriteNPYs(DSname+'_'+'expectations'+notes, graph, trainDat, maxWorkers=nCores, order = intOrder, estimator = estimator, nResamps=nResamps)

    elif estimationMethod == 'probabilities':
        estimator = calcInteraction_binTrick
        calcInteractionsAndWriteNPYs(DSname+'_'+estimationMethod+notes, graph, trainDat, maxWorkers=nCores, order = intOrder, estimator = estimator, nResamps=nResamps)

    elif estimationMethod == 'expectations':
        estimator = calcInteraction_expectations_numba
        calcInteractionsAndWriteNPYs(DSname+'_'+estimationMethod+notes, graph, trainDat, maxWorkers=nCores, order = intOrder, estimator = estimator, nResamps=nResamps)
    else:
        print('Invalid estimation method -- terminating...')        
        return 1
        
    
      
    
    print('***********DONE***********')

if __name__=="__main__":
    main()
    
    
    



