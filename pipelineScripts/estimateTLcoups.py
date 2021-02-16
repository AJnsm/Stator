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

if PrintBool: print('Modules imported \n')
    
dataPath = sys.argv[1]
trainDat = pd.read_csv(dataPath)


graphPath = sys.argv[2]
DSname = graphPath.split('.')[0]
adjMat = pd.read_csv(graphPath, index_col=0)
  
    
graph = ig.Graph.Adjacency(adjMat.values.tolist()) 


# Creating empty control graph
graph_ctrl = graph.copy()

for e in graph_ctrl.es():
    graph_ctrl.es.delete(e) 

    
    
if PrintBool: print('data import and graph construction done')

def findMarkovBlanket(v, g):
    '''
    Markov blanket is all parents, children, and spouses.
    i.e. the parents of your chidren, except yourself
    '''
    parents = g.neighbors(v, mode='IN')
    children = g.neighbors(v, mode='OUT')
    spouses = [spouse for s in [g.neighbors(child, mode='IN') for child in children] 
               for spouse in s if (spouse != v)] #Nested for loops to flatten list
    
    return list(set(parents + children + spouses)) #Sets to keep uniques
   
    
def conditionOnMB(genes, PCgraph, dataSet, mode='0'):
    '''
    Calculate the MB for each gene in genes, and set all to zero. 
    mode=='Min' uses the smallest blanket, 
    while an integer specifies a particular gene's MB
    '''
    MBs = [findMarkovBlanket(gene, PCgraph) for gene in genes]
    
    if (mode=='Min'):
        MB=min(MBs, key=len)        
    else:
        try:
            MB = MBs[int(mode)]
        except:
            print('Invalid mode')
    
    MB = list(set(MB) - set(genes)) #Remove the interacting genes from the markov Blanket.     
    data_conditioned = dataSet[(dataSet.iloc[:, MB]==0).all(axis=1)] #Set whole MB to zero.     
    return data_conditioned.iloc[:, genes]

def calcInteraction_expectations(conditionedGenes):
    '''
    Calc the interactions from the conditioned data using expectation values. 
    Could see if using probabilities is faster. 
    Currently, orders are implemented separately.

    NOTE: previously, each interactions was corrected by a factor to match the {-1, 1} basis of the Ising model.
    I now just take logs and leave them in the {0, 1} basis, which makes more sense for gene expression.
    '''
    
    order = len(conditionedGenes.columns)
    
    if(order==1):
        
        E = conditionedGenes.iloc[:].mean()[0]
        if E==1:
            return np.nan
        return np.log(E/(1-E))


    if(order==2):
        E1 = conditionedGenes[conditionedGenes.iloc[:, 1]==1].iloc[:, 0].mean()
        E0 = conditionedGenes[conditionedGenes.iloc[:, 1]==0].iloc[:, 0].mean()
        if (min(E1, E0)==0) | (max(E1, E0)==1):
            return np.nan
        return np.log(E1*(1-E0)/(E0*(1-E1)))
        
    elif(order==3):
        E11 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2]]==1).all(axis=1)].iloc[:, 0].mean()
        E00 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2]]==0).all(axis=1)].iloc[:, 0].mean()
        
        E10 = conditionedGenes[(conditionedGenes.iloc[:, 1]==1) & (conditionedGenes.iloc[:, 2]==0)].iloc[:, 0].mean()
        E01 = conditionedGenes[(conditionedGenes.iloc[:, 1]==0) & (conditionedGenes.iloc[:, 2]==1)].iloc[:, 0].mean()
    
        if (min(E11, E00, E10, E01)==0) | (max(E11, E00, E10, E01)==1):
            return np.nan
        
        else:
            return np.log(E11*(1-E01)*E00*(1-E10)/(E01*(1-E11)*E10*(1-E00)))
    else:
        print('Order not yet implemented')
        return np.nan
    
def calcInteraction_binTrick(conditionedGenes):
    order = len(conditionedGenes.columns)
    nStates = 2**order
    if order==1:
        binCs = np.bincount(conditionedGenes.values.flatten(), minlength=nStates)
        if sum(binCs==0)>0:
            return np.nan
        else:
            return np.log(binCs[1]/binCs[0])
    
    elif order==2:
        binCs = np.bincount(2 * conditionedGenes.iloc[:, 0] +  conditionedGenes.iloc[:, 1], minlength=nStates)
        if sum(binCs==0)>0:
            return np.nan
        else:
            return np.log(binCs[0]*binCs[3]/(binCs[1]*binCs[2]))
    
    elif order==3:
        binCs = np.bincount(4 * conditionedGenes.iloc[:, 0] + 2 * conditionedGenes.iloc[:, 1] +  conditionedGenes.iloc[:, 2], minlength=nStates)
        if sum(binCs==0)>0:
            return np.nan
        else:
            return np.log(np.prod(binCs[[1, 2, 4, 7]])/np.prod(binCs[[0, 3, 5, 6]]))

    else:
        print('Order not implemented, use order-agnostic version!')
        return np.nan

def calcInteraction_binTrick_allOrders(conditionedGenes):
    
    order = len(conditionedGenes.columns)
    nStates = 2**order
    powers = 2*np.array([np.base_repr(i).count('1')%2==order%2 for i in range(2**order)]).astype(float)-1
    
    f = lambda x: ''.join(map(str, x))
    binCounts = np.bincount(list(map(lambda x: int(x, 2), list(map(f, conditionedGenes.values)))), minlength=nStates)
        
    return np.log(np.prod(np.array([x**p for (x, p) in zip(binCounts, powers)])))
    
    
def calcInteraction_binTrick_withCI(genes, PCgraph, dataSet, nResamps=1000):
    '''
    Add 95% confidence interval bounds from bootstrap resamples,
    and the F value: the proportion of resamples with a different sign.
    '''
    
    
    conditionedGenes = conditionOnMB(genes, PCgraph, dataSet, mode='0')
        
    val0 = calcInteraction_binTrick(conditionedGenes)
    vals = np.zeros(nResamps)
    if np.isnan(val0):
        return [np.nan, np.nan, np.nan, np.nan, genes]
    
    for i in range(nResamps):
        genes_resampled = conditionedGenes.sample(frac=1, replace=True)
        vals[i] = calcInteraction_binTrick(genes_resampled)
    
    vals.sort()
    vals_noNan = vals[~np.isnan(vals)]
    CI = (vals_noNan[int(np.around(len(vals_noNan)/40))], vals_noNan[int(np.floor(len(vals_noNan)*39/40))])

    propDifSign = sum(np.sign(vals)==-np.sign(val0))/nResamps
    
    if(len(vals_noNan) >= 0.999*nResamps): # Threshold to set allowed nans in BS distribution. 
        return [val0, CI[0], CI[1], propDifSign, genes]
    else:
        return [np.nan, np.nan, np.nan, np.nan, genes]      
    
def calcInteraction_binTrick_withCI_parallel(args, nResamps=1000):
    '''
    wrapper to unpack function arguments so that it can be mapped over process pool with one arg.
    (I actually think there is something like executor.starmap that could do this for us)
    '''
    genes, PCgraph, dataSet, nResamps = args
    
    return calcInteraction_binTrick_withCI(genes, PCgraph, dataSet, nResamps=nResamps)       
                  
        
        
def calcInteractionsAndWriteNPYs(ID, graph, trainDat, maxWorkers, order=2, nResamps=1000):
    
    # if PrintBool: print(f'Starting with {ID}...')
    n = len(trainDat.columns)


    #We're now doing every interaction twice, but that's ok since they come with different markov blankets
    #(As long as mode is not set to 'Min')
    
    if (order==1):
        args = [([x], graph, trainDat, nResamps) for x in range(n)]

    if (order==2):
        args = [([x, y], graph, trainDat, nResamps) for x in range(n) for y in range(n)]
    
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
        
        args = [(triplet, graph, trainDat, nResamps) for triplet in trips]
    
    start = time.perf_counter()
    with concurrent.futures.ProcessPoolExecutor(max_workers=maxWorkers) as executor:
        results = executor.map(calcInteraction_binTrick_withCI_parallel, args)  
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

    if (order==2):
        TLcoups = resultArr[:, 0].reshape([n for i in range(order)])
        TLcoups_LB = resultArr[:, 1].reshape([n for i in range(order)])
        TLcoups_UB = resultArr[:, 2].reshape([n for i in range(order)])
        TLcoups_nonZero = resultArr[:, 3].reshape([n for i in range(order)])

    elif (order==3):
        TLcoups, TLcoups_LB, TLcoups_UB, TLcoups_nonZero = np.empty((n, n, n)), np.empty((n, n, n)), np.empty((n, n, n)), np.empty((n, n, n))
        TLcoups[:], TLcoups_LB[:], TLcoups_UB[:], TLcoups_nonZero[:] = np.nan, np.nan, np.nan, np.nan
        for r in resultArr:
            TLcoups[tuple(r[-1])] = r[0]
            TLcoups_LB[tuple(r[-1])] = r[1]
            TLcoups_UB[tuple(r[-1])] = r[2]
            TLcoups_nonZero[tuple(r[-1])] = r[3]
            
    np.save(f'interactions_order{order}_{ID}', TLcoups)
    np.save(f'interactions_order{order}_{ID}_CI_LB', TLcoups_LB)
    np.save(f'interactions_order{order}_{ID}_CI_UB', TLcoups_UB)
    np.save(f'interactions_order{order}_{ID}_CI_F', TLcoups_nonZero)

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
    
    
    notes = ''
    
    
    if(len(sys.argv)<5):
        print('Not enough arguments, terminating...')
    else:        
        intOrder = int(sys.argv[3])
        nResamps = int(sys.argv[4])
        nCores = int(sys.argv[5])
        print('Starting calculation on ' + DSname)
#         print('NOTE: low resample rate (100), should perhaps be larger.\n')                                         
        
        print(f'Calculating interactions at order {intOrder}')
        print(f'With {nResamps} bootstrap resamples')
        print(f'Parallelised over {nCores} cores. ')
        calcInteractionsAndWriteNPYs(DSname+notes, graph, trainDat, maxWorkers=nCores, order = intOrder, nResamps=nResamps)
          
        
        print('***********DONE***********')

if __name__=="__main__":
    main()
    
    
    



