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

import scipy
from scipy.stats import kstest
from scipy.interpolate import interp1d

# import hartiganDip

if PrintBool: print('Modules imported \n')

    
if(len(sys.argv)<8):
    print('Not enough arguments -- terminating...')
dataPath = sys.argv[1]
graphPath = sys.argv[2]
intOrder = int(sys.argv[3])
nResamps = int(sys.argv[4])
nCores = int(sys.argv[5])   
pValPath = sys.argv[6]
estimationMethod = sys.argv[7]
edgeListAlpha = float(sys.argv[8])


trainDat = pd.read_csv(dataPath)
pVals = pd.read_csv(pValPath, index_col=0)
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
   
    
def conditionOnMB(genes, graph, dataSet, mode='0'):
    '''
    Calculate the MB for each gene in genes, and set all to zero. 
    mode=='Min' uses the smallest blanket, 
    while an integer specifies a particular gene's MB
    '''
    MBs = [findMarkovBlanket(gene, graph) for gene in genes]
    
    if (mode == 'Min'):
        MB = min(MBs, key=len) 
    elif (mode == 'All'):
        MB = [gene for MB in MBs for gene in MB]      
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
        if (E==1) | (E==0):
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
        print('Order not yet implemented, change estimation method to probabilities.')
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
        print('Order not implemented, using slower order-agnostic version!')

        return calcInteraction_binTrick_allOrders(conditionedGenes)

def calcInteraction_binTrick_allOrders(conditionedGenes):
    
    order = len(conditionedGenes.columns)
    nStates = 2**order
    powers = 2*np.array([np.base_repr(i).count('1')%2==order%2 for i in range(2**order)]).astype(float)-1
    
    f = lambda x: ''.join(map(str, x))
    binCounts = np.bincount(list(map(lambda x: int(x, 2), list(map(f, conditionedGenes.values)))), minlength=nStates)
        
    return np.log(np.prod(np.array([x**p for (x, p) in zip(binCounts, powers)])))
    
    
def calcInteraction_withCI(genes, graph, dataSet, estimator, nResamps=1000):
    '''
    Add 95% confidence interval bounds from bootstrap resamples,
    and the F value: the proportion of resamples with a different sign.
    '''
    
    if estimator is calcInteraction_expectations:
        MBmode = '0' # Use first gene to get MB
    else:
        MBmode = 'All' # Use MB of all genes -- safer, so used as else statement. 

    conditionedGenes = conditionOnMB(genes, graph, dataSet, mode=MBmode)
        
    val0 = estimator(conditionedGenes)
    vals = np.zeros(nResamps)
    if np.isnan(val0):
        return [np.nan, np.nan, np.nan, np.nan, genes]
    
    for i in range(nResamps):
        genes_resampled = conditionedGenes.sample(frac=1, replace=True)
        vals[i] = estimator(genes_resampled)
    
    vals.sort()
    vals_noNan = vals[~np.isnan(vals)]
    # dip = hartiganDip.diptst(vals_noNan)[0]
    # dipPval = fromD2P(dip, len(vals_noNan))

    ksStat = kstest(vals_noNan, lambda x: scipy.stats.norm.cdf(x, loc=vals_noNan.mean(), scale=vals_noNan.std()))[1]

    CI = (vals_noNan[int(np.around(len(vals_noNan)/40))], vals_noNan[int(np.floor(len(vals_noNan)*39/40))])

    propDifSign = sum(np.sign(vals)==-np.sign(val0))/nResamps
    
    # # If it's *really* close to a unimodal distribution according to Dip or KS test, or doesn't have undef. resamples:
    # if((len(vals_noNan) == nResamps) | (dipPval>=0.99) | (ksStat>0.01)) : 

    if((len(vals_noNan) == nResamps) | (ksStat>0.01)) : 

        return [val0, CI[0], CI[1], propDifSign, genes]
    else:
        return [np.nan, np.nan, np.nan, np.nan, genes]      
    
def calcInteraction_withCI_parallel(args, nResamps=1000):
    '''
    wrapper to unpack function arguments so that it can be mapped over process pool with one arg.
    (I actually think there is something like executor.starmap that could do this for us)
    '''
    genes, graph, dataSet, estimator, nResamps = args
    
    return calcInteraction_withCI(genes, graph, dataSet, estimator, nResamps=nResamps)       
                  
        
        
def calcInteractionsAndWriteNPYs(ID, graph, trainDat, maxWorkers, order, estimator, nResamps=1000):
    
    # if PrintBool: print(f'Starting with {ID}...')
    genes = trainDat.columns
    n = len(genes)


    #We're now doing every interaction multiple times, but that's ok since they come with different markov blankets
    #(As long as mode is not set to 'Min')
    
    if (order==1):
        args = [([x], graph, trainDat, estimator, nResamps) for x in range(n)]

    if (order==2):
        args = [([x, y], graph, trainDat, estimator, nResamps) for x in range(n) for y in range(n)]
    
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
        
        args = [(triplet, graph, trainDat, estimator, nResamps) for triplet in trips]
    
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
        df.reset_index(drop=True)
        return df


    if (order==2):
        with open("edgeList_interactions_order{order}_{ID}.csv", 'w', encoding = 'utf-8') as f:
            f.write('G1,G2,coup,1-F\n')
            for i, row in onlyUniques_mostSig(arr2SIF(TLcoups, TLcoups_nonZero, alpha = edgeListAlpha)).iterrows():
                s = f"{genes[row['genes'][0]]},{genes[row['genes'][1]]},"
                f.write(s)
                f.write(str(round(row['coup'], 5)) + ',')
                f.write(str(round(1-row['F'], 5)))
                f.write('\n')

    if (order==3):
        with open("edgeList_interactions_order{order}_{ID}.csv", 'w', encoding = 'utf-8') as f:
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


        with open("edgeList_interactions_order{order}_collapsed_{ID}.csv", 'w', encoding = 'utf-8') as f:
            f.write('S1,C1,S2,C2\n')
            arr = onlyUniques_mostSig(arr2SIF(TLcoups, TLcoups_nonZero, alpha = edgeListAlpha))
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
                        f.write(str(round(arr.iloc[i]['coup'], 5)))
                        f.write('\n')



    if PrintBool: print(f'DONE with {ID}...\n')

def fromD2P(D, n):
    if np.isnan(D):
        return np.nan
    if n <=3:
        p=1
    else:
        Ps = pVals.columns.values.astype('float')
        nn = pVals.index.values
        max_n = max(nn)
        
        if n>max_n:
            n1, n0 = max_n, max_n
            i_2, i_n = len(nn), len(nn)
            f_n = 0
            
        else:
            i_n = np.argmin(list(map(lambda x: n-x if (n-x)>0 else np.nan, nn))) - 1
            i_2 = i_n + 1 
            n_0 = nn[i_n]
            n_1 = nn[i_2]
            f_n = (n-n_0)/(n_1-n_0)

        y_0 = np.sqrt(n_0) * pVals.iloc[i_n, :]
        y_1 = np.sqrt(n_1) * pVals.iloc[i_2, :]
        sD = np.sqrt(n) * D
        f = interp1d((y_0 + f_n * (y_1 - y_0)).values, Ps, kind='linear', fill_value=(0, 1), bounds_error=False) # fill_value='extrapolate')

        p = 1-f(sD)

    return p
    
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
        estimator = calcInteraction_expectations
        calcInteractionsAndWriteNPYs(DSname+'_'+'expectations'+notes, graph, trainDat, maxWorkers=nCores, order = intOrder, estimator = estimator, nResamps=nResamps)

    elif estimationMethod == 'probabilities':
        estimator = calcInteraction_binTrick
        calcInteractionsAndWriteNPYs(DSname+'_'+estimationMethod+notes, graph, trainDat, maxWorkers=nCores, order = intOrder, estimator = estimator, nResamps=nResamps)

    elif estimationMethod == 'expectations':
        estimator = calcInteraction_expectations
        calcInteractionsAndWriteNPYs(DSname+'_'+estimationMethod+notes, graph, trainDat, maxWorkers=nCores, order = intOrder, estimator = estimator, nResamps=nResamps)
    else:
        print('Invalid estimation method -- terminating...')        
        return 1
        
    
      
    
    print('***********DONE***********')

if __name__=="__main__":
    main()
    
    
    



