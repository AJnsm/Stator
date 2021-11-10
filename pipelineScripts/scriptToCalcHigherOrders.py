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

dataPath = sys.argv[1]
graphPath = sys.argv[2]
coups_3pts_F_Path = sys.argv[3]

nRands = sys.argv[4]

nResamps = int(sys.argv[5])
nCores = int(sys.argv[6])   



trainDat = pd.read_csv(dataPath)
adjMat = pd.read_csv(graphPath, index_col=0)
graph = ig.Graph.Adjacency(adjMat.values.tolist()) 

coups_3pts_F = np.load(coups_3pts_F_Path, allow_pickle=True)


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
        
    elif(order==4):
        E111 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==1).all(axis=1)].iloc[:, 0].mean()
        E000 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==0).all(axis=1)].iloc[:, 0].mean()
        
        E001 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==[0, 0, 1]).all(axis=1)].iloc[:, 0].mean()
        E010 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==[0, 1, 0]).all(axis=1)].iloc[:, 0].mean()
        E100 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==[1, 0, 0]).all(axis=1)].iloc[:, 0].mean()
        
        E011 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==[0, 1, 1]).all(axis=1)].iloc[:, 0].mean()
        E101 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==[1, 0, 1]).all(axis=1)].iloc[:, 0].mean()
        E110 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==[1, 1, 0]).all(axis=1)].iloc[:, 0].mean()
        
       
        if (min(E111, E000, E001, E010, E100, E011, E101, E110)==0) | (max(E111, E000, E001, E010, E100, E011, E101, E110)==1):
            return np.nan
        
        else:
            return np.log(E111*(1-E011)*(1-E101)*E001*E010*(1-E110)*E100*(1-E000)/((1-E111)*E011*E101*(1-E001)*(1-E010)*E110*(1-E100)*E000))
        
        
    elif(order==5):
        E1111 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==1).all(axis=1)].iloc[:, 0].mean()
        E0000 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==0).all(axis=1)].iloc[:, 0].mean()
        
        E0001 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[0, 0, 0, 1]).all(axis=1)].iloc[:, 0].mean()
        E0010 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[0, 0, 1, 0]).all(axis=1)].iloc[:, 0].mean()
        E0100 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[0, 1, 0, 0]).all(axis=1)].iloc[:, 0].mean()
        E1000 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[1, 0, 0, 0]).all(axis=1)].iloc[:, 0].mean()
        
        E0011 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[0, 0, 1, 1]).all(axis=1)].iloc[:, 0].mean()
        E0101 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[0, 1, 0, 1]).all(axis=1)].iloc[:, 0].mean()
        E0110 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[0, 1, 1, 0]).all(axis=1)].iloc[:, 0].mean()
        E1010 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[1, 0, 1, 0]).all(axis=1)].iloc[:, 0].mean()
        E1100 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[1, 1, 0, 0]).all(axis=1)].iloc[:, 0].mean()
        E1001 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[1, 0, 0, 1]).all(axis=1)].iloc[:, 0].mean()
        
        E1110 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[1, 1, 1, 0]).all(axis=1)].iloc[:, 0].mean()
        E1101 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[1, 1, 0, 1]).all(axis=1)].iloc[:, 0].mean()
        E1011 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[1, 0, 1, 1]).all(axis=1)].iloc[:, 0].mean()
        E0111 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3, 4]]==[0, 1, 1, 1]).all(axis=1)].iloc[:, 0].mean()
        
        
        allEs = [E1111, E0000, E0001, E0010, E0100, E1000, E0011, E0101, E0110, E1010, E1100, E1001, E1110, E1101, E1011, E011]
        if (min(allEs)==0) | (max(allEs)==1):
            return np.nan
        else:
            return np.log(E1111*E1100*E1010*E0110*E0101*E0011*E1001*E0000/((1-E1111)*(1-E1100)*(1-E1010)*(1-E0110)*(1-E0101)*(1-E0011)*(1-E1001)*(1-E0000)) \
                          * (1-E0111)*(1-E1011)*(1-E1101)*(1-E1110)*(1-E0001)*(1-E0010)*(1-E0100)*(1-E1000)/(E0111*E1011*E1101*E1110*E0001*E0010*E0100*E1000))

    else:
        print('Order not yet implemented, change estimation method to probabilities.')
        return np.nan


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
                  
        
        
def calcInteractionsAndWriteNPYs(ID, graph, trainDat, maxWorkers, estimator, nResamps=1000):
    
    # if PrintBool: print(f'Starting with {ID}...')
    genes = trainDat.columns
    n = len(genes)


    print('Generating random triplets...')
    randTrips = np.array([np.random.choice(np.arange(500), 3, replace=False) for i in range(nRands)]).astype(int)
    args_randTrips = [(triplet, graph, trainDat, estimator, nResamps) for triplet in randTrips]

    print('Generating random quads...')
    randQuads = np.array([np.random.choice(np.arange(500), 4, replace=False) for i in range(nRands)]).astype(int)
    args_randQuads = [(quad, graph, trainDat, estimator, nResamps) for quad in randQuads]

    print('Generating random quads...')
    randPents = np.array([np.random.choice(np.arange(500), 4, replace=False) for i in range(nRands)]).astype(int)
    args_randPents = [(pent, graph, trainDat, estimator, nResamps) for pent in randPents]


    sigTrips = np.array(np.where(coups_3pts_F<0.01)).T

    quads = []
    pents = []
    for trip in sigTrips:

        nbs = list(set(graph.neighbors(trip[0])) - set(trip))
        newInd = np.random.choice(nbs, 1)
        newInd2 = np.random.choice(nbs, 2, replace=False)
        # Good to keep trip on the left, since that apparently made it estimable
        quads.append(np.hstack([trip, newInd]))
        pents.append(np.hstack([newQuad, newInd2]))


        args_Quads = [(quad, graph, trainDat, estimator, nResamps) for quad in quads]
        args_Pents = [(pent, graph, trainDat, estimator, nResamps) for pent in pents]




    for ID, args in [['random_Trips', args_randTrips], ['random_Quads', args_randQuads], ['random_Pents', args_randPents],
                        ['connected_Quads', args_Quads], ['connected_Pents', args_Pents]]:

        start = time.perf_counter()
        with concurrent.futures.ProcessPoolExecutor(max_workers=maxWorkers) as executor:
            results = executor.map(calcInteraction_withCI_parallel, args)  
        finish = time.perf_counter()
        if PrintBool: print(f'Time elapsed: {round(finish-start, 2)} secs')
        if PrintBool: print('calculation done, storing results...')

        resultArr = np.array(list(results), dtype=object)
        if PrintBool: print('writing files...')
        



        np.save(f'interactions_{n}point_{ID}', resultArr, allow_pickle=True)

        if PrintBool: print(f'********** DONE with {ID} **********\n')




def main():
    np.random.seed(0)

    print('Starting calculation on ' )
    print('Using estimation method:  ', estimationMethod)

    print(f'Calculating interactions at order {intOrder}')
    print(f'With {nResamps} bootstrap resamples')
    print(f'Parallelised over {nCores} cores. ')



    notes = ''
    
    estimationMethod = 'expectations'
    estimator = calcInteraction_expectations
    calcInteractionsAndWriteNPYs(estimationMethod+notes, graph, trainDat, maxWorkers=nCores, estimator = estimator, nResamps=nResamps)
    


        
    
      
    
    print('***********DONE***********')

if __name__=="__main__":
    main()
    

