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
# from scipy.stats import kstest
# from scipy.interpolate import interp1d

# import hartiganDip

if PrintBool: print('Modules imported \n')

parser = argparse.ArgumentParser(description='Args for coupling estimation')
parser.add_argument("--dataPath", type=str, nargs='?', help="Path to training data")
parser.add_argument("--graphPath", type=str, nargs=1, help="Path to graph file")
parser.add_argument("--intOrder", type=int, nargs='?', help="order of interaction")
parser.add_argument("--nResamps", type=int, nargs=1, help="Number of BS resamples")
parser.add_argument("--nCores", type=int, nargs='?', help="Number of cores")
parser.add_argument("--estimationMethod", type=str, nargs=1, help="Estimation method to use")
parser.add_argument("--edgeListAlpha", type=str, nargs='?', help="Significance threshold for edge list inclusion")
parser.add_argument("--genesToOne", type=str, nargs='?', help="Path to list of genes that should be set to 1")

args = parser.parse_args()

dataPath = args.dataPath
graphPath = args.graphPath[0]
intOrder = args.intOrder
nResamps = args.nResamps[0]
nCores = args.nCores
estimationMethod = args.estimationMethod[0]
edgeListAlpha = args.edgeListAlpha
genesToOnePath = args.genesToOne



trainDat = pd.read_csv(dataPath)
DSname = graphPath.split('.')[0]
adjMat = pd.read_csv(graphPath, index_col=0)
graph = ig.Graph.Adjacency(adjMat.values.tolist()) 

try:
    if PrintBool: print('Loading genes to condition on 1')
    genesToOne = pd.read_csv(args.genesToOnePath).columns.values
    genesToOneIndices = np.where([gene in genesToOne for gene in geneNames])[0]

except:
    if PrintBool: print('NOTE: all genes conditioned on 0s.')
    genesToOneIndices = []




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
   
    
def conditionOnMB_new(genes, graph, dataSet, mode='0', genesToOne=genesToOneIndices):
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
    
    condState = [1 if gene in genesToOne else 0 for gene in MB]
    data_conditioned = dataSet[(dataSet.iloc[:, MB]==condState).all(axis=1)] #Set whole MB to conditioned state.     
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
        num = E
        denom = 1-E
        
    elif(order==2):
        E1 = conditionedGenes[conditionedGenes.iloc[:, 1]==1].iloc[:, 0].mean()
        E0 = conditionedGenes[conditionedGenes.iloc[:, 1]==0].iloc[:, 0].mean()

        num = E1*(1-E0)
        denom = E0*(1-E1)
        
    elif(order==3):
        E11 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2]]==1).all(axis=1)].iloc[:, 0].mean()
        E00 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2]]==0).all(axis=1)].iloc[:, 0].mean()
        
        E10 = conditionedGenes[(conditionedGenes.iloc[:, 1]==1) & (conditionedGenes.iloc[:, 2]==0)].iloc[:, 0].mean()
        E01 = conditionedGenes[(conditionedGenes.iloc[:, 1]==0) & (conditionedGenes.iloc[:, 2]==1)].iloc[:, 0].mean()
        
        num = E11*(1-E01)*E00*(1-E10)
        denom = E01*(1-E11)*E10*(1-E00)
        
    elif(order==4):
        E111 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==1).all(axis=1)].iloc[:, 0].mean()
        E000 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==0).all(axis=1)].iloc[:, 0].mean()
        
        E001 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==[0, 0, 1]).all(axis=1)].iloc[:, 0].mean()
        E010 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==[0, 1, 0]).all(axis=1)].iloc[:, 0].mean()
        E100 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==[1, 0, 0]).all(axis=1)].iloc[:, 0].mean()
        
        E011 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==[0, 1, 1]).all(axis=1)].iloc[:, 0].mean()
        E101 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==[1, 0, 1]).all(axis=1)].iloc[:, 0].mean()
        E110 = conditionedGenes[(conditionedGenes.iloc[:, [1, 2, 3]]==[1, 1, 0]).all(axis=1)].iloc[:, 0].mean()
        
       
        num = E111*(1-E011)*(1-E101)*E001*E010*(1-E110)*E100*(1-E000)
        denom = (1-E111)*E011*E101*(1-E001)*(1-E010)*E110*(1-E100)*E000
        
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
        
        num = E1111*E1100*E1010*E0110*E0101*E0011*E1001*E0000 * (1-E0111)*(1-E1011)*(1-E1101)*(1-E1110)*(1-E0001)*(1-E0010)*(1-E0100)*(1-E1000)
        denom = (1-E1111)*(1-E1100)*(1-E1010)*(1-E0110)*(1-E0101)*(1-E0011)*(1-E1001)*(1-E0000)*(E0111*E1011*E1101*E1110*E0001*E0010*E0100*E1000)

    else:
        print('Order not yet implemented, change estimation method to probabilities.')
        return np.nan

    if ((num==0) & (denom==0)):
            return np.nan 
    elif num==0:
        return -np.inf
    elif denom==0:
        return np.inf
    else:
        return np.log(num/denom)
    
def calcInteraction_expectations_np(conditionedGenes):
    '''
    Same as calcInteraction_expectations, but might be 10x faster?!
    '''
    
    order = len(conditionedGenes.columns)
    
    conditionedGenes_np = conditionedGenes.values
    
    
    if(order==1):
        E = conditionedGenes_np[:].mean()
        num = E
        denom = 1-E
        
    elif(order==2):
        E1 = conditionedGenes_np[conditionedGenes_np[:, 1]==1][:, 0].mean()
        E0 = conditionedGenes_np[conditionedGenes_np[:, 1]==0][:, 0].mean()

        num = E1*(1-E0)
        denom = E0*(1-E1)
        
    elif(order==3):
        E11 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2]]==1).all(axis=1)][:, 0].mean()
        E00 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2]]==0).all(axis=1)][:, 0].mean()
        
        E10 = conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==0)][:, 0].mean()
        E01 = conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==1)][:, 0].mean()
        
        num = E11*(1-E01)*E00*(1-E10)
        denom = E01*(1-E11)*E10*(1-E00)
        
    elif(order==4):
        E111 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3]]==1).all(axis=1)][:, 0].mean()
        E000 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3]]==0).all(axis=1)][:, 0].mean()
        
        E001 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3]]==[0, 0, 1]).all(axis=1)][:, 0].mean()
        E010 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3]]==[0, 1, 0]).all(axis=1)][:, 0].mean()
        E100 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3]]==[1, 0, 0]).all(axis=1)][:, 0].mean()
        
        E011 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3]]==[0, 1, 1]).all(axis=1)][:, 0].mean()
        E101 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3]]==[1, 0, 1]).all(axis=1)][:, 0].mean()
        E110 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3]]==[1, 1, 0]).all(axis=1)][:, 0].mean()
        
       
        num = E111*(1-E011)*(1-E101)*E001*E010*(1-E110)*E100*(1-E000)
        denom = (1-E111)*E011*E101*(1-E001)*(1-E010)*E110*(1-E100)*E000
        
    elif(order==5):
        E1111 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==1).all(axis=1)][:, 0].mean()
        E0000 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==0).all(axis=1)][:, 0].mean()
        
        E0001 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[0, 0, 0, 1]).all(axis=1)][:, 0].mean()
        E0010 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[0, 0, 1, 0]).all(axis=1)][:, 0].mean()
        E0100 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[0, 1, 0, 0]).all(axis=1)][:, 0].mean()
        E1000 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[1, 0, 0, 0]).all(axis=1)][:, 0].mean()
        
        E0011 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[0, 0, 1, 1]).all(axis=1)][:, 0].mean()
        E0101 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[0, 1, 0, 1]).all(axis=1)][:, 0].mean()
        E0110 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[0, 1, 1, 0]).all(axis=1)][:, 0].mean()
        E1010 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[1, 0, 1, 0]).all(axis=1)][:, 0].mean()
        E1100 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[1, 1, 0, 0]).all(axis=1)][:, 0].mean()
        E1001 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[1, 0, 0, 1]).all(axis=1)][:, 0].mean()
        
        E1110 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[1, 1, 1, 0]).all(axis=1)][:, 0].mean()
        E1101 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[1, 1, 0, 1]).all(axis=1)][:, 0].mean()
        E1011 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[1, 0, 1, 1]).all(axis=1)][:, 0].mean()
        E0111 = conditionedGenes_np[(conditionedGenes_np[:, [1, 2, 3, 4]]==[0, 1, 1, 1]).all(axis=1)][:, 0].mean()
        
        num = E1111*E1100*E1010*E0110*E0101*E0011*E1001*E0000 * (1-E0111)*(1-E1011)*(1-E1101)*(1-E1110)*(1-E0001)*(1-E0010)*(1-E0100)*(1-E1000)
        denom = (1-E1111)*(1-E1100)*(1-E1010)*(1-E0110)*(1-E0101)*(1-E0011)*(1-E1001)*(1-E0000)*(E0111*E1011*E1101*E1110*E0001*E0010*E0100*E1000)

    else:
        print('Order not yet implemented, change estimation method to probabilities.')
        return np.nan

    if ((num==0) & (denom==0)):
            return np.nan 
    elif num==0:
        return -np.inf
    elif denom==0:
        return np.inf
    else:
        return np.log(num/denom)
    
@njit
def safeMean(a):
    if len(a)>0:
        return a.mean()
    else:
        return 0

@njit
def calcInteraction_expectations_numba(conditionedGenes_np):
    '''
    Calc the interactions from the conditioned data using expectation values. 
    Could see if using probabilities is faster. 
    Currently, orders are implemented separately.

    NOTE: previously, each interactions was corrected by a factor to match the {-1, 1} basis of the Ising model.
    I now just take logs and leave them in the {0, 1} basis, which makes more sense for gene expression.
    '''
    if (len(conditionedGenes_np)==0):
        return np.nan
    
    order = conditionedGenes_np.shape[-1]
    
    
    
    if(order==1):
        E = conditionedGenes_np[:].mean()
        num = E
        denom = 1-E
        
    elif(order==2):
        E1 = safeMean(conditionedGenes_np[conditionedGenes_np[:, 1]==1][:, 0])
        E0 = safeMean(conditionedGenes_np[conditionedGenes_np[:, 1]==0][:, 0])

        num = E1*(1-E0)
        denom = E0*(1-E1)
        
    elif(order==3):
        E11 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==1)][:, 0])
        E00 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==0)][:, 0])
        E10 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==0)][:, 0])
        E01 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==1)][:, 0])
        num = E11*(1-E01)*E00*(1-E10)
        denom = E01*(1-E11)*E10*(1-E00)
        
    elif(order==4):
        E111 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==1) & (conditionedGenes_np[:, 3]==1)][:, 0])
        E000 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==0) & (conditionedGenes_np[:, 3]==0)][:, 0])
        
        E001 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==0) & (conditionedGenes_np[:, 3]==1)][:, 0])
        E010 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==1) & (conditionedGenes_np[:, 3]==0)][:, 0])
        E100 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==0) & (conditionedGenes_np[:, 3]==0)][:, 0])
        
        E011 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==1) & (conditionedGenes_np[:, 3]==1)][:, 0])
        E101 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==0) & (conditionedGenes_np[:, 3]==1)][:, 0])
        E110 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==1) & (conditionedGenes_np[:, 3]==0)][:, 0])
        
       
        num = E111*(1-E011)*(1-E101)*E001*E010*(1-E110)*E100*(1-E000)
        denom = (1-E111)*E011*E101*(1-E001)*(1-E010)*E110*(1-E100)*E000
        
    elif(order==5):
        E1111 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==1) & (conditionedGenes_np[:, 3]==1) & (conditionedGenes_np[:, 4]==1)][:, 0])  
        E0000 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==0) & (conditionedGenes_np[:, 3]==0) & (conditionedGenes_np[:, 4]==0)][:, 0])  
        
        E0001 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==0) & (conditionedGenes_np[:, 3]==0) & (conditionedGenes_np[:, 4]==1)][:, 0])  
        E0010 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==0) & (conditionedGenes_np[:, 3]==1) & (conditionedGenes_np[:, 4]==0)][:, 0])  
        E0100 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==1) & (conditionedGenes_np[:, 3]==0) & (conditionedGenes_np[:, 4]==0)][:, 0])  
        E1000 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==0) & (conditionedGenes_np[:, 3]==0) & (conditionedGenes_np[:, 4]==0)][:, 0])  
        
        E0011 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==0) & (conditionedGenes_np[:, 3]==1) & (conditionedGenes_np[:, 4]==1)][:, 0])  
        E0101 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==1) & (conditionedGenes_np[:, 3]==0) & (conditionedGenes_np[:, 4]==1)][:, 0])  
        E0110 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==1) & (conditionedGenes_np[:, 3]==1) & (conditionedGenes_np[:, 4]==0)][:, 0])  
        E1010 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==0) & (conditionedGenes_np[:, 3]==1) & (conditionedGenes_np[:, 4]==0)][:, 0])  
        E1100 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==1) & (conditionedGenes_np[:, 3]==0) & (conditionedGenes_np[:, 4]==0)][:, 0])  
        E1001 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==0) & (conditionedGenes_np[:, 3]==0) & (conditionedGenes_np[:, 4]==1)][:, 0])  
        
        E1110 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==1) & (conditionedGenes_np[:, 3]==1) & (conditionedGenes_np[:, 4]==0)][:, 0])  
        E1101 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==1) & (conditionedGenes_np[:, 3]==0) & (conditionedGenes_np[:, 4]==1)][:, 0])  
        E1011 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==1) & (conditionedGenes_np[:, 2]==0) & (conditionedGenes_np[:, 3]==1) & (conditionedGenes_np[:, 4]==1)][:, 0])  
        E0111 = safeMean(conditionedGenes_np[(conditionedGenes_np[:, 1]==0) & (conditionedGenes_np[:, 2]==1) & (conditionedGenes_np[:, 3]==1) & (conditionedGenes_np[:, 4]==1)][:, 0])  
        
        num = E1111*E1100*E1010*E0110*E0101*E0011*E1001*E0000 * (1-E0111)*(1-E1011)*(1-E1101)*(1-E1110)*(1-E0001)*(1-E0010)*(1-E0100)*(1-E1000)
        denom = (1-E1111)*(1-E1100)*(1-E1010)*(1-E0110)*(1-E0101)*(1-E0011)*(1-E1001)*(1-E0000)*(E0111*E1011*E1101*E1110*E0001*E0010*E0100*E1000)

    else:
        print('Order not yet implemented, change estimation method to probabilities.')
        return np.nan
    
    if ((num==0) & (denom==0)):
            return np.nan 
    elif num==0:
        return -np.inf
    elif denom==0:
        return np.inf
    else:
        return np.log(num/denom)


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
    
    # This function is relatively slow compared to using expectation values, but will work for any order of interactions. 

    order = len(conditionedGenes.columns)
    nStates = 2**order

    # This assigns every state to numerator or denominator, depending on the number of 1s:
    # In numerator: states where the number of ones has same parity as order itself. 
    # In denom: When this is not the case.
    powers = 2*np.array([np.base_repr(i).count('1')%2==order%2 for i in range(2**order)]).astype(float)-1
    
    f = lambda x: ''.join(map(str, x))
    binCounts = np.bincount(list(map(lambda x: int(x, 2), list(map(f, conditionedGenes.values)))), minlength=nStates)
        
    return np.log(np.prod(np.array([x**p for (x, p) in zip(binCounts, powers)])))  


def calcInteraction_withCI_andBounds(genes, graph, dataSet, estimator, nResamps=1000):
    '''
    Add 95% confidence interval bounds from bootstrap resamples,
    and the F value: the proportion of resamples with a different sign.

    Note that to check for function equality, you need to use bytecode
    '''
    
    if estimator.__code__.co_code == calcInteraction_expectations.__code__.co_code:
        # if PrintBool: print('Detected old')
        MBmode = '0' # Use first gene to get MB
    elif estimator.__code__.co_code == calcInteraction_expectations_np.__code__.co_code:
        # if PrintBool: print('Detected numpy')
        MBmode = '0' # Use first gene to get MB
    elif estimator.__code__.co_code == calcInteraction_expectations_numba.__code__.co_code:
        # if PrintBool: print('Detected numba')
        MBmode = '0' # Use first gene to get MB

    else:
        if PrintBool: print('Warning, using ALL mode for MB selection, make sure you want this. ')
        MBmode = 'All' # Use MB of all genes -- safer, so used as else statement. 

    conditionedGenes = conditionOnMB(genes, graph, dataSet, mode=MBmode)
    
    

    if estimator.__code__.co_code == calcInteraction_expectations_numba.__code__.co_code:
        val0 = estimator(conditionedGenes.values)

    else:
        val0 = estimator(conditionedGenes)
    vals = np.zeros(nResamps)
    
    # Stores if estimate is real val: 0, or UB/LB: 1/-1
    boundVal = 0

    if np.isnan(val0):
        # Then both num and denom are zero and we can't do anything
        return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, genes]


    if np.isinf(val0):

        # if empirical is +/- inf, then we can add artificial cells to put bounds:

        order = len(genes)
        # This assigns every state to numerator or denominator, depending on the number of 1s:
        # In numerator: states where the number of ones has same parity as order itself. 
        # In denom: When this is not the case.
        # powers = 2*np.array([np.base_repr(i).count('1')%2==order%2 for i in range(2**order)]).astype(float)-1
        nStates = 2**order
        f = lambda x: ''.join(map(str, x))
        binCounts = np.bincount(list(map(lambda x: int(x, 2), list(map(f, conditionedGenes.values)))), minlength=nStates)

        # find the binary rep of the state that was missing, and see if we can put upper/lower bound
        boundVal = -1*int(np.sign(val0))
        statesToAdd = np.array([np.array(list(np.binary_repr(i, order))).astype(int) for i in range(2**order)])[np.where(binCounts==0)[0]]
        conditionedGenes = conditionedGenes.append(pd.DataFrame(statesToAdd, columns=conditionedGenes.columns), ignore_index=True)
        
    
    if estimator.__code__.co_code == calcInteraction_expectations_numba.__code__.co_code:
        rng = np.random.default_rng()
        conditionedGenes_np = conditionedGenes.values
        for i in range(nResamps):
            resampled = conditionedGenes_np[rng.choice(len(conditionedGenes_np), len(conditionedGenes_np), replace=True)]
            vals[i] = estimator(resampled)
        vals.sort()

    else:
        for i in range(nResamps):
            genes_resampled = conditionedGenes.sample(frac=1, replace=True)
            vals[i] = estimator(genes_resampled)
        vals.sort()
    vals_noNan = vals[~np.isnan(vals)]

    # ksStat = kstest(vals_noNan, lambda x: scipy.stats.norm.cdf(x, loc=vals_noNan.mean(), scale=vals_noNan.std()))[1]

    CI = (vals_noNan[int(np.around(len(vals_noNan)/40))], vals_noNan[int(np.floor(len(vals_noNan)*39/40))])

    propDifSign = sum(np.sign(vals_noNan)==-np.sign(val0))/len(vals_noNan)
    propUndefined = sum(np.isnan(vals))/len(vals)
    propInfinite = sum(np.isinf(vals))/len(vals)



    return [val0, CI[0], CI[1], propDifSign, propUndefined, propInfinite, boundVal, genes]


    # # If it's *really* close to a unimodal distribution according to KS test, or doesn't have undef. resamples:

    # # if((len(vals_noNan) == nResamps) | (ksStat>0.01)) : 
    # if(len(vals_noNan) >= nResamps*0.95):
    #     return [val0, CI[0], CI[1], propDifSign, propUndefined, propInfinite, genes, boundVal]
    # else:
    #     return [np.nan, np.nan, np.nan, np.nan, genes, boundVal]      


    
def calcInteraction_withCI_parallel(args, nResamps=1000):
    '''
    wrapper to unpack function arguments so that it can be mapped over process pool with one arg.
    (I actually think there is something like executor.starmap that could do this for us)
    '''
    genes, graph, dataSet, estimator, nResamps = args
    
    return calcInteraction_withCI_andBounds(genes, graph, dataSet, estimator, nResamps=nResamps)       
                  
        
        
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
    
    
    



