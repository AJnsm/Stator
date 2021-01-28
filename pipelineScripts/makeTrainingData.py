import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scrublet as scr
import scanpy as sc
sc.settings.figdir = ''
import pandas as pd
import scipy.io
import sys
import os
import seaborn
np.random.seed(0)
print('libs loaded')
# print('If figures come out blank, change the matplotlib backend. ')

print("loading data")
scObj = sc.read_10x_h5(sys.argv[1])
scObj.var_names_make_unique()
scObj.obs['index'] = np.arange(len(scObj))
print("data loaded!")



# ---------- Standard QC -----------
print('adding doublet and cluster data')
clusters = pd.read_csv(sys.argv[2], index_col=0)
nGenes = int(sys.argv[3])
nCells = int(sys.argv[4])
cl = int(sys.argv[5])

scObj.obs['cluster'] = clusters

try:
    print('loading doublet data')
    doubs = pd.read_csv(sys.argv[6], index_col=0)
    scObj.obs['doublet'] = doubs
except:
    print('WARNING: continuing without doublet annotation.')
    scObj.obs['doublet'] = False
    
clObj = scObj[(scObj.obs['doublet']==False) & (scObj.obs['cluster']==cl)]

def filterForMouseB(bcs):
    libs_mouseB = np.array(list(map(lambda x: (x.split('-'))[1], bcs))).astype('int')>=70
    return bcs[libs_mouseB]

clObj = clObj[filterForMouseB(clObj.obs.index)]


# Remove cells with high mito, low n_genes.
sc.pp.filter_cells(clObj, min_genes=600)

mito_genes = clObj.var_names.str.startswith('mt-')
clObj.obs['percent_mito'] = np.sum(
    clObj[:, mito_genes].X, axis=1).A1 / np.sum(clObj.X, axis=1).A1

clObj = clObj[clObj.obs['percent_mito']<.12]

sc.pl.violin(clObj, ['n_genes'],
             jitter=0.4, multi_panel=True, save='QC_n_genes'+'{:0>2}'.format(cl) + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.png')

sc.pl.violin(clObj, ['percent_mito'],
             jitter=0.4, multi_panel=True, save='QC_percent_mito'+'{:0>2}'.format(cl) + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.png')


# ------------ Making data with only HVG for each cluster -------------


sc.pp.normalize_total(clObj, target_sum=1e4)
sc.pp.log1p(clObj)
clObj.raw = clObj
sc.pp.highly_variable_genes(clObj, min_mean=0.0125, max_mean=7, min_disp=0.2)
sc.pl.highly_variable_genes(clObj, save=f'QC_HVG_selection_CL'+'{:0>2}'.format(cl) + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.png')
print('selected genes: ', sum(clObj.var['highly_variable']))
hvgObj = clObj[:,clObj.var['highly_variable'].values]
sorted_HVG = hvgObj.var.sort_values('dispersions_norm', ascending=False).index
genes = sorted_HVG[:nGenes]
selected_genes = np.unique(genes)
print('Final selected genes:   ', selected_genes.shape)
clObjBin = clObj.copy()
clObjBin.X = (clObjBin.X>0)*1

selectedCellsAndGenes = clObjBin[:,clObjBin.var.index.isin(selected_genes)]
clDF = pd.DataFrame(selectedCellsAndGenes.X.todense())
print('Final data set size: ', clDF.shape)
clDF.columns = selectedCellsAndGenes.var.index
clDF = clDF.sample(frac=1).reset_index(drop=True) #Shufle full cluster so that any selection is randomised. 

clDF.iloc[:nCells].to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS1_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.csv', index=False)
clDF.iloc[nCells:2*nCells].to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS2_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.csv', index=False)
  
print('****DONE****')





