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


clusters = pd.read_csv(sys.argv[2], index_col=0)
nGenes = int(sys.argv[3])
nCells = int(sys.argv[4])
cellType = int(sys.argv[5])


try:
    doubs = pd.read_csv(sys.argv[6], index_col=0)
except:
    doubs=None
    
print("loading data")
scObj = sc.read_10x_h5(sys.argv[1])
scObj.var_names_make_unique()
scObj.obs['index'] = np.arange(len(scObj))
print("data loaded!")
scObj.obs['doublet'] = False


# THIS IS A TEMPORARY CLAUSE, NO DOUBLETS ARE DETECTED
# SHOULD BE REPLACED BY MEMORY FRIENDLY DOUB DETECTION.
print('WARNING: aborting doublet detection, insuff. memory. Different solution should be found.')
if False:

    # ---------- Calculate doublets -----------
    if ((len(doubs)<len(scObj)) | (doubs is None)):
        

        print('constructing count matrix')
        counts_matrix = scObj.X
        print('starting doublet detection')
        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.076)
        print('calling doublets')
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                                  min_cells=3,
                                                                  min_gene_variability_pctl=85,
                                                                  n_prin_comps=40)

        # NOTE: could choose to do in log-space by adding log_transform=True

        print('No. of doublets detected: ', sum(scrub.call_doublets()))
        scObj.obs['doublet'] = scrub.call_doublets()
        doubs = pd.DataFrame(scObj.obs[['doublet']])
        doubs.to_csv("bcDoublets.csv", index=0)
        print('DONE') 

        fig, ax = scrub.plot_histogram()
        fig.savefig('doubletScore_hist.png')

        print('Running UMAP...')
        scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
        print('Done with UMAP.')


        x = scrub._embeddings['UMAP'][:,0]
        y = scrub._embeddings['UMAP'][:,1]

        fig = plt.figure(figsize = (5, 5))
        coldat = scrub.predicted_doublets_
        o = np.argsort(coldat)
        plt.scatter(x[o], y[o], c = coldat[o], cmap=scr.custom_cmap([[.7,.7,.7], [0,0,0]]), s = 2)
        plt.xticks([]), plt.yticks([])
        plt.title('Doublets: UMAP embedding')
        plt.xlabel('UMAP x')
        plt.ylabel('UMAP y')
        plt.tight_layout()
        fig.savefig('doublets.png', dpi=None)



# ---------- Standard QC -----------

# First: remove cells with high mito, low counts, and doublets.
sc.pp.filter_cells(scObj, min_genes=400)
sc.pp.filter_genes(scObj, min_cells=0.02*len(scObj))

mito_genes = scObj.var_names.str.startswith('mt-')
scObj.obs['percent_mito'] = np.sum(
    scObj[:, mito_genes].X, axis=1).A1 / np.sum(scObj.X, axis=1).A1
scObj.obs['n_counts'] = scObj.X.sum(axis=1).A1


sc.pl.violin(scObj, ['n_genes'],
             jitter=0.4, multi_panel=True, save='QC_n_genes.png')

sc.pl.violin(scObj, ['n_counts'],
             jitter=0.4, multi_panel=True, save='QC_n_counts.png')

sc.pl.violin(scObj, ['percent_mito'],
             jitter=0.4, multi_panel=True, save='QC_percent_mito.png')


# ------------ Making data with only HVG for each cluster -------------

cl = cellType
print(cellType)
print('\n')
print(scObj.shape)
print('\n')
print('Making data from cluster ', cl)
#     clObj = scObj[scObj.obs['doublet']==False]
clusters_noDoubs = clusters[doubs['doublet']==False]
clObj = scObj[clusters_noDoubs.index[clusters_noDoubs['Cluster']==cl]]

sc.pp.normalize_total(clObj, target_sum=1e4)
sc.pp.log1p(clObj)
clObj.raw = clObj
sc.pp.highly_variable_genes(clObj, min_mean=0.0125, max_mean=7, min_disp=0.2)
sc.pl.highly_variable_genes(clObj, save=f'QC_HVG_selection_CL{cl}.png')
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
clDF.columns = selectedCellsAndGenes.var.index
clDF = clDF.sample(frac=1).reset_index(drop=True) #Shufle full cluster so that any selection is randomised. 
#         clDF.to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.csv', index=False)
clDF.iloc[:nCells].to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS1_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.csv', index=False)
clDF.iloc[nCells:2*nCells].to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS2_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.csv', index=False)

#     clDF.iloc[:nCells].apply(np.random.permutation, axis=1, result_type='broadcast').to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS1_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_shuffled.csv', index=False)
#     clDF.iloc[nCells:2*nCells].apply(np.random.permutation, axis=1, result_type='broadcast').to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS2_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_shuffled.csv', index=False)
    
print('****DONE****')





