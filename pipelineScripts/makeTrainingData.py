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
import argparse
np.random.seed(0)
print('libs loaded')
# print('If figures come out blank, change the matplotlib backend. ')

def none_or_str(value):
    if value == 'None':
        return None
    return value


parser = argparse.ArgumentParser(description='Prepare the training data')

parser.add_argument("--dataType", type=str, nargs='?', help="From which experiment does the data come (10X/Zeisel)")
parser.add_argument("--rawData", type=str, nargs=1, help="Path to the raw data file")
parser.add_argument("--nGenes", type=int, nargs=1, help="Number of genes to keep")
parser.add_argument("--userGenes", type=str, nargs='?', help="List of genes to always include")
parser.add_argument("--nCells", type=int, nargs=1, help="Number of cells to keep")
parser.add_argument("--bcDoublets", type=str, nargs='?', help="Path to file with booleans for doublets")
parser.add_argument("--fracMito", type=float, nargs='?', help="Max percentage of mitochondrial transcripts")
parser.add_argument("--minGenes", type=int, nargs='?', help="Min number of genes expressed in a cell")
parser.add_argument("--minCells", type=int, nargs='?', help="Min number of cells in which a gene is present")

args = parser.parse_args()

rawData = args.rawData[0]
nGenes = int(args.nGenes[0])
nCells = int(args.nCells[0])



# ****************** data-agnostic: minimal manipulation  ******************
if args.dataType=='agnostic':

    print("loading data")
    scObj = sc.read_csv(rawData)
    print("data loaded!")

    scObj.var_names_make_unique()
    scObj.obs['index'] = np.arange(len(scObj))


    try:
        print('loading user-defined genes')
        print('file with user genes:  ', args.userGenes)
        userGenes = pd.read_csv(args.userGenes).columns.values
        userGenes = [g for g in userGenes if g in scObj.var.index.values]
        nGenes = max(nGenes, len(userGenes))
    except Exception as e:
        print(e)
        print('NOTE: continuing without user-defined genes')
        userGenes = np.array([])


    print('adding doublet data')

    try:
        print('loading doublet data')
        doubs = pd.read_csv(args.bcDoublets[0])
        scObj.obs['doublet'] = doubs
    except:
        print('NOTE: continuing without doublet annotation.')
        scObj.obs['doublet'] = False
        
    scObj = scObj[scObj.obs['doublet']==False]


# ------------ add embedding coords -------------
    print('Running PCA...')
    sc.tl.pca(scObj, n_comps=2)

    print('Running UMAP...')
    sc.pp.neighbors(scObj)
    sc.tl.umap(scObj)
    print('DONE with embeddings.')




    selected_genes = np.hstack([userGenes, [s for s in scObj.var.index.values if not s in userGenes]])
    selected_genes = selected_genes[:nGenes]
    print('Number of genes selected:   ', selected_genes.shape)

    sc.pp.subsample(scObj, fraction=1.) #Shufle full cluster so that any selection is randomised. 

    scObj[:nCells].write('unbinarised_cell_data.h5ad')
    
    scObjBin = scObj.copy()
    scObjBin.X = (scObjBin.X>0)*1

    selectedCellsAndGenes = scObjBin[:,scObjBin.var.index.isin(selected_genes)]

    clDF = pd.DataFrame(selectedCellsAndGenes.X.toarray())
    clDF.columns = selectedCellsAndGenes.var.index
    print('Final QCd data set size: ', clDF.shape)

    clDF.iloc[:nCells].to_csv('trainingData_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.csv', index=False)
    pd.DataFrame(scObjBin.obsm['X_pca'][:nCells]).to_csv('trainingData_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_PCAcoords.csv', index=False)
    pd.DataFrame(scObjBin.obsm['X_umap'][:nCells]).to_csv('trainingData_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_UMAPcoords.csv', index=False)

      
    print('****DONE****')




# ****************** General expression-data QC analysis and selection: ******************


elif args.dataType=='expression':

    print("loading data")
    scObj = sc.read_csv(rawData)
    print("data loaded!")

    scObj.var_names_make_unique()
    scObj.obs['index'] = np.arange(len(scObj))


    try:
        print('loading user-defined genes')
        print('file with user genes:  ', args.userGenes)
        userGenes = pd.read_csv(args.userGenes).columns.values
        userGenes = [g for g in userGenes if g in scObj.var.index]

        if len(userGenes)>nGenes:
            print('WARNING: more user-defined genes than expected, I am including them all!')

        nGenes = max(nGenes, len(userGenes))


    except Exception as e:
        print(e)
        print('NOTE: continuing without user-defined genes')
        userGenes = np.array([])

    print('adding doublet data')

    try:
        print('loading doublet data')
        doubs = pd.read_csv(args.bcDoublets[0])
        scObj.obs['doublet'] = doubs
    except:
        print('NOTE: continuing without doublet annotation.')
        scObj.obs['doublet'] = False

    scObj = scObj[scObj.obs['doublet']==False]

    print('Starting QC')
    print(f'Filtering -- Keeping cells with less than {args.fracMito} mitochondrial transcripts...')
    print(f'Filtering -- Keeping cells with at least {args.minGenes} genes expressed...')
    print(f'Filtering -- Keeping genes expressed in at least {args.minCells} cells...')

    sc.pp.filter_cells(scObj, min_genes=args.minGenes)
    sc.pp.filter_genes(scObj, min_cells=args.minCells)

    mito_genes = scObj.var_names.str.startswith('mt-')
    scObj.obs['percent_mito'] = np.sum(
        scObj[:, mito_genes].X, axis=1) / np.sum(scObj.X, axis=1)

    scObj = scObj[scObj.obs['percent_mito']<(args.fracMito)]

    sc.pl.violin(scObj, ['n_genes'],
                 jitter=0.4, multi_panel=True, save='QC_n_genes_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.png')

    sc.pl.violin(scObj, ['percent_mito'],
                 jitter=0.4, multi_panel=True, save='QC_percent_mito_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.png')



    print('Normalise & log-transform...')
    sc.pp.normalize_total(scObj, target_sum=1e6)
    sc.pp.log1p(scObj)
    scObj.raw = scObj

    print('Selecting highly variable genes...')
#     Min mean is low enough to include all genes
#     Max mean 14, because that is ln(10^6)
#     Min dispersion doesn't matter because only the top N are taken anyway
    sc.pp.highly_variable_genes(scObj, min_mean=0.01, max_mean=14, min_disp=0.0)
    sc.pl.highly_variable_genes(scObj, save=f'QC_HVG_selection_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.png')
    print('selected genes: ', sum(scObj.var['highly_variable']))

    # ------------ add embedding coords -------------
    print('Running PCA...')
    sc.tl.pca(scObj, n_comps=2)

    print('Running UMAP...')
    sc.pp.neighbors(scObj)
    sc.tl.umap(scObj)
    print('DONE with embeddings.')

    hvgObj = scObj[:,scObj.var['highly_variable'].values]
    sorted_HVG = hvgObj.var.sort_values('dispersions_norm', ascending=False).index
    selected_genes = np.hstack([userGenes, [s for s in sorted_HVG if not s in userGenes]])
    selected_genes = selected_genes[:nGenes]
    print('Number of genes selected:   ', selected_genes.shape)

    print('Shuffling data...')
    sc.pp.subsample(scObj, fraction=1.) #Shufle data so that any selection is randomised. 

    print('Saving unbinarised expression of selected cells...')
    scObj[:nCells].write('unbinarised_cell_data.h5ad')
    
    print('Binarising data...')
    scObjBin = scObj.copy()
    scObjBin.X = (scObjBin.X>0)*1

    selectedCellsAndGenes = scObjBin[:,scObjBin.var.index.isin(selected_genes)]

    print('Saving final binarised count matrix of selected cells and genes...')
    clDF = pd.DataFrame(selectedCellsAndGenes.X.toarray())
    clDF.columns = selectedCellsAndGenes.var.index
    print('Final QCd data set size: ', clDF.shape)

    # Output the selected cells and genes:
    clDF.iloc[:nCells].to_csv('trainingData_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.csv', index=False)

    print('Saving embedding coordinates...')
    # Output the embedding coordinates of the selected cells:
    pd.DataFrame(scObjBin.obsm['X_pca'][:nCells]).to_csv('trainingData_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_PCAcoords.csv', index=False)
    pd.DataFrame(scObjBin.obsm['X_umap'][:nCells]).to_csv('trainingData_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_UMAPcoords.csv', index=False)

    print('****DONE****')
    
else:
    print('ERROR: invalid dataType, choose agnostic or expression mode.')
    sys.exit()




