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
parser.add_argument("--clusters", type=str, nargs='?', help="Path to file with cluster annotation")
parser.add_argument("--nGenes", type=int, nargs=1, help="Number of genes to keep")
parser.add_argument("--userGenes", type=str, nargs='?', help="List of genes to always include")
parser.add_argument("--nCells", type=int, nargs=1, help="Number of cells to keep")
parser.add_argument("--cluster", type=int, nargs='?', help="Which cluster/cell Type to use")
parser.add_argument("--bcDoublets", type=str, nargs='?', help="Path to file with booleans for doublets")
parser.add_argument("--twoReplicates", type=int, nargs='?', help="Boolean for constructing two replicates")

args = parser.parse_args()


rawData = args.rawData[0]
nGenes = int(args.nGenes[0])
nCells = int(args.nCells[0])
twoReps = bool(args.twoReplicates)

try:
    print('loading user-defined genes')
    userGenes = pd.read_csv(args.userGenes).columns.values
except:
    print('NOTE: continuing without user-defined genes')
    userGenes = np.array([])


# ****************** data-agnostic QC analysis and selection: ******************


if args.dataType=='agnostic':


    print("loading data")
    scObj = sc.read_csv(rawData)
    print("data loaded!")

    scObj.var_names_make_unique()
    scObj.obs['index'] = np.arange(len(scObj))


    try:
        print('loading user-defined genes')
        userGenes = pd.read_csv(args.userGenes).columns.values
        userGenes = [g for g in userGenes if g in scObj.var.index]
        nGenes = max(nGenes, len(userGenes))
    except:
        userGenes = np.array([])


    print('adding doublet and cluster data')


    try:
        print('loading cluster data')
        clusters = pd.read_csv(args.clusters, index_col=0)
        scObj.obs['cluster'] = clusters.values
        cl = int(args.cluster)

    except:
        print('NOTE: continuing without cluster file.')
        scObj.obs['doublet'] = False
        scObj.obs['cluster'] = 1
        cl = 1

    try:
        print('loading doublet data')
        doubs = pd.read_csv(args.bcDoublets[0], index_col=0)
        scObj.obs['doublet'] = doubs
    except:
        print('NOTE: continuing without doublet annotation.')
        scObj.obs['doublet'] = False
        
    scObj = scObj[(scObj.obs['doublet']==False) & (scObj.obs['cluster']==cl)]


    sc.pp.normalize_total(scObj, target_sum=1e6)
    sc.pp.log1p(scObj)
    scObj.raw = scObj
    sc.pp.highly_variable_genes(scObj, min_mean=0.0125, max_mean=7, min_disp=0.2)
    sc.pl.highly_variable_genes(scObj, save=f'QC_HVG_selection_CL'+'{:0>2}'.format(cl) + '_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.png')
    print('selected genes: ', sum(scObj.var['highly_variable']))

    # # ------------ add embedding coords -------------
    # print('Running PCA...')
    # sc.tl.pca(scObj, n_comps=2)

    # print('Running UMAP...')
    # sc.pp.neighbors(scObj)
    # sc.tl.umap(scObj)
    # print('DONE with embeddings.')

    hvgObj = scObj[:,scObj.var['highly_variable'].values]
    sorted_HVG = hvgObj.var.sort_values('dispersions_norm', ascending=False).index
    selected_genes = np.hstack([userGenes, [s for s in sorted_HVG if not s in userGenes]])
    selected_genes = selected_genes[:nGenes]
    print('Number of genes selected:   ', selected_genes.shape)
    scObjBin = scObj.copy()
    scObjBin.X = (scObjBin.X>0)*1

    sc.pp.subsample(scObjBin, fraction=1.) #Shufle full cluster so that any selection is randomised. 


    selectedCellsAndGenes = scObjBin[:,scObjBin.var.index.isin(selected_genes)]
    clDF = pd.DataFrame(selectedCellsAndGenes.X.toarray())
    clDF.columns = selectedCellsAndGenes.var.index
    print('Final QCd data set size: ', clDF.shape)

    print('Two replicates?  ', twoReps)
    if twoReps:
        clDF.iloc[:nCells].to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS1_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.csv', index=False)
        clDF.iloc[nCells:2*nCells].to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS2_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.csv', index=False)

    else:
        clDF.iloc[:nCells].to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS1_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.csv', index=False)


    # pd.DataFrame(scObjBin.obsm['X_pca'][:nCells]).to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS1_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_PCAcoords.csv', index=False)
    # pd.DataFrame(scObjBin.obsm['X_pca'][nCells:2*nCells]).to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS2_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_PCAcoords.csv', index=False)

    # pd.DataFrame(scObjBin.obsm['X_umap'][:nCells]).to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS1_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_UMAPcoords.csv', index=False)
    # pd.DataFrame(scObjBin.obsm['X_umap'][nCells:2*nCells]).to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS2_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_UMAPcoords.csv', index=False)


      
    print('****DONE****')











# ****************** 10X QC analysis and selection: ******************
elif args.dataType=='10X':
    rawData = args.rawData[0]
    nGenes = args.nGenes
    nCells = args.nCells
    cluster = args.cluster
    bcDoublets = args.bcDoublets[0]
    print("loading data")
    scObj = sc.read_10x_h5(rawData)
    scObj.var_names_make_unique()
    scObj.obs['index'] = np.arange(len(scObj))
    print("data loaded!")


    # ---------- Standard QC -----------
    print('adding doublet and cluster data')
    clusters = pd.read_csv(args.clusters, index_col=0)
    nGenes = int(args.nGenes)
    nCells = int(args.nCells)
    cl = int(args.cluster)

    scObj.obs['cluster'] = clusters

    try:
        print('loading doublet data')
        doubs = pd.read_csv(bcDoublets, index_col=0)
        scObj.obs['doublet'] = doubs
    except:
        print('NOTE: continuing without doublet annotation.')
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
                 jitter=0.4, multi_panel=True, save='QC_n_genes_CL'+'{:0>2}'.format(cl) + '_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.png')

    sc.pl.violin(clObj, ['percent_mito'],
                 jitter=0.4, multi_panel=True, save='QC_percent_mito_CL'+'{:0>2}'.format(cl) + '_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.png')


    # ------------ Making data with only HVG for each cluster -------------


    sc.pp.normalize_total(clObj, target_sum=1e4)
    sc.pp.log1p(clObj)
    clObj.raw = clObj
    sc.pp.highly_variable_genes(clObj, min_mean=0.0125, max_mean=7, min_disp=0.2)
    sc.pl.highly_variable_genes(clObj, save=f'QC_HVG_selection_CL'+'{:0>2}'.format(cl) + '_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.png')
    print('selected genes: ', sum(clObj.var['highly_variable']))

    # ------------ add embedding coords -------------
    print('Running PCA...')
    sc.tl.pca(clObj, n_comps=2)

    print('Running UMAP...')
    sc.pp.neighbors(clObj)
    sc.tl.umap(clObj)
    print('DONE with embeddings.')

    hvgObj = clObj[:,clObj.var['highly_variable'].values]
    sorted_HVG = hvgObj.var.sort_values('dispersions_norm', ascending=False).index
    genes = sorted_HVG[:nGenes]
    selected_genes = np.unique(genes)
    print('Final selected genes:   ', selected_genes.shape)
    clObjBin = clObj.copy()
    clObjBin.X = (clObjBin.X>0)*1

    sc.pp.subsample(clObjBin, fraction=1.) #Shufle full cluster so that any selection is randomised. 


    selectedCellsAndGenes = clObjBin[:,clObjBin.var.index.isin(selected_genes)]
    clDF = pd.DataFrame(selectedCellsAndGenes.X.todense())
    print('Final data set size: ', clDF.shape)
    clDF.columns = selectedCellsAndGenes.var.index

    clDF.iloc[:nCells].to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS1_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.csv', index=False)
    clDF.iloc[nCells:2*nCells].to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS2_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes.csv', index=False)

    pd.DataFrame(clObjBin.obsm['X_pca'][:nCells]).to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS1_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_PCAcoords.csv', index=False)
    pd.DataFrame(clObjBin.obsm['X_pca'][nCells:2*nCells]).to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS2_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_PCAcoords.csv', index=False)

    pd.DataFrame(clObjBin.obsm['X_umap'][:nCells]).to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS1_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_UMAPcoords.csv', index=False)
    pd.DataFrame(clObjBin.obsm['X_umap'][nCells:2*nCells]).to_csv('trainingData_CL'+'{:0>2}'.format(cl)+ '_DS2_' + '{:0>5}'.format(nCells) + 'Cells_'+'{:0>4}'.format(nGenes) + 'Genes_UMAPcoords.csv', index=False)


      
    print('****DONE****')







# ****************** ZEISEL selection: ******************

elif args.dataType=='Zeisel':
    nCells = args.nCells
    nGenes = args.nGenes
    neurons = sc.read_loom(args.rawData[0]+"/ZEISEL_neuronsCNS.loom")
    astrocytes = sc.read_loom(args.rawData[0]+"/ZEISEL_astrocytes.loom")
    # neuroblasts_G = sc.read_loom(args.rawData+"/ZEISEL_neuroblastsGlut.loom")
    # neuroblasts_NG = sc.read_loom(args.rawData+"/ZEISEL_neuroblastsNon-Glut.loom")

    neurons.obs_names_make_unique()
    neurons.var_names_make_unique()

    astrocytes.obs_names_make_unique()
    astrocytes.var_names_make_unique()

    # neuroblasts_G.obs_names_make_unique()
    # neuroblasts_G.var_names_make_unique()

    # neuroblasts_NG.obs_names_make_unique()
    # neuroblasts_NG.var_names_make_unique()


    nG = '{:0>4}'.format(str(nGenes))
    nC = '{:0>5}'.format(str(nCells))
    genesNeurons = np.array(pd.read_csv(args.rawData[0]+f'/trainingData_CL07_DS1_{nC}Cells_{nG}Genes.csv').columns)
    genesAstros = np.array(pd.read_csv(args.rawData[0]+f'/trainingData_CL10_DS1_{nC}Cells_{nG}Genes.csv').columns)


    def f(x):
        if x=='':
            return 99999
        else: return int(x)

    neuronsYoungerThanP22 = neurons[neurons.obs['Age'].str[1:3].apply(f)<22]

    neurons500Genes = pd.DataFrame(data = neuronsYoungerThanP22.X[:, neuronsYoungerThanP22.var.index.isin(genesNeurons)].todense(), 
                            columns = neuronsYoungerThanP22.var.index[neuronsYoungerThanP22.var.index.isin(genesNeurons)])
    neurons500Genes = neurons500Genes[genesNeurons]

    astros500Genes = pd.DataFrame(data = astrocytes.X[:, astrocytes.var.index.isin(genesAstros)].todense(), 
                            columns = astrocytes.var.index[astrocytes.var.index.isin(genesAstros)])
    astros500Genes = astros500Genes[genesAstros]


    neurons_combined = neurons500Genes.sample(20000, replace=False).reset_index(drop=True)
    neurons_combined_bin = (neurons_combined>0)*1

    astros_combined_bin = (astros500Genes>0)*1

    N = args.nCells
    neurons_combined_bin.iloc[:N].to_csv(f'trainingData_ZEISEL_neurons_DS1_{nC}Cells_{nG}Genes.csv', index=False)
    neurons_combined_bin.iloc[N:2*N].to_csv(f'trainingData_ZEISEL_neurons_DS2_{nC}Cells_{nG}Genes.csv', index=False)

    N = min(args.nCells, int(len(astros_combined_bin)/2))

    print(f'Keeping {N} astrocytes')    
    astros_combined_bin.iloc[:N].to_csv(f'trainingData_ZEISEL_astrocytes_DS1_{nC}Cells_{nG}Genes.csv', index=False)
    astros_combined_bin.iloc[N:2*N].to_csv(f'trainingData_ZEISEL_astrocytes_DS2_{nC}Cells_{nG}Genes.csv', index=False)

else:
    print('ERROR: invalid data type, choose 10X or Zeisel.')
    sys.exit()




