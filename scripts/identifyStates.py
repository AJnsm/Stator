import pandas as pd
import numpy as np
import igraph as ig
import argparse
import sys
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cut_tree
from sklearn.metrics import pairwise_distances
from utilities import *

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.offsetbox import OffsetImage,AnnotationBbox
import io
from PIL import Image
import seaborn as sns
sns.set(style='white')
sns.set_palette('colorblind')

print('Modules imported \n')

parser = argparse.ArgumentParser(description='Args for coupling estimation')

parser.add_argument("--dataPath", type=str, help="Path to training data")
parser.add_argument("--PCApath", type=str, help="Path to PCA embedding of training data")
parser.add_argument("--devStates", type=str, help="Path to list of most deviating states")
parser.add_argument("--diffCutoff", type=float, help="Dice distance to slice dendrogram")
parser.add_argument("--bsResamps", type=int, help="Number of bootstrap resamples")
parser.add_argument("--auThreshold", type=float, help="AU threshold for state identification")


args = parser.parse_args()

bsResamps = args.bsResamps
auThreshold = args.auThreshold
diffCutoff = args.diffCutoff
trainDat = pd.read_csv(args.dataPath)
pcaCoords= pd.read_csv(args.PCApath)
devStates = pd.read_csv(args.devStates, dtype=str, index_col=0)

# Max number of genes in a label of the dendrogram:
maxInLab = 5

# If there are no deviating states, terminate
if len(devStates)==0:
    print('no deviating states, terminating...')
    sys.exit()

if len(devStates)==1:
    print('Only one deviation state, terminating...')
    sys.exit()


devStates.columns = ['genes', 'state', 'enrichment', 'pval', 'pval_corrected', 'CellIDs']

# Add the gene states to the devStates dataframe
def makeGeneState(dtupleDFrow):
    genes = dtupleDFrow.genes.split('_')
    state = list(dtupleDFrow.state)
    zipped =  list(zip(genes, ['+' if int(s) else '-' for s in state]))
    return [a + b for (a, b) in zipped]
    
devStates['geneState'] = devStates.apply(makeGeneState, axis=1)


# Binreps is the binary represenations of the interactions: binReps[i] is 1 if cell i is in the deviating state, 0 otherwise.  
binReps = np.array(devStates.apply(lambda x: (trainDat[x['genes'].rsplit('_')]==[int(g) for g in list(str(x['state']))]).all(axis=1), axis=1))*1
n = len(binReps)


# Labels that combine the genes and their states---once as list, once as string with newlines
labsWithStates = devStates.apply(lambda x: [''.join(g) for g in list(zip(x['genes'].split('_'), ['+' if int(s)==1 else '-' for s in x['state']]))], axis=1)
labsWithStates_str = labsWithStates.apply(lambda x: '\n'.join(x)).values

# add the cell counts to the labels:
cellCounts = binReps.sum(axis=1)
labsWithStates_str = list(map(lambda x: x[0]+f'\n({x[1]} cells)', list(zip(labsWithStates_str, cellCounts))))


print('Calculating linkage matrix...')
# linkage defines the distances between the binReps, using the Dice-distance: https://en.wikipedia.org/wiki/Sørensen–Dice_coefficient
linked_full = linkage(binReps, 'average', metric='dice')
print('Linkage matrix calculated \n')

# pd.DataFrame(linked_full).to_csv('fullLinkageMatrix.csv')
# pd.DataFrame(labsWithStates_str).to_csv('fullLinkageMatrix_labels.csv')


def fromImToArr(img):
    '''
    Takes an ImageFile and converts it to a a numpy array of floats so that imshow() can plot it.  
    '''
    return np.array(img.getdata()).reshape(img.height, img.width, -1)

def add_imgLab(xCoord, yCoord, img, ax, zoom = 0.25):
    '''
    Adds an image as a label to the ax.
    '''
    # Crop the image a little bit
    img = img[5:-5, 5:-5, :]
    im = OffsetImage(img, zoom=zoom)
    im.image.axes = ax

    ab = AnnotationBbox(im, (xCoord, yCoord),  xybox=(0., 30.), frameon=False,
                        xycoords='data',  boxcoords="offset points", pad=0)

    ax.add_artist(ab)

def extract_levels(row_clusters, labels=np.arange(len(binReps))):
    '''
    Extracts the individual states in each branch of the dendrogram.
    '''
    labels=np.arange(len(row_clusters)+1)
    clusters = {}
    for row in range(row_clusters.shape[0]):
        cluster_n = row + len(labels)
        # which clusters / labels are present in this row
        glob1, glob2 = row_clusters[row, 0], row_clusters[row, 1]

        # if this is a cluster, pull the cluster
        this_clust = []
        for glob in [glob1, glob2]:
            if glob > (len(labels)-1):
                this_clust += clusters[glob]
            # if it isn't, add the label to this cluster
            else:
                this_clust.append(glob)

        clusters[cluster_n] = this_clust
    return clusters


# Create image labels for the full dendrogram, stores these in a buffer
alph=0.8
sns.set_style('white')
statePlots = {}
d = dendrogram(linked_full, labels = labsWithStates_str, distance_sort='descending', no_plot=True)
index = 0

# Try to plot the dendrogram, though this can fail due to memory issues if the dendrogram is too large
try:
    for img, geneStr in enumerate(d['ivl']):
        cl = []
        plt.figure(figsize=[6, 3])

        geneState = geneStr.rsplit('\n')[:-1]
        genes = [x[:-1] for x in geneState]
        state = [1 if x[-1]=='+' else 0 for x in geneState]

        # Plot all cells/observations in grey
        plt.plot(pcaCoords.values[:, 0], pcaCoords.values[:, 1], '.', color=(0.5, 0.5, 0.5, 0.05))
        charCells = (trainDat[genes]==state).all(axis=1)
        plt.plot(pcaCoords.values[charCells, 0], pcaCoords.values[charCells, 1], 'o', alpha=alph)

        plt.xticks([])
        plt.yticks([])
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight',pad_inches = 0)
        buf.seek(0)
        statePlots[img] = fromImToArr(Image.open(buf))
        plt.close()
except:
    print('Error in plotting the dendrogram, most likely caused by too many states in the dendrogram.')


try:
    # Only calculate the bootstrap statistics if there are not too many d-tuples:
    if len(labsWithStates_str) < 250:
        # Run bootstrap resampling to get the bootstrap statistics of the dendrogram.
        X = pd.DataFrame(binReps.T)
        pv = PvClust(X, method="average", metric="dice", nboot=bsResamps, parallel=True)
        pvalues = pv._result[['AU', 'BP']].values
        pv.result.to_csv('bootstrapStats.csv')
        print('Done with first clustering')

        # clusterDict keeps track of which original states are contained in which cluster index:
        originalclusterDict = {i:[i] for i in range(n)}
        newClusterDict = extract_levels(linked_full)
        clusterDict = {**originalclusterDict, **newClusterDict}

        # Export the bootstrap statistics to a csv file
        stateDF = []
        for i, merge in enumerate(linked_full):
                stateDF.append([labsWithStates[clusterDict[n+i]].values, pvalues[i][0], pvalues[i][1]])
        stateDF = pd.DataFrame(stateDF)
        stateDF.columns = ['State', 'AU', 'BP']
        stateDF.to_csv('statesWithBootstrapStats.csv')

except:
    print('Error in bootstrapping the dendrogram, most likely caused by too many states in the dendrogram.')

xcoord = d["icoord"]
ycoord = d["dcoord"]
# Obtaining the coordinates of all nodes above leaves
x = {i: (j[1]+j[2])/2 for i, j in enumerate(xcoord)}
y = {i: j[1] for i, j in enumerate(ycoord)}
pos = node_positions(y, x)

# Try to plot the dendrogram, though this can fail due to memory issues if the dendrogram is too large
try:
    # Plot the full dendrogram with the bootstrap values:
    plt.figure(figsize=(len(labsWithStates_str)+20, 10))
    plt.tight_layout()
    set_link_color_palette(['c', 'g'])
    d = dendrogram(linked_full, labels = labsWithStates_str, distance_sort='descending', above_threshold_color='c',
                color_threshold=1.0)


    ax = plt.gca()
    # Adding the bootstrap stats to the dendrogram plot
    for node, (x, y) in pos.items():

        if node == (len(pos.items())-1):
            ax.text(x-6, y, 'AU', fontsize=14, fontweight='bold',
                    color='black')
            ax.text(x+1, y, 'BP', fontsize=14, fontweight='bold',
                    color='black')
        else:
            if pvalues[node][0] >= auThreshold:
                ax.text(x-5, y, f' {pvalues[node][0]*100:.0f}', fontsize=10,
                        color='green', fontweight='bold')
                ax.text(x+1, y, f'{pvalues[node][1]*100:.0f}', fontsize=10,
                        color='green', fontweight='bold')
            else:
                ax.text(x-5, y, f' {pvalues[node][0]*100:.0f}', fontsize=10,
                        color='black')
                ax.text(x+1, y, f'{pvalues[node][1]*100:.0f}', fontsize=10,
                        color='black')

    ax = plt.gca()
    x_labLocs = [x.get_position()[0] for x in ax.get_xmajorticklabels()]

    # Adding the PCA embeddings as labels to the dendrogram plot
    for i, xL in enumerate(x_labLocs):
        add_imgLab(xL, -0.03, statePlots[i], ax, zoom=0.18)
        
    plt.xticks(x_labLocs, labels = d['ivl'], rotation=0, fontsize=10)    
    sns.despine(left=True, top=True, right=True, bottom=True)
    plt.ylabel('Dice-distance')
    plt.ylim(-0.03, 1)
    plt.savefig('dendrogram_dTuples.png', bbox_inches='tight')
    plt.close()

except:
    print('Error in plotting the dendrogram, most likely caused by too many states in the dendrogram.')


# Finally, create the states that result from a simple cutoff:



def modularity_score(adjMat, cluster_labels, verbose=False):
    '''
    Calculate the modularity score for a given clustering of a graph.
    '''
    n = len(adjMat)
    m = np.sum(adjMat) / 2.0
    k = np.sum(adjMat, axis=1)
    q = 0.0
    for i in range(n):
        for j in range(n):
            if (cluster_labels[i] == cluster_labels[j]) & (i!=j):
                q += (adjMat[i][j] - k[i]*k[j]/(2.0*m))       
    if verbose:
        print(f'n:{n}')
        print(f'm:{m}')
        print(f'k:{k}')
        print(f'q:{q}')
    return q/(2.0*m)

# function to calculate the cluster labels for a given cutoff
cutAt = lambda x: fcluster(linked_full, x, criterion = 'distance')


pairwiseDists = pairwise_distances(binReps, metric='dice')

# Modularity calculation for a range of N cutoffs (set to 50 for round numbers, but also forced to be round to 2 decimal places)
N = 50
cutoffs = np.linspace(0.01, 0.99, N)
cutoffs = np.round(cutoffs, 2)

modScores = [modularity_score((1-pairwiseDists), cutAt(d)) for d in cutoffs]


# The cutoff is either what the user specified, or what maximises the modularity score:
if diffCutoff == -1:
    diffCutoff = cutoffs[np.argmax(modScores)]
    print(f'Using cutoff of {diffCutoff} to maximise modularity score')

print(f'Max modularity score: {max(modScores)}')
pd.DataFrame(zip(cutoffs, modScores), columns=['Cutoff', 'Modularity score']).to_csv('modularity_scores.csv')


# Each interaction is put in a cluster by cutting the dendrogram at a threshold
devStates['cluster'] = fcluster(linked_full, diffCutoff, criterion = 'distance')

devStates.to_csv(f'top_DTuples_withStatesFromCut.csv')

# Try to plot the dendrogram, though this can fail due to memory issues if the dendrogram is too large
try:

    # Create image labels 
    statePlots = {}
    truncatedClusters = max(devStates['cluster']) # The total number of clusters after truncation


    R_full = dendrogram(linked_full, labels = labsWithStates_str, distance_sort='descending', no_plot=True)
    R_trunc =dendrogram(linked_full, labels = labsWithStates_str, distance_sort='descending', no_plot=True, p=truncatedClusters, truncate_mode='lastp')

    intsInClusters = []
    cellsInClusters = []
    index = 0
    for img, geneStr in enumerate(R_trunc['ivl']):
        cl = []
        cellNumber = 0
        plt.figure(figsize=[6, 3])

        # Truncation renames a cluster to '(n)' if it is the result of merging n clusters
        # To get proper labelling, I check if we are plotting a singleton cluster, or a bigger one
        # The index variable keeps track of how many clusters are plotted so the correct states are added to the correct embeddings.
        # This is a pretty confusing process: after plotting a cluster with n interactions, we increment the index by n. 
        if geneStr in labsWithStates_str:
            # Remove the last entry, which is the cell count of that state. 
            geneState = geneStr.rsplit('\n')[:-1]
            genes = [x[:-1] for x in geneState]
            state = [1 if x[-1]=='+' else 0 for x in geneState]
            
            # Plot all cells/observations in grey
            plt.plot(pcaCoords.values[:, 0], pcaCoords.values[:, 1], '.', color=(0.5, 0.5, 0.5, 0.05))
            charCells = (trainDat[genes]==state).all(axis=1)
            plt.plot(pcaCoords.values[charCells, 0], pcaCoords.values[charCells, 1], 'o', alpha=alph)
            cl.append([''.join(s) for s in list(zip(genes, ['+' if x==1 else '-' for x in state]))])
            cellNumber += charCells.sum()
            index+=1

        else:
            # Extract the number of interactings in this cluster/branch from the label it has
            toAdd = int(geneStr.replace('(', '').replace(')', ''))
            
            # Plot all cells/observations in grey    
            plt.plot(pcaCoords.values[:, 0], pcaCoords.values[:, 1], '.', color=(0.5, 0.5, 0.5, 0.05))

            # Loop over all interactions in this cluster/branch
            charCellArray = np.zeros(trainDat.shape[0], dtype=bool)
            for i in range(index, index+toAdd):
                geneState = R_full['ivl'][i].rsplit('\n')[:-1]
                genes = [x[:-1] for x in geneState]
                state = [1 if x[-1]=='+' else 0 for x in geneState]
                
                charCells = (trainDat[genes]==state).all(axis=1)
                plt.plot(pcaCoords.values[charCells, 0], pcaCoords.values[charCells, 1], 'o', alpha=alph)
                cl.append([''.join(s) for s in list(zip(genes, ['+' if x==1 else '-' for x in state]))])
                # The cells corresponding to each sub branch has to be ORred together to get the cells in the full branch.
                # Ignoring this will result in double counting of cells.
                charCellArray = charCellArray | charCells

            cellNumber += charCellArray.sum()
            index += toAdd
        
        
        intsInClusters.append(cl)
        cellsInClusters.append(cellNumber)
        plt.xticks([])
        plt.yticks([])
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight',pad_inches = 0)
        buf.seek(0)
        statePlots[img] = fromImToArr(Image.open(buf))
        plt.close()


    def add_imgLab(xCoord, yCoord, img, ax):
        '''
        Adds an image from statePlots to the ax.
        '''
        
        # Crop the image a little bit
        img = img[5:-5, 5:-5, :]
        im = OffsetImage(img, zoom=0.25)
        im.image.axes = ax

        ab = AnnotationBbox(im, (xCoord, yCoord),  xybox=(0., 30.), frameon=False,
                            xycoords='data',  boxcoords="offset points", pad=0)

        ax.add_artist(ab)


    sns.set_style('white')

    # create list of top occuring genes across interactions in a cluster to add to the label

    labs = []
    for (cellNr, intList) in list(zip(cellsInClusters, intsInClusters)):
        gs = np.array([a for b in intList for a in b])
        sortedGenes = list(np.unique(gs)[np.argsort(-np.unique(gs, return_counts=True)[1])])
        if len(sortedGenes) > maxInLab:
            sortedGenes = sortedGenes[:maxInLab] + [f'...']
        sortedGenes = sortedGenes + [f'({len(intList)} d-tuples)'] + [f'({cellNr} cells)']
        labs.append('\n'.join(sortedGenes))

    plt.figure(figsize=(truncatedClusters*2, 5))
    R = dendrogram(linked_full,
                orientation='top',
                truncate_mode='lastp',
                p = truncatedClusters,
                distance_sort='descending',
                show_leaf_counts=True)
    plt.ylim(diffCutoff-0.03, 1)

    ax = plt.gca()
    x_labLocs = [x.get_position()[0] for x in ax.get_xmajorticklabels()]

    # Adding the PCA embeddings as labels to the dendrogram plot
    for i, xL in enumerate(x_labLocs):
        add_imgLab(xL, diffCutoff-0.03, statePlots[i], ax)
        
    plt.xticks(x_labLocs, labels = labs, rotation=0)    

    sns.despine(left=True, top=True, right=True, bottom=True)
    plt.yticks(np.linspace(diffCutoff, 1.0, 4))
    plt.ylabel('Dice-distance')
    plt.savefig('dendrogram_dTuples_cut.png', bbox_inches='tight')
    plt.close()

except:
    print('Error in plotting the dendrogram, most likely caused by too many states in the dendrogram.')

