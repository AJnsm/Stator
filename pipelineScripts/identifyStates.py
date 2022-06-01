import pandas as pd
import numpy as np
import igraph as ig
import argparse
import sys
import scipy
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cut_tree

import matplotlib.pyplot as plt
from matplotlib import cm
import io
from PIL import Image
import seaborn as sns
sns.set(style='darkgrid')
sns.set_palette('colorblind')

print('Modules imported \n')

parser = argparse.ArgumentParser(description='Args for coupling estimation')

parser.add_argument("--dataPath", type=str, help="Path to training data")
parser.add_argument("--PCApath", type=str, help="Path to PCA embedding of training data")
parser.add_argument("--devStates", type=str, help="Path to list of most deviating states")
parser.add_argument("--diffCutoff", type=float, help="Dice distance to slice dendrogram")

args = parser.parse_args()
diffCutoff = args.diffCutoff
trainDat = pd.read_csv(args.dataPath)
pcaCoords= pd.read_csv(args.PCApath)

devStates = pd.read_csv(args.devStates, dtype=str)
if len(devStates)==0:
    print('no deviating states, terminating...')
    sys.exit()

devStates.columns = ['genes', 'state', 'dev']

# Binreps is the binary represenations of the interactions: binReps[i] is 1 if cell i is in the maxDevState, 0 otherwise.  
binReps = np.array(devStates.apply(lambda x: (trainDat[x['genes'].rsplit('_')]==[int(x) for x in list(str(x['state']))]).all(axis=1), axis=1))*1

# linkage defines the distances between the binReps, using the Dice-distance: https://en.wikipedia.org/wiki/Sørensen–Dice_coefficient
linked = linkage(binReps, 'average', metric='dice')

# Each interaction is put in a cluster by cutting the dendrogram at a threshold
devStates['cluster'] = fcluster(linked, diffCutoff, criterion = 'distance')


# Plotting interaction clusters separately:

alph=0.8 # plotting transparancy
nCol=4 # number of subfigure columns
nCl = max(devStates['cluster'])-1
fig, ax = plt.subplots(int(np.ceil(nCl/nCol)), nCol, figsize=[7*nCol, 1.3*nCol*int(np.ceil(nCl/nCol))])
r=0
for c in range(1, max(devStates['cluster'])+1):
    
    if len(ax.flatten())<6:
        a = ax[(c-1)%nCol]
    else:
        a = ax[r, (c-1)%nCol]

    a.plot(pcaCoords.values[:, 0], pcaCoords.values[:, 1], '.', color=(0.5, 0.5, 0.5, 0.05))
    for i, row in devStates[devStates['cluster']==c].iterrows():
        genes = row['genes'].rsplit('_')
        state = [int(x) for x in list(str(row['state']))]
        labels = [''.join(x) for x in list(zip(genes, ['+' if s else '-' for s in state]))]
        charCells = (trainDat[genes]==state).all(axis=1)
        a.plot(pcaCoords.values[charCells, 0], pcaCoords.values[charCells, 1], 'o', alpha=alph, label = ", ".join(labels))
    if c%nCol==0:
        r+=1

for a in ax.flatten():            
    a.set_xticks([])
    a.set_yticks([])
    a.legend()

plt.savefig('distinctDeviatingStates.png')
plt.close(fig)


# Plotting interactions in dendrogram representation:
# I plot the embedding of each cluster
# These embeddings are saved in a buffer and stored in a dict, and then used as axis-labels when plotting the dendrogram. 

def fromImToArr(img):
    '''
    Takes an ImageFile and converts it to a a numpy array of floats so that imshow() can plot it.  
    '''
    return np.array(img.getdata()).reshape(img.height, img.width, -1)

# Create image labels 
sns.set_style('white')
statePlots = {}
intsInClusters = []

truncatedClusters = max(devStates['cluster']) # The total number of clusters after truncation

R_full = dendrogram(linked, labels = devStates['genes'].values, distance_sort='descending', no_plot=True)
R_trunc =dendrogram(linked, labels = devStates['genes'].values, distance_sort='descending', no_plot=True, p=truncatedClusters, truncate_mode='lastp')

index = 0
for img, geneStr in enumerate(R_trunc['ivl']):
    cl = []
    plt.figure()

    # Truncation renames a cluster to '(n)' if it is the result of merging n clusters
    # To get proper labelling, I check if we are plotting a singleton cluster, or a bigger one
    # The index variable keeps track of how many clusters are plotted so the correct states are added to the correct embeddings.
    # This is a pretty confusing process: after plotting a cluster with n interactions, we increment the index by n. 
    if geneStr in devStates['genes'].values:
        genes = geneStr.rsplit('_')

        # Plot all cells/observations in grey
        plt.plot(pcaCoords.values[:, 0], pcaCoords.values[:, 1], '.', color=(0.5, 0.5, 0.5, 0.05))
        state = [int(x) for x in devStates[devStates['genes']==geneStr]['state'].values[0]]
        charCells = (trainDat[genes]==state).all(axis=1)
        plt.plot(pcaCoords.values[charCells, 0], pcaCoords.values[charCells, 1], 'o', alpha=alph, label = ", ".join(labels))
        cl.append([''.join(s) for s in list(zip(genes, ['+' if x==1 else '-' for x in state]))])
        index+=1

    else:
        # Extract the number of interactings in this cluster/branch from the label it has
        toAdd = int(geneStr.replace('(', '').replace(')', ''))
        
        # Plot all cells/observations in grey    
        plt.plot(pcaCoords.values[:, 0], pcaCoords.values[:, 1], '.', color=(0.5, 0.5, 0.5, 0.05))

        # Loop over all interactions in this cluster/branch
        for i in range(index, index+toAdd):
            genes = R_full['ivl'][i].rsplit('_')
            state = [int(x) for x in devStates[devStates['genes']==R_full['ivl'][i]]['state'].values[0]]
            charCells = (trainDat[genes]==state).all(axis=1)
            plt.plot(pcaCoords.values[charCells, 0], pcaCoords.values[charCells, 1], 'o', alpha=alph, label = ", ".join(labels))
            cl.append([''.join(s) for s in list(zip(genes, ['+' if x==1 else '-' for x in state]))])
        index += toAdd
    
    
    intsInClusters.append(cl)
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
    im = OffsetImage(img, zoom=0.34)
    im.image.axes = ax

    ab = AnnotationBbox(im, (xCoord, yCoord),  xybox=(0., 30.), frameon=False,
                        xycoords='data',  boxcoords="offset points", pad=0)

    ax.add_artist(ab)


sns.set_style('white')

# create list of top 4 occuring genes across interactions in a cluster to add to the label
labs = []
for c in intsInClusters:
    gs = np.array([a for b in c for a in b])
    labs.append('\n'.join(np.unique(gs)[np.argsort(-np.unique(gs, return_counts=True)[1])][:4]))


plt.figure(figsize=(truncatedClusters*2, 5))
R = dendrogram(linked,
            orientation='top',
            truncate_mode='lastp',
            p = truncatedClusters,
            distance_sort='descending',
            show_leaf_counts=True)
plt.ylim(diffCutoff-0.02, 1)

ax = plt.gca()
x_labLocs = [x.get_position()[0] for x in ax.get_xmajorticklabels()]

# Adding the PCA embeddings as labels to the dendrogram plot
for i, xL in enumerate(x_labLocs):
    add_imgLab(xL, diffCutoff-0.02, statePlots[i], ax)
    
plt.xticks(x_labLocs, labels = labs, rotation=0)    

sns.despine(left=True, top=True, right=True, bottom=True)
plt.yticks(np.linspace(diffCutoff, 1.0, 4))
plt.ylabel('Dice-distance')
plt.savefig('distinctDeviatingStates_dendrogram.png')
plt.close(fig)






















