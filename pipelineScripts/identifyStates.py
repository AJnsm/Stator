import pandas as pd
import numpy as np
import igraph as ig
import argparse
import sys
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cut_tree
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

if len(devStates)==0:
    print('no deviating states, terminating...')
    sys.exit()

if len(devStates)==1:
    print('Only one deviation state, terminating...')
    sys.exit()

devStates.columns = ['genes', 'state', 'enrichment', 'pval', 'pval_corrected', 'CellIDs']

# Binreps is the binary represenations of the interactions: binReps[i] is 1 if cell i is in the deviating state, 0 otherwise.  
binReps = np.array(devStates.apply(lambda x: (trainDat[x['genes'].rsplit('_')]==[int(g) for g in list(str(x['state']))]).all(axis=1), axis=1))*1
n = len(binReps)


# Labels that combine the genes and their states---once as list, once as string with newlines
labsWithStates = devStates.apply(lambda x: [''.join(g) for g in list(zip(x['genes'].split('_'), ['+' if int(s)==1 else '-' for s in x['state']]))], axis=1)
labsWithStates_str = labsWithStates.apply(lambda x: '\n'.join(x)).values

# linkage defines the distances between the binReps, using the Dice-distance: https://en.wikipedia.org/wiki/Sørensen–Dice_coefficient
linked_full = linkage(binReps, 'average', metric='dice')


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
for img, geneStr in enumerate(d['ivl']):
    cl = []
    plt.figure(figsize=[6, 3])

    geneState = geneStr.rsplit('\n')
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


d = dendrogram(linked_full, labels = labsWithStates_str, distance_sort='descending', no_plot=True)
xcoord = d["icoord"]
ycoord = d["dcoord"]
# Obtaining the coordinates of all nodes above leaves
x = {i: (j[1]+j[2])/2 for i, j in enumerate(xcoord)}
y = {i: j[1] for i, j in enumerate(ycoord)}
pos = node_positions(y, x)


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
plt.savefig('dendrogram_all_dTuples.png', bbox_inches='tight')
plt.close()



# Keep all states that have significant bootstrap stats
AU = pvalues[:, 0]
states = []
for i, merge in enumerate(linked_full):
    if AU[i] >=auThreshold:
        states.append(clusterDict[n+i])
        
#       if one of the two merged clusters is original and get merged with a signficiant cluster, add it separately as well
        if ((merge[0]<n) & (merge[1]>=n)):
            if AU[int(merge[1]-(n))]>=auThreshold:
                states.append([int(merge[0])])
            
        elif ((merge[1]<n) & (merge[0]>=n)):
            if AU[int(merge[0]-(n))]>=auThreshold:
                states.append([int(merge[1])])

states = [[int(x) for x in state] for state in states ]
singletonStates = [[i] for i in range(n) if not i in set.union(*[set(state) for state in states])]
for s in singletonStates:
    states.append(s)

# Create binary representation of the significant states:
fullStateBinReps = []
for stateList in states:
    stateList_binRep = [binReps[state] for state in stateList]
    fullStateBinRep = np.array(stateList_binRep).any(axis=0)*1
    fullStateBinReps.append(fullStateBinRep)
fullStateBinReps = np.array(fullStateBinReps)


if len(fullStateBinReps)==0:
    print('No robust combined states, terminating...')
    sys.exit()

if len(fullStateBinReps)==1:
    print('Only one robust combined state, terminating...')
    sys.exit()


# Redo the bootstrap analysis with the significant states:
X = pd.DataFrame(fullStateBinReps.T)
pv = PvClust(X, method="average", metric="dice", nboot=bsResamps, parallel=True)
print('Done with second clustering')

# Create the labels for the significant states from the top 6 occuring genes
labs = []
for state in states:
    genes = list(map(lambda x: [''.join(g) for g in list(zip(devStates.iloc[x]['genes'].split('_'), ['+' if int(s)==1 else '-' for s in devStates.iloc[x]['state']]))], state))
    gs = np.array([a for b in genes for a in b])
    labs.append('\n'.join(np.unique(gs)[np.argsort(-np.unique(gs, return_counts=True)[1])][:6]))

# Create the linkage matrix and dendrogram for the significant states only:
linked_sig = linkage(fullStateBinReps, 'average', metric='dice')
d = dendrogram(linked_sig, labels = labs, distance_sort='descending', no_plot=True)

# Create image labels for the significant states:
sns.set_style('white')
statePlots = {}
alph=0.8
for img, leaf in enumerate(d['leaves']):    
    plt.figure(figsize=[6, 3])
    plt.plot(pcaCoords.values[:, 0], pcaCoords.values[:, 1], '.', color=(0.5, 0.5, 0.5, 0.05))
    
    for stateIndex in states[leaf]:
        genes = devStates.iloc[stateIndex]['genes'].rsplit('_')
        state = [int(x) for x in devStates.iloc[stateIndex]['state']]
        charCells = (trainDat[genes]==state).all(axis=1)
        plt.plot(pcaCoords.values[charCells.astype(bool), 0], pcaCoords.values[charCells.astype(bool), 1], 'o', alpha=alph)
        
    plt.xticks([])
    plt.yticks([])
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight',pad_inches = 0)
    buf.seek(0)
    statePlots[img] = fromImToArr(Image.open(buf))
    plt.close()

xcoord = d["icoord"]
ycoord = d["dcoord"]
# Obtaining the coordinates of all nodes above leaves
x = {i: (j[1]+j[2])/2 for i, j in enumerate(xcoord)}
y = {i: j[1] for i, j in enumerate(ycoord)}
pos = node_positions(y, x)

plt.figure(figsize=[len(fullStateBinReps)+20, 20])     
d = dendrogram(linked_sig, labels=labs, above_threshold_color='c',
               color_threshold=0.0, leaf_font_size=15)

pvalues = pv._result[['AU', 'BP']].values
ax = plt.gca()
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
plt.savefig('dendrogram_reclustered_robust_states.png', bbox_inches='tight')
plt.close()




# Finally, create the states that result from a simple cutoff:

# Each interaction is put in a cluster by cutting the dendrogram at a threshold
devStates['state'] = fcluster(linked_full, diffCutoff, criterion = 'distance')

devStates.to_csv(f'top_DTuples_withStatesFromCut.csv')

# Create image labels 
statePlots = {}
truncatedClusters = max(devStates['state']) # The total number of clusters after truncation

# Labels that combine the genes and their states---once as list, once as string with underscores

R_full = dendrogram(linked_full, labels = labsWithStates_str, distance_sort='descending', no_plot=True)
R_trunc =dendrogram(linked_full, labels = labsWithStates_str, distance_sort='descending', no_plot=True, p=truncatedClusters, truncate_mode='lastp')

intsInClusters = []
index = 0
for img, geneStr in enumerate(R_trunc['ivl']):
    cl = []
    plt.figure(figsize=[6, 3])

    # Truncation renames a cluster to '(n)' if it is the result of merging n clusters
    # To get proper labelling, I check if we are plotting a singleton cluster, or a bigger one
    # The index variable keeps track of how many clusters are plotted so the correct states are added to the correct embeddings.
    # This is a pretty confusing process: after plotting a cluster with n interactions, we increment the index by n. 
    if geneStr in labsWithStates_str:
        geneState = geneStr.rsplit('\n')
        genes = [x[:-1] for x in geneState]
        state = [1 if x[-1]=='+' else 0 for x in geneState]
        
        # Plot all cells/observations in grey
        plt.plot(pcaCoords.values[:, 0], pcaCoords.values[:, 1], '.', color=(0.5, 0.5, 0.5, 0.05))
        charCells = (trainDat[genes]==state).all(axis=1)
        plt.plot(pcaCoords.values[charCells, 0], pcaCoords.values[charCells, 1], 'o', alpha=alph)
        cl.append([''.join(s) for s in list(zip(genes, ['+' if x==1 else '-' for x in state]))])
        index+=1

    else:
        # Extract the number of interactings in this cluster/branch from the label it has
        toAdd = int(geneStr.replace('(', '').replace(')', ''))
        
        # Plot all cells/observations in grey    
        plt.plot(pcaCoords.values[:, 0], pcaCoords.values[:, 1], '.', color=(0.5, 0.5, 0.5, 0.05))

        # Loop over all interactions in this cluster/branch
        for i in range(index, index+toAdd):
            geneState = R_full['ivl'][i].rsplit('\n')
            genes = [x[:-1] for x in geneState]
            state = [1 if x[-1]=='+' else 0 for x in geneState]
            
            charCells = (trainDat[genes]==state).all(axis=1)
            plt.plot(pcaCoords.values[charCells, 0], pcaCoords.values[charCells, 1], 'o', alpha=alph)
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
    im = OffsetImage(img, zoom=0.25)
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
plt.savefig('dendrogram_all_dTuples_cut.png', bbox_inches='tight')
plt.close()
