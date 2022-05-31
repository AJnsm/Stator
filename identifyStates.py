import pandas as pd
import numpy as np
import igraph as ig
import argparse

import scipy
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cut_tree

import matplotlib.pyplot as plt
from matplotlib import cm

import seaborn as sns
sns.set(style='darkgrid')
sns.set_palette('colorblind')



# input: data, PCA embeddings, deviating states. 

print('Modules imported \n')

parser = argparse.ArgumentParser(description='Args for coupling estimation')

parser.add_argument("--dataPath", type=str, nargs=1, help="Path to training data")
parser.add_argument("--PCApath", type=str, nargs=1, help="Path to PCA embedding of training data")
parser.add_argument("--devStates", type=str, nargs=1, help="Path to list of most deviating states")
parser.add_argument("--diffCutoff", type=str, nargs=1, help="Dice distance to slice dendrogram")

args = parser.parse_args()
diffCutoff = args.diffCutoff
trainDat = pd.read_csv(dataPath)
pcaCoords= pd.read_csv(PCApath)

devStates = pd.read_csv(devStates, dtype=str)
devStates.columns = ['genes', 'state', 'dev']


binReps = np.array(devStates.apply(lambda x: (trainDat[x['genes'].rsplit('_')]==[int(x) for x in list(str(x['state']))]).all(axis=1), axis=1))*1
linked = linkage(binReps, 'average', metric='dice')
devStates['cluster'] = fcluster(linked, diffCutoff, criterion = 'distance')
alph=0.8
nCol=4
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



# Create image labels 

sns.set_style('white')


statePlots = {}
intsInClusters = []

truncatedClusters = max(devStates['cluster'])
R_full = dendrogram(linked, labels = devStates[ds]['genes'].values, distance_sort='descending', no_plot=True)
R_trunc =dendrogram(linked, labels = devStates[ds]['genes'].values, distance_sort='descending', no_plot=True, p=truncatedClusters, truncate_mode='lastp')

index = 0
for img, geneStr in enumerate(R_trunc['ivl']):
    cl = []
    plt.figure()
    if geneStr in devStates[ds]['genes'].values:
        genes = geneStr.rsplit('_')
        plt.plot(pcaCoords.values[:, 0], pcaCoords.values[:, 1], '.', color=(0.5, 0.5, 0.5, 0.05))
        state = [int(x) for x in devStates[ds][devStates[ds]['genes']==geneStr]['state'].values[0]]
        charCells = (data[ds][genes]==state).all(axis=1)
        plt.plot(pcaCoords.values[charCells, 0], pcaCoords.values[charCells, 1], 'o', alpha=alph, label = ", ".join(labels))
        cl.append([''.join(s) for s in list(zip(genes, ['+' if x==1 else '-' for x in state]))])
        index+=1

    else:
        toAdd = int(geneStr.replace('(', '').replace(')', ''))
        plt.plot(pcaCoords.values[:, 0], pcaCoords.values[:, 1], '.', color=(0.5, 0.5, 0.5, 0.05))

        for i in range(index, index+toAdd):
            genes = R_full['ivl'][i].rsplit('_')
            state = [int(x) for x in devStates[ds][devStates[ds]['genes']==R_full['ivl'][i]]['state'].values[0]]
            charCells = (data[ds][genes]==state).all(axis=1)
            plt.plot(pcaCoords.values[charCells, 0], pcaCoords.values[charCells, 1], 'o', alpha=alph, label = ", ".join(labels))
            cl.append([''.join(s) for s in list(zip(genes, ['+' if x==1 else '-' for x in state]))])
        index += toAdd
    
    
    intsInClusters.append(cl)
    plt.xticks([])
    plt.yticks([])
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight',pad_inches = 0)
    buf.seek(0)
    statePlots[img] = Image.open(buf)
    plt.close()


def add_imgLab(xCoord, yCoord, ID, ax, ds):
    img = statePlots[ID]
    
    img = np.array(img.getdata()).reshape(img.height, img.width, -1)
    img = img[5:-5, 5:-5, :]
    im = OffsetImage(img, zoom=0.34)
    im.image.axes = ax

    ab = AnnotationBbox(im, (xCoord, yCoord),  xybox=(0., 30.), frameon=False,
                        xycoords='data',  boxcoords="offset points", pad=0)

    ax.add_artist(ab)


sns.set_style('white')

# create list of top 4 occuring genes across interactions in a cluster
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
for i, xL in enumerate(x_labLocs):
    add_imgLab(xL, diffCutoff-0.02, i, ax, ds)
    
plt.xticks(x_labLocs, labels = labs, rotation=0)    

sns.despine(left=True, top=True, right=True, bottom=True)
plt.yticks(np.linspace(diffCutoff, 1.0, 4))
plt.ylabel('Dice-distance')
plt.savefig('distinctDeviatingStates_dendrogram.png')
plt.close(fig)






















