import requests 
import os
import pandas as pd
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import seaborn as sns
viridis = cm.get_cmap('viridis', 12)
sns.set_style("dark")

sns.set_palette('colorblind')
from utilities import *
import collections
from collections import OrderedDict
import itertools
import random

import hypernetx as hnx
import scanpy as sc

import upsetplot as up
from upsetplot import generate_counts
from upsetplot import from_indicators


if PrintBool: print('Modules imported \n')

parser = argparse.ArgumentParser(description='Args for coupling estimation')

parser.add_argument("--dataPath", type=str, nargs=1, help="Path to training data")
parser.add_argument("--PCApath", type=str, nargs=1, help="Path to PCA embedding of training data")
parser.add_argument("--CPDAGgraphPath", type=str, nargs=1, help="Path to CPDAG graph file")
parser.add_argument("--MCMCgraphPath", type=str, nargs=1, help="Path to MCMC graph file")
parser.add_argument("--pathTo2pts", type=str, nargs='?', help="Path to calculated 2-point interactions")
parser.add_argument("--pathTo2pts_CI_F", type=str, nargs='?', help="Path to calculated 2-point interactions: significance")
parser.add_argument("--pathTo2pts_undef", type=str, nargs='?', help="Path to calculated 2-point interactions: number of undef. resamples")
parser.add_argument("--pathTo2pts_inf", type=str, nargs='?', help="Path to calculated 2-point interactions: number of inf. resamples")
parser.add_argument("--pathTo3pts", type=str, nargs='?', help="Path to calculated 3-point interactions")
parser.add_argument("--pathTo4pts", type=str, nargs='?', help="Path to calculated 4-point interactions")
parser.add_argument("--pathTo5pts", type=str, nargs='?', help="Path to calculated 5-point interactions")

args = parser.parse_args()

trainDat = pd.read_csv(dataPath)
pcaCoords= pd.read_csv(PCApath)
DSname = MCMCgraphPath.split('.')[0]
MCMCadjMat = pd.read_csv(MCMCgraphPath, index_col=0)
MCMCgraph = ig.Graph.Adjacency(MCMCadjMat.values.tolist())

CPDAGadjMat = pd.read_csv(CPDAGgraphPath, index_col=0)
CPDAGgraph = ig.Graph.Adjacency(CPDAGadjMat.values.tolist())

genes = trainDat.columns.values
scObj = sc.AnnData(trainDat)
scObj.obsm['X_pca'] = pcaCoords.values


coups_2pts = np.load(pathTo2pts, allow_pickle=True).astype(np.float32)
coups_2pts_CI_F = np.load(pathTo2pts_CI_F, allow_pickle=True).astype(np.float32)
coups_2pts_undef = np.load(pathTo2pts_undef, allow_pickle=True).astype(np.float32)
coups_2pts_inf = np.load(pathTo2pts_inf, allow_pickle=True).astype(np.float32)

coups_2pts[(coups_2pts_undef>0) | (coups_2pts_inf>0)] = np.nan 
coups_2pts_CI_F[(coups_2pts_undef>0) | (coups_2pts_inf>0)] = np.nan 


HHOIs = {}
alpha=0.05
for order, intPath in enumerate([pathTo3pts, pathTo4pts, pathTo5pts]):
	ints = np.load(intPath, allow_pickle=True)
	perfectSigEsts = list(map(lambda x: (((x[[4, 5, 6]]==0).all()) & (x[3]<alpha)), ints))
	HHOIs[f'n{order+3}'] = ints[perfectSigEsts][:, [0, -1]]
			
pairs = np.array(np.where(coups_2pts_CI_F<alpha)).T
vals = coups_2pts[coups_2pts_CI_F<alpha]
HHOIs[f'n2'] = np.array([list(x) for x in list(zip(vals, pairs))])




def findsubsets(s, n):
	return list(itertools.combinations(s, n))

def plotUpsetPlot(d, fig, title='', legend=False, save=False, filename=''):
	genes = d.columns.values
	
	f = lambda x: ''.join(map(str, x))
	
	order = len(d.columns)
	nStates = 2**order

	binCounts = np.bincount(list(map(lambda x: int(x, 2), list(map(f, d.values)))), minlength=nStates)
	upSetObj = generate_counts(1, n_samples=100000, n_categories=order)
	upSetObj.loc[:] = binCounts
	upSetObj.index.names = genes

	means = d.mean(axis=0)
	expected = np.array([np.prod([m if state[i] else 1-m for i, m in enumerate(means)]) for state in upSetObj.index])*len(d)

	upSetObj2 = pd.DataFrame()
	upSetObj2['count'] = upSetObj
	upSetObj2['deviation'] = (upSetObj - expected)/(expected)*100
	upSetObj2['expected'] = expected

	upset = up.UpSet(upSetObj2, subset_size='sum', sum_over='count', intersection_plot_elements=5, min_subset_size=0)

	upset.add_stacked_bars(by='deviation', sum_over='deviation', colors='coolwarm', elements=5)
	upset.add_catplot(value='expected', kind='strip', elements=1)

	for gs in powerset(genes):
		if (len(genes) - len(gs))%2==0:
			upset.style_subsets(present=gs, absent=[g for g in genes if not (g in gs)],
						facecolor="cornflowerblue",
						label="numerator")
		else:
			upset.style_subsets(present=gs, absent=[g for g in genes if not (g in gs)],
						facecolor="indianred",
						label="denominator")


	ax = upset.plot(fig=fig)
	ax['extra1'].set_ylabel('% deviation')
	ax['extra1'].get_legend().remove()
	
	if not legend:
		ax['intersections'].get_legend().remove()
	
	ax['intersections'].set_ylabel('No. of cells')
	ax['intersections'].set_yscale('log')
	ax['intersections'].plot([x.get_offsets()[0][1] for x in ax['extra2'].get_children()[:nStates]], 'k+', zorder=10, label='expected')
	ax['extra2'].remove()

	fig.suptitle(f'{title},   I = ' + 
			str(round(calcInteraction_binTrick_allOrders(d), 3)), fontsize=20)

	if(save):
		plt.savefig(filename)
		plt.close(fig)
	else:
		return




f=2
kwargs = {'with_node_counts': True, 'with_node_labels':True, 'with_edge_labels':False}
concatInts = lambda x: ''.join(map(str, x))
deviations = {}

# dicts to store plots in
plotHypergraph = {}
plotCPDAG = {}
plotPCA = {}
plotupset = {}
plotUpset_cond = {}
plotUpset_uncond = {}
plotMaxDev = {}


for order in [3, 4, 5]:
	nStates = 2**order
	binStates = [np.array(list(format(x, f"0{order}b"))).astype(bool) for x in range(nStates)]
	devs = []

	for w, geneTuple in HHOIs[f'n{order}'][:, [0, -1]]:
		ID = '_'.join(genes[geneTuple])	  
		g = findLocalGraph(geneTuple, CPDAGgraph, order=0)
		layout_c = g.layout('circle')

		tmp = dict(zip(g.vs['label'], np.array(layout_c)*np.array([1, -1])))
		edges = {'maxOrder': genes[geneTuple]}
		weights = [w]

		for i in range(2, order):
			for subset in findsubsets(geneTuple, i):
				for tup in HHOIs[f'n{i}']:
					if sorted(tup[-1]) == sorted(subset):
						edges[len(edges)+1] =  [genes[g] for g in subset]
						weights.append(tup[0])
						break

		#  ************************ Interaction Hypergraph  ************************ 
		H = hnx.Hypergraph(edges)
		cList = ['red' if w<0 else 'green' for w in weights]
		layout_fn = lambda x: tmp
		plt.figure(figsize=[10, 10])

		hnx.draw(H,  layout = layout_fn,
					 label_alpha=0,
					 node_labels_kwargs={
						'fontsize': 24
					},
					 edges_kwargs={
						 'edgecolors': cList,
						 'linewidths': 3,
						 'dr': 0.05
					 },
					 **kwargs)  
		plt.savefig(f'{ID}_HOIs.png')
		buf = io.BytesIO()
		plt.savefig(buf, format='png')
		buf.seek(0)
		plotHypergraph[ID] = Image.open(buf)
		plt.close()

		#  ************************ CPDAG ************************ 

		# ig.plot(g, f'{ID}_CPDAG.png', layout=layout_c, bbox=(600/f, 600/f), vertex_size=120/f, vertex_color='white', margin=100/f)

		fig, ax = plt.subplots(figsize=[6, 6])
		ig.plot(g, layout=layout_c, bbox=(600/f, 600/f), vertex_size=120/f, vertex_color='white', margin=100/f, vertex_label_color='black', target=ax)
		buf = io.BytesIO()
		plt.savefig(buf, format='png')
		buf.seek(0)
		plotCPDAG[ID] = Image.open(buf)
		plt.close()

		#  ************************ PCA embeddings ************************ 
		fig, ax = plt.subplots(1, len(fsplit[0].split('_')), figsize=[20, 4])
		for g in geneTuple:
			sc.pl.embedding(scObj,'pca', color=g, 
								size=30, color_map="viridis", add_outline=True, ncols=len(fsplit[0]), show=False, frameon=False, ax = ax[i])
		fig = plt.gcf()
		fig.axes[-1].remove()
		# plt.savefig(f'{ID}_Expression.png')
		buf = io.BytesIO()
		plt.savefig(buf, format='png')
		buf.seek(0)
		plotPCA[ID] = Image.open(buf)
		plt.close() 

		#  ************************ Upset plots ************************ 

		unConditionedGenes = trainDat.iloc[:, geneTuple]
		conditionedGenes = conditionOnMB(geneTuple, MCMCgraph, trainDat, mode='Min')

		fig.figure(figsize=[10, 10])
		plotUpsetPlot(d = conditionedGenes,fig=fig, legend=False, title = 'Conditioned on MB', filename=ID + '_Upset_conditioned.png', save=True)
		buf = io.BytesIO()
		plt.savefig(buf, format='png')
		buf.seek(0)
		plotUpset_cond[ID] = Image.open(buf)
		plt.close()

		fig.figure(figsize=[10, 10])
		plotUpsetPlot(d = unConditionedGenes,fig=fig, legend=False, title = 'Unonditioned', filename=ID + '_Upset_unconditioned.png', save=True)
		buf = io.BytesIO()
		plt.savefig(buf, format='png')
		buf.seek(0)
		plotUpset_uncond[ID] = Image.open(buf)
		plt.close()
		

		#  ************************ Calculate deviations ************************ 

		conditionedGenes = conditionOnMB(interactors, graph, dat, mode='Min')
		binCounts = np.bincount(list(map(lambda x: int(x, 2), list(map(concatInts, conditionedGenes.values)))), minlength=nStates)
		means = conditionedGenes.mean(axis=0)
		expected = np.array([np.prod([m if state[i] else 1-m for i, m in enumerate(means)]) for state in binStates])*len(conditionedGenes)
		
		deviation = (binCounts - expected)/(expected)
		devs.append([deviation, interactors])
	
	devs  = np.array(devs, dtype=object)
	devs = devs[(-np.array(list(map(np.max, devs[:, 0])))).argsort()]
	deviations[f'n{order}'] = devs

		
#  ************************ PCA embedding on max deviating state ************************ 
for order in [3, 4, 5]:
	for devs, interactors in deviations[f'n{order}']:
		ID = '_'.join(genes[interactors])
		
		maxDevState = format(np.argmax(devs), f"0{order}b")
		maxDevState_embedded = scObj.obsm['X_pca'][(scObj.X[:, interactors] == np.array(list(maxDevState)).astype(float)).all(axis=1)]
		xs = scObj.obsm['X_pca'][:, 0]
		ys = scObj.obsm['X_pca'][:, 1]
		plt.plot(xs, ys, 'o', color = viridis(0), alpha=0.1)
		plt.plot(maxDevState_embedded[:, 0], maxDevState_embedded[:, 1], 'o', color = viridis(0.99), alpha=0.9)
		plt.title(', '.join(genes[interactors]) + ' = ' + ', '.join(maxDevState), fontsize=8)
		plt.xticks([])
		plt.yticks([])

		buf = io.BytesIO()
		plt.savefig(buf, format='png')
		buf.seek(0)
		plotMaxDev[ID] = Image.open(buf)
		plt.close()
		# plt.savefig(f'{ID}_Expression_maxDevState.png')
		plt.close('all') 


#  ************************ Plot summary figures ************************ 

sns.set_style("white")

for order in [3, 4, 5]:
	for w, geneTuple in HHOIs[f'n{order}'][:, [0, -1]]:
		ID = '_'.join(genes[geneTuple])	  
		fig = plt.figure(figsize=(15, 10))

		axCPDAG = fig.add_axes([0, 0.66, 0.33, 0.33])
		axPC = fig.add_axes([0.33, 0.66, 0.33, 0.33])
		axHOI = fig.add_axes([0.66, 0.66, 0.33, 0.33])
		axEXP = fig.add_axes([0., 0.33, 0.99, 0.33])

		axUPS1 = fig.add_axes([0, 0, 0.33, 0.33])
		axUPS2 = fig.add_axes([0.33, 0, 0.33, 0.33])
		axEXP_maxDev = fig.add_axes([0.66, 0.0, 0.33, 0.33])

		axes = [axCPDAG, axPC, axHOI, axEXP, axEXP_maxDev, axUPS1, axUPS2]

		for a in axes:
			a.axis('off')

		axCPDAG.imshow(plotCPDAG)
		axHOI.imshow(plotHypergraph[100:-100, 100:-50])
		axEXP.imshow(plotPCA[:, 400:, :])
		axEXP_maxDev.imshow(plotMaxDev[:, :, :])

		axUPS1.imshow(plotUpset_cond[:, 20:])
		axUPS2.imshow(plotUpset_uncond[:, 20:])

		plt.savefig(f'{ID}_summary.png')
		plt.close(fig) 











