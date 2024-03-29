import pandas as pd
import numpy as np
import igraph as ig
import matplotlib
# No interactive backend necessary
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import argparse
import seaborn as sns
viridis = cm.get_cmap('viridis', 12)
sns.set_style("white")

sns.set_palette('colorblind')
from utilities import *
import collections
from collections import OrderedDict
import itertools
import random

import hypernetx as hnx
import scanpy as sc

from scipy.stats import binom_test
from statsmodels.stats.multitest import multipletests

import upsetplot as up
from upsetplot import generate_counts
from upsetplot import from_indicators

import io
from PIL import Image




print('Modules imported \n')

# DEBUGGING
print('Before loading data:')


parser = argparse.ArgumentParser(description='Args for coupling estimation')

parser.add_argument("--dataPath", type=str, help="Path to training data")
parser.add_argument("--PCApath", type=str, help="Path to PCA embedding of training data")
parser.add_argument("--CPDAGgraphPath", type=str, help="Path to CPDAG graph file")
parser.add_argument("--MCMCgraphPath", type=str, help="Path to MCMC graph file")
parser.add_argument("--estimationMode", type=str, help="Estimation mode used for interaction estimation (LOR or MFI)")
parser.add_argument("--pathTo2pts", type=str, help="Path to calculated 2-point interactions")
parser.add_argument("--pathTo3pts", type=str, help="Path to calculated 3-point interactions")
parser.add_argument("--pathTo4pts", type=str, help="Path to calculated 4-point interactions")
parser.add_argument("--pathTo5pts", type=str, help="Path to calculated 5-point interactions")
parser.add_argument("--minStateDeviation", type=float, help="Minimum enrichment factor of a particular state.")
parser.add_argument("--stateDevAlpha", type=float, help="significance threshold to call a state deviating")
parser.add_argument("--plotPairwiseUpsets", type=int, help="Boolean int to decide whether to plot pairwise upset plots.")
parser.add_argument("--sigHOIthreshold", type=float, help="Threshold for significance of HOIs")


args = parser.parse_args()

trainDat = pd.read_csv(args.dataPath)
pcaCoords= pd.read_csv(args.PCApath)
DSname = args.MCMCgraphPath.split('.')[0]
MCMCadjMat = pd.read_csv(args.MCMCgraphPath, index_col=0)
MCMCgraph = ig.Graph.Adjacency(MCMCadjMat.values.tolist())
MCMCgraph.vs['label'] = trainDat.columns.values

CPDAGadjMat = pd.read_csv(args.CPDAGgraphPath, index_col=0)
CPDAGgraph = ig.Graph.Adjacency(CPDAGadjMat.values.tolist())
CPDAGgraph.vs['label'] = trainDat.columns.values

genes = trainDat.columns.values

# Data is loaded as scanpy Anndata objects for easy embedding plots. 
scObj = sc.AnnData(trainDat)
scObj.obsm['X_pca'] = pcaCoords.values


# Load all the 2-point information so that we can select only perfect and signigicant estimations. 
# coups_2pts = np.load(args.pathTo2pts, allow_pickle=True).astype(np.float32)
# coups_2pts_CI_F = np.load(args.pathTo2pts_CI_F, allow_pickle=True).astype(np.float32)
# coups_2pts_undef = np.load(args.pathTo2pts_undef, allow_pickle=True).astype(np.float32)
# coups_2pts_inf = np.load(args.pathTo2pts_inf, allow_pickle=True).astype(np.float32)

# All imperfect estiamtions get set to NaN
# coups_2pts[(coups_2pts_undef>0) | (coups_2pts_inf>0)] = np.nan 
# coups_2pts_CI_F[(coups_2pts_undef>0) | (coups_2pts_inf>0)] = np.nan 


# HHOIs stores all the significant n-point interactions present among interacting pairs, triplets, quadruplets, and pentuplets
HHOIs = {}
alpha = args.sigHOIthreshold
for order, intPath in enumerate([args.pathTo2pts, args.pathTo3pts, args.pathTo4pts, args.pathTo5pts]):
	ints = np.load(intPath, allow_pickle=True)
	if len(ints)>0:
		perfectSigEsts = list(map(lambda x: (((x[[4, 5, 6]]==0).all()) & (x[3]<alpha)), ints))
		HHOIs[f'n{order+2}'] = ints[perfectSigEsts][:, [0, -1]]
	else: HHOIs[f'n{order+2}'] = []

# pairs = np.array(np.where(coups_2pts_CI_F<alpha)).T
# vals = coups_2pts[coups_2pts_CI_F<alpha]
# HHOIs[f'n2'] = np.array([list(x) for x in list(zip(vals, pairs))])

# DEBUGGING
print('After loading data:')


def findsubsets(s, n):
	return list(itertools.combinations(s, n))


def plotUpsetPlot(d, fig, title='', legend=False, save=False, filename=''):
	'''
	Using the upsetplot module to generate the upset plots. 
	The upset library allows you to add extra plots to the axis dictionary, which I use to plot the deviations etc. 
	'''
	genes = d.columns.values
	
	f = lambda x: ''.join(map(str, x))
	
	order = len(d.columns)
	nStates = 2**order

	binCounts = np.bincount(list(map(lambda x: int(x, 2), list(map(f, d.values)))), minlength=nStates)

	# This is an ugly work-around. I need the whole Upset object to be populated with counts, so I just generate a lot of random counts.
	# This results in each possible state occuring at least once, which allows me to reassign each count to its correct value.
	# If this ever fails, it would raise an error when reassigning the counts with upSetObj.loc[:] = binCounts.
	# It would be better to immediately instantiate the Upset object with the correct counts, but I haven't been able to. 
	upSetObj = generate_counts(1, n_samples=100000, n_categories=order)
	upSetObj.loc[:] = binCounts
	upSetObj.index.names = genes

	means = d.mean(axis=0)
	# The expected count is generated from the null hypothesis of independent Bernoulli processes.
	# The probability of success (a 1) for each variable is its mean expression in the conditioned data.
	expected = np.array([np.prod([m if state[i] else 1-m for i, m in enumerate(means)]) for state in upSetObj.index])*len(d)

	upSetObj2 = pd.DataFrame()
	upSetObj2['count'] = upSetObj
	upSetObj2['deviation'] = (upSetObj - expected)/(expected)*100
	upSetObj2['expected'] = expected

	upset = up.UpSet(upSetObj2, subset_size='sum', sum_over='count', intersection_plot_elements=5, min_subset_size=0)

	# Plot the deviations, and instantiate the 'expected' plot, which will be removed later again. 
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

	# Extract the values from the 'expected' plot, and then remove it. 
	ax['intersections'].plot([x.get_offsets()[0][1] for x in ax['extra2'].get_children()[:nStates]], 'k+', zorder=10, label='expected')
	ax['extra2'].remove()

	fig.suptitle(f'{title},   I = ' + 
			str(round(calcInteraction_binTrick_allOrders(d), 3)), fontsize=20)

	if(save):
		plt.savefig(filename)
		plt.close(fig)
	else:
		return

def fromImToArr(img):
	'''
	Takes an ImageFile and converts it to a a numpy array of floats so that imshow() can plot it.  
	'''
	return np.array(img.getdata()).reshape(img.height, img.width, -1)

f=2
kwargs = {'with_node_counts': False, 'with_node_labels':False, 'with_edge_labels':False}
concatInts = lambda x: ''.join(map(str, x))


# I write the plots to a buffer, and store them in dictionaries where the entries are the names of the interacting variables
# Once they are stored like this, I can compose them more easily into a summary figure. 

# dicts to store plots in
plotHypergraph = {}
plotCPDAG = {}
plotPCA = {}
plotupset = {}
plotUpset_cond = {}
plotUpset_uncond = {}
plotMaxDev = {}

if (args.plotPairwiseUpsets):
	ordersToPlot = [2, 3, 4, 5]
else:
	ordersToPlot = [3, 4, 5]

enrichments = {}

for order in ordersToPlot:
	nStates = 2**order

	# Binary representation of each possible state:
	binStates = [np.array(list(format(x, f"0{order}b"))).astype(bool) for x in range(nStates)]
	
	# Stores the log 2-fold change of each d-tuple
	log2Fs = []

	if len(HHOIs[f'n{order}'])>0:
		# loop over the HOIs, use their value (weight w) and gene tuple:
		for w, geneTuple in HHOIs[f'n{order}'][:, [0, -1]]:
			ID = '_'.join(genes[geneTuple])

			# The local CPDAG structure will be plotted, and its layout will be used for the hypergraph of interactions. 
			g = findLocalGraph(geneTuple, CPDAGgraph, order=0)
			layout_c = g.layout('circle')

			# Layout needs to be manipulated a bit so that hypernetx can use it for the hypergraph:
			# tmp = dict(zip([x.replace('.', '-') for x in g.vs['label']], np.array(layout_c)*np.array([1, -1])))
			# NOTE: I have removed the string replacement, as it no longer seems to be necessary and was causing problems. 
			tmp = dict(zip([x for x in g.vs['label']], np.array(layout_c)*np.array([1, -1])))


			# Constructing the hypergraph
			edges = {'maxOrder': genes[geneTuple]}
			weights = [w] # Stores the actual interactions -- used to colour the edges in the hypergraph. 
			for i in range(2, order):
				for subset in findsubsets(geneTuple, i):
					for tup in HHOIs[f'n{i}']:
						if sorted(tup[-1]) == sorted(subset):
							edges[len(edges)+1] =  [genes[g] for g in subset]
							weights.append(tup[0])
							break

			#  ************************ Interaction Hypergraph  ************************ 
			# DEBUGGING
			print('Hypergraph:')
			
			
			H = hnx.Hypergraph(edges)
			cList = ['red' if w<0 else 'green' for w in weights]
			layout_fn = lambda x: tmp
			# print(CPDAGadjMat.columns.values)
			# print(CPDAGgraph.vs['label'])
			# print(list(g.vs()))
			# print(edges)
			# print(list(H.edges()))
			# print(list(H.nodes()))
			fig, ax = plt.subplots(figsize=[10, 10])
			hnx.draw(H, ax=ax, layout = layout_fn,
						 label_alpha=0,
						 node_labels_kwargs={
							'fontsize': 34
						},
						 edges_kwargs={
							 'edgecolors': cList,
							 'linewidths': 3,
							 'dr': 0.05
						 },
						 **kwargs)  
			ly = np.array([x for x in tmp.values()])
			for i, v in enumerate(g.vs['label']):
				plt.text(ly[i][0] - 0.0, ly[i][1], v, fontsize=40)
			plt.show()
			buf = io.BytesIO()
			plt.savefig(buf)
			buf.seek(0)
			plotHypergraph[ID] = fromImToArr(Image.open(buf))
			buf.seek(0)
			buf.truncate(0)
			plt.clf()
			plt.close()

			#  ************************ CPDAG ************************ 
			# DEBUGGING
			print('CPDAG:')
			

			fig, ax = plt.subplots(figsize=[6, 6])
			ig.plot(g, layout=layout_c, edge_arrow_size=20, edge_arrow_width=10,
					vertex_size=80/f, vertex_color='lightgrey',
					target=ax)
			# Have to add labels manually -- igraph does not work properly with matplotlib backend
			for i, v in enumerate(g.vs['label']):
				plt.text(layout_c[i][0] - 0.0, layout_c[i][1]-0.03, v, fontsize=40)
			plt.xlim(-1.2, 1.2)
			plt.ylim(-1.2, 1.2)
			ax.axis('off')
			buf = io.BytesIO()
			plt.savefig(buf)
			buf.seek(0)
			plotCPDAG[ID] = fromImToArr(Image.open(buf))
			buf.seek(0)
			buf.truncate(0)
			plt.clf()
			plt.close()

			#  ************************ PCA embeddings ************************ 

			# DEBUGGING
			print('PCA:')
			

			fig, ax = plt.subplots(1, len(geneTuple), figsize=[20, 4])
			for i, g in enumerate(genes[geneTuple]):
				sc.pl.embedding(scObj,'pca', color=g, 
									size=30, color_map="viridis", add_outline=True, show=False, frameon=False, ax = ax[i])
				tmpFig = plt.gcf()
				tmpFig.axes[-1].remove()
			buf = io.BytesIO()
			plt.savefig(buf)
			buf.seek(0)
			plotPCA[ID] = fromImToArr(Image.open(buf))
			buf.seek(0)
			buf.truncate(0)
			plt.clf()
			plt.close() 

			#  ************************ Upset plots ************************ 
			
			# DEBUGGING
			print('Upset:')
			
			
			unConditionedGenes = trainDat.iloc[:, geneTuple]


			if args.estimationMode=='MFI':
				conditionedGenes = conditionOnMB(geneTuple, MCMCgraph, trainDat, mode='Min')
				fig = plt.figure(figsize=[10, 10])
				buf = io.BytesIO()
				plotUpsetPlot(d = conditionedGenes,fig=fig, legend=False, title = 'Conditioned on MB', filename=buf, save=True)
				buf.seek(0)
				plotUpset_cond[ID] = fromImToArr(Image.open(buf))
				buf.seek(0)
				buf.truncate(0)
				plt.clf()
				plt.close()
			
			else:
				# create empty plot for conditioned genes when not using MFI estimation, since MFI is likely not estimable.
				conditionedGenes = unConditionedGenes
				fig = plt.figure(figsize=[10, 10])
				buf = io.BytesIO()
				plt.savefig(buf)
				buf.seek(0)
				plotUpset_cond[ID] = fromImToArr(Image.open(buf))
				buf.seek(0)
				buf.truncate(0)
				plt.clf()
				plt.close()
				

			fig = plt.figure(figsize=[10, 10])
			buf = io.BytesIO()
			plotUpsetPlot(d = unConditionedGenes,fig=fig, legend=False, title = 'Unconditioned', filename=buf, save=True)
			buf.seek(0)
			plotUpset_uncond[ID] = fromImToArr(Image.open(buf))
			buf.seek(0)
			buf.truncate(0)
			plt.clf()
			plt.close()
			

			#  ************************ Calculate deviations ************************ 

			# Count the occurence of each of the possible binary states in the conditioned data:
			binCounts = np.bincount(list(map(lambda x: int(x, 2), list(map(concatInts, conditionedGenes.values)))), minlength=nStates)

			# The means determine the null hypothesis
			means = conditionedGenes.mean(axis=0)
			pStates = np.array([np.prod([m if state[i] else 1-m for i, m in enumerate(means)]) for state in binStates])
			expected = pStates*len(conditionedGenes)
			
			# Calculate the log 2-fold enrichment, and p-values based on the binomial null
			# Vectorised over the 2^n states
			log2FoldEnrichment = np.log2(binCounts/expected)
			enrichment_pval = np.array([binom_test(binCounts[i], p = pStates[i], n = len(conditionedGenes), alternative='greater') for i in range(len(binCounts))])
			
			log2Fs.append([log2FoldEnrichment, enrichment_pval, geneTuple])
		
		log2Fs  = np.array(log2Fs, dtype=object)
		enrichments[f'n{order}'] = log2Fs
	else: enrichments[f'n{order}'] = []
		
#  ************************ PCA embedding on max deviating state ************************ 

# DEBUGGING
print('Before embedding:')


for order in ordersToPlot:
	for devs, pvals, interactors in enrichments[f'n{order}']:
		ID = '_'.join(genes[interactors])
		
		# maxDevState is the most deviating state
		maxDevState = format(np.argmax(devs), f"0{order}b")

		# The find the PCA coordinates of cells/observations that are in this maxDevState
		maxDevState_embedded = scObj.obsm['X_pca'][(scObj.X[:, interactors] == np.array(list(maxDevState)).astype(float)).all(axis=1)]

		xs = scObj.obsm['X_pca'][:, 0]
		ys = scObj.obsm['X_pca'][:, 1]
		plt.plot(xs, ys, 'o', color = viridis(0), alpha=0.1)
		plt.plot(maxDevState_embedded[:, 0], maxDevState_embedded[:, 1], 'o', color = viridis(0.99), alpha=0.9)
		plt.title(', '.join(genes[interactors]) + ' = ' + ', '.join(maxDevState), fontsize=20)
		plt.xticks([])
		plt.yticks([])

		buf = io.BytesIO()
		plt.savefig(buf)
		buf.seek(0)
		plotMaxDev[ID] = fromImToArr(Image.open(buf))
		buf.seek(0)
		buf.truncate(0)
		plt.clf()
		plt.close()
		# plt.savefig(f'{ID}_Expression_maxDevState.png')
		plt.close('all') 


#  ************************ Plot summary figures ************************ 
# DEBUGGING
print('Before summary figs:')


sns.set_style("white")

for order in ordersToPlot:
	# DEBUGGING
	print(f'Order: {order}')
	

	if len(HHOIs[f'n{order}'])>0:
		for w, geneTuple in HHOIs[f'n{order}'][:, [0, -1]]:
			ID = '_'.join(genes[geneTuple])	  
			fig = plt.figure(figsize=(15, 10))

			# Defining custom axes to position individial plots. 
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

			axCPDAG.imshow(plotCPDAG[ID])
			axHOI.imshow(plotHypergraph[ID][100:-100, 100:-50])
			axEXP.imshow(plotPCA[ID][:, :, :])
			axEXP_maxDev.imshow(plotMaxDev[ID][:, :, :])

			axUPS1.imshow(plotUpset_cond[ID][:, :])
			axUPS2.imshow(plotUpset_uncond[ID][:, :])

			plt.savefig(f'{ID}_summary.png')
			plt.close(fig) 

# DEBUGGING
print('After summary figs:')


# Construct a dataframe with the deviations of each positively enriched state:
devDict = []
for order in [3, 4, 5]:
	for devs, pvals, interactors in enrichments[f'n{order}']:

		ID = '_'.join(genes[interactors])
		for devStateInd in np.where((devs>=0))[0]:
			# Format the binary state
			maxDevState = format(devStateInd, f"0{order}b")
			devDict.append([ID, maxDevState, devs[devStateInd], pvals[devStateInd]])

# DEBUGGING
print('After dev dicts:')


if len(devDict)>0:
	deviators = pd.DataFrame(devDict, columns=['genes', 'state', 'enrichment', 'pval'])
	deviators = deviators.sort_values(by='pval', ascending=True)
	# Add the Benjamini-Yekutieli corrected p-values:
	deviators['pval_corrected'] = multipletests(deviators['pval'], method='fdr_by')[1]

	# Binreps is the binary represenation of the interactions: binReps[i] is 1 if cell i is in the deviating state, 0 otherwise.  
	binReps = np.array(deviators.apply(lambda x: (trainDat[x['genes'].rsplit('_')]==[int(g) for g in list(str(x['state']))]).all(axis=1), axis=1))*1

	# Add the cell IDs to each of the d-tuples:
	# (has to be made into a list to avoid the ellipsis in the output)
	deviators['cellIDs'] = [np.where(binRep)[0].astype(int).tolist() for binRep in binReps]

	# Filter out the d-tuples that are significantly enriched:
	strongDeviators = deviators.query(f'pval_corrected <= {args.stateDevAlpha} and enrichment >= {args.minStateDeviation}').copy()

	deviators.to_csv(f'all_DTuples.csv')
	strongDeviators.to_csv(f'top_DTuples.csv')

	pd.DataFrame(binReps).to_csv(f'DTuples_binaryReps.csv')

else: 
	pd.DataFrame(data=[]).to_csv(f'all_DTuples.csv')
	pd.DataFrame(data=[]).to_csv(f'top_DTuples.csv')

