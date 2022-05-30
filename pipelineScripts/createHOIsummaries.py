import requests 
import os
import pandas as pd
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
from matplotlib import cm
%config InlineBackend.figure_format = 'retina'
import math
import seaborn as sns
sns.set(style='darkgrid')
sns.set_palette('colorblind')
from utilities import *
import collections
from collections import OrderedDict
import itertools
import random

from gprofiler import GProfiler
import hypernetx as hnx
import scanpy as sc

import upsetplot as up
from upsetplot import generate_counts
from upsetplot import from_indicators


if PrintBool: print('Modules imported \n')

parser = argparse.ArgumentParser(description='Args for coupling estimation')

parser.add_argument("--dataPath", type=str, nargs=1, help="Path to training data")
parser.add_argument("--graphPath", type=str, nargs=1, help="Path to graph file")
parser.add_argument("--pathToNpts", type=str, nargs='?', help="Path to calculated 5-point interactions")

args = parser.parse_args()

pathToNpts = args.pathToNpts

trainDat = pd.read_csv(dataPath)
DSname = graphPath.split('.')[0]
adjMat = pd.read_csv(graphPath, index_col=0)
graph = ig.Graph.Adjacency(adjMat.values.tolist())


for order in [3, 4, 5]:
    ints = np.load(f'{ds}_{ct}/higherOrderOutput/interactions_MB_{order}pts_expectations.npy', allow_pickle=True)
    perfectSigEsts = list(map(lambda x: (((x[[4, 5, 6]]==0).all()) & (x[3]<alpha)), ints))
    HHOIs[f'{ds}_{ct}_n{order}'] = ints[perfectSigEsts][:, [0, -1]]
            
        pairs = np.array(np.where(coups_perfectEstimations_MCMC[f'{ds}_{ct}_n2_CI_F']<alpha)).T
        vals = coups_perfectEstimations_MCMC[f'{ds}_{ct}_n2'][coups_perfectEstimations_MCMC[f'{ds}_{ct}_n2_CI_F']<alpha]
        HHOIs[f'{ds}_{ct}_n2'] = np.array([list(x) for x in list(zip(vals, pairs))])