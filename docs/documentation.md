
# Table of contents
* [Requirements](#requirements)
* [Input files](#input-files)
* [Output files](#output-files)
* [Parameters](#parameters)
* [Usage](#usage)
* [Unit tests](#unit-tests)
* [Integration tests](#integration-tests)



## Pipeline design

![Pipeline flow](/docs/nextflowDAG.png)

## Requirements

Nextflow needs to be installed on the machine you want to run the pipeline on, which is as simple as running:

```bash
curl -s https://get.nextflow.io | bash
```

You must have access to either Docker or Singularity. Most HPC clusters use Singularity for security reasons, but both will automatically pull the right container from DockerHub. (Conda environments will be supported, but are yet not recommended or guaranteed to work.)

## Input files

* `rawDataPath`: A count matrix in csv format, where the rows are cells, and the columns genes. The first row should contain the gene names (does not matter in which format). 
* `clusterFile` (optional): A list of integer cluster annotations per cell in csv format: This should be in the same order as the cells in the count matrix. 
* `userGenes` (optional): A list of genes that should be included in the final analysis, irrespective of their variability.  
* `genesToOne` (optional): A list of genes that should be conditioned to a 1 instead of a 0.  
* `doubletFile` (optional): A list of booleans that indicate whether a cell should be included (e.g. on the basis of being a suspected doublet), in the same order as the cells in the count matrix. 


## Output files
* `output/`
    * `trainingData_CL{$cluster}_{$nCells}Cells_{$nGenes}Genes.csv`: The nCells x nGenes count matrices after filtering.
    * `unbinarised_cell_data.h5ad`: The unbinarised expression of all genes in the selected cells. 
    * `{$graphType}graph_CL{$cluster}_{$nCells}Cells_{$nGenes}Genes.csv`: The nGenes x nGenes adjacency matrices for the graphs used in the estimation steps.
    * `trainingData_CL01_02000Cells_0020Genes_{$embedding}coords.csv`: Embedding coordinates of the selected cells. 
    * `.png`: Figures with basic QC metrics.
* `coupling_output/`
    * `interactions_order{$order}[...]_coup.npy`: The interaction point estimates at order 1 & 2.
    * `interactions_order{$order}[...]_CI_U(L)B.npy`: The upper (lower) bound of the 95% confidence interval.
    * `interactions_order{$order}[...]_CI_F.npy`: The fraction of resamples with a different sign than the point estimate (F-value).
    * `interactions_order{$order}[...]_inf(undef).npy`: The fraction of resamples that were infinite (undefined).
    * `interactions_withinMB_{$order}pts[...].npy`: The interaction estimates and statistics at order 3-7 (within Markov blanket and `nRandomHOIs` random ones).
* `HOIsummaries/`
    * `{$genes}_summary.png`: Figures that summarise the significant 3-, 4-, and 5-point interactions. 
    * `all_DTuples.csv`: A list of the positively enriched d-tuples. 
    * `top_DTuples.csv`: A list of the positively and significantly enriched d-tuples. 
    * `DTuples_binaryReps.csv`: Binary representations of all positively enriched d-tuples.
    * `distinctDeviatingStates_dendrogram.png`: The characteristic states embedded in PCA coordinates, in a dendrogram. 
* `states_output/`
    * `modularity_scores.csv`: Cutoffs and modularity scores of the hierarchical clustering of the binary representation.
    * `top_DTuples_withStatesFromCut.csv`: A list of the positively and significantly enriched d-tuples, with the associated cluster resulting from the cut.
    * `bootstrapStats.csv`: The bootstrap statistics of the dendrogram branches.
    * `statesWithBootstrapStats.csv`: The cell states with their associated bootstrap statistics. 
    * `dendrogram_all_dTuples.png`: A figure of the full dendrogram of all d-tuples.
    * `dendrogram_all_dTuples_cut.png`: A figure of the dendrogram after the cut.
* `reports/`
    * Some reports from Nextflow on resource usage etc.
* `work/`
    * The working directory of Nextflow, useful for debugging. 


## Parameters
The pipeline takes in a number of parameters, which can be set by the user in the params.json file. 

These affect the calculation and the results:
| Parameter | Default | Description | Required? | 
| :----- | :----- | :----- | :-- |
| dataType | 'agnostic' | type of input data: determines preprocessing | Yes |
| rawDataPath | ' ' | absolute path to count matrix (`.csv`) | Yes |
| nGenes | ' ' | Number of genes to keep | Yes |
| nCells | ' ' | Number of cells to keep | Yes |
| doubletFile | ' ' | absolute path to doublet annotation (`.csv`) | No |
| userGenes | ' ' | absolute path to list of required genes (`.csv`) | No |
| fracMito | 1 | cells with more than `fracMito` mitochondrial reads get dismissed | Only when `dataType=='expression'` |
| minGenes | 0 | cells with fewer than `minGenes` expressed get dismissed | Only when `dataType=='expression'` |
| minCells | 0 | genes expressed in fewer than `minCells` get dismissed | Only when `dataType=='expression'` |
| PCalpha | 0.05 | Significance threshold to use for the PC-algorithm | Yes |
| asympBool | 0 | Boolean that determines if the variance is estimated from bootstrap resamples (0) or an asymptotic approximation (1) | Yes |
| boundBool | 0 |  Boolean that determines if inestimable interactions be bounded | Yes |
| bsResamps | 1000 | Number of bootstrap resamples to use when calculating confidence intervals on interactions | Only when `asympBool==0` |
| calcAll2pts | 0 |  Boolean that determines if all 2-points should be calculated (1) or only the Markov-connected ones (0, default) | Yes |
| estimationMode | 'MFI' | Setting this to `MFI` (default) yields estimates of model-free interactions by conditioning on the Markov-blanket, setting it to `LOR` yields unconditioned log-odds ratios. | Yes |
| nRandomHOIs | 1000 | How many random 6 & 7-point interactions to calculate | Yes |
| plotPairwiseUpsets | 0 | Boolean to determine if pairwise upset plots should be generated | Yes |
| sigHOIthreshold | 0.05 | Significance threshold on F-value to decide which HOIs get summarised and used for states | Yes |
| minStateDeviation | 5 | Min. enrichment factor for characteristic states | Yes |
| stateDevAlpha | 0.05 | Min. enrichment significance for characteristic states | Yes |
| dendCutoff | 0.88 | Dice distance at which the dendrogram gets cut | Yes |
| bsResamps_HC | 100 | Number of bootstrap resampled state dendrograms to generate | Yes |
| auThreshold | 0.95 | AU threshold that determines if a dendrogram branch is significant | Yes |



And these affect the resources accessible to each of the processes, but shouldn't affect the results:
| Parameter | Default | Description | Required? | 
| :----- | :----- | :----- | :-- |
| executor | ' ' | How to execute the pipeline ('sge' for cluster usage) | Yes |
| maxQueueSize | ' ' | How many jobs are allowed to be scheduled at the same time | Only when `executor=='sge'` |
| [cores, mem, time]\_makeData | ' ' | How many [cores, memory, hours] are available when preparing the data | Only when `executor=='sge'` |
| [cores, mem, time]\_PC | ' ' | How many [cores, memory, hours] are available to the PC-algorithm | Only when `executor=='sge'` |
| [cores, mem, time]\_MCMC | ' ' | How many [cores, memory, hours] are available to the MCMC scheme | Only when `executor=='sge'` |
| [cores, mem, time]\_1pt | ' ' | How many [cores, memory, hours] are available to the 1-point estimation | Only when `executor=='sge'` |
| [cores, mem, time]\_2pt | ' ' | How many [cores, memory, hours] are available to the 2-point estimation | Only when `executor=='sge'` |
| [cores, mem, time]\_3pt | ' ' | How many [cores, memory, hours] are available to the 3-point estimation | Only when `executor=='sge'` |
| [cores, mem, time]\_HOIs_MB | ' ' | How many [cores, memory, hours] are available to the 3-5-point estimations within Markov blankets | Only when `executor=='sge'` |
| [cores, mem, time]\_HOIs_6n7 | ' ' | How many [cores, memory, hours] are available to the 6- and 7-point estimations | Only when `executor=='sge'` |
| [cores, mem, time]\_HOIs_plots | ' ' | How many [cores, memory, hours] are available for plotting the HOI summaries and characteristic states | Only when `executor=='sge'` |


## Usage

First, you need to pull the latest version of the pipeline from Github. This is done by:


```bash
nextflow pull AJnsm/Stator -r branch
```
where branch is set to either `main` (most stable) or `develop`. On a cluster, you need to load Singularity. On Eddie this is done with 


```bash
module load singularity
```

 Then, from the directory where you want the output directories to be generated, the pipeline can be run with the command:

```bash
nextflow run AJnsm/Stator -r main -profile eddie_singularity -params-file params.json
```

Where ```-r main ``` specifies the branch/revision, ```-profile eddie_singularity``` selects the right profile for the Eddie environment, and ```-params-file params.json``` specifies a JSON file with the necessary parameters. An example JSON file is provided in this repository.


**NOTE: On a cluster, you need to make sure you are on a node that allows automated job submission. On Eddie these are known as the wild-west nodes.**

## Unit tests
The `pipelineScripts` directory contains a file `pipelineScripts/unit_tests.py`. This is a collection of `pytest` unit tests to assert that various parts of the estimation procedure function properly. As it is not necessary to run the pipeline, `pytest` is not part of the containerised `Python` environment. To run the unit test, install `pytest` using e.g. `pip` or `conda`, and run the following command from the `pipelineScripts` directory:
```bash
pytest unit_tests.py
```
This should run over 20 unit tests. Compiling the `numba` estimation function can take a few minutes. Running all tests takes around 5 minutes on a 2018 Macbook Pro (2.9 GHz 6-Core Intel Core i9, 16GB memory). 

## Integration tests
**These are tests I used in the development, but require some dependencies.**
The `integration tests` directory contains a bash script `runTests.py` that generates 100k states from a nearest neighbour Ising model, using the `magneto` simulation software, and estimates the 1st and 2nd order interactions. This is already done, and all generated files are located in the same directory. Running 

```bash
sh runTests.py
pytest -rP integration_tests.py
```
verifies that the inferred interactions are approximately correct. This done once using the full pipeline---PC algorithm, MCMC optimisation, and MFI estimation---and once using the true underlying conditional dependency graph---the nearest-neighbour strucure. 

**Results:**
- Using the MCMC graph, the 67% of the estimable 1-point interactions were significant at alpha=0.1 (in theory all should be significant), and the 2-point interactions identified the nearest-neighbour structure with an F1-score of  0.53, vs. 0.37 for randomly shuffled 2-point interactions. 
- Using the known nearest-neighbour structure as the causal graph, the 100% of the estimable 1-point interactions were significant at alpha=0.1, and the 2-point interactions identified the nearest-neighbour structure with an F1-score of 0.92, vs. 0.54 for randomly shuffled 2-point interactions. 

