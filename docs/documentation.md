
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

Nextflow needs to be installed on the machine you want to run the pipeline on, which can be as simple as running:

```bash
curl -s https://get.nextflow.io | bash
```
For more information on setting up Nextflow, see [their official guide](https://www.nextflow.io/docs/latest/getstarted.html).

Since runs typically take a few days (for, say, 1000 genes and 20,000 cells), it is recommended to run the pipeline on a cluster, but local execution is supported as well. 

You must have access to either Docker or Singularity. Most HPC clusters use Singularity for security reasons, but both will automatically pull the right container from DockerHub. 

To schedule job submission on a cluster, this repo comes with a profile for Sun Grid Engine (SGE) compatible platforms, which is used on the Edinburgh University compute cluster Eddie. This can be used as a template for other cluster schedulers, but is not guaranteed to work, so contact your local Nextflow users for help. 


## Input files

* `rawDataPath`: A count matrix in csv format, where the rows are cells, and the columns genes. The first row should contain the gene names (does not matter in which format). 
* `userGenes` (optional): A list of genes that should be included in the final analysis, irrespective of their variability. WARNING: If more than `nGenes` genes are provided, `nGenes` will be overwritten and set to the number of genes in the list. The gene names should be spelled exactly as they are in the count matrix header. 


## Output files
* `output/`
    * `trainingData_{$nCells}Cells_{$nGenes}Genes.csv`: The nCells x nGenes count matrices after filtering.
    * `unbinarised_cell_data.h5ad`: The unbinarised expression of all genes in the selected cells. 
    * `{$graphType}graph_{$nCells}Cells_{$nGenes}Genes.csv`: The nGenes x nGenes adjacency matrices for the graphs used in the estimation steps.
    * `trainingData_{$nCells}Cells_{$nGenes}Genes_{$embedding}coords.csv`: Embedding coordinates of the selected cells. 
    * `.png`: Figures with basic QC metrics (only if running in `expression` mode).
* `coupling_output/`
    * `interactions_withinMB_{$order}pts_[...].npy`: 
    The interaction estimates and statistics at order 2-7. Each interaction is a tuple of length 8. The values are, in order.
    1. The point estimate
    2. Lower bound of 95% CI
    3. Upper bound of 95% CI
    4. F value
    5. Proportion of undefined bootstrap resamples, 
    6. Proportion of divergent bootstrap resamples
    7. Bound indicator (0 by default)
    8. Tuple of gene indices
    * `interactions_random_{$order}pts_[...].npy`: The same as above, but for `nRandomHOIs` random interactions.
* `dtuples_output/`
    * `{$genes}_summary.png`: Figures that summarise the significant 2-, 3-, 4-, and 5-point interactions and contain upset plots for the d-tuples.
    * `all_DTuples.csv`: A list of the positively enriched d-tuples. 
    * `top_DTuples.csv`: A list of the positively and significantly enriched d-tuples. 
    * `DTuples_binaryReps.csv`: Binary representations of all positively enriched d-tuples.
* `states_output/`
    * `modularity_scores.csv`: Cutoffs and modularity scores of the hierarchical clustering of the binary representation.
    * `top_DTuples_withStatesFromCut.csv`: A list of the positively and significantly enriched d-tuples, with the associated cluster resulting from the cut that maximised the modularity score.
    * `bootstrapStats.csv`: The bootstrap statistics of the dendrogram branches.
    * `dendrogram_all_dTuples.png`: A figure of the full dendrogram of all d-tuples.
* `reports/`
    * Some reports from Nextflow on resource usage etc.
* `work/`
    * The working directory of Nextflow, useful for debugging, and resuming previous runs. 


## Parameters
The pipeline takes in a number of parameters, which can be set by the user in the params.json file. An example file called `params_example.json` is included in this directory. 

The `dataType` parameter determines preprocessing of input data, and can take two values:
* `dataType = "expression"`: Standard scRNA-seq QC is performed (further specified by `fracMito`, `minGenes`, `minCells`), the count values are normalised and log-transformed, a total of `(nGenes - |userGenes|)` of the most highly variable genes are included beyond the optionally provided `userGenes`. 
* `dataType = "agnostic"` (default): No QC is performed, and only the `userGenes` and first `(nGenes - |userGenes|)` genes are included. 

In both cases, if more than `nCells` are provided in the count matrix, a random selection of `nCells` is included in the analysis.

Here follows an overview of the other parameters that can affect the results:
| Parameter | Default | Description | Required? | 
| :----- | :----- | :----- | :-- |
| dataType | 'agnostic' | See explanation above this table | Yes |
| rawDataPath | ' ' | absolute path to count matrix (`.csv`) | Yes |
| nGenes | ' ' | Number of genes to keep | Yes |
| nCells | ' ' | Number of cells to keep | Yes |
| userGenes | ' ' | absolute path to list of required genes (`.csv`) | No |
| fracMito | 1 | cells with more than `fracMito` mitochondrial reads get dismissed | Only when `dataType=='expression'` |
| minGenes | 0 | cells with fewer than `minGenes` expressed get dismissed | Only when `dataType=='expression'` |
| minCells | 0 | genes expressed in fewer than `minCells` get dismissed | Only when `dataType=='expression'` |
| PCalpha | 0.05 | Significance threshold to use for the PC-algorithm | Yes |
| asympBool | 0 | Boolean that determines if the variance is estimated from bootstrap resamples (0) or an asymptotic approximation (1) | Yes |
| bsResamps | 1000 | Number of bootstrap resamples to use when calculating confidence intervals on interactions | Only when `asympBool==0` |
| estimationMode | 'MFI' | Setting this to `MFI` (default) yields estimates of model-free interactions by conditioning on the Markov-blanket, setting it to `LOR` yields unconditioned log-odds ratios. | Yes |
| nRandomHOIs | 1000 | How many random 6 & 7-point interactions to calculate | Yes |
| plotPairwiseUpsets | 0 | Boolean to determine if pairwise upset plots should be generated | Yes |
| sigHOIthreshold | 0.05 | Significance threshold on F-value to decide which HOIs get summarised and used for states **(will soon be replaced by CI)** | Yes |




And these affect the resources accessible to each of the processes, but shouldn't affect the results:
| Parameter | Default | Description | Required? | 
| :----- | :----- | :----- | :-- |
| executor | ' ' | How to execute the pipeline ('sge' for cluster usage, 'local' for local usage) | Yes |
| maxQueueSize | ' ' | How many jobs are allowed to be scheduled at the same time | Only when `executor=='sge'` |
| [cores, mem, time]\_makeData | ' ' | How many [cores, memory, hours] are available when preparing the data | Only when `executor=='sge'` |
| [cores, mem, time]\_PC | ' ' | How many [cores, memory, hours] are available to the PC-algorithm | Only when `executor=='sge'` |
| [cores, mem, time]\_MCMC | ' ' | How many [cores, memory, hours] are available to the MCMC scheme | Only when `executor=='sge'` |
| [cores, mem, time]\_1pt | ' ' | How many [cores, memory, hours] are available to the 1-point estimation | Only when `executor=='sge'` |
| [cores, mem, time]\_HOIs_MB | ' ' | How many [cores, memory, hours] are available to the 3-5-point estimations within Markov blankets | Only when `executor=='sge'` |
| [cores, mem, time]\_HOIs_6n7 | ' ' | How many [cores, memory, hours] are available to the 6- and 7-point estimations | Only when `executor=='sge'` |
| [cores, mem, time]\_HOIs_plots | ' ' | How many [cores, memory, hours] are available for plotting the HOI summaries and d-tuples | Only when `executor=='sge'` |


## Usage

First, you need to pull the latest version of the pipeline from Github. This is done by:


```bash
nextflow pull AJnsm/Stator -r branch
```
where branch is set to either `main` (most stable) or `develop`. 

To run `Stator` locally, you need to have Docker installed. Then, from the directory where you want the output directories to be generated, the pipeline can be run with the command:

```bash
nextflow run AJnsm/Stator -r main -profile docker -params-file params.json
```
where `-r main ` specifies the branch/revision, `-profile docker` selects the Docker profile that links to Dockerhub, and `-params-file params.json` specifies the JSON file with the necessary parameters. An example JSON file is provided in this repository. 

On a cluster, you need to load either Docker or Singularity. On the Edinburgh University compute cluster (Eddie) this is done with 
```bash
module load singularity
```
Then, from the directory where you want the output directories to be generated, the pipeline can be run with the command:
```bash
nextflow run AJnsm/Stator -r main -profile eddie_singularity -params-file params.json
```
Clusters can have different scheduling software, so you might need to create a custom profile. The one included here--[eddie_singularity.config](../configs/eddie_singularity.config)--is for the Sun Grid Engine (SGE) scheduler installed on the Edinburgh compute cluster. It's probably easiest to contact your local Nextflow users for help if you need to create a custom profile.

**NOTE: On a cluster, you need to make sure you are on a node that allows automated job submission. On Eddie these are known as the wild-west nodes.**


## Unit tests
The `scripts` directory contains a file `scripts/unit_tests.py`. This is a collection of `pytest` unit tests to assert that various parts of the estimation procedure function properly. As it is not necessary to run the pipeline, `pytest` is not part of the containerised `Python` environment. To run the unit test, install `pytest` using e.g. `pip` or `conda`, and run the following command from the `scripts` directory:
```bash
pytest unit_tests.py
```
This should run over 20 unit tests. Compiling the `numba` estimation function can take a few minutes. Running all tests takes around 5 minutes on a 2018 Macbook Pro (2.9 GHz 6-Core Intel Core i9, 16GB memory). 

## Integration tests
**These are tests I use during development, but require some dependencies.**
The `integration tests` directory contains a bash script `runTests.py` that generates 100k states from a nearest neighbour Ising model, using the `magneto` simulation software, and estimates the 1st and 2nd order interactions. This is already done, and all generated files are located in the same directory. Running 

```bash
sh runTests.py
pytest -rP integration_tests.py
```
verifies that the inferred interactions are (approximately) correct. This done once using the full pipeline---PC algorithm, MCMC optimisation, and MFI estimation---and once using the true underlying conditional dependency graph---the nearest-neighbour strucure. 

**Results:**
- Using the MCMC graph, the 67% of the estimable 1-point interactions were significant at alpha=0.1 (in theory all should be significant), and the 2-point interactions identified the nearest-neighbour structure with an F1-score of  0.53, vs. 0.37 for randomly shuffled 2-point interactions. 
- Using the known nearest-neighbour structure as the causal graph, 100% of the estimable 1-point interactions were significant at alpha=0.1, and the 2-point interactions identified the nearest-neighbour structure with an F1-score of 0.92, vs. 0.54 for randomly shuffled 2-point interactions. 

