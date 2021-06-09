# CoComix
**Co**nditional **Co**uplings in (transcript-)**omics** data. 


## Introduction
This pipeline takes in single cell RNA-seq count matrices, and estimates gene-interactions at first, second, and third order. 
The pipeline can be run directly from the command line. It pulls all code from Github, and the container from Dockerhub. It can run on your local machine, or on a cluster with the Sun Grid Engine scheduling system (like Eddie). 


## Design/DAG

![Pipeline flow](diagram.pdf)

## Requirements

Nextflow needs to be installed on the machine you want to run the pipeline on, which is as simple as running:

```bash
curl -s https://get.nextflow.io | bash
```

Secondly, you must have access to either Docker or Singularity. Most clusters (like Eddie) use Singularity for security reasons, but both will automatically pull the right container from DockerHub. (Conda environments are supported, but not recommended or guaranteed to work.)

## Input files

* A count matrix in csv format, where the rows are genes, and the columns cells. The index column should contain the gene names (does not matter in which format), and the first row/header should contain cell identifiers (either barcodes or an index). 
* A list of integer cluster annotations per cell in csv format: This should be in the same order as the cells in the count matrix. 
* (optional) A list of genes that should be included in the final analysis, irrespective of their variability.  
* (optional) A list of booleans on whether or not a cell should be included (e.g. on the basis of being a suspected doublet), in the same order as the cells in the count matrix. 


## Output files
* .png plot with the mean-variance relationship among the genes (in `plots/`)
* .csv files with the nGenes x nCells count matrices after filtering (in `output/`). 
* .csv files with the nGenes x nGenes adjacency matrices for the graphs used in the estimation steps (in `output/`). 
* .npy files with the interactions at each order for each of the data sets (in `coupling_output/`). 
* .npy files with the bounds on the 95% confidence interval of the error on the interactions (in `coupling_output/`). 
* .npy files with the proportion of bootstrap resamples that fall on the other side of zero (in `coupling_output/`). 
* some reports from Nextflow on resource usage etc. (in `reports/`).


## Parameters
The pipeline takes in a number of parameters, which can be set by the user in the params.json file. 

These affect the calculation and the results:
| Parameter | Default | Description | Required? | 
| :----- | :----- | :----- | :-- |
| rawDataPath | ' ' | absolute path to count matrix .csv | Yes |
| clusterFile | ' ' | absolute path to cluster annotation .csv | No |
| userGenes | ' ' | absolute path to list of required genes .csv | No |
| nGenes | 5 | Number of genes to keep | Yes |
| nCells | 20 | Number of cells to keep | Yes |
| clusterArray | 0 | List of which clusters to keep | No |
| pcAlpha | 0.05 | Significance threshold to use for the PC-algorithm | Yes |
| bsResamps | 1000 | Number of bootstrap resamples to use when calculating confidence intervals on interactions | Yes |
| twoReplicates | false | Boolean to decide on constructing two replicates from the data | Yes |
| nCells | 20 | Number of cells to keep | Yes |


And these only affect the resources accessible to each of the processes:
| Parameter | Default | Description | Required? | 
| :----- | :----- | :----- | :-- |
| executor | 'local' | How to execute the pipeline ('sge' for cluster usage) | Yes |
| maxQueueSize | 25 | How many processes are allowed to be scheduled at the same time | Only with 'sge' |
| [cores, mem, time]\_makeData | [1, '4G', '1h'] | How many [cores, memory, hours] are available when preparing the data | Only with 'sge' |
| [cores, mem, time]\_PC | [1, '4G', '1h'] | How many [cores, memory, hours] are available to the PC-algorithm | Only with 'sge' |
| [cores, mem, time]\_makeData | [1, '4G', '1h'] | How many [cores, memory, hours] are available to the MCMC scheme | Only with 'sge' |
| [cores, mem, time]\_makeData | [1, '4G', '1h'] | How many [cores, memory, hours] are available to the 1-point estimation | Only with 'sge' |
| [cores, mem, time]\_makeData | [1, '4G', '1h'] | How many [cores, memory, hours] are available to the 2-point estimation | Only with 'sge' |
| [cores, mem, time]\_makeData | [1, '4G', '1h'] | How many [cores, memory, hours] are available to the 3-point estimation | Only with 'sge' |


## Usage

First, you need to pull the latest version of the pipeline from Github. This is done by:


```bash
nextflow pull AJnsm/NF_TL_pipeline
```
Currently, the branch to use is called `develop`. On a cluster, you need to load Singularity. On Eddie this is done with 


```bash
module load singularity
```

 Then, from the directory where you want the output directories to be generated, the pipeline can be run with the command:

```bash
nextflow run AJnsm/NF_TL_pipeline -r develop -profile eddie_singularity -params-file params.json
```

Where ```-r develop ``` specifies the branch/revision, ```-profile eddie_singularity``` selects the right profile for the Eddie environment, and ```-params-file params.json``` specifies a JSON file with the necessary parameters. An example JSON file is provided in this repository.


**NOTE: On a cluster, you need to make sure you are on a node that allows automated job submission. On Eddie these are known as the wild-west nodes.**

## To do

* make sure running locally is working, has a good profile, and is explained
* show a full run on some example data. 













