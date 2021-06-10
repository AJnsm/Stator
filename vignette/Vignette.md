# Minimal working example

We will run the pipeline on some public data from 10X. We need two files:
* `10X_sampleData.csv`: the count matrices where each row is a cell, and each column a gene
* `10X_sampleClusters.csv`: a cluster annotation for each cell, in the same order as the cells in `10X_sampleData.csv`

Let's have a look at the first 10 rows and columns of `10X_sampleData.csv`:
```bash
$ head -n 10 10X_sampleData.csv | cut -d ',' -f1-10
> ,Xkr4,Gm1992,Gm37381,Rp1,Rp1-1,Sox17,Gm37323,Mrpl15,Lypla1
> AACCATGAGCTGAAAT-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
> CACCAGGGTAGCACGA-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
> GATCGATCAGGACGTA-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0
> TAGTTGGCATGCGCAC-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
> TTAGTTCCACCCAGTG-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
> AGAGCGACAGCTCGAC-2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
> CACCAGGTCGCCATAA-2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
> CGTAGGCAGGATGTAT-2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
> CTACGTCAGTTAAGTG-2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0

```

Looking at `10X_sampleClusters.csv` shows that both files are indeed in the same order:

```bash
$ head -n 10 10X_sampleData.csv | cut -d ',' -f1-10
> ,cluster
> AACCATGAGCTGAAAT-1,10
> CACCAGGGTAGCACGA-1,10
> GATCGATCAGGACGTA-1,10
> TAGTTGGCATGCGCAC-1,10
> TTAGTTCCACCCAGTG-1,7
> AGAGCGACAGCTCGAC-2,10
> CACCAGGTCGCCATAA-2,10
> CGTAGGCAGGATGTAT-2,10
> CTACGTCAGTTAAGTG-2,10

```

We might be interested in some particular genes that should definitely be included. We've listed some genes in `10X_userGenes.csv`:

```bash
$ cat 10X_userGenes.csv
> Jun,Atf3,Fos,Igf2,Btg2,Vim,Apoe
```

These are all the input files we need. The only other thing we need to do is to configure the pipeline parameters in a JSON file. There is an example file in this repo with the following settings:

```bash
$ cat 10X_sampleParams.json
> {
>       "dataType"    : "agnostic",
>       "rawDataPath" : "/absolute/path/to/this/vignette/10X_sampleData.csv",
>       "clusterFile" : "/absolute/path/to/this/vignette/10X_sampleClusters.csv",
>       "userGenes"   : "/absolute/path/to/this/vignette/10X_userGenes.csv",
>       "nGenes"      : 10,
>       "nCells"      : 20,
>       "clusterArray": [7, 10],
>       "PCalpha"     : 0.05,
>       "bsResamps"   : 500,
>       "twoReplicates": false,
> 
>       "executor"    : "sge",
>       "maxQueueSize"  : 25,
> 
>       "cores_makeData": 1,
>       "cores_PC"      : 2,
>       "cores_MCMC"    : 2,
>       "cores_1pt"     : 2,
>       "cores_2pt"     : 2,
>       "cores_3pt"     : 2,
> 
>       "mem_makeData"  : "4G",
>       "mem_PC"        : "4G",
>       "mem_MCMC"      : "4G",
>       "mem_1pt"       : "4G",
>       "mem_2pt"       : "4G",
>       "mem_3pt"       : "4G",
> 
>       "time_makeData" : "1h",
>       "time_PC"       : "1h",
>       "time_MCMC"     : "1h",
>       "time_1pt"      : "1h",
>       "time_2pt"      : "1h",
>       "time_3pt"      : "1h"
>  }
```
Some explanation:
* `"dataType"    : "agnostic"` just means that you're putting in .csv's, rather than some special data-types. I will probably remove this option soon. 
* Then we sepcify the absolute paths to the files we've defined above. 
* We will keep `nGenes = 10` genes (that is including the ones defined in `10X_userGenes.csv`), and `nCells = 20` cells. In reality, you would probably be interested in way more cells and genes, but to keep things light, I've put in some very small numbers here.
* `"clusterArray": [7, 10]`: As seen in `10X_sampleClusters.csv`, the cells all belong to clusters 7 and 10 (what these clusters are is not relevant, they come from 10X's clustering on a much larger data set), so let's calculate interactions in both of these. That means we'll calculate interactions on 20 cells from cluster 7, and 20 cells from cluster 10. 
* `"PCalpha"     : 0.05`: When running the PC-algorithm, we will use this as the significance threshold to conclude a conditional dependency. 
* `"twoReplicates": false`: We don't want to create two replicates per cluster. This is nice if you want to validate your results in disjoint sets of cells from the same cluster/cell type. Make sure that there are at least `2 x nCells` cells in each cluster when this is set to true!
* `"executor"    : "sge"`: We want to run this on a cluster with the sun grid engine scheduler.
* `"maxQueueSize"  : 25`: And we don't want to be impolite and submit more than 25 jobs at the same time. 
* `"cores_*` specifies the number of cores/cpus each of the processes should request. The data preparation process can not parallelised so one is fine, but the others will be significantly faster on more cores. I've used up to 32 on Eddie.
* The `mem_*` and `time_*` specify how much memory and time each process gets. If these are too low, the whole pipeline terminates when the scheduler starts to complain. What I do: set them relatively high the first time, and then look at the reports from nextflow to see the peak usage of each of the resources. 




We are now ready to run the pipeline. Make sure you are on a node of the scheduler that allows for automated job submission (e.g. the wild-west nodes on Eddie), and that you are in the right directory:

```bash
[you@wild-west-node]$ ls
> 10X_sampleClusters.csv  10X_sampleData.csv  10X_userGenes.csv  params_test.json  

```

Then pull the latest version of the pipeline:

```bash
$ nextflow pull AJnsm/NF_TL_pipeline -r develop
> Checking AJnsm/NF_TL_pipeline ...
>  Already-up-to-date - revision: 86158df0d8 [develop]
```

And run the pipeline:

```bash
$ nextflow run AJnsm/NF_TL_pipeline -r develop -profile eddie_singularity -params-file 10X_sampleParams.json
> N E X T F L O W  ~  version 21.04.1
> Launching `AJnsm/NF_TL_pipeline` [grave_becquerel] - revision: 86158df0d8 [develop]
> [06/50026e] process > makeData (2)       [  0%] 0 of 2
> [-        ] process > estimatePCgraph    -
> [-        ] process > iterMCMCscheme     -
> [-        ] process > estimateCoups_1pts -
> [-        ] process > estimateCoups_2pts -
> [-        ] process > estimateCoups_3pts -
> Pulling Singularity image docker://ajnsm/py_nf_container_new [cache /gpfs/igmmfs01/eddie/ponting-lab/abel/10X_cells/NF_TL/vignette/work/singularity/ajnsm-py_nf_container_new.img]
```

It should automatically start pulling the singularity image from docker. This could take a few minutes the first time. Once that is done, nextflow should automatically start submitting jobs to the scheduler. 

Once the pipeline has finished, the output should look something like this: 


```bash
> N E X T F L O W  ~  version 21.04.1
> Launching `AJnsm/NF_TL_pipeline` [grave_becquerel] - revision: 86158df0d8 [develop]
[06/50026e] process > makeData (2)           [100%] 2 of 2 ✔
[0a/51146a] process > estimatePCgraph (2)    [100%] 2 of 2 ✔
[22/7e55a2] process > iterMCMCscheme (2)     [100%] 2 of 2 ✔
[f7/faa0b1] process > estimateCoups_1pts (4) [100%] 6 of 6 ✔
[44/a10772] process > estimateCoups_2pts (3) [100%] 6 of 6 ✔
[98/517d23] process > estimateCoups_3pts (3) [100%] 6 of 6 ✔
Completed at: 10-Jun-2021 14:33:19
Duration    : 5m 47s
CPU hours   : (a few seconds)
Succeeded   : 24

```

And we should see that the output directories have been generated:

```bash
$ ls -l
> Jun 10 13:40 10X_sampleClusters.csv
> Jun 10 13:34 10X_sampleData.csv
> Jun 10 14:07 10X_sampleParams.json
> Jun 10 14:01 10X_userGenes.csv
> Jun 10 14:01 10X_userGenes.csv~
> Jun 10 14:33 coupling_output/
> Jun 10 14:28 embeddings/
> Jun 10 14:31 output/
> Jun 10 14:05 params_test.json~
> Jun 10 14:28 plots/
> Jun 10 14:33 reports/
> Jun 10 14:31 work/

```


