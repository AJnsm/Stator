# Minimal working example: Finding Blood contamination among astrocytes

We will run Stator on some gene expression data from the mouse brain, publicly available from 10X. This directory contains the necessary files:
* `astrocytesVignette.csv.csv`: the raw count matrix (5000 rows × 100 columns) where each row is a cell, and each column a gene, with a header of gene names. Previous analysis has identified these cells as astrocyte-like. 
* `userGenes.csv`: A list of genes that should be included in the analysis. 

Let's have a look at the first 10 rows and columns of `astrocytesVignette.csv`:
```bash
$ head -n 10 astrocytesVignette.csv | cut -d ',' -f1-10
> ,Sox17,Gm37323,Mrpl15,Lypla1,Gm37988,Tcea1,Rgs20,Gm16041,Atp6v1h
> 0,0,0,2,0,0,1,1,0,0
> 1,0,0,1,0,0,0,0,0,1
> 2,0,0,0,0,0,0,1,0,0
> 3,0,0,1,0,0,2,2,0,2
> 4,0,0,1,0,0,2,0,0,0
> 5,0,0,0,0,0,0,1,0,0
> 6,0,0,1,0,0,1,0,0,1
> 7,0,0,1,0,0,2,0,0,0
> 8,0,0,0,1,0,1,1,0,0
```

For illustrative purposes, the file `userGenes.csv` contains the names of 19 genes--arbitrarily chosen from the columns of `astrocytesVignette.csv`--that should be included in the analysis, irrespective of any QC procedures that might follow:

```bash
$ cat userGenes.csv
> Serpine2,Ptgds,Cldn11,Neurog2,Isg15,Hbb-bt,Hbb-y,Crym,Ifitm3,Bst2,Fgf17,Eomes,Hba-a1,Hba-a2,Olig1,Hspa1a,Ttr,Iigp1,Ifit1
```

Next, we can set the pipeline parameters in a JSON file. There is an example file in this repo with the following settings:

```bash
$ cat params.json.json
>{
>      "dataType"    : "expression",
>      "rawDataPath" : "/exports/igmm/eddie/ponting-lab/abel/STATOR_vignette_run/astrocytesVignette.csv",
>      "userGenes"   : "/exports/igmm/eddie/ponting-lab/abel/STATOR_vignette_run/userGenes.csv",
>      "nGenes"      : 50,
>      "fracMito"    : 0.1,
>      "minGenes"    : 1,
>      "minCells"    : 1,
>      "nCells"      : 4000,
>      "PCalpha"     : 0.1,
>      "bsResamps"   : 1000,
>
>      "nRandomHOIs" : 10,
>      "sigHOIthreshold":0.1,
>      "minStateDeviation": 0,
>      "stateDevAlpha": 0.1,
>
>      "dendCutoff"  : -1,
>      "asympBool"   : 0,
>      "calcAll2pts" : 0,
>      "estimationMode" : "MFI",
>
>      "executor"    : "sge",
>      "maxQueueSize"  : 25,
>      "cores_makeData": 1,
>      "cores_PC"      : 6,
>      "cores_MCMC"    : 2,
>      "cores_1pt"     : 4,
>      "cores_2pt"     : 6,
>      "cores_HOIs_MB" : 6,
>      "cores_HOIs_6n7" : 6,
>      "cores_HOIs_plots": 6,
>
>      "mem_makeData"  : "16G",
>      "mem_PC"        : "32G",
>      "mem_MCMC"      : "16G",
>      "mem_1pt"       : "16G",
>      "mem_2pt"       : "16G",
>      "mem_HOIs_MB"   : "16G",
>      "mem_HOIs_6n7"   : "16G",
>      "mem_HOIs_plots" : "16G",
>
>      "time_makeData" : "1h",
>      "time_PC"       : "1h",
>      "time_MCMC"     : "1h",
>      "time_1pt"      : "1h",
>      "time_2pt"      : "1h",
>      "time_HOIs_MB"  : "1h",
>      "time_HOIs_6n7"  : "1h",
>      "time_HOIs_plots"  : "1h"
>  }
```
Some explanation:
* `"dataType" : "expression"` indicates that the pipeline should assume that the data is expression data, and perform some scRNA-seq specific QC, like filtering out low-read cells. More info on this is available from the documentation. 
* `nGenes : 50` indicates that a total of 50 genes should be included, which includes the 19 genes specified in `userGenes.csv`.
* `nCells : 4000`: The total number of cells used in the estimation. In practice, you would probably be interested in more cells and genes, but to keep things light, I've put in some very small numbers here.
* `"PCalpha" : 0.1`: When running the PC-algorithm, this is used as the significance threshold to conclude a conditional dependency. Setting this to a higher value makes it easier to reject the null hypothesis of no dependence, which leaves more causal edges present, which is more conservative.
* `"nRandomHOIs" : 10`: Beyond the Markov-connected tuples, Stator calculates this many interactions among random tuples at orders 3-7.
* `"sigHOIthreshold" : 0.1`: The significance threshold at which a higher-order interaction (HOI) is deemed significant. In practice, you might want to set this to 0.05 to be more conservative.
* `"minStateDeviation" : 0`: The minimum log-2 fold deviation from the expected number of occurrences of a particular gene state (d-tuple). Setting this to 0, as is done here, merely requires the deviation to be positive. 
* `"stateDevAlpha" : 0.1`: The significance threshold (after Benjamini-Yekutieli correction) to determine that the occurrence of a particular d-tuple deviates significantly. In practice, you might want to set this to 0.05 to be more conservative. 
* `"dendCutoff"  : -1`: The value of -1 tells Stator to cut the dendrogram of d-tuples at whichever height maximises the modularity score. Setting this to a value between 0 and 1 cuts at that particular Dice-coefficient. 
* `"asympBool" : 0`: This boolean 0 indicates that the confidence intervals and the uncertainty in the interaction estimation should be calculated using bootstrap resampling, not using an asymptotic estimation (`"asympBool" : 1`). 
* `"executor" : "sge"`: Run Stator on a cluster with the Sun Grid Engine scheduler.
* `"maxQueueSize"  : 25`: We don't want to be impolite and submit more than 25 jobs at the same time. 
* `"cores_*` specifies the number of cores/cpus each of the processes should request. The data preparation process can not parallelised, so 1 core is fine, but the others will be significantly faster on more cores. I've used up to 32 on Eddie.
* The `mem_*` and `time_*` specify how much memory and time each process gets. If these are too low, the whole pipeline terminates when the scheduler starts to complain. In practice, it can be useful to set these parameters relatively high the first time, and then look at the reports from Nextflow to see the peak usage of each of the resources. 


We are now ready to run the pipeline. Make sure you are on a node of the scheduler that allows for automated job submission (e.g. the wild-west nodes on Eddie), and that you are in the right directory:

```bash
[you@wild-west-node]$ pwd
> /absolute/path/to/this/vignette/
[you@wild-west-node]$ ls
> params.json  userGenes.csv  astrocytesVignette.csv  Vignette.md

```

Then pull the latest version of the pipeline:

```bash
$ nextflow pull AJnsm/Stator -r main
> Checking AJnsm/Stator ...
> Already-up-to-date - revision: 3ecffdbbdd [main]
```

And run the pipeline:

```bash
$ NXF_VER=23.04.4 nextflow run AJnsm/Stator -r main -profile eddie_singularity -params-file params.json

> N E X T F L O W  ~  version 23.04.4
> Launching `https://github.com/AJnsm/Stator` [modest_mccarthy] DSL2 - revision: 3ecffdbbdd [main]
executor >  sge (1)
[28/5d1f37] process > makeData                       [  0%] 0 of 1
[-        ] process > estimatePCgraph                -
[-        ] process > iterMCMCscheme                 -
[-        ] process > estimateCoups_2345pts_WithinMB -
[-        ] process > estimateCoups_6n7pts           -
[-        ] process > identifyDTuples                -
[-        ] process > identifyStates                 -
> Pulling Singularity image docker://ajnsm/py_nf_container_new
```

It should automatically start pulling the singularity image from docker. This could take a few minutes the first time. Once that is done, Nextflow will automatically start submitting jobs to the scheduler. 

Once the pipeline has finished, the output should look something like this: 


```bash
> N E X T F L O W  ~  version 22.04.3
> Launching `https://github.com/AJnsm/Stator` [jolly_turing] DSL1 - revision: 4902322ba5 [develop]
executor >  sge (11)
[28/5d1f37] process > makeData                       [100%] 1 of 1 ✔
[4c/882e6e] process > estimatePCgraph                [100%] 1 of 1 ✔
[73/532a35] process > iterMCMCscheme                 [100%] 1 of 1 ✔
[1b/e9e1e9] process > estimateCoups_2345pts_WithinMB [100%] 1 of 1 ✔
[96/b78172] process > estimateCoups_6n7pts           [100%] 1 of 1 ✔
[eb/b71633] process > identifyDTuples                [100%] 1 of 1 ✔
[fd/85ec7a] process > identifyStates                 [100%] 1 of 1 ✔
Completed at: 10-Jun-2024 11:32:49
Duration    : 4m 7s
CPU hours   : 0.1
Succeeded   : 7
```

It shows that in total, 0.1 CPU hours were used, which amounted to about 4 minutes on the Eddie cluster. Parallellisation happens both within in a process--where individual calculations get distributed over cores--and throughout the pipeline, where multiple processes can run in parallel if they don't depend on each other's output. 

After finishing, the following output directories should have been generated:

```bash
$ ls -l
> 10 jun 11:15 astrocytesVignette.csv
> 10 jun 12:00 coupling_output
> 10 jun 11:16 createVignetteData.ipynb
> 10 jun 12:00 dtuples_output
> 10 jun 11:58 output
> 10 jun 11:28 params.json
> 10 jun 12:00 reports
> 10 jun 12:00 states_output
> 10 jun 10:29 userGenes.csv
> 10 jun 12:00 work
```
The file `dtuples_output/top_DTuples.csv` shows the d-tuples that deviate and passed our significance threshold:

```bash
$ cat dtuples_output/top_DTuples.csv | cut -d ',' -f1-6 
> ,genes,state,enrichment,pval,pval_corrected
> 7,Hbb-bt_Hba-a1_Hba-a2,111,3.165433984247194,8.85974904654049e-45,1.926362578404946e-43
> 4,Adhfe1_Crym_Ifitm3,111,3.9210589776345297,1.7501130095553385e-16,1.902622857530875e-15
> 5,Hbb-bt_Hba-a1_Hba-a2,000,0.0746580402787477,2.0103515310271956e-06,1.4570262048682817e-05
> 2,Adhfe1_Crym_Ifitm3,101,0.5718584068512097,0.001955427518307478,0.01062914529594279
> 0,Adhfe1_Crym_Ifitm3,000,0.04538528788780138,0.0029326464768470203,0.012752822679317613
> 6,Hbb-bt_Hba-a1_Hba-a2,011,0.31033257847745366,0.017996142320096557,0.06521459193139752
> 3,Adhfe1_Crym_Ifitm3,110,0.8788154265835793,0.023879948910279794,0.0741740453907058
```
where the final column, that contains the cell IDs of the cells that are in that d-tuple state, is not shown. Looking at the rightmost value--the BY-corrected p-value--it seems that only the first 3 d-tuples here are actually significant in our data.

Most deviating and striking is the Globin-positive state where the genes *Hbb-bt*, *Hba-a1*, and *Hba-a2* are all expressed. This state is over 9 times enriched in the data set (2^3.17), which is very significant even after the BY-correction. In fact, the 3-point interaction among these genes led to a second significantly enriched d-tuple where all three genes are not expressed, though this state is only around 5% enriched in the data (2^0.07).

The d-tuples are combined to form cell states by cutting a hierarchical clustering at the Dice-coefficient that maximised the modularity score, which was 0.29 in this case (as can be seen in the file `states_output/modularity_scores.csv`). A figure showing this cut dendrogram is saved in `states_output/dendrogram_all_dTuples_cut.png`:

![an example dendrogram](https://github.com/AJnsm/Stator/blob/develop/vignette/dendrogram_vignette_example.png)

(Note that for this low number of genes and cells, the PCA embedding of the cells is a bit useless.)

To find exactly which d-tuples make up the *Hb+* state, the file `states_output/` lists exactly which cluster each d-tuple ended up in:

```bash
$ cat states_output/top_DTuples_withStatesFromCut.csv | awk -F',' '{print $2","$3","$NF}'
> ,genes,state,enrichment,pval,pval_corrected,CellIDs,geneState,cluster
> 7,Hbb-bt_Hba-a1_Hba-a2,111,3.165433984247194,8.85974904654049e-45,1.926362578404946e-43,4
> 4,Adhfe1_Crym_Ifitm3,111,3.9210589776345297,1.7501130095553385e-16,1.902622857530875e-15,6
> 5,Hbb-bt_Hba-a1_Hba-a2,000,0.0746580402787477,2.0103515310271956e-06,1.4570262048682817e-05,1
> 2,Adhfe1_Crym_Ifitm3,101,0.5718584068512097,0.001955427518307478,0.01062914529594279,2
> 0,Adhfe1_Crym_Ifitm3,000,0.04538528788780138,0.0029326464768470203,0.012752822679317613,1
> 6,Hbb-bt_Hba-a1_Hba-a2,011,0.31033257847745366,0.017996142320096557,0.06521459193139752,3
> 3,Adhfe1_Crym_Ifitm3,110,0.8788154265835793,0.023879948910279794,0.0741740453907058,5
```
which shows that the *Hb*+ d-tuple forms a state by itself (namely, the one indicated as cluster 4). That this population of astrocyte-like cells contain blood contamination has been confirmed in a larger data set involving a thousand genes in our study [here](https://www.biorxiv.org/content/10.1101/2023.12.18.572232v1).

As the interaction estimates depend on the Markov blankets, there should sufficient genes present in the analysis to get a good estimate of the causal relationships. We have verified in multiple gene expression data sets that the Markov blankets only stabilise once at least a few hundred of the most highly-variable genes are included in the analysis, so the estimates from this 50-gene data set cannot be taken seriously, but serve only as an illustrative example of the workflow. 


The output of Stator can be further analysed with the provided Shiny App, available at https://shiny.igc.ed.ac.uk/MFIs/. An interactive overview of the significant d-tuples and states at different Dice similarities can is provided when the appropriate files are uploaded:


![shiny app overview](https://github.com/AJnsm/Stator/blob/develop/vignette/shiny_app_table.png)

For illustrative purposes, we have created two artificial external cell type annotations that are just copies of the two most significant Stator states: the Globin expressing 'blood' state, and the *Adhfe1*+ *Crym*+ *Ifitm3*+ state. This annotation is provided in a csv file `annotation_vignette.csv`. With this file uploaded to the Shiny app, an overview of the enrichment of different Stator states in the external annotation is shown:

![shiny app overview](https://github.com/AJnsm/Stator/blob/develop/vignette/shiny_app_heatmap.png)

It can be seen that the app correctly identifies that the first two states 2 and 4 are significantly enriched in the external annotations for blood and *Adhfe1*+ *Crym*+ *Ifitm3*+, while the other states are not.



For more information, see [our study](https://www.biorxiv.org/content/10.1101/2023.12.18.572232v1) or [my thesis](https://abeljansma.nl/assets/JansmaAbel_PhDThesis_corrected.pdf), or raise an issue on this repository. 
