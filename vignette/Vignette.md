# Minimal working example

NOTE: this is slightly outdated since the latest version only calculates Markov-connected 2-point interactions by default, and has some file name changes. 

We will run Stator on some gene expression data from astrocytes, publicly available from 10X. This directory contains the necessary files:
* `vignetteAstrocyteData.csv`: the raw count matrix (23900 rows × 35 columns) where each row is a cell, and each column a gene, with a header of gene names.
* `vignetteDoublets.csv`: a doublet annotation for each cell, in the same order as the cells in `vignetteAstrocyteData.csv`. In this case, it is just a list of zeros, which indicates none of the listed cells are annotated as doublets, so it might as well be excluded.
* `userGenes.csv`: A list of genes that should be included in the analysis. 

Let's have a look at the first 10 rows and columns of `vignetteAstrocyteData.csv`:
```bash
$ head -n 10 vignetteAstrocyteData.csv | cut -d ',' -f1-10
> Serpine2,Ptgds,Cldn11,Neurog2,Isg15,Hbb-bt,Hbb-y,Crym,Ifitm3,Bst2
> 0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0
> 0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0
> 0.0,0.0,0.0,0.0,0.0,2.0,0.0,0.0,0.0,0.0
> 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
> 4.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
> 0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
> 0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0
> 6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
> 0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0
```

For illustrative purposes, the file `userGenes.csv` contains the names of 19 genes--taken from the columns of `vignetteAstrocyteData.csv`--that should be included in the analysis, irrespective of any QC procedures that might follow:

```bash
$ cat userGenes.csv
> Serpine2,Ptgds,Cldn11,Neurog2,Isg15,Hbb-bt,Hbb-y,Crym,Ifitm3,Bst2,Fgf17,Eomes,Hba-a1,Hba-a2,Olig1,Hspa1a,Ttr,Iigp1,Ifit1
```

Next, we can set the pipeline parameters in a JSON file. There is an example file in this repo with the following settings:

```bash
$ cat params.json.json
> {
>    "dataType"    : "agnostic",
>    "rawDataPath" : "/exports/igmm/eddie/ponting-lab/abel/10X_cells/NF_TL/paperRuns/vignette_test/vignetteAstrocyteData.csv",
>    "doubletFile" : "/exports/igmm/eddie/ponting-lab/abel/10X_cells/NF_TL/paperRuns/vignette_test/vignette_doublets.csv",
>    "userGenes"   : "/exports/igmm/eddie/ponting-lab/abel/10X_cells/NF_TL/paperRuns/vignette_test/userGenes.csv",
>    "nGenes"      : 19,
>    "nCells"      : 2000,
>    "PCalpha"     : 0.05,
>
>    "nRandomHOIs" : 10,
>    "sigHOIthreshold":0.4,
>    "minStateDeviation": 0,
>    "stateDevAlpha": 0.9,
>
>    "dendCutoff"  : -1,
>    "asympBool"   : 1,
>
>    "executor"    : "sge",
>    "maxQueueSize"  : 25,
>    "cores_makeData": 1,
>    "cores_PC"      : 6,
>    "cores_MCMC"    : 2,
>    "cores_1pt"     : 4,
>    "cores_2pt"     : 12,
>    "cores_HOIs_MB" : 6,
>    "cores_HOIs_6n7" : 6,
>    "cores_HOIs_plots": 6,
>
>    "mem_makeData"  : "32G",
>    "mem_PC"        : "16G",
>    "mem_MCMC"      : "16G",
>    "mem_1pt"       : "32G",
>    "mem_2pt"       : "64G",
>    "mem_HOIs_MB"   : "64G",
>    "mem_HOIs_6n7"   : "16G",
>    "mem_HOIs_plots" : "64G",
>
>    "time_makeData" : "1h",
>    "time_PC"       : "1h",
>    "time_MCMC"     : "1h",
>    "time_1pt"      : "1h",
>    "time_2pt"      : "1h",
>    "time_HOIs_MB"  : "1h",
>    "time_HOIs_6n7"  : "1h",
>    "time_HOIs_plots"  : "1h"
>  }
```
Some explanation:
* `"dataType" : "agnostic"` indicates that the pipeline should be agnostic with respect to the data type, and not perform any scRNA-seq specific QC. If the data set was larger than this simple test set, this could be set to `expression`, and additional QC parameters could be set to e.g. filter out low-read cells. More info on this in the documentation. 
* `nGenes : 25` indicates that 25 genes should be included, which in this case includes the 19 genes specified in `userGenes.csv`.
* `nCells : 2000`: The total number of cells used in the estimation. In reality, you would probably be interested in more cells and genes, but to keep things light, I've put in some very small numbers here.
* `"PCalpha" : 0.05`: When running the PC-algorithm, this is used as the significance threshold to conclude a conditional dependency. Setting this to a higher value makes it easier to reject the null hypothesis of no dependence, which leaves more causal edges present, which is more conservative.
* `"nRandomHOIs" : 10`: Beyond the Markov-connected tuples, Stator calculates this many interactions among random tuples at orders 3-7.
* `"sigHOIthreshold" : 0.4`: The significance threshold at which a higher-order interaction (HOI) is deemed significant. In real data sets, this should be set to a lower value like 0.05. 
* `"minStateDeviation" : 0`: The minimum log-2 fold deviation from the expected number of occurrences of a particular gene state (d-tuple). Setting this to 0, as is done here, merely requires the deviation to be positive. 
* `"stateDevAlpha" : 0.9`: The significance threshold (after Benjamini-Yekutieli correction) to determine that the occurrence of a particular d-tuple deviates significantly. Should be set to a lower value like 0.05 in practice. 
* `"dendCutoff"  : -1`: The value of -1 tells Stator to cut the dendrogram of d-tuples at whichever height maximises the modularity score. Setting this to a value between 0 and 1 cuts at that particular Dice-coefficient. 
* `"asympBool" : 1`: This boolean 1 indicates that the confidence intervals and the uncertainty in the interaction estimation should be calculated using an asymptotic estimation. The alternative (`"asympBool" : 0`) would estimate these statistics using bootstrap resampling.

* `"executor" : "sge"`: Run Stator on a cluster with the Sun Grid Engine scheduler.
* `"maxQueueSize"  : 25`: We don't want to be impolite and submit more than 25 jobs at the same time. 
* `"cores_*` specifies the number of cores/cpus each of the processes should request. The data preparation process can not parallelised, so 1 core is fine, but the others will be significantly faster on more cores. I've used up to 32 on Eddie.
* The `mem_*` and `time_*` specify how much memory and time each process gets. If these are too low, the whole pipeline terminates when the scheduler starts to complain. What I do: set them relatively high the first time, and then look at the reports from Nextflow to see the peak usage of each of the resources. 


We are now ready to run the pipeline. Make sure you are on a node of the scheduler that allows for automated job submission (e.g. the wild-west nodes on Eddie), and that you are in the right directory:

```bash
[you@wild-west-node]$ pwd
> /absolute/path/to/this/vignette/
[you@wild-west-node]$ ls
> params.json  userGenes.csv  vignetteAstrocyteData.csv  vignetteDoublets.csv Vignette.md

```

Then pull the latest version of the pipeline:

```bash
$ nextflow pull AJnsm/Stator -r develop
> Checking AJnsm/Stator ...
>  Already-up-to-date - revision: 520967e399 [develop]
```

And run the pipeline:

```bash
$ nextflow run AJnsm/Stator -r develop -profile eddie_singularity -params-file params.json -resume

> N E X T F L O W  ~  version 22.04.3
> Launching `https://github.com/AJnsm/Stator` [elegant_rosalind] DSL1 - revision: 520967e399 [develop]
executor >  sge (1)
[0e/6c9394] process > makeData (1)                  [0%] 0 of 1
[-        ] process > estimatePCgraph               -
[-        ] process > iterMCMCscheme                -
[-        ] process > estimateCoups_1pts            -
[-        ] process > estimateCoups_2pts            -
[-        ] process > estimateCoups_345pts_WithinMB -
[-        ] process > estimateCoups_6n7pts          -
[-        ] process > createHOIsummaries            -
[-        ] process > identifyStates                -
> Pulling Singularity image docker://ajnsm/py_nf_container_new
```

It should automatically start pulling the singularity image from docker. This could take a few minutes the first time. Once that is done, Nextflow should automatically start submitting jobs to the scheduler. 

Once the pipeline has finished, the output should look something like this: 


```bash
> N E X T F L O W  ~  version 22.04.3
> Launching `https://github.com/AJnsm/Stator` [jolly_turing] DSL1 - revision: 4902322ba5 [develop]
executor >  sge (11)
[99/96d88e] process > makeData                          [100%] 1 of 1 ✔
[cd/06c207] process > estimatePCgraph (1)               [100%] 1 of 1 ✔
[de/86eeb2] process > iterMCMCscheme (1)                [100%] 1 of 1 ✔
[cf/47e197] process > estimateCoups_1pts (2)            [100%] 2 of 2 ✔
[2b/d30c10] process > estimateCoups_2pts (2)            [100%] 2 of 2 ✔
[98/c3dc1b] process > estimateCoups_345pts_WithinMB (1) [100%] 1 of 1 ✔
[be/871de0] process > estimateCoups_6n7pts (1)          [100%] 1 of 1 ✔
[a2/2cde80] process > createHOIsummaries (1)            [100%] 1 of 1 ✔
[ee/51b66d] process > identifyStates (1)                [100%] 1 of 1 ✔
Completed at: 06-Mar-2023 14:39:32
Duration    : 20m 37s
CPU hours   : 5.2
Succeeded   : 11
```

It shows that in total, 5.2 CPU hours were used, which amounted to about 20 minutes on the Eddie cluster. The parallellisation happens both within in a process--where individual calculations get distributed over cores--and throughout the pipeline, where multiple processes can run in parallel if they don't depend on each other's output. 

After finishing, the following output directories should have been generated:

```bash
$ ls -l
> Mar  6 14:37 HOIsummaries/
> Mar  6 14:39 coupling_output/
> Mar  6 14:23 output/
> Mar  6 14:13 params.json
> Mar  6 14:39 reports/
> Mar  6 14:38 states_output/
> Mar  6 13:20 userGenes.csv
> Mar  6 13:50 vignetteAstrocyteData.csv
> Mar  6 14:14 vignetteDoublets.csv
> Mar  6 14:37 work/
```

The file `HOIsummaries/top_DTuples.csv` shows the d-tuples that deviate and passed our significance threshold:

```bash
$ head -n 10 HOIsummaries/top_DTuples.csv | cut -d ',' -f1-6
> ,genes,state,enrichment,pval,pval_corrected
> 41,Hba-a1_Hbb-bt_Hba-a2,111,3.2907442436316536,5.042731009532022e-34,1.> 5021325049175918e-31
> 37,Neurog2_Ifitm3_Eomes,101,2.156699225328587,1.1010156294693394e-17,1.> 6398568179841119e-15
> 32,Neurog2_Crym_Ifitm3,011,1.879413199401285,1.870648016253587e-07,1.> 8574334589482004e-05
> 39,Hba-a1_Hbb-bt_Hba-a2,000,0.10297011997190542,8.219424660847979e-06,0.> 0006121020996084621
> 23,Hspa1a_Crym_Ifitm3,011,1.4106091093783257,3.5257908890442744e-05,0.> 0019585150895565476
> 28,Lypla1_Crym_Ifitm3,011,1.5227306274895436,3.944897567589089e-05,0.> 0019585150895565476
> 5,Bst2_Crym_Ifitm3,011,1.63104410329908,5.8684719756655876e-05,0.0024972926354453187
> 18,Isg15_Bst2_Ifitm3,011,1.3791641581601983,0.0013545459037925451,0.050436644043764516
> 54,Isg15_Bst2_Hba-a2_Ifitm3,0101,1.3237674496085021,0.001924503628319066,0.06369695421137356

```
where I've cut off the final column, which contains the cell IDs of the cells that are in that d-tuple state. Looking at the rightmost value--the BY-corrected p-value--I would argue that the first 7 d-tuples here are actually significant in our data.

Most deviating and striking is the Globin-positive state where the genes *Hbb-bt*, *Hba-a1*, and *Hba-a2* are all expressed. This state is almost 10 times enriched in the data set (2^3.3), which is very significant even after the BY-correction. In fact, the 3-point interaction among these genes led to a second significantly enriched d-tuple where all three genes are not expressed, though this state is only around 7% enriched in the data (2^0.1).

The d-tuples are combined to form cell states by cutting a hierarchical clustering at the Dice-coefficient that maximised the modularity score, which was 0.83 in this case (as can be seen in the file `states_output/modularity_scores.csv`). A figure showing this cut dendrogram is saved in `states_output/dendrogram_all_dTuples_cut.png`:

![an example dendrogram](https://github.com/AJnsm/Stator/blob/develop/vignette/dendrogram_vignette_example.png)

(Note that for this crude selection of genes and cells, the PCA embedding of the cells is a bit useless. In practice, gene selection would at least in part be based QC metrics like high variance *etc.*, which Stator can automatically apply when running in `expression` mode.)

To find exactly which d-tuples make up the *Hb+* state, the file `states_output/` lists exactly which cluster each d-tuple ended up in:

```bash
$ cat states_output/top_DTuples_withStatesFromCut.csv | cut -d '[' -f3
> ,genes,state,enrichment,pval,pval_corrected,CellIDs,geneState,cluster
> 'Hba-a1+', 'Hbb-bt+', 'Hba-a2+']",5
> 'Neurog2+', 'Ifitm3-', 'Eomes+']",3
> 'Neurog2-', 'Crym+', 'Ifitm3+']",6
> 'Hba-a1-', 'Hbb-bt-', 'Hba-a2-']",4
> 'Hspa1a-', 'Crym+', 'Ifitm3+']",6
> 'Lypla1-', 'Crym+', 'Ifitm3+']",6
> 'Bst2-', 'Crym+', 'Ifitm3+']",6
> 'Isg15-', 'Bst2+', 'Ifitm3+']",9
> 'Isg15-', 'Bst2+', 'Hba-a2-', 'Ifitm3+']",9
> 'Neurog2+', 'Ifitm3+', 'Eomes+']",1
> 'Lypla1+', 'Crym+', 'Ifitm3+']",6
> 'Bst2+', 'Hba-a2-', 'Ifitm3+']",9
> 'Isg15+', 'Bst2+', 'Hba-a2+', 'Ifitm3+']",10
> 'Isg15+', 'Hba-a2-', 'Ifitm3+']",8
> 'Isg15+', 'Bst2-', 'Ifitm3+']",8
> 'Bst2+', 'Crym-', 'Ifitm3+']",9
> 'Atp6v1h+', 'Hba-a2-', 'Eomes+']",3
> 'Isg15+', 'Bst2+', 'Hba-a2+']",10
> 'Hspa1a+', 'Crym+', 'Ifitm3-']",2
> 'Lypla1+', 'Crym-', 'Ifitm3+']",5
> 'Isg15+', 'Bst2-', 'Hba-a2-', 'Ifitm3+']",8
> 'Isg15+', 'Bst2+', 'Ifitm3+']",10
> 'Neurog2-', 'Ifitm3-', 'Eomes-']",4
> 'Hspa1a+', 'Crym+', 'Ifitm3+']",7
> 'Isg15+', 'Bst2+', 'Hba-a2-', 'Ifitm3+']",10
> 'Atp6v1h-', 'Hba-a2+', 'Eomes-']",5
> 'Atp6v1h+', 'Hba-a2+', 'Eomes+']",1
> 'Hspa1a-', 'Crym-', 'Ifitm3-']",4
> 'Isg15+', 'Bst2+', 'Hba-a2+', 'Ifitm3-']",11
> 'Hspa1a+', 'Neurog2-', 'Ifitm3+']",7
> 'Bst2+', 'Hba-a2+', 'Ifitm3+']",10
> 'Lypla1-', 'Crym-', 'Ifitm3-']",4
> 'Bst2-', 'Crym-', 'Ifitm3-']",4
> 'Isg15-', 'Bst2-', 'Hba-a2-', 'Ifitm3-']",4
> 'Isg15-', 'Hba-a2+', 'Ifitm3+']",5
```
which shows that the *Hb*-genes ended up in cluster 5 with three other d-tuples. However, none of these other d-tuples are actually significant to any reasonable level, so the *Hb*-gene d-tuple seems to form a state by itself. (This is something I have confirmed in a larger data set involving a thousand genes)

As the interaction estimates depend on the Markov blankets, there should sufficient genes present in the analysis to get a good estimate of the causal relationships. I have verified in multiple gene expression data sets that the Markov blankets only stabilise once at least a few hundred of the most highly-variable genes are included in the analysis, so the estimates from this 25-gene data set cannot be taken seriously, but serve only as an illustrative example of the workflow. 

For more information, see my thesis [link]. 