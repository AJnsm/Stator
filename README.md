# Stator
Preprint now available at https://www.biorxiv.org/content/10.1101/2023.12.18.572232v1

NOTE: this repo will soon move to the account of the [Edinburgh Biomedical AI Lab](https://edbiomed.ai)

## Table of contents
* [Introduction](#introduction)
* [Docs](#docs)
* [To do](#to-do)


## Introduction
The Stator pipeline takes in single cell RNA-seq count matrices, and estimates gene-gene interactions at up to seventh order. The 3-, 4-, and 5-point interactions among Markov-connected genes are used to find characteristic, multi-type states present in the cell population. 

The pipeline can be run directly from the command line. It pulls all code from Github, and the required containers from Dockerhub. It can run on your local machine (though this is not tested thoroughly or recommended), or on a Sun Grid Engine (SGE) compatible platform, like the Edinburgh University compute cluster Eddie.

Subsequent analysis can be done with our bespoke Stator Shiny app, available from https://shiny.igc.ed.ac.uk/MFIs/

## Docs
Documentation on installation and usage are available [here](/docs)

## To do
- [ ] switch to DSL2
- [X] update example params file
- [X] change `pipelineScripts` to `scripts`
- [X] put all configs in a single directory
- [X] put all environment data in a single directory
- [ ] update vignette
- [ ] remove unnecessary dendrogram outputs and bootstrapping
- [ ] remove boundval from output
- [ ] output order 1 interactions in same format as the rest.
- [ ] remove references to Eddie and change to generic "SGE"
- [ ] remove unnecessary configs
- [ ] make sure running locally is working, has a good profile, and is explained
- [ ] add and check conda environment yamls
- [X] add option to just calculate log-odds ratios
- [ ] improve documentation and tutorial
- [ ] Update vignette
- [ ] Create more unit tests for state inference
- [ ] add option to also use 2-points for state inference













