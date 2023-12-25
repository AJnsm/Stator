# Stator
Preprint now available at https://www.biorxiv.org/content/10.1101/2023.12.18.572232v1

NOTE: this repo will soon move to the account of the [Edinburgh Biomedical AI Lab](https://edbiomed.ai)

## Table of contents
* [Introduction](#introduction)
* [Docs](#docs)
* [To do](#to-do)


## Introduction
The Stator pipeline takes in single cell RNA-seq count matrices, and estimates gene-gene interactions at up to seventh order. The higher-order interactions are then used to find characteristic, multi-type states present in the cell population. 

In [our research](https://www.biorxiv.org/content/10.1101/2023.12.18.572232v1), we used Stator to find cell identities that were invisible to clustering and NMF methods. We found sub-phases of the cell cycle, future neuronal fate states, and cell states predictive of liver cancer survival. 

The pipeline can be run directly from the command line, both locally or on an HPC cluster. It pulls all code from Github, and the required containers from Dockerhub. 

Subsequent analysis can be done with our bespoke Stator Shiny app, available at https://shiny.igc.ed.ac.uk/MFIs/

## Docs
Documentation on installation and usage are available [here](/docs)

## Changes in this version (V1.1)
- [X] Switched to Nextflow DSL2 and the latest Nextflow version (23.10)
- [X] Renamed multiple scripts, files, and directories
- [X] Updated and tested Docker profile for local runs
- [X] Removed 1-point calculation
- [X] Simplified Nextflow config files
- [X] Removed conda yamls
- [X] Switched back to numpy (was numba in V1.0) for interaction estimation (might be reverted in the future, depending on performance)

## To do
- [ ] switch to using CIs for significance estimation, abandon F-value. 
- [ ] update vignette
- [ ] make sure running locally is working, has a good profile, and is explained
- [ ] improve documentation and tutorial
- [ ] Create more unit tests for state inference
- [ ] add option to also use pairwise interactions for state inference













