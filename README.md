# Stator
Preprint now available at https://www.biorxiv.org/content/10.1101/2023.12.18.572232v1

NOTE: this repo will soon move to the account of the [Edinburgh Biomedical AI Lab](https://edbiomed.ai)


NOTE: this is the main branch, but the most recent version is in the `develop` branch, which is much more streamlined and up-to-date. As soon as the develop branch is thoroughly tested, it will be merged into main.

## Table of contents
* [Introduction](#introduction)
* [Docs](#docs)
* [Changes in V1.1](#changes-in-this-version-v11)
* [To do](#to-do)


## Introduction
The Stator pipeline takes in single cell RNA-seq count matrices, and estimates gene-gene interactions at up to seventh order. Up to fifth-order interactions are then used to find characteristic, multi-type states present in the cell population. 

In [our research](https://www.biorxiv.org/content/10.1101/2023.12.18.572232v1), Stator found cell identities that were invisible to clustering and NMF methods. Stator found sub-phases of the cell cycle, future neuronal fate states, and liver cancer states predictive of patient survival. 

The pipeline can be run directly from the command line, both locally or on an HPC cluster. It pulls all code from Github, and the required containers from Dockerhub. 

Subsequent analysis can be done with our bespoke Stator Shiny app, available at https://shiny.igc.ed.ac.uk/MFIs/

## Docs
Documentation on installation and usage are available [here](/docs).

A small tutorial is available [here](/vignette/Vignette.md).

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
- [ ] improve documentation and tutorial
- [ ] Create more unit tests for state inference
- [ ] add option to also use pairwise interactions for state inference

