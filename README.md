# Stator
Preprint now available at https://www.biorxiv.org/content/10.1101/2023.12.18.572232v1

NOTE: this repo will soon move to the account of the [Edinburgh Biomedical AI Lab](https://edbiomed.ai)


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
Documentation on installation and usage are available [here](/docs/documentation.md).

A small tutorial/vignette is available [here](/vignette/Vignette.md).

## Changes in this version (V1.2)
- [X] Removed the identifyStates process--this is replaced by the Shiny R app.
- [X] Fixed a bug where the QC would filter out unexpressed genes, but not update nGenes in the output filenames. 
- [X] Default settings changed to output *all* positively enriched d-tuples (further selection should be done with the Shiny app).
- [X] Removed documentation on the `doubletFile` and `genesToOne` parameters (these will be deprecated), but functionality remains available for now. 

## To do
- [ ] switch to using CIs for significance estimation, abandon F-value. 
- [ ] Create more unit tests for state inference
- [ ] add option to also use pairwise interactions for state inference
- [ ] if pairwise not used, *optionally* skip their computation. 