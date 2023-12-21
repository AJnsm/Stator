# Stator
Preprint now available at https://www.biorxiv.org/content/10.1101/2023.12.18.572232v1

NOTE: this repo will soon move to the account of the [Edinburgh Biomedical AI Lab](https://edbiomed.ai)

## Table of contents
* [Introduction](#introduction)
* [Docs](#docs)
* [To do](#to-do)


## Introduction
The Stator pipeline takes in single cell RNA-seq count matrices, and estimates gene-gene interactions at up to seventh order. The 3-, 4-, and 5-point interactions among Markov-connected genes are used to find characteristic, multi-type states present in the cell population. 

The pipeline can be run directly from the command line. It pulls all code from Github, and the required containers from Dockerhub. It can run on your local machine (though I have not tested this thoroughly), or on a cluster with the Sun Grid Engine scheduling system (like Eddie). 

Subsequent analysis can be done with our bespoke R Shiny app, available from https://shiny.igc.ed.ac.uk/MFIs/

## To do
- [ ] remove references to eddie and change to generic "SGE"
- [ ] remove unnecessary configs
- [ ] make sure running locally is working, has a good profile, and is explained
- [ ] add and check conda environment yamls
- [ ] add option to just calculate log-odds ratios
- [ ] create better documentation
- [ ] show a full run on some example data. 
- [ ] create better tests













