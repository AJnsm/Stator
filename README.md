# CoComics
**Co**nditional **Co**uplings in (transcript-)**omics** data. 


## Introduction
This pipeline takes in single cell RNA-seq count matrices, and estimates interactions between the expressed genes. 


Once nextflow is installed on your machine, the pipeline can be pulled from Github, and runs in a Docker/Singularity container on the data you provide it with. It can run on your local machine, or on a cluster with the Sun Grid Engine scheduling system (like Eddie). 


## Requirements

Nextflow needs to be installed on the machine you want to run the pipeline on, which is as simple as running:

```bash
curl -s https://get.nextflow.io | bash
```

Secondly, you must have access to either Docker or Singularity. Most clusters (like Eddie) use Singularity for security reasons, but both will automatically pull the right container from DockerHub. (Conda environments are supported, but not recommended or guaranteed to work.)

