if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("graph")

BiocManager::install("RBGL")

BiocManager::install("graph")

install.packages("ggm")

install.packages("pcalg")

BiocManager::install("Rgraphviz")

install.packages("BiDAG")


## Dependencies: conda create -n rEnv r-base=4.0.3  bioconductor-graph bioconductor-rbgl r-ggm r-pcalg bioconductor-rgraphviz r-stringr
## Then add BiDAG from within R, make sure .libPaths() is pointing to the rEnv directory. 