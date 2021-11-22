.libPaths(R.home("library"))
library(graph)
cat('graph loaded')
library(pcalg)
cat('pcalg loaded')
library(igraph)
cat('igraph loaded')
library(stringr)
cat('stringr loaded')
# library(devtools)
# cat('devtools loaded')
# if (!require("BiDAG")) install.packages("BiDAG", repos='https://cloud.r-project.org', INSTALL_opts = '--no-lock')

if (!require("BiDAG")){
	install_github("cran/BiDAG")
}
library(BiDAG)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("3 arguments must be supplied. \n")
}
cat(args)
nGenes = as.numeric(args[3])

PCgraph_csv <- args[1]
PCadjMat <- read.csv(PCgraph_csv, row.names=1)
# PCgraph  <- graph_from_adjacency_matrix(as(PCadjMat , 'matrix'), mode="directed")



dataPath <- args[2]
dataSet <- read.csv(dataPath, header=T)
    
# This is safer, since sometimes there are more userGenes than nGenes:

nGenes <- dim(dataSet)[2]

DSname<-str_sub(str_split(dataPath, "_", n=2)[[1]][2], 1, -5)

print(paste('Starting MCMC MAP estimation on cluster', DSname))
cat(dim(dataSet))
cat('\n')
cat(nGenes)
sp <- scoreparameters(n = nGenes, scoretype = "bde", dataSet)
system.time(result_iterMCMC <- iterativeMCMC(sp, MAP=TRUE, verbose=TRUE, chainout=FALSE,
                                         startspace=PCadjMat, cpdag=TRUE))




write.csv(result_iterMCMC$endspace, paste('MCMCgraph_', DSname, '.csv', sep=''), row.names=TRUE)

# To read out and plot: tmp <- read.csv('iterMCMC...'), then set row.names(tmp) <- colnames(tmp), and run
# plotFn(graph_from_adjacency_matrix(as(tmp, 'matrix')), ...)
