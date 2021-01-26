# .libPaths(c(R.home("library"), .libPaths()))
.libPaths(R.home("library"))
library('pcalg')
# library('Rgraphviz')
library('igraph')
library('parallel')
library('stringr')
args = commandArgs(trailingOnly=TRUE)



### Please cite the following paper when using the code
### Paper: A fast PC algorithm for high dimensional causal discovery with multi-core PCs
skeleton_parallel_custom <- function(suffStat, indepTest, alpha, labels, p,
  method = c("parallel"),  mem.efficient=FALSE, workers, num_workers,
  m.max = Inf, fixedGaps = NULL, fixedEdges = NULL,
  NAdelete = TRUE, verbose = FALSE, identifier='UNIDENTIFIED', dirName="")
{
  
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
      } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
      else
      message("No need to specify 'p', when 'labels' is given")
    }
    seq_p <- seq_len(p)
    method <- match.arg(method)
    
  ## G := !fixedGaps, i.e. G[i,j] is true  iff  i--j  will be investigated
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
    diag(G) <- FALSE
    } else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
    else if (!identical(fixedGaps, t(fixedGaps)) )
    stop("fixedGaps must be symmetric")
    else
    G <- !fixedGaps
    
  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
    } else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
    else if (fixedEdges != t(fixedEdges))
    stop("fixedEdges must be symmetric")
    
    pval <- NULL
  sepset <- lapply(seq_p, function(.) vector("list",p))# a list of lists [p x p]
  ## save maximal p value
  pMax <- matrix(-Inf, nrow = p, ncol = p)
  diag(pMax) <- 1
  done <- FALSE
  ord <- 0
  n.edgetests <- numeric(1)# final length = max { ord}
  
  # edge test function, conditioning on x's neighbours
  edge_test_xy <- function(x, y) {
    G_xy <- TRUE
    num_tests_xy <- 0
    pMax_xy <- pMax[x, y]
    sepset_xy <- NULL
    done_xy <- TRUE
    if (G_xy && !fixedEdges[y, x]) {
      nbrsBool <- G.l[[x]]
      nbrsBool[y] <- FALSE
      nbrs <- seq_p[nbrsBool]
      #rm(nbrsBool)
      length_nbrs <- length(nbrs)
      if (length_nbrs >= ord) {
        if (length_nbrs > ord) done_xy <- FALSE
        S <- seq_len(ord)
        repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
        num_tests_xy <- num_tests_xy + 1
        pval <- indepTest(x, y, nbrs[S], suffStat)
          if(is.na(pval)) pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
          if (pMax_xy < pval) pMax_xy <- pval
          if(pval >= alpha) { # independent
            G_xy <- FALSE
            sepset_xy <- nbrs[S]
            break
            } else {
              nextSet <- getNextSet(length_nbrs, ord, S)
              if (nextSet$wasLast)
              break
              S <- nextSet$nextSet
            #rm(nextSet)
          } ## if (pval >= alpha)
        } ## repeat
        #rm(S)
      } ## if (length_nbrs >= ord)
    } ## if(!done)
    list(G_xy, sepset_xy, num_tests_xy, pMax_xy, done_xy)
  }
  
  # edge test function
  edge_test <- function(i) {
    x <- ind[i, 1]
    y <- ind[i, 2]
    num_tests_i <- 0
    G_i <- TRUE
    pMax_xy <- pMax[x, y]
    pMax_yx <- pMax[y, x]
    sepset_xy <- NULL
    sepset_yx <- NULL
    done_i <- TRUE
    
    # conditioning on neighbors of x
    res_x <- edge_test_xy(x, y)
    G_i <- res_x[[1]]
    sepset_xy <- res_x[[2]]
    num_tests_i <- num_tests_i + res_x[[3]]
    pMax_xy <- res_x[[4]]
    done_i <- done_i & res_x[[5]]
    
    if (G_i) {
      if (ord == 0) {
        num_tests_i <- num_tests_i + 1
        } else {
        # conditioning on neighbors of y
        res_y <- edge_test_xy(y, x)
        G_i <- res_y[[1]]
        sepset_yx <- res_y[[2]]
        num_tests_i <- num_tests_i + res_y[[3]]
        pMax_yx <- res_y[[4]]
        done_i <- done_i & res_y[[5]]
      }
    }
    
    # cleanup
    #rm(x)
    #rm(y)
    #rm(res_x)
    #rm(res_y)
    
    list(i, G_i, sepset_xy, sepset_yx, num_tests_i, pMax_xy, pMax_yx, done_i)
  }
  
  edge_tests <- function(l) {
    res <- vector("list",length(l))
    for (k in 1:length(l)) {
      res[[k]] <- edge_test(l[[k]])
    }
    res
  }

total_mem <- function() {
  tryCatch({
    if (Sys.info()[["sysname"]] == "Linux") {
      total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^MemFree:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Cached:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Inactive:' /proc/meminfo", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("egrep '^Buffers:' /proc/meminfo", intern = TRUE))))/1000
      return(total)
      } else if (Sys.info()[["sysname"]] == "Windows") {
      #total <- as.numeric(memory.limit())
      total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("wmic OS get FreePhysicalMemory /Value", intern=TRUE))[3]))/1000
      return(total)
    } else if (Sys.info()[["sysname"]] == "Darwin") { # Mac OS X
      total <- 4096*(as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages free'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages inactive'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages speculative'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vm_stat | grep 'Pages purgeable'", intern = TRUE))))/1000000
      return(total)
    }
    else { # other OS, i.e. Solaris
      total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vmstat -s | grep 'free memory'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vmstat -s | grep 'inactive memory'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vmstat -s | grep 'buffer memory'", intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("vmstat -s | grep 'swap cache'", intern = TRUE))))/1000
      return(total)
    }
    }, error=function(e){
      return(1024)
      }, warning=function(e){
        return(1024)
        })
}

parallel_threshold <- 100
if (mem.efficient) {
    mem_per_test <- 2 #MB
    tests_per_batch <- as.integer(total_mem() / mem_per_test)
  }
  
  start_total <- proc.time()
  
  while (!done && any(G) && ord <= m.max) {
    n.edgetests[ord1 <- ord+1L] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    ## For comparison with C++ sort according to first row
    ind <- ind[order(ind[, 1]), ]
    ## Consider only unique edge
    ind <- subset(ind, ind[, 1] < ind[, 2])
    remainingEdgeTests <- nrow(ind)
    ## Order-independent version: Compute the adjacency sets for any vertex
    ## Then don't update when edges are deleted
    G.l <- split(G, gl(p,p))
    
    if (!mem.efficient) {
      tests_per_batch <- remainingEdgeTests
    }
    
    for (j in seq(1, remainingEdgeTests, by=tests_per_batch)) {
      l <- min(remainingEdgeTests, j + tests_per_batch - 1)
      if (l - j + 1 < num_workers) {
        num_workers <- l - j + 1
      }
      res <- NULL
      if (l - j + 1 < parallel_threshold) {
        res <- lapply(j:l, edge_test)
        } else if (Sys.info()[['sysname']] == 'Windows') {
          res <- do.call("c", clusterApply(workers, clusterSplit(workers, j:l), edge_tests))
          } else {
            res <- mclapply(j:l, edge_test, mc.cores=num_workers, mc.set.seed=FALSE, mc.cleanup=TRUE, mc.allow.recursive=FALSE)
          }
          
      # synchronize
      for (p_obj in res) {
        i <- p_obj[[1]]
        x <- ind[i, 1]
        y <- ind[i, 2]
        n.edgetests[ord1] <- n.edgetests[ord1] + p_obj[[5]]
        pMax[x, y] <- p_obj[[6]]
        pMax[y, x] <- p_obj[[7]]
        G[x, y] <- G[y, x] <- p_obj[[2]]
        if (!p_obj[[2]]) {
          if (!is.null(p_obj[[3]])) sepset[[x]][[y]] <- p_obj[[3]]
          if (!is.null(p_obj[[4]])) sepset[[y]][[x]] <- p_obj[[4]]
        }
        done <- done & p_obj[[8]]
      }
    }
    
        # increase the nbrs size
        ord <- ord + 1

      # Write out partially finished graph
      ## transform matrix to graph object :
      Gobject <-
      if (sum(G) == 0) {
        new("graphNEL", nodes = labels)
        } else {
          colnames(G) <- rownames(G) <- labels
          as(G,"graphNEL")
        }
        tmp = new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
            max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
            sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
      
        write.csv(as(tmp,"matrix"), paste('SKEL_', identifier, '_alpha', alpha, '_level', ord-1, '.csv', sep=""))
      
        total_t = proc.time()-start_total
        cat('Num CI Tests=', n.edgetests, ',Total CI Tests=', sum(unlist(n.edgetests)), ',Total Time=', total_t[3], '\n', sep=" ")
        flush.console()

  } ## while()
  
  total_t = proc.time()-start_total
  
  # write results
  cat('n=', suffStat[[2]], ',p=', p, '\n', sep="")
  cat('Num CI Tests=', n.edgetests, ',Total CI Tests=', sum(unlist(n.edgetests)), ',Total Time=', total_t[3], '\n', sep=" ")
  
  for (i in 1:(p - 1)) {
    for (j in 2:p)
    pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j,i])
  }
  
  ## transform matrix to graph object :
  Gobject <-
  if (sum(G) == 0) {
    new("graphNEL", nodes = labels)
    } else {
      colnames(G) <- rownames(G) <- labels
      as(G,"graphNEL")
    }
    
  ## final object
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
    max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
    sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
  
}## end{ skeleton }
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
            
pc_parallel_custom <- function(suffStat, indepTest, alpha, labels, p,
      method = c("parallel"),  mem.efficient=FALSE, workers, num_workers,
      m.max = Inf, fixedGaps = NULL, fixedEdges = NULL,
      NAdelete = TRUE, verbose = FALSE, identifier='UNIDENTIFIED', dirName="",
      u2pd = c("relaxed", "rand", "retry"),
      conservative = FALSE, maj.rule = FALSE, solve.confl = FALSE)
{

  ## Initial Checks
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }

  u2pd <- match.arg(u2pd)
  if(u2pd != "relaxed") {
    if (conservative || maj.rule)
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")

    if (solve.confl)
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }

  if (conservative && maj.rule) stop("Choose either conservative PC or majority rule PC!")

  ## Skeleton
  skel <- skeleton_parallel_custom(suffStat = suffStat, indepTest = indepTest, alpha=alpha, labels = labels, method = method,
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                   NAdelete=NAdelete, m.max=m.max, verbose=verbose,
                   workers=NULL, num_workers=num_workers, identifier=identifier)
  skel@call <- cl # so that makes it into result

  ## Orient edges
  if (!conservative && !maj.rule) {
    switch (u2pd,
            "rand" = udag2pdag(skel),
            "retry" = udag2pdagSpecial(skel)$pcObj,
            "relaxed" = udag2pdagRelaxed(skel, verbose = verbose, solve.confl = solve.confl))
  }
  else { ## u2pd "relaxed" : conservative _or_ maj.rule

    ## version.unf defined per default
    ## Tetrad CPC works with version.unf=c(2,1)
    ## see comment on pc.cons.intern for description of version.unf
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                          version.unf = c(2,1), maj.rule = maj.rule, verbose = verbose)
      
    
    udag2pdagRelaxed(pc.$sk, verbose = verbose,
                     unfVect = pc.$unfTripl, solve.confl = solve.confl)
  }
} ## {pc}     
          
                   
                   
                   
                   
                   
                   
cat(args)
cat('\n')
cat(args[1])
cat('\n')                
dataPath = args[1]                   
nC = args[2]                   
DSname<-str_sub(str_split(dataPath, "_", n=2)[[1]][2], 1, -5)
alpha_<-args[3]

           
data = read.csv(dataPath)
cat('\n')
cat('dimension of used data:  ')   
cat(dim(data))
cat('Starting parallel PC...')
pcPar <- pc_parallel_custom(suffStat = list(dm = data, adaptDF = FALSE), binCItest, 
         method='parallel', labels=colnames(data), alpha=alpha_, workers=NULL, num_workers=nC, ##num_workers=detectCores(), 
         identifier=DSname, conservative = FALSE, maj.rule = TRUE, solve.confl = TRUE)

# Conservative = TRUE results in almost only bidirectionals

graphMat = as(pcPar@graph, 'matrix')                   
write.csv(graphMat, paste('PCgraph_', DSname, '.csv', sep=""))   

# Delete all edges, and write out empty control graph
graphMat[,]=0
write.csv(graphMat, paste('CTRLgraph_', DSname, '.csv', sep=""))
                   
cat('*****DONE*****')

                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   

