library(Seurat)
library(PISCES)
library(apcluster)
library(aricode)
library(reshape2)
library(ggplot2)
library(mcclust)
library(philentropy)

#' A function for an exhaustive query of Affinity Clustering results when used with
#' four possible distance metrics: The Pearson & Spearman correlations between
#' samples, as well as Manhattan & Euclidean distances.  The function will cluster
#' using Affinity Propagation clustering, with each sample's own similarity being
#' set as the 0th (min), 25th, 50th (median), 75th, and 100th (max) percentiles
#' of all data in the similarity matrix.  The fewest clusters should be generated
#' when the diagonal is set to the minimum.  This function can also be applied
#' only on 'npcs' Principal Components of the data, and an iteration option is
#' available to cluster on 2,3,...,npcs principal components using the distance
#' metric 'iteratemethod' iteratively, to view how the number of PCs impacts
#' the number of clusters.  data.obj should be transformed with SCTransform
#' 
#' @param data.obj a Seurat object. Cells are columns and features are rows.
#' @param skipnonpca set TRUE if you want to skip any calculations not on PCs
#' @param spearman test the Spearman distance metric (both non-PCA and PCA)
#' @param pearson test the Pearson distance metric (both non-PCA and PCA)
#' @param negmanhattan test the Manhattan distance metric (negative) (both non-PCA and PCA)
#' @param negeuclidean test the Euclidean distance metric (negative) (both non-PCA and PCA)
#' @param PCA set TRUE (default) to perform tests on 'npcs' many PCs of the data
#' @param npcs number of PCs to test on.  Iterating also only reaches this many PCs
#' @param iterate set TRUE if you want clusters vs. npcs for 'iteratemethod' metric
#' @param iteratemethod the distance metric used for PC iteration (must choose only one)
#' @param step the step size when iterating over PCs
#' @return A named list of lists. APCluster results
#' @export
APTest <- function(data.obj, skipnonpca = FALSE, spearman = TRUE, pearson = TRUE, negmanhattan = TRUE, 
                   negeuclidean = TRUE, PCA = TRUE, npcs = 30, iterate = FALSE, iteratemethod = 'manhattan', step = 1) { 
  
  ##### NON-PCA PHASE
  if(!skipnonpca) {
    print("BEGIN NON-PCA PHASE")
    nonpca <- list()
    
    cpm.mat <- as.matrix(data.obj@assays$SCT@data)
    sct.mat <- as.matrix(data.obj@assays$SCT@scale.data) #sct.mat less features than cpm
    
    if (spearman) {
      print("BEGIN NON-PCA SPEARMAN")
      #Initially, do spearman correlation on both matrices
      cpm.temp <- cor(cpm.mat, method = 'spearman')
      sct.temp <- cor(sct.mat, method = 'spearman')
      
      #cluster all quantiles and put into list
      list <- append(ClusterSimQuantiles(cpm.temp, "cpm.spearman"), ClusterSimQuantiles(sct.temp, "sct.spearman")) 
      nonpca <- append(nonpca, list)
    }
    
    if (pearson) {
      print("BEGIN NON-PCA PEARSON")
      #Initially, do pearson correlation on both matrices
      cpm.temp <- cor(cpm.mat, method = 'pearson')
      sct.temp <- cor(sct.mat, method = 'pearson')
      
      #cluster all quantiles and add into list
      list <- append(ClusterSimQuantiles(cpm.temp, "cpm.pearson"), ClusterSimQuantiles(sct.temp, "sct.pearson")) 
      nonpca <- append(nonpca, list)
    }
    
    #Negative manhattan distance
    if (negmanhattan) {
      print("BEGIN NON-PCA NEGMANHATTAN")
      # Must transpose because *.mat is featsxcells
      cpm.temp <- as.matrix(dist(t(cpm.mat), method = 'manhattan'))
      sct.temp <- as.matrix(dist(t(sct.mat), method = 'manhattan'))
      # Since we've measured distance, and we want similarity, we negate absolute val.
      # High manhattan distance (neg or pos) will now be absolutely negative.
      cpm.temp <- -abs(cpm.temp)
      sct.temp <- -abs(sct.temp)
      
      #cluster all quantiles and add into list
      list <- append(ClusterSimQuantiles(cpm.temp, "cpm.manhattan"), ClusterSimQuantiles(sct.temp, "sct.manhattan")) 
      nonpca <- append(nonpca, list)
    }
    
    if (negeuclidean) {
      print("BEGIN NON-PCA NEGEUCLIDEAN")
      # Must transpose because *.mat is featsxcells
      cpm.temp <- as.matrix(dist(t(cpm.mat), method = 'euclidean'))
      sct.temp <- as.matrix(dist(t(sct.mat), method = 'euclidean'))
      # Since we've measured distance, and we want similarity, we negate absolute val.
      # High manhattan distance (neg or pos) will now be absolutely negative.
      cpm.temp <- -abs(cpm.temp)
      sct.temp <- -abs(sct.temp)
      
      #cluster all quantiles and add into list
      list <- append(ClusterSimQuantiles(cpm.temp, "cpm.euclidean"), ClusterSimQuantiles(sct.temp, "sct.euclidean")) 
      nonpca <- append(nonpca, list)
    }
  }
  
  
  ##### PCA PHASE
  if (is.null(data.obj@reductions[["pca"]])) {
    PCA = FALSE
    print("Warning: Tried to cluster PCs on non-PCA transformed Seurat obj")
  }
  if (PCA) {
    print("BEGIN PCA PHASE")
    if (npcs > length(data.obj@reductions[["pca"]]) || npcs < 2) {
      npcs <- length(data.obj@reductions[["pca"]])
      print("Warning: Tried to cluster invalid # of PCs")
    }
    pca.mat <- as.matrix(data.obj@reductions[["pca"]]@cell.embeddings)[,1:npcs]
    pca <- list()
    
    if (spearman) {
      print("BEGIN PCA SPEARMAN")
      #Initially, do spearman correlation on both matrices.  Must transpose mat
      #as *.mat is cellsxpcs (needs pcsxcells)
      pca.temp <- cor(t(pca.mat), method = 'spearman')
      #cluster all quantiles and put into list
      list <- ClusterSimQuantiles(pca.temp, "pca.spearman")
      pca <- append(pca, list)
    }
    
    if (pearson) {
      print("BEGIN PCA PEARSON")
      #Initially, do spearman correlation on both matrices.  Must transpose mat
      #as *.mat is cellsxpcs (needs pcsxcells)
      pca.temp <- cor(t(pca.mat), method = 'pearson')
      #cluster all quantiles and put into list
      list <- ClusterSimQuantiles(pca.temp, "pca.pearson")
      pca <- append(pca, list)
    }
    
    if (negmanhattan) {
      print("BEGIN PCA NEGMANHATTAN")
      # DON'T transpose PCA because *.mat is cellsxpcs
      pca.temp <- as.matrix(dist(pca.mat, method = 'manhattan'))
      # Since we've measured distance, and we want similarity, we negate absolute val.
      # High manhattan distance (neg or pos) will now be absolutely negative.
      pca.temp <- -abs(pca.temp)
      
      #cluster all quantiles and add into list
      list <- ClusterSimQuantiles(pca.temp, "pca.manhattan")
      pca <- append(pca, list)
    }
    
    if (negeuclidean) {
      print("BEGIN PCA NEGEUCLIDEAN")
      # DON'T transpose PCA because *.mat is cellsxpcs
      pca.temp <- as.matrix(dist(pca.mat, method = 'euclidean'))
      # Since we've measured distance, and we want similarity, we negate absolute val.
      # High manhattan distance (neg or pos) will now be absolutely negative.
      pca.temp <- -abs(pca.temp)
      
      #cluster all quantiles and add into list
      list <- ClusterSimQuantiles(pca.temp, "pca.euclidean")
      pca <- append(pca, list)
    }
    
    if (iterate) {
      print("BEGIN PCA ITERATION")
      #GENERATE PC SETTINGS TO ITERATE OVER
      pcasteps <- list()
      steps <- c()
      i <- 2
      while (i <= npcs) {
        steps <- append(steps, i)
        i <- step + i
      }
      
      #Just in case user selected something other than the 4 options
      if (iteratemethod != 'spearman' && iteratemethod != 'pearson' &&
          iteratemethod != 'manhattan' && iteratemethod != 'euclidean') {
        print("Warning: Method not defined for PCA iteration; choosing manhattan")
        iteratemethod <- 'manhattan'
      }
      
      if (iteratemethod == 'spearman' || iteratemethod == 'pearson') {
        pcasteps<- c()
        for (i in steps) {
          print(paste("ITERATION WITH PCs: ", i, sep=''))
          pca.temp <- cor(t(pca.mat[,1:i]), method = iteratemethod)
          list <- list(ClusterSimQuantiles(pca.temp, paste(iteratemethod, "_", i, "PCs", sep="")))
          pcasteps <- append(pcasteps, list)
        }
      }
      
      if (iteratemethod == 'manhattan' || iteratemethod == 'euclidean') {
        pcasteps<- c()
        for (i in steps) {
          print(paste("ITERATION WITH PCs: ", i, sep=''))
          pca.temp <- as.matrix(dist(pca.mat[,1:i], method = iteratemethod))
          pca.temp <- -abs(pca.temp)
          list <- list(ClusterSimQuantiles(pca.temp, paste(iteratemethod, "_", i, "PCs", sep="")))
          pcasteps <- append(pcasteps, list)
        }
      }
      pcastepslist <- list(pcasteps)
      names(pcastepslist) <- (paste(npcs, "PCit_", iteratemethod, "_", step, "=step",sep = ''))
    }
  }
  
  ##### EXPORT RESULTS PHASE
  print("EXPORTING RESULTS")
  apresults <- list()
  if (!skipnonpca && PCA) {
    apresults <- list(nonpca, pca)
    names(apresults) <- c("NON-PCA", paste(npcs, "-PCS", sep="")) 
  } else if (!skipnonpca) {
    apresults <- list(nonpca)
    names(apresults) <- c("NON-PCA")
  } else if (PCA) {
    apresults <- list(pca)
    names(apresults) <- c(paste(npcs, "-PCS", sep=""))
  }
  if (iterate == TRUE) {
    apresults <- append(apresults, pcastepslist)
  }
  return(apresults)
}

#' Takes an apresult returned by APTest and makes a data table that summarizes
#' the number of exemplars produced by each distance metric, using each quantile
#' of the data, on the CPM matrix returned from Seurat::SCTransform.  
#' All similarity metrics must be used in this APTest, and there must
#' be a NON-PCA analysis (skipnonpca = FALSE).  
#'
#' @param apresult an AP clustering result returned by APTest()
#' @return a table that summarizes the NON-PCA CPM experiment
#' @export
SummarizeAPTest_NONPCA_CPM <- function(apresult) {
  apresult <- apresult[["NON-PCA"]]
  # Normalized Counts (SCT Pre-Scale)
  cpm_mat <- matrix(, nrow=4, ncol=5)
  for (i in 1:4) {
    for (j in 1:5) {
      k <- j+(i-1)*10
      cpm_mat[i,j] <- length(aptest1[[k]]@exemplars)
    }
  }
  rownames(cpm_mat) <- c("Spearman", "Pearson", "Manhattan", "Euclidean")
  colnames(cpm_mat) <- c("Min", "25th", "Median", "75th", "Max")
  return(cpm_mat)
}

#' Takes an apresult returned by APTest and makes a data table that summarizes
#' the number of exemplars produced by each distance metric, using each quantile
#' of the data, on the SCT matrix returned from Seurat::SCTransform.  
#' All similarity metrics must be used in this APTest, and there must
#' be a NON-PCA analysis (skipnonpca = FALSE).  
#'
#' @param apresult an AP clustering result returned by APTest()
#' @return a table that summarizes the NON-PCA SCT experiment
#' @export
SummarizeAPTest_NONPCA_SCT <- function(apresult) {
  apresult <- apresult[["NON-PCA"]]
  # SCT Scaled Data
  sct_mat <- matrix(, nrow=4, ncol=5)
  for (i in 1:4) {
    for (j in 1:5) {
      k <- j+(i-1)*10+5
      sct_mat[i,j] <- length(aptest1[[k]]@exemplars)
    }
  }
  rownames(sct_mat) <- c("Spearman", "Pearson", "Manhattan", "Euclidean")
  colnames(sct_mat) <- c("Min", "25th", "Median", "75th", "Max")
  return(sct_mat)
}



#' Helper function to run APCluster on a given similarity matrix with diagonals
#' set to different quantiles of the similarity matrix.  Diagonals are not counted
#' when calculating quantiles.  Quantiles used include 0th (min), 25th, 50th (med),
#' 75th, 100th(max) percentiles  
#' Note that on a 5000x5000 similarity matrix, the total size of the output list
#' is typically around 2.5 Megabytes
#' 
#' @param sim.mat A similarity matrix; diagonals will be ignored
#' @param name The name of the data whose quantiles are being clustered. Required.
#' @return A named list of APCluster results 
#' @export
ClusterSimQuantiles <- function(sim.mat, name) {
  min <- q2 <- median <- q3 <- max <- sim.mat
  min <- apcluster(min, q=0)
  q2 <- apcluster(q2, q=0.25)
  median <- apcluster(median, q=0.5)
  q3 <- apcluster(q3, q=0.75)
  max <- apcluster(max, q=1)
  list <- list(min, q2, median, q3, max)
  names(list) <- c(paste(name, "_Min",sep = ""), paste(name, "_25thQ",sep = ""), 
                   paste(name, "_Median",sep = ""), paste(name, "_75thQ",sep = ""), 
                   paste(name, "_Max",sep = ""))
  return(list)
}


#' Helper function to run APCluster on a given similarity matrix with diagonals
#' set to different quantiles of the similarity matrix.  Diagonals are not counted
#' when calculating quantiles.  Quantiles used include 0th (min), 25th, 50th (med),
#' 75th, 100th(max) percentiles  
#' Note that on a 5000x5000 similarity matrix, the total size of the output list
#' is typically around 2.5 Megabytes
#' 
#' @param sim.mat A similarity matrix; diagonals will be ignored
#' @param name The name of the data whose quantiles are being clustered. Required.
#' @return A named list of APCluster results 
#' @export
ClusterSimQuantiles <- function(sim.mat, name) {
  min <- q2 <- median <- q3 <- max <- sim.mat
  min <- apcluster(min, q=0)
  q2 <- apcluster(q2, q=0.25)
  median <- apcluster(median, q=0.5)
  q3 <- apcluster(q3, q=0.75)
  max <- apcluster(max, q=1)
  list <- list(min, q2, median, q3, max)
  names(list) <- c(paste(name, "_Min",sep = ""), paste(name, "_25thQ",sep = ""), 
                   paste(name, "_Median",sep = ""), paste(name, "_75thQ",sep = ""), 
                   paste(name, "_Max",sep = ""))
  return(list)
}

#' Finds the APResult in a given Seurat object and returns a cluster vector
#' recognizable by PISCES.  The Seurat object must have been transformed using
#' the apcluster library.
#' 
#' @param apresult.specific A specific APCluster result, not the list returned by APTest
#' @return A named int vector with sample names and respective cluster #
#' @export
APResultToVec <- function(apresult.specific) {
  affinity.cluster.vect <- c()
  for (i in 1:length(apresult.specific@clusters)) {
    clust.vect <- rep(i, each=length(apresult.specific@clusters[[i]]))
    names(clust.vect) <- names(apresult.specific@clusters[[i]])
    affinity.cluster.vect <- append(affinity.cluster.vect, clust.vect)
  }
  return(affinity.cluster.vect)
}

#' A function for generating a 2D plot of the variance explained by each PC of a
#' Seurat 'data.obj' transformed with RunPCA.  Will plot PCs on X-axis and variance
#' explained on Y-axis.  Function does not return a result, only prints.Graph 
#' should be saved by user with ggsave().
#' 
#' @param data.obj a Seurat object transformed with PCA
#' @export
graphPCvariance <- function(data.obj) {
  pca.var <- data.obj@reductions$pca@stdev **2
  # pct var explained by pcs. inflect pts appear 9, 16, where best numpcs should be
  pca.var <- pca.var/sum(pca.var)
  pca.var.df <- data.frame(pctvar = pca.var)
  pca.var.scatter <- ggplot(pca.var.df, aes(x=1:50, y=pctvar)) + geom_line() + labs(title = "PCA Variance (lung)",
                                               y = 'Percent', x='PC') +
    theme(plot.title = element_text(hjust = 0.5, colour = "black", size = 25), 
          axis.title = element_text(hjust = 0.5, colour = "black", size = 16), 
          axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12), 
          axis.text.y = element_text(hjust = 1, colour = "black", size = 12), 
          legend.title = element_text(hjust = 0, colour = "black", size = 16), 
          legend.text = element_text(hjust = 0, colour = "black", size = 12), legend.position = "right")
  plot(pca.var.scatter)
}

#' A function for testing a list of cluster vectors on 5 different metrics.
#' The 5 metrics used for comparison are the Variational Information (vi.dist),
#' the Rand Index (RI), the Normalized Information Distance (NID), the Adjusted
#' Rand Index (ARI), and the Modified Adjusted Rand Index (MARI).  All of these
#' are metrics defined in Chiquet et al. 2020, and are various computations of
#' cluster comparisons.  The function returns a table with the metrics for all
#' vectors compared to the first vector, so the first vector should be the
#' assignments vector if the computation is intended to compare to true labels.
#' 
#' @param clustvect.list A named list of cluster vectors, all of same cardinality
#' @return A table, with a row for each vector in the list, and a column for each comparison metric
#' @export
compareClusters <- function(clustvect.list) {
  assignments <- clustvect.list[[1]]
  results.table <- matrix(0, nrow = length(clustvect.list), ncol = 5)
  rownames(results.table) <- names(clustvect.list)
  colnames(results.table) <- c("VI", "RI", "NID", "ARI", "MARI")
  sort(assignments)
  for (i in 1:length(clustvect.list)) {
    vect <- compare[[i]]
    sort(vect)
    results.table[i,1] <- vi.dist(assignments, vect)
    results.table[i,2] <- RI(assignments, vect)
    results.table[i,3] <- NID(assignments, vect)
    results.table[i,4] <- ARI(assignments, vect)
    results.table[i,5] <- MARI(assignments, vect)
  }
  results.table <- signif(results.table, digits=3)
  return(results.table)
}

#' A wrapper function for ggplot2's heatmap and stats::heatmap, with an option
#' for hierarchical clustering.  This function turns a table of numerical values
#' into a heatmap, with color intensity based on the value in each cell.
#' 
#' @param table A table with numerical values in each cell
#' @param dodge Allow multiple overlaps for the labels
#' @param hcluster Perform default hierarchical clustering on heatmap
#' @export
plotHeatmap <- function(table, dodge = 1, hcluster = TRUE) {
  if (hcluster) {
    stats::heatmap(t(as.matrix(table)), keep.dendro = FALSE)
  } else {
    m <- data.frame(table)
    ggplot(m, aes(x = Var1, y = Var2, fill = Freq)) + 
      scale_x_discrete(guide = guide_axis(n.dodge = dodge))+
      geom_tile() + 
      scale_fill_distiller(palette = "Reds", direction = 1)
  }
}

#' Takes 2 cluster vectors, which should be named, and makes a matrix with all
#' incident datapoints.  For example, if data point i is incident on cluster
#' 2 of vector 1 and cluster 5 of vector 2, it will add one count to the cell
#' (2,5).  
#' 
#' @param rowVect the cluster vector to index the rows
#' @param colVect the cluster vector to index the cols
#' @param rowNames a prefix for each rowVect cluster name (1 -> CLUSTER.1)
#' @param colNames a prefix for each colVect cluster name (T-cell -> CLUSTER.T-cell)
#' @return a matrix with all datapoints incident on cluster assignments from rowVect & colVect
#' @export
clusterMatrix <- function(rowVect, colVect, rowNames = 'CLUSTER', colNames = 'CLUSTER') {
  row.vs.col.table <- table(rowVect, colVect)
  rownames(row.vs.col.table) <- paste(rowNames, rownames(row.vs.col.table), sep = '.')
  colnames(row.vs.col.table) <- paste(colNames, colnames(row.vs.col.table), sep = '.')
  return(row.vs.col.table)
}