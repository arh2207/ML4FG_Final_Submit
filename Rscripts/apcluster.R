#################
### Install all required packages
#################
install.packages(c("cluster", "ggplot2", "devtools", "Seurat", 
                   "pheatmap", "BiocManager", "RColorBrewer", "stringr",
                   "SingleR", "apcluster", "ClusterR", "mcclust", "philentropy",
                   "Matrix", "reshape2", "scales"))
BiocManager::install("viper")
BiocManager::install("biomaRt")
devtools::install_github("JEFworks/MUDAN")
devtools::install_github(repo = "califano-lab/PISCES", force = TRUE, build_vignettes = TRUE)
devtools::install_github("https://github.com/jchiquet/aricode")

#################
### Load packages, set working directory
#################
set.seed(343)
#setwd('C:/Users/Charlie/Documents/COMS4762_Project')
#setwd('/Users/andrew/Library/CloudStorage/OneDrive-Personal/Columbia/2022-1 Spring/COMSE4762/APCluster_ML4FG')
setwd('SET YOUR WORKING DIRECTORY TO THE BASE GITHUB DIRECTORY')
library(data.table)
library(stringr)
library(Seurat)
library(PISCES)
library(SingleR)
library(apcluster)
library(ClusterR)
library(mcclust)
library(philentropy)
library(Matrix)
library(aricode)
library(scales)

source('Rscripts/pipeline.R')

#################
### Preparing PBMC Data >> THIS IS ONE-OFF! >> DO NOT RUN
#################
#matrix_dir = "/Users/andrew/Library/CloudStorage/OneDrive-Personal/Columbia/2022-1 Spring/COMSE4762/APCluster_ML4FG/local/pbmc/GSE112845_RAW/"
#barcode.path <- paste(matrix_dir, "barcodes.tsv.gz", sep = "")
#features.path <- paste(matrix_dir, "features.tsv.gz", sep = "")
#matrix.path <- paste(matrix_dir, "matrix.mtx.gz", sep = "")
#mat <- readMM(file = matrix.path)
#feature.names = read.delim(features.path,
#                           header = FALSE,
#                           stringsAsFactors = FALSE)
#barcode.names = read.delim(barcode.path,
#                           header = FALSE,
#                           stringsAsFactors = FALSE)
#colnames(mat) = barcode.names$V1
#rownames(mat) = feature.names$V1
#counts.mat <- mat
#colnames(counts.mat) <- paste("pbmc", 1:length(colnames(counts.mat)), sep='.')
#data.obj <- CreateSeuratObject(counts = counts.mat, project = paste('hpa_pbmc', sep = ''),
#                               min.cells = 200,
#                               min.features = 3)
#export.mat <- GetAssayData(data.obj, assay = "RNA")
#export.mat <- as.matrix(export.mat)
#saveRDS(export.mat, paste('tissue_counts/pbmc_counts.rds', sep = ''))


#################
### Initialize Directories and Load Data
#################
tissue.name <- 'pbmc'
analysis.directory <- paste('tissue_analysis/', tissue.name, '/', sep = '')
figures.directory <- paste('tissue_analysis/', tissue.name, '/figures/', sep = '')
tables.directory <- paste('tissue_analysis/', tissue.name, '/tables/', sep = '')

dir.create(paste('tissue_analysis/', sep=''))
dir.create(analysis.directory)
dir.create(figures.directory)
dir.create(tables.directory)
dir.create('local/')
dir.create(paste('local/', tissue.name, '/', sep=''))

#################
### Data filtration and Pre-Processing, load HPA clusters
#################
# Loads counts matrix
counts.mat <- readRDS(paste('tissue_counts/', tissue.name, '_counts.rds', sep = ''))

hpca.se <- HumanPrimaryCellAtlasData()
meta.mat <- readRDS(file = 'hpa_meta-data.rds') 

data.obj <- CreateSeuratObject(counts = counts.mat, project = paste('hpa_', tissue.name, sep = ''))
hpa.samples <- which(meta.mat$Tissue == tissue.name)
hpa.clust <- meta.mat$Cluster[hpa.samples]
names(hpa.clust) <- paste(tissue.name, meta.mat$Cell[hpa.samples], sep = '.')
data.obj[["hpa.clusters"]] <- hpa.clust
mt.features <- intersect(mt.genes$hum.ensg, rownames(data.obj))
data.obj[["percent.mt"]] <- PercentageFeatureSet(object = data.obj, features = mt.features)

jpeg(paste(figures.directory, tissue.name, '_qc.jpg', sep = ''),
     width = 8, height = 4, units = "in", res = 2000)
qc.plot <- QCPlots(data.obj)
print(qc.plot)
dev.off()
print(qc.plot)


#*This If condition is here because every dataset should have a different
#*cutoff for nCount_RNA and percent.mt, and so we would like to preserve the
#*analysis by making a condition every time we use a new tissue type.
if (tissue.name == 'lung') {
  data.obj <- subset(data.obj, subset = nCount_RNA < 6000 & percent.mt < 8)
} else if (tissue.name == 'pbmc') {
  data.obj <- subset(data.obj, subset = nCount_RNA < 7500 & percent.mt < 10)
}

data.obj <- SCTransform(data.obj, vars.to.regress = 'percent.mt', verbose = FALSE)
hugo.mat <- GeneNameConvert(data.obj@assays$RNA@data, species = 'human', from = 'ensg', to = 'gn')


#################
### Run SingleR, Louvain method (LouvainKRange)
#################
# Generate SingleR Labels and Clusters
singleR.obj <- SingleR(hugo.mat, ref = hpca.se, labels = hpca.se$label.main, de.method = "wilcox")
data.obj@misc[['singleR.obj']] <- singleR.obj
singleR.flabels <- singleR.obj$first.labels; names(singleR.flabels) <- colnames(hugo.mat)
data.obj@misc[['singleR.flabels']] <- singleR.flabels

## Seurat default Louvain clustering
data.obj <- RunPCA(data.obj, verbose = FALSE)
data.obj <- RunUMAP(data.obj, dims = 1:30, verbose = FALSE)
data.obj <- FindNeighbors(data.obj, dims = 1:30, verbose = FALSE)
data.obj <- FindClusters(data.obj, verbose = FALSE)

# Modified Louvain clustering (PISCES pipeline; LouvainKRange vs. Louvain)
data.obj <- CorDist(data.obj, pca.feats = 10)
data.obj <- LouvainKRange(data.obj, kmax = 100)
saveRDS(data.obj, file = paste('local/', tissue.name, '/', tissue.name, '_pre-pisces.rds', sep = ''))

#################
### Running APTest twice.  Once on original SCT and CPM matrices. 
### Then, do PC iteration.  For the PC-iteration step, we graph those results.
#################
# APTest was only run on the Drop-seq data because of long runtime.  Results
# are saved already.  This analysis can easily be done on any other datatype.
if (tissue.name == 'lung') {
  data.obj <- readRDS(file = paste('local/', tissue.name, '/', tissue.name, '_pre-pisces.rds', sep = ''))
  ##  THESE ARE COMMENTED OUT DUE TO LONG RUNTIME.  TO RUN THE APTEST FUNCTION,
  ##  UNCOMMENT BELOW.
  #print("Generating aptest using all distance metrics, no PC iteration...")
  #aptest.apvect <- APTest(data.obj, spearman = TRUE, pearson = TRUE, negmanhattan = TRUE, negeuclidean = TRUE,
  #                         skipnonpca = FALSE, PCA=FALSE, iterate=FALSE)
  #saveRDS(aptest.apvect, file = paste('tissue_analysis/', tissue.name, '/', tissue.name, '_aptest.rds', sep = ''))
  #
  #print("Doing APTest PC Iteration using negative manhattan distance...")
  #manhattan14pcit.apvect <- APTest(data.obj, iterate=FALSE, iteratemethod = 'manhattan', npcs = 15, step=1)
  #saveRDS(manhattan14pcit.apvect, file = paste('tissue_analysis/', tissue.name, '/', tissue.name, '_manhattan14pcit.rds', sep = ''))
  #
  #print("Doing APTest PC Iteration using negative euclidean distance...")
  #euclidean14pcit.apvect <- APTest(data.obj, iterate=FALSE, iteratemethod = 'euclidean', npcs = 15, step=1)
  #saveRDS(euclidean14pcit.apvect, file = paste('tissue_analysis/', tissue.name, '/', tissue.name, '_euclidean14pcit.rds', sep = ''))
  #
  #print("Doing APTest PC Iteration using pearson correlation...")
  #pearson14pcit.apvect <- APTest(data.obj, iterate=FALSE, iteratemethod = 'pearson', npcs = 15, step=1)
  #saveRDS(pearson14pcit.apvect, file = paste('tissue_analysis/', tissue.name, '/', tissue.name, '_pearson14pcit.rds', sep = ''))
  #
  #print("Doing APTest PC Iteration using spearman correlation...")
  #spearman14pcit.apvect <- APTest(data.obj, iterate=FALSE, iteratemethod = 'spearman', npcs = 15, step=1)
  #saveRDS(spearman14pcit.apvect, file = paste('tissue_analysis/', tissue.name, '/', tissue.name, '_spearman14pcit.rds', sep = ''))
  
  manhattan14pcit <- readRDS(file = paste('tissue_analysis/', tissue.name, '/', tissue.name, '_manhattan14pcit.rds', sep = ''))
  euclidean14pcit <- readRDS(file = paste('tissue_analysis/', tissue.name, '/', tissue.name, '_euclidean14pcit.rds', sep = ''))
  pearson14pcit <- readRDS(file = paste('tissue_analysis/', tissue.name, '/', tissue.name, '_pearson14pcit.rds', sep = ''))
  spearman14pcit <- readRDS(file = paste('tissue_analysis/', tissue.name, '/', tissue.name, '_spearman14pcit.rds', sep = ''))
  
  numclust_l1_min <- c()
  numclust_l2_min <- c()
  pearson14pcit_min <- c()
  spearman14pcit_min <- c()
  
  for (i in 1:14) {
    numclust_l1_min <- append(numclust_l1_min, length(manhattan14pcit[[2]][[i]][[1]]@exemplars))
    numclust_l2_min <- append(numclust_l2_min, length(euclidean14pcit[[2]][[i]][[1]]@exemplars))
    pearson14pcit_min <- append(pearson14pcit_min, length(pearson14pcit[[2]][[i]][[1]]@exemplars))
    spearman14pcit_min <- append(spearman14pcit_min, length(spearman14pcit[[2]][[i]][[1]]@exemplars))
  }
  colors <- c("L1" = "darkred", "L2" = "blueviolet", "r" = "aquamarine3", "s" = "darkorange")
  pca.minit.df <- data.frame(l1 = numclust_l1_min, l2 = numclust_l2_min, p = pearson14pcit_min, s = spearman14pcit_min)
  pca.minit.scatter <- ggplot(pca.minit.df, aes(x=2:15)) + 
    geom_smooth(aes(y = l1, color = "L1"), se=F) + 
    geom_smooth(aes(y = l2, color="L2"), se=F) + 
    geom_smooth(aes(y = p, color = "r"), se=F) +
    geom_smooth(aes(y = s, color = "s"), se=F) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    labs(title = "Using Diagonal = Min", y = 'Number of Clusters', x='Principal Components', color = "Legend") +
    theme(plot.title = element_text(hjust = 0.5, colour = "black", size = 25), 
          axis.title = element_text(hjust = 0.5, colour = "black", size = 16), 
          axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12), 
          axis.text.y = element_text(hjust = 1, colour = "black", size = 12), 
          legend.title = element_text(hjust = 0, colour = "black", size = 16), 
          legend.text = element_text(hjust = 0, colour = "black", size = 12)) +
    scale_color_manual(values = colors)
  print(pca.minit.scatter)
  ggsave(filename = paste(figures.directory, 'pca_minit.jpg', sep=''), dpi = 500, height = 7, width = 10)
  
  numclust_l1_25th <- c()
  numclust_l2_25th <- c()
  pearson14pcit_25th <- c()
  spearman14pcit_25th <- c()
  for (i in 1:14) {
    numclust_l1_25th <- append(numclust_l1_25th, length(manhattan14pcit[[2]][[i]][[2]]@exemplars))
    numclust_l2_25th <- append(numclust_l2_25th, length(euclidean14pcit[[2]][[i]][[2]]@exemplars))
    pearson14pcit_25th <- append(pearson14pcit_25th, length(pearson14pcit[[2]][[i]][[2]]@exemplars))
    spearman14pcit_25th <- append(spearman14pcit_25th, length(spearman14pcit[[2]][[i]][[2]]@exemplars))
  }
  colors <- c("L1" = "darkred", "L2" = "blueviolet", "r" = "aquamarine3", "s" = "darkorange")
  pca.25thit.df <- data.frame(l1 = numclust_l1_25th, l2 = numclust_l2_25th, p = pearson14pcit_25th, s = spearman14pcit_25th)
  pca.25thit.scatter <- ggplot(pca.25thit.df, aes(x=2:15)) + 
    geom_smooth(aes(y = l1, color = "L1"), se=F) + 
    geom_smooth(aes(y = l2, color="L2"), se=F) + 
    geom_smooth(aes(y = p, color = "r"), se=F) +
    geom_smooth(aes(y = s, color = "s"), se=F) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    labs(title = "Using Diagonal = 25th", y = 'Number of Clusters', x='Principal Components', color = "Legend") +
    theme(plot.title = element_text(hjust = 0.5, colour = "black", size = 25), 
          axis.title = element_text(hjust = 0.5, colour = "black", size = 16), 
          axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12), 
          axis.text.y = element_text(hjust = 1, colour = "black", size = 12), 
          legend.title = element_text(hjust = 0, colour = "black", size = 16), 
          legend.text = element_text(hjust = 0, colour = "black", size = 12)) +
    scale_color_manual(values = colors)
  print(pca.25thit.scatter)
  ggsave(filename = paste(figures.directory, 'pca_25thit.jpg'), dpi = 500, height = 7, width = 10)
  
  numclust_l1_med <- c()
  numclust_l2_med <- c()
  pearson14pcit_med <- c()
  spearman14pcit_med <- c()
  for (i in 1:14) {
    numclust_l1_med <- append(numclust_l1_med, length(manhattan14pcit[[2]][[i]][[3]]@exemplars))
    numclust_l2_med <- append(numclust_l2_med, length(euclidean14pcit[[2]][[i]][[3]]@exemplars))
    pearson14pcit_med <- append(pearson14pcit_med, length(pearson14pcit[[2]][[i]][[3]]@exemplars))
    spearman14pcit_med <- append(spearman14pcit_med, length(spearman14pcit[[2]][[i]][[3]]@exemplars))
  }
  colors <- c("L1" = "darkred", "L2" = "blueviolet", "r" = "aquamarine3", "s" = "darkorange")
  pca.medit.df <- data.frame(l1 = numclust_l1_med, l2 = numclust_l2_med, p = pearson14pcit_med, s = spearman14pcit_med)
  pca.medit.scatter <- ggplot(pca.medit.df, aes(x=2:15)) + 
    geom_smooth(aes(y = l1, color = "L1"), se=F) + 
    geom_smooth(aes(y = l2, color="L2"), se=F) + 
    geom_smooth(aes(y = p, color = "r"), se=F) +
    geom_smooth(aes(y = s, color = "s"), se=F) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    labs(title = "Using Diagonal = med", y = 'Number of Clusters', x='Principal Components', color = "Legend") +
    theme(plot.title = element_text(hjust = 0.5, colour = "black", size = 25), 
          axis.title = element_text(hjust = 0.5, colour = "black", size = 16), 
          axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12), 
          axis.text.y = element_text(hjust = 1, colour = "black", size = 12), 
          legend.title = element_text(hjust = 0, colour = "black", size = 16), 
          legend.text = element_text(hjust = 0, colour = "black", size = 12)) +
    scale_color_manual(values = colors)
  print(pca.medit.scatter)
  ggsave(filename = paste(figures.directory, 'pca_medit.jpg'), dpi = 500, height = 7, width = 10)
  
  numclust_l1_75th <- c()
  numclust_l2_75th <- c()
  pearson14pcit_75th <- c()
  spearman14pcit_75th <- c()
  for (i in 1:14) {
    numclust_l1_75th <- append(numclust_l1_75th, length(manhattan14pcit[[2]][[i]][[4]]@exemplars))
    numclust_l2_75th <- append(numclust_l2_75th, length(euclidean14pcit[[2]][[i]][[4]]@exemplars))
    pearson14pcit_75th <- append(pearson14pcit_75th, length(pearson14pcit[[2]][[i]][[4]]@exemplars))
    spearman14pcit_75th <- append(spearman14pcit_75th, length(spearman14pcit[[2]][[i]][[4]]@exemplars))
  }
  colors <- c("L1" = "darkred", "L2" = "blueviolet", "r" = "aquamarine3", "s" = "darkorange")
  pca.75thit.df <- data.frame(l1 = numclust_l1_75th, l2 = numclust_l2_75th, p = pearson14pcit_75th, s = spearman14pcit_75th)
  pca.75thit.scatter <- ggplot(pca.75thit.df, aes(x=2:15)) + 
    geom_smooth(aes(y = l1, color = "L1"), se=F) + 
    geom_smooth(aes(y = l2, color="L2"), se=F) + 
    geom_smooth(aes(y = p, color = "r"), se=F) +
    geom_smooth(aes(y = s, color = "s"), se=F) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    labs(title = "Using Diagonal = 75th", y = 'Number of Clusters', x='Principal Components', color = "Legend") +
    theme(plot.title = element_text(hjust = 0.5, colour = "black", size = 25), 
          axis.title = element_text(hjust = 0.5, colour = "black", size = 16), 
          axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12), 
          axis.text.y = element_text(hjust = 1, colour = "black", size = 12), 
          legend.title = element_text(hjust = 0, colour = "black", size = 16), 
          legend.text = element_text(hjust = 0, colour = "black", size = 12)) +
    scale_color_manual(values = colors)
  print(pca.75thit.scatter)
  ggsave(filename = paste(figures.directory, 'pca_75thit.jpg', sep=''), dpi = 500, height = 7, width = 10)
  
  numclust_l1_Max <- c()
  numclust_l2_Max <- c()
  pearson14pcit_Max <- c()
  spearman14pcit_Max <- c()
  for (i in 1:14) {
    numclust_l1_Max <- append(numclust_l1_Max, length(manhattan14pcit[[2]][[i]][[5]]@exemplars))
    numclust_l2_Max <- append(numclust_l2_Max, length(euclidean14pcit[[2]][[i]][[5]]@exemplars))
    pearson14pcit_Max <- append(pearson14pcit_Max, length(pearson14pcit[[2]][[i]][[5]]@exemplars))
    spearman14pcit_Max <- append(spearman14pcit_Max, length(spearman14pcit[[2]][[i]][[5]]@exemplars))
  }
  colors <- c("L1" = "darkred", "L2" = "blueviolet", "r" = "aquamarine3", "s" = "darkorange")
  pca.Maxit.df <- data.frame(l1 = numclust_l1_Max, l2 = numclust_l2_Max, p = pearson14pcit_Max, s = spearman14pcit_Max)
  pca.Maxit.scatter <- ggplot(pca.Maxit.df, aes(x=2:15)) + 
    geom_smooth(aes(y = l1, color = "L1"), se=F) + 
    geom_smooth(aes(y = l2, color="L2"), se=F) + 
    geom_smooth(aes(y = p, color = "r"), se=F) +
    geom_smooth(aes(y = s, color = "s"), se=F) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    labs(title = "Using Diagonal = Max", y = 'Number of Clusters', x='Principal Components', color = "Legend") +
    theme(plot.title = element_text(hjust = 0.5, colour = "black", size = 25), 
          axis.title = element_text(hjust = 0.5, colour = "black", size = 16), 
          axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12), 
          axis.text.y = element_text(hjust = 1, colour = "black", size = 12), 
          legend.title = element_text(hjust = 0, colour = "black", size = 16), 
          legend.text = element_text(hjust = 0, colour = "black", size = 12)) +
    scale_color_manual(values = colors)
  print(pca.Maxit.scatter)
  ggsave(filename = paste(figures.directory, 'pca_Maxit.jpg', sep=''), dpi = 500, height = 7, width = 10)
}

#################
### Generating APCluster cluster vector (Pipeline or non-pipeline)
#################
data.obj <- readRDS(file = paste('local/', tissue.name, '/', tissue.name, '_pre-pisces.rds', sep = ''))

if (tissue.name == 'lung') {
  ## We only ran aptest on lung!
  ## aptest is a bigresult returned from function APTest(), where we choose the best cluster
  aptest <- readRDS(file = paste('tissue_analysis/', tissue.name, '/', tissue.name, '_aptest.rds', sep = ''))
  ## Chosen cluster is L1 Min on CPM because of number of cluster centers.  This analysis
  ## is done manually from the figures.
  affinity.cluster <- aptest[["cpm.manhattan_Min"]]
  saveRDS(affinity.cluster, file=paste(analysis.directory, tissue.name, "_al1min-cpm.rds", sep=''))
} else {
  ## In the else case ('pbmc'), we do apclustering on the best generated vector
  # from the APTest run in the lung data
  pca.data <- Embeddings(object = data.obj, reduction = "pca")[,1:15]
  sct.sim <- -distance(pca.data, method="manhattan", est.prob="empirical", use.row.names = TRUE)
  affinity.cluster <- apcluster(sct.sim, q=0.01)
  saveRDS(affinity.cluster, file=paste(analysis.directory, tissue.name, "_al1min-cpm.rds", sep=''))
}

# We want to transform this into a named vector and store it in the object
affinity.cluster <- readRDS(file=paste(analysis.directory, tissue.name, "_al1min-cpm.rds", sep=''))
affinity.cluster.vect <- APResultToVec(affinity.cluster)
data.obj@assays[[data.obj@active.assay]]@misc[['affinity.cluster.vect']] <- affinity.cluster.vect
saveRDS(data.obj, file = paste('local/', tissue.name, '/', tissue.name, '_pre-pisces.rds', sep = ''))

#################
### Analyzing NON-PCA, PCA, and PCA-ITERATE results from APTest
#################
## We only ran aptest on lung!
if (tissue.name == 'lung') {
  cpm_mat <- SummarizeAPTest_NONPCA_CPM(aptest)
  sct_mat <- SummarizeAPTest_NONPCA_SCT(aptest)
  write.table(cpm_mat, sep = ',', quote = FALSE, col.names=NA, 
              file = paste(tables.directory, "cpm_table.csv", sep=''))
  write.table(sct_mat, sep = ',', quote = FALSE, col.names=NA, 
              file = paste(tables.directory, "sct_table.csv", sep=''))
}

#################
### Comparing Best AP Clustering result -L1 (Min) to hpa, louvain, singleR
#################
data.obj <- readRDS(file = paste('local/', tissue.name, '/', tissue.name, '_pre-pisces.rds', sep = ''))

# Now we compare to the HPA labels
hpa.affinity.table <- clusterMatrix(data.obj$hpa.clusters[names(data.obj@assays$SCT@misc$affinity.cluster.vect)],
                                        data.obj@assays$SCT@misc$affinity.cluster.vect,
                                        rowNames = 'HPA',
                                        colNames = 'AFFINITY')
write.table(hpa.affinity.table, sep = ',', quote = FALSE, col.names=NA,
            file = paste(tables.directory, 'hpa-al1min.csv', sep=''))
jpeg(file = paste(figures.directory, 'hpa-al1min.jpg', sep=''),
     width = 1300,
     height = 700,
     units = "px",
     quality = 200)
plotHeatmap(hpa.affinity.table, hcluster = TRUE)
dev.off()
# Comparing to Louvain method
louvain.affinity.table <- clusterMatrix(data.obj@assays$SCT@misc$pisces.cluster,
                                  data.obj@assays$SCT@misc$affinity.cluster.vect,
                                  rowNames = 'PISCES',
                                  colNames = 'AFFINITY')
write.table(louvain.affinity.table, sep = ',', quote = FALSE, col.names=NA,
            file = paste(tables.directory, 'louvain-al1min.csv', sep=''))
jpeg(file = paste(figures.directory, 'louvain-al1min.jpg', sep=''),
     width = 1300,
     height = 700,
     units = "px",
     quality = 200)
plotHeatmap(louvain.affinity.table, hcluster = TRUE)
dev.off()
# Comparing to SingleR
singler.affinity.table <- clusterMatrix(data.obj@assays$SCT@misc$affinity.cluster.vect[names(data.obj@misc$singleR.flabels)],
                                   data.obj@misc$singleR.flabels,
                                   colNames = 'sR',
                                   rowNames = 'AFFINITY')
write.table(singler.affinity.table, sep = ',', quote = FALSE, col.names=NA,
            file = paste(tables.directory, 'singler-al1min.csv', sep=''))
jpeg(file = paste(figures.directory, 'singler-al1min.jpg', sep=''),
     width = 1300,
     height = 700,
     units = "px",
     quality = 100)
plotHeatmap(singler.affinity.table, hcluster = TRUE)
dev.off()

#################
### Performing K-means on Single Cell data
#################
## perform k means on single cell data
pre.sct.data1 <- GetAssayData(object = data.obj, slot = "data", assay = 'RNA')
pre.sct.data <- GetAssayData(object = data.obj, slot = "data", assay = 'SCT')
post.sct.data <- GetAssayData(object = data.obj, slot = "scale.data", assay = 'SCT')
pca.data <- Embeddings(object = data.obj, reduction = "pca")[,1:15]

pca.kmeans <- kmeans(pca.data, centers = 3, nstart = 20)
kmeans.clust.vect <- pca.kmeans$cluster
saveRDS(kmeans.clust.vect, paste("tissue_analysis/", tissue.name, "/", tissue.name, "_kmeans.rds", sep = ''))

#post.sct.kmeans <- kmeans(post.sct.data, centers = 3, nstart = 20)
#saveRDS(data.obj, file = paste('local/', tissue.name, '/', tissue.name, '_clustered.rds', sep = ''))

#################
### Making APCluster vectors on KL-divergence, cosine similarity
#################
data.obj <- readRDS(file = paste('local/', tissue.name, '/', tissue.name, '_pre-pisces.rds', sep = ''))
pca.data <- Embeddings(object = data.obj, reduction = "pca")[,1:15]
pca.data.kl.sim <- -distance(pca.data, method="kullback-leibler", est.prob="empirical", use.row.names = TRUE)
pca.data.kl.sim[is.nan(pca.data.kl.sim)] <- -100
pca.data.cos.sim <- distance(pca.data, method="cosine", est.prob="empirical", use.row.names = TRUE)
apkl.clust.vect <- APResultToVec(apcluster(pca.data.kl.sim, q=0))
saveRDS(apkl.clust.vect, file=paste(analysis.directory, tissue.name, '_ap-KL.rds', sep=''))
apcos.clust.vect <- APResultToVec(apcluster(pca.data.cos.sim, q=0))
saveRDS(apcos.clust.vect, file=paste(analysis.directory, tissue.name, '_ap-cos.rds', sep=''))

#################
### Comparing all clusters to the original hpa clusters
#################
data.obj <- readRDS(file = paste('local/', tissue.name, '/', tissue.name, '_pre-pisces.rds', sep = ''))
hpa.clust.vect <- data.obj@meta.data$hpa.clusters
names(hpa.clust.vect) <- colnames(data.obj)
# Some values in hpa.clust.vect are NA
hpa.clust.vect <- replace(hpa.clust.vect, is.na(hpa.clust.vect), 0) 
louv.clust.vect <- data.obj@assays$SCT@misc$pisces.cluster
kmeans.clust.vect <- readRDS(paste("tissue_analysis/", tissue.name, "/", tissue.name, "_kmeans.rds", sep = ''))
apbest.clust.vect <- data.obj@assays$SCT@misc$affinity.cluster.vect
apkl.clust.vect <- readRDS(file=paste(analysis.directory, tissue.name, '_ap-KL.rds', sep=''))
apcos.clust.vect <- readRDS(file=paste(analysis.directory, tissue.name, '_ap-cos.rds', sep=''))
# Using variation of information (lower is better)
assignments <- hpa.clust.vect
compare <- list(assignments, louv.clust.vect, kmeans.clust.vect,
                apbest.clust.vect, apkl.clust.vect, apcos.clust.vect)
names(compare) <- c("true", "louv", "kmeans", "apbest", "apkl", "apcos")
cluster.comparison.table <- compareClusters(compare)
write.table(cluster.comparison.table, sep = ',', quote = FALSE, col.names=NA,
            file = paste(tables.directory, 'cluster_comparison.csv', sep=''))

