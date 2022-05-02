A Pipeline for Single-Cell Cluster Analysis

--------------------------------------------------------------------------------

The data is also available on Google Drive here: https://drive.google.com/drive/folders/1UgMCUQZf2arqn1XnzBTR0x-3Uf2AHYX7?usp=sharing

Note to the grader:
Due to the long runtime of some algorithms in the pipeline, we have also uploaded object outputs from the pipeline, as well as have commented out code where the functions with long runtimes are run.  The rest of the experiment is intended to be generated on your machine de novo.  In addition to being uploaded on Google Drive, the data are already located in the directory structure for the Github repo, so manually downloading from the Google Drive is redundant.  Important directory structures and R files:

Working directory ‘wd’ should be the base Github directory
Put any ‘\*.R’ files in wd/Rscripts/.
Put any tissue counts ‘\*.rds’ files in wd/tissue_counts/.
Keep the ‘hpa_meta-data.rds’ in wd/.

apcluster.R

This R file does all package downloading, loads packages, runs the pipeline, and generates all data figures and tables in the wd/tissue_analysis/. directory.  
*Note: to run this script on the PBMC vs lung set, simply set ‘tissue.name’ on line 65 to ‘lung’ from ‘pbmc’.  
Note: You must set the working directory with setwd() on line 20.*

pipeline.R

This contains the pipeline with various functions called in apcluster.R.  This file does not need to be touched or modified in any way, but it should be kept in wd/Rscripts/.

--------------------------------------------------------------------------------

Data From:
https://www.proteinatlas.org/about/assays+annotation

Tissue Cell Type: Using GTEx bulk RNAseq data to profile gene cell type specificity

GTEx data was used in an integrative network analysis to determine the cell type specificity of all protein coding genes within a given tissue type. For more details on this analysis and the classifications, see the Tissue Cell Type section Methods Summary. 

scRNA-seq data

--------------------------------------------------------------------------------

Inclusion criteria

The single cell RNA sequencing dataset is based on meta-analysis of literature on single cell RNA sequencing and single cell databases that include healthy human tissue. To avoid technical bias and to ensure that the single cell dataset can best represent the corresponding tissue, the following data selection criteria were applied: (1) Single cell transcriptomic datasets were limited to those based on the Chromium single cell gene expression platform from 10X Genomics (version 2 or 3); (2) Single cell RNA sequencing was performed on single cell suspension from tissues without pre-enrichment of cell types; (3) Only studies with >4,000 cells and 20 million read counts were included, (4) Only dataset whose pseudo-bulk transcriptomic expression profile is highly correlated with the transcriptomic expression profile of the corresponding HPA tissue bulk sample were included. It should be noted that exceptions were made for lung (~7.3 million reads), pancreas (3,719 cells) and rectum (3,898 cells) to include various cell types in the analysis.


--------------------------------------------------------------------------------

Single cell transcriptomics datasets

In total, single cell transcriptomics data for 25 tissues and peripheral blood mononuclear cells (PBMCs) were analyzed. These datasets were respectively retrieved from the Single Cell Expression Atlas, the Human Cell Atlas, the Gene Expression Omnibus, the Allen Brain Map, and the European Genome-phenome Archive. The complete list of references is shown in the table below. 

Tissue	Data source	No. of M reads	No. of cells	Correlation with
HPA bulk RNA	Reference
Colon	GSE116222	12.9	11167	0.811	Parikh K et al. (2019)
Eye	GSE137537	22.6	20091		Menon M et al. (2019)
Heart muscle	GSE109816	396.7	9182	0.797	Wang L et al. (2020)
Small intestine	GSE125970	59.6	6167		Wang Y et al. (2020)
Kidney	GSE131685	56.2	25279	0.867	Liao J et al. (2020)
Liver	GSE115469	42	8439	0.837	MacParland SA et al. (2018)
Lung	GSE130148	6.9	4599	0.863	Vieira Braga FA et al. (2019)
Placenta	E-MTAB-6701	347	18547	0.879	Vento-Tormo R et al. (2018)
Prostate	GSE117403	177.4	35862	0.756	Henry GH et al. (2018)
Rectum	GSE125970	60.4	3898	0.756	Wang Y et al. (2020)
PBMC	GSE112845	19.4	4972	0.756	Chen J et al. (2018)
Testis	GSE120508	70.5	6490	0.756	Guo J et al. (2018)
Pancreas	GSE131886	92.3	3719	0.829	Qadir MMF et al. (2020)
Skin	GSE130973	56.4	15798	0.756	Solé-Boldo L et al. (2020)
Brain		1357.2	76533	0.661	Allen brain map
Bronchus		85.3	17521		Lukassen S et al. (2020)
Endometrium	GSE111976	624.1	71032	0.807	Wang W et al. (2020)
Skeletal muscle	GSE143704	77.8	22030	0.697	De Micheli AJ et al. (2020)
Ovary	GSE146512	259.2	43636	0.808	Man L et al. (2020)
Adipose tissue	GSE155960	418.5	83536	0.813	Hildreth AD et al. (2021)
Esophagus	159929-GSM4850580	31.5	9117	0.840	He S et al. (2020)
Lymph node	GSE159929-GSM4850583	14.2	7771	0.849	He S et al. (2020)
Bone marrow	GSE159929-GSM4850584	8.6	3230	0.818	He S et al. (2020)
Spleen	GSE159929-GSM4850589	14.8	4512	0.804	He S et al. (2020)
Stomach	GSE159929-GSM4850590	18.9	5318	0.814	He S et al. (2020)
Breast	GSE164898	342.3	47662	0.839	Bhat-Nakshatri P et al. (2021)

--------------------------------------------------------------------------------

Clustering of single cell transcriptomics data

For each of the single cell transcriptomics datasets, the quantified raw sequencing data were downloaded from the corresponding depository database based on the accession number provided by the corresponding study in the available format (total cells, read, and feature counts, or count tables). Unfiltered data were used as input for downstream analysis with an in-house pipeline using Scanpy (version 1.4.4.post1) in Python 3.7.3 for the 13 tissues and PBMC published in HPA v20 and Scanpy (version 1.7.2) in Python 3.8.5 for the 12 tissues published in HPA v21. In the pipeline, the data were filtered using two criteria: a cell is considered as valid if at least 200 genes are detected and a gene is considered as valid if it is expressed in at least 10% of the cells. Specially, in the HPA v21, tissues which contain more than 10,000 cells used 1000 cells as their cutoff. Subsequently, the cell counts were normalized to have a total count per cell of 10000. The valid cells were then clustered using Louvain clustering function within Single-Cell Analysis in Python (Scanpy). Default values of parameters were used in clustering. More in detail, the features of cells were projected into a PCA space with 50 components using UMAP, and a k-nearest neighbours (KNN) graph was generated. 15 neighbours were used in the network for Louvain, while the resolution of clustering was set as 1.0. The total read counts for all genes in each cluster was calculated by adding up the read counts of each gene in all cells belonging to the corresponding cluster. Finally, the read counts were normalized to transcripts per million protein coding genes (pTPM) for each of the single cell clusters. When calculating the expression profile for pseudo-bulk samples based on single cell transcriptomics, we added the read counts for all genes from all cells of the sample, and normalized it to pTPM in the same way as for the cluster ones. 

--------------------------------------------------------------------------------

Defining cell types

Each of the 444 different cell type clusters were manually annotated based on an extensive survey of >500 well-known tissue and cell type-specific markers, including both markers from the original publications, and additional markers used in pathology diagnostics. For each cluster, one main cell type was chosen by taking into consideration the expression of different markers. For a few clusters, no main cell type could be selected, and these clusters were not used for gene classification. The most relevant markers are presented in a heatmap on the Cell Type Atlas, in order to clarify cluster annotation to visitors. 

Cell type dendrogram

The cell type dendrogram presented on the Single Cell Type section shows the relationship between the single cell types based on genome-wide expression. The dendrogram is based on agglomerative clustering of 1 - Spearman's rho between cell types using Ward's criterion. The dendrogram was then transformed into a hierarchical graph, and link distances were normalized to emphasize graph connections rather than link distances. Link width is proportional to the distance from the root, and links are colored according to cell type group if only one cell type group is present among connected leaves.
