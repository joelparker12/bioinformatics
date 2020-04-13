library(scRNAseq)
library(SingleCellExperiment)
#### Simple Workflow. 
#Bring in data
sce <- MacoskoRetinaData()


#############################Quality Control####################################
library(scater)

is.mito <- grepl("MT^", rownames(sce))

qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))

filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent") ## filter
sce <- sce[,!filtered$discard] ## takes out the low quality cells. 



############################### Normalization####################################
sce <- logNormCounts(sce)



##############################Feature Selection ################################
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=.1) ## Gets high value genes.


#############################Diminsionality Reduction###########################
set.seed(1234)
sce <- runPCA(sce, ncomponents=25, subset_row = hvg)
sce < runUMAP(sce, dimred = 'PCA', external_neighbors =TRUE )

############################ Clustering #########################################

g <- buildSNNGraph(sce, use.dimred = 'PCA')
colLabels(sce) = factor(igraph::cluster_louvain(g)$membership)


########################### Visualization ######################################
plotUMAP(sce, colour_by = "label")




























