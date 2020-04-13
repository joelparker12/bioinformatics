library(SingleCellExpirement)
BiocManager::install(c('scater', 'scran', 'uwot'))
BiocManager::install('SingleCellExperiment')
library(SingleCellExperiment)


counts_matrix <- data.frame(cell_1 = rpois(10,10),
                            cell_2= rpois(10,10),
                            cell_3 = rpois(10,30))
rownames(counts_matrix) <- paste0("gene_",1:10)

counts_matrix <- as.matrix(counts_matrix) ### must be a matrix

#### Now we will make a single cell expirement. 
sce <- SingleCellExperiment(assays = list(counts=counts_matrix))


#### To view our counts 
assay(sce,"counts")

#### Now we will add another assay looking at the log counts
sce <- scater::logNormCounts(sce)
sce

logcounts(sce) ## viewing logcounts


assays(sce) ## looking at the available assays

counts_100 <- counts(sce)+100 ## adds 100 to our count data
assay(sce,"counts_100") <- counts_100 ## creates new assay


#### Looking at metadata, cells 1 and 2 in batch 1, and cell 3 in batch 2
cell_metadata <- data.frame(batch = c(1,1,2))
rownames(cell_metadata) <- paste0("cell_", 1:3) 

sce <- SingleCellExperiment(assays= list(counts=counts_matrix), colData =cell_metadata)
sce
colData(sce) # looking at coldata


sce <- scater::addPerCellQC(sce)
colData(sce)[, 1:5] #### adds percell quality score. 

##### TO add more meta data
sce$more_stuff <- runif(ncol(sce))
colnames(colData(sce))

#### filtering expirement
sce[,sce$batch==1]


#### Rowdata slots contains contains genes and stores information on the gene.
rowRanges(sce) ## its empty because we havent added anything. 

#### we are going to add per feature QC for our row data
sce <- scater::addPerFeatureQC(sce)
rowData(sce)


library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]] # Human, Ensembl v97.


#### We can look at specific genes like regular. 
sce[c("gene_1","gene_4"), ]


##### We can store anything we want in our metadata slot. 
my_genes <- c("gene_1", "gene_2")
metadata(sce) <- list(fav_genes = my_genes)
metadata(sce)

#### to append genes 
your_genes <- c('gene_4', 'gene_8')
metadata(sce)$your_genes <- your_genes
metadata(sce)


##### Dimentionality reduction
sce <- scater::logNormCounts(sce)
sce <- scater::runPCA(sce)
reducedDim(sce,"PCA") ## view data for our PCA


#### We can also run TSNE
sce <- scater::runTSNE(sce,perplexity = 0.1)
reducedDim(sce,"TSNE")
reducedDim(sce)


## generating your own dim reduction
u <- uwot::umap(t(logcounts(sce)), n_neighbors = 2)
reducedDim(sce, "UMAP_uwot") <- u
reducedDims(sce) #

#### Alternative Expirenments for "good for spike counts"
spike_counts <- cbind(cell_1 <- rpois(5,10),
                      cell_2 <- rpois(5,10),
                      cell_3 <- rpois(5,30))
rownames(spike_counts) <- paste0("spike_",1:5)
spike_se <- SummarizedExperiment(list(counts= spike_counts))
spike_se


#### now store this new summarized expirement in old expiremement. 
altExp(sce, "spike") <- spike_se
altExp(sce)




###Size factors allows us to use or get a set a numeric value for each cell for normalization reasons. 
sce <- scran::computeSumFactors(sce)
sizeFactors(sce)
sizeFactors(sce) <- scater::librarySizeFactors(sce)
sizeFactors(sce)

####### COllabels allows us to set labels to the columns for clusering
colLables(sce) <- LETTERS[1:3]
colLabels(sce)

  










