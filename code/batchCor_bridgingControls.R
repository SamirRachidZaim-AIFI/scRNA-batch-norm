################################################################################################
################################################################################################

# Title: driver_scMerge 
# Creator: Monica Chaudhari
# Last Modified: Samir Rachid Zaim 
# Date Modified: 4/19/21

# Desc: 
#     This file is the "driver" script in scMerge normalization in which a set of datasets,
#     a bridging control, and parameters are provided in order to "merge" these datasets in 
#     order to analyze them together. 

#     Technical Details of scMerge are found here: https://www.pnas.org/content/116/20/9775
#     The bioconductor R package is found here: https://bioconductor.org/packages/release/bioc/html/scMerge.html

# Future Work:
#     The future work is to parameterize this script and modularize its components
#     in order to eventually be able to package it into the scRNA normalization pipeline.

################################################################################################
################################################################################################

# Load libraries

source('scMerge_helperFunctions.R')
library(SingleCellExperiment)
require(scMerge)
require(BiocParallel)

#############

# setwd('/home/')
# sce <- readRDS('sceFltrd_kpGns.rds')

# sceLP<-sce[,sce$subjid=="IMM19"] # dim(sceLP):  22758 x 9699; table(sceLP$batch, sceLP$cellType)
# assay(sceLP,1)<-as.matrix(assay(sceLP,1))
# assay(sceLP,2)<-as.matrix(assay(sceLP,2))
# sceLP<-sceLP[which(rowSums(assay(sceLP,2)) != 0), which(colSums(assay(sceLP,2)) != 0)] #dim(sceLP): 19524  9699
# table(sceLP$cellType,sceLP$batch); dim(sceLP)

correct_bridging_controls <- function(tmp){
    ## remove cols/rows of all zeros
    tmp<-tmp[which(rowSums(assay(tmp,2)) != 0), which(colSums(assay(tmp,2)) != 0)] 
    n.ct<-length(unique(colData(tmp)$cellType))
    n.batch<-length(unique(colData(tmp)$batch))
    
    data("segList_ensemblGeneID", package = "scMerge")
    seg_index <- segList_ensemblGeneID$human$human_scSEG
    cmSEGs <- intersect(row.names(tmp), seg_index)

    # subsample.tmp <- tmp[,tmp$barcodes %in% sample(tmp$barcodes, 8000)]
    repMat <- scMerge::scReplicate(sce_combine = tmp, 
                                     batch = tmp$batch, 
                                     kmeansK = rep(n.ct,n.batch),
                                     exprs = 'logcpm', 
                                     hvg_exprs = 'counts', 
                                     marker = NULL, 
                                     marker_list = NULL, 
                                     replicate_prop = 1, 
                                     cell_type = tmp$cell_type, 
                                     cell_type_match = FALSE, 
                                     cell_type_inc = NULL, 
                                     dist = 'cor', 
                                     WV = NULL, 
                                     WV_marker = NULL, 
                                     verbose = TRUE)

    exprs_mat <- SummarizedExperiment::assay(tmp, 'logcpm')
    ctl <- which(row.names(tmp) %in% cmSEGs)
    ruvK <- 1:5

    message("Step 2: Performing RUV normalisation. This will take minutes to hours. \n")

    ## if you subsample too small a number.. a few rows will be zeros
    ## which will lead to 0 denoms and NAs 
    ruv3res <- scMerge::scRUVIII(Y = exprs_mat, 
                                   M = repMat, 
                                   ctl = ctl, 
                                   batch=tmp$batch,                          
                                   k = ruvK,
                                   cell_type = tmp$cellType)

    k<-ruv3res$optimal_ruvK
    assay(tmp,"normalized")<-t(ruv3res[[k]]$newY) #newY_mc
    # saveRDS(object=list(sce_object =tmp, scNormalization=ruv3res), '../data/tmp_scMerge.RDS')
    return(tmp)

}