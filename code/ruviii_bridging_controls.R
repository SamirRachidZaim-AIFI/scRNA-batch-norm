################################################################################################
################################################################################################

# Title: ruviii_bridging_controls 
# Creator: Monica Chaudhari
# Last Modified: Samir Rachid Zaim 
# Date Modified: 4/19/21

# Desc: 
#     This file provides the script to do an ruv-iii correction on the bridging controls. 

#     Technical Details of scMerge are found here: https://www.pnas.org/content/116/20/9775
#     The bioconductor R package is found here: https://bioconductor.org/packages/release/bioc/html/scMerge.html

# Inputs
#     tmp = sce object containing bridging controls to normalize
#     exprsname= string containing assayname to normalize
#     cellType = string containing column name with cell type information
#     batchColname= string containing column name with batch info

# Output
#     2-dimensional list
#     - normed_sce = an sce object returned with a 3rd assay containing the normalized expressions  
#     - ruv3res= an scRUVIII object containing results of scMerge normalization

################################################################################################
################################################################################################

# Load libraries

source('/home/jupyter/scRNA-batch-norm/code/scMerge_helperFunctions.R')
require('DelayedArray')
library(SingleCellExperiment)
require(scMerge)
require(BiocParallel)

ruviii_bridging_controls <- function(tmp, 
                                      exprsname='logcounts', 
                                      cellType='seurat_pbmc_type',
                                      batchColname='batch_id'
                                     ){
    
    require("org.Hs.eg.db")
    ## remove cols/rows of all zeros
    tmp<-tmp[which(rowSums(assay(tmp,2)) != 0), which(colSums(assay(tmp,2)) != 0)] 
    n.ct<-length(unique(colData(tmp)[,cellType]))
    n.batch<-length(unique(colData(tmp)[,batchColname]))
    
    data("segList_ensemblGeneID", package = "scMerge")
    seg_index <- segList_ensemblGeneID$human$human_scSEG
    symbols <- mapIds(org.Hs.eg.db, keys = seg_index, keytype = "ENSEMBL", column="SYMBOL")

    cmSEGs <- intersect(row.names(tmp), symbols)
    

    # subsample.tmp <- tmp[,tmp$barcodes %in% sample(tmp$barcodes, 8000)]
    repMat <- scMerge::scReplicate(sce_combine = tmp, 
                                     batch = colData(tmp)[,batchColname], 
                                     kmeansK = rep(n.ct,n.batch),
                                     exprs = exprsname, 
                                     hvg_exprs = 'counts', 
                                     marker = NULL, 
                                     marker_list = NULL, 
                                     replicate_prop = 1, 
                                     cell_type = colData(tmp)[,cellType], 
                                     cell_type_match = FALSE, 
                                     cell_type_inc = NULL, 
                                     dist = 'cor', 
                                     WV = NULL, 
                                     WV_marker = NULL, 
                                     verbose = TRUE)

    exprs_mat <- SummarizedExperiment::assay(tmp, exprsname)
    ctl <- which(row.names(tmp) %in% cmSEGs)
    ruvK <- 1:5

    message("Step 2: Performing RUV normalisation. This will take minutes to hours. \n")

    ## if you subsample too small a number.. a few rows will be zeros
    ## which will lead to 0 denoms and NAs when you standardize
    ## genes element-wise 
    ruv3res <- scMerge::scRUVIII(Y = exprs_mat, 
                                   M = repMat, 
                                   ctl = ctl, 
                                   batch=colData(tmp)[,batchColname],                          
                                   k = ruvK,
                                   cell_type = colData(tmp)[,cellType])

    k<-ruv3res$optimal_ruvK
    assay(tmp,"normalized")<-t(ruv3res[[k]]$newY) #newY_mc
    # saveRDS(object=list(sce_object =tmp, scNormalization=ruv3res), '../data/tmp_scMerge.RDS')
    return(list(normed_sce =tmp, ruv3res=ruv3res))

}