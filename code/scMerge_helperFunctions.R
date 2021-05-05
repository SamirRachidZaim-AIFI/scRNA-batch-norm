################################################################################################
################################################################################################

# Title: scMerge_helperFunctions 
# Author: Monica Chaudhari
# Last Modified: Samir Rachid Zaim 
# Date Modified: 4/19/21

# Desc: 
#     This file contain a list of auxiliary helper functions that serve as a companion for 
#     for driver_scMerge.R script that executes single-cell RNA batch correction and 
#     data set merging via the use of our bridge controls. 

################################################################################################
################################################################################################

# Func: solve_axb_mod
# Inputs:
#     a = m x n matrix 
#     b = n x p matrix 

# output:
#     x = resulting vector

solve_axb_mod = function(a, b){
  x = solve(DelayedArray::t(a) %*% a) %*% DelayedArray::t(a) %*% b
  return(x)
}

################################################################################################
################################################################################################


################################################################################################
################################################################################################

# Func: standardize2_mod
# Inputs:
#     Y = gene expression matrix (m x n matrix)
#     batch = m x 1 vector indicating the batch for sample

# output2:
#     stand_Y = standardized gene expression matrix, 
#     stand_mean = standardized means from gene expression matrix, 
#     stand_var = standardized variances from gene expression matrix

standardize2_mod <- function(Y, batch) {
  num_cell <- ncol(Y)
  num_batch <- length(unique(batch))
  batch <- as.factor(batch)
  stand_mean <- DelayedArray::rowMeans(Y)
  if(num_batch>1){
    design <- stats::model.matrix(~-1 + batch)    
  }else{
    design<-matrix(as.numeric(batch),nrow=ncol(Y))
  }
  B.hat = solve_axb_mod(a = DelayedArray::t(design) %*% design, b = DelayedArray::t(Y %*% design))
  stand_var <- DelayedArray::rowSums(((Y - DelayedArray::t(B.hat) %*% DelayedArray::t(design))^2))/(num_cell - num_batch)
  stand_Y <- (Y-stand_mean)/sqrt(stand_var)
  return(res = list(stand_Y = stand_Y, stand_mean = stand_mean, stand_var = stand_var))
}

################################################################################################

visualize_normalization<- function(sub.sce, numComponentsPCA = NULL){
    
        #################### LP: Run PCA - PRIOR/POST ##################################
    

        ## Calculate principal components for T-SNE & UMAP Dimension Reduction 
        BC.prior.smpl = scater::runPCA(sub.sce, exprs_values = "logcpm",ncomponents=5) 
        PCs.prior.smpl<- data.frame(cbind(reducedDim(BC.prior.smpl),
                        "cellType"=colData(sub.sce)$cellType,"batch"=colData(sub.sce)$batch))
        PCs.prior.smpl[,1:5]<-apply(PCs.prior.smpl[,1:5],2,
                                    function(x){as.numeric(as.character(x))})

        ## Run Multi-variate anova to evaluate batch and batch X celltype effects
        op <- options(contrasts = c("contr.helmert", "contr.poly"))
        maov.prior.smpl<-manova(cbind(PC1,PC2,PC3,PC4,PC5) ~
                                as.factor(cellType)*as.factor(batch),data=PCs.prior.smpl)   

        summary(maov.prior.smpl)

        fname <- paste('../../plots/clinical_PCAs.pdf')
        pdf(fname) 
        p1<-scater::plotPCA(BC.prior.smpl, colour_by = "cellType", shape_by = "batch")
        p2<-scater::plotPCA(BC.prior.smpl, colour_by = "batch")

        ggarrange(p1,p2, ncol=1)

        BC.post.smpl = scater::runPCA(sub.sce, exprs_values = "normalized",ncomponents=5); 
        PCs.post.smpl<-data.frame(cbind(reducedDim(BC.post.smpl),
                                            "cellType"=colData(sub.sce)$cellType,
                                            "batch"=colData(sub.sce)$batch))
        PCs.post.smpl[,1:5]<-apply(PCs.post.smpl[,1:5],2,
                                   function(x){as.numeric(as.character(x))})
        op <- options(contrasts = c("contr.helmert", "contr.poly"))
        maov.post.smpl<-manova(cbind(PC1,PC2,PC3,PC4,PC5) ~ 
                               as.factor(cellType)*as.factor(batch),
                                   data=PCs.post.smpl)
        summary(maov.post.smpl)

        p1<-scater::plotPCA(BC.post.smpl, colour_by = "cellType", shape_by = "batch")
        p2<-scater::plotPCA(BC.post.smpl, colour_by = "batch")

        ggarrange(p1,p2, ncol=1) 
        dev.off()

        ##################### LP: Run TSNE - PRIOR/POST #############################
        fname <- paste('../../plots/clinical_tsne.pdf')
        pdf(fname) 
        set.seed(1000)
        BC.prior.smpl <- scater::runTSNE(BC.prior.smpl, perplexity=50, 
                                         dimred="PCA", exprs_values ="logcpm",  n_dimred=5)

        p1<-scater::plotTSNE(BC.prior.smpl, colour_by = "cellType", shape_by = "batch")
        p2<-scater::plotTSNE(BC.prior.smpl, colour_by = "batch")
        ggarrange(p1,p2, ncol=1) 


        BC.post.smpl <- scater::runTSNE(BC.post.smpl, 
                                        perplexity=50, dimred="PCA", 
                                        exprs_values = "normalized",  n_dimred=5)
        p3<-scater::plotTSNE(BC.post.smpl, colour_by = "cellType", shape_by = "batch")
        p4<-scater::plotTSNE(BC.post.smpl, colour_by = "batch")

        ggarrange(p3,p4, ncol=1) 
        dev.off()

        ##################### LP: Run UMAP - PRIOR/POST #####################################
        BC.prior.smpl <- scater::runUMAP(BC.prior.smpl,  exprs_values = "logcpm",  
                                             n_dimred=5)
        BC.post.smpl <- scater::runUMAP(BC.post.smpl, exprs_values = "normalized", 
                                             n_dimred=5)    

        outfile = "../../plots/clinical_UMAPs.pdf" 
        pdf(file = outfile, width = 7, height = 5)
        p1<-scater::plotUMAP(BC.prior.smpl, colour_by = "cellType", shape_by = "batch")
        p2<-scater::plotUMAP(BC.prior.smpl, colour_by = "batch")
        ggarrange(p1,p2, ncol=1) 

        p3<-scater::plotUMAP(BC.post.smpl, colour_by = "cellType", shape_by = "batch")
        p4<-scater::plotUMAP(BC.post.smpl, colour_by = "batch")
        ggarrange(p3,p4, ncol=1) 

        dev.off()

}

quantify_variance_variancePartition <- function(tmp){
    require(variancePartition)
    # head(tmp)
    geneExpr_norm <-   as.matrix(tmp@assays@data$normalized)
    geneExpr_logcpm <- as.matrix(tmp@assays@data$logcpm)
    info <- data.frame(tmp@colData)
    info$cellType <- as.factor(info$cellType)
    info$batch <- as.factor(info$batch)
    varPartfrmla <- ~ cellType*batch
    varPart_pre <- fitExtractVarPartModel( geneExpr_norm, varPartfrmla, info )
    varPart_post <- fitExtractVarPartModel( geneExpr_logcpm, varPartfrmla, info )
    
    vp1 <- sortCols( varPart_pre )
    vp2 <- sortCols( varPart_post )
    
    outfile = "../../plots/clinical_variancePartitions.pdf" 
    pdf(file = outfile, width = 7, height = 5)

    par(mfrow=c(1,2))
    plotVarPart( vp1 )
    plotVarPart( vp2 )
    dev.off()
    
}

# quantify_variance_variancePartition(sub.sce)

quantify_variance_pvca <- function(tmp){
    source('../code/pvca.R')
    
    gdata <- as.matrix(assay(tmp,'normalized'))
    gdata <- data.frame(gdata)
    mdata <- info <- data.frame(tmp@colData)

    names(gdata) <- gsub('X','', names(gdata))
    
    sid <- 'barcodes'
    batch.factors <- c('cellType','batch')
    threshold <- 0.8
    interaction=T
   
    pvcaObject<-runPVCA(gdata=gdata, mdata=mdata, sid=sid, factors = batch.factors, 0.8, T)

    outfile = "../../plots/clinical_pvca.pdf" 
    pdf(file = outfile, width = 7, height = 5)

    PlotPVCA(pvcaObject, "PVCA") 
    dev.off()
    
}




################################################################################################
################################################################################################

# Do not use! (Not stable - still in dev.) 


# Func: scMerge.mod
# Inputs:
#     - TBD 

# Outputs:
#     - TBD 

# scMerge.mod <- function (sce_combine, ctl = NULL, kmeansK = NULL, exprs = "logcounts", 
#                         hvg_exprs = "counts", marker = NULL, marker_list = NULL, 
#                         ruvK = 20, replicate_prop = 0.5, cell_type = NULL, cell_type_match = FALSE, 
#                         cell_type_inc = NULL, 
# #                       fast_svd = FALSE,
#                         rsvd_prop = 0.1, 
#                         dist = "cor", WV = NULL, WV_marker = NULL, 
#                         parallel = FALSE, 
#                         parallelParam = NULL, 
#                          return_all_RUV = FALSE, assay_name = NULL, 
#                         verbose = FALSE) 
# {
#   cellNames = colnames(sce_combine)
#   if (length(cellNames) != length(unique(cellNames))) {
#     stop("Please make sure column names are unique.")
#   }
#   if (is.null(assay_name)) {
#     stop("assay_name is NULL, please provide a name to store the results under")
#   }
#   if (length(ruvK) > 1) {
#     message("You chose more than one ruvK. The argument return_all_RUV is forced to be TRUE.")
#     return_all_RUV = TRUE
#   }
#   if (return_all_RUV) {
#     message("You chose return_all_RUV = TRUE. The result will contain all RUV computations. This could be a very large object.")
#     if (length(assay_name) != length(ruvK)) {
#       stop("You chose return_all_RUV = TRUE. In this case, the length of assay_name must be equal to the length of ruvK")
#     }
#   }
#   if (is.null(exprs)) {
#     stop("exprs is NULL.")
#   }
#   if (!exprs %in% SummarizedExperiment::assayNames(sce_combine)) {
#     stop(paste("No assay named", exprs))
#   }
#   exprs_mat <- SummarizedExperiment::assay(sce_combine, exprs)
#   if (!is.matrix(exprs_mat)) {
#     stop(paste0("The assay named '", exprs, "' must be of class 'matrix', please convert this."))
#   }
#   sce_rownames <- rownames(sce_combine)
#   hvg_exprs_mat <- SummarizedExperiment::assay(sce_combine, 
#                                                hvg_exprs)
#   if (!is.matrix(hvg_exprs_mat)) {
#     stop(paste0("The assay named '", hvg_exprs, "' must be of class 'matrix', please convert this."))
#   }
#   if (any(base::rowSums(exprs_mat) == 0) | any(base::colSums(exprs_mat) == 
#                                                0)) {
#     stop("There are rows or columns that are all zeros in the expression matrix. Please remove these rows/columns.")
#   }
#   if (is.null(ctl)) {
#     stop("Negative control genes are needed. \n")
#   }
#   else {
#     if (is.character(ctl)) {
#       ctl <- which(sce_rownames %in% ctl)
#     }
#     if (length(ctl) == 0) {
#       stop("Could not find any negative control genes in the row names of the expression matrix", 
#            call. = FALSE)
#     }
#   }
#   if (is.null(sce_combine$batch)) {
#     stop("Could not find a 'batch' column in colData(sce_combine)", 
#          call. = FALSE)
#   }
#   if (is.factor(sce_combine$batch)) {
#     batch <- droplevels(sce_combine$batch)
#   }
#   else {
#     batch <- sce_combine$batch
#   }
#   if (!is.null(parallelParam)) {
#     message("Step 1: Computation will run in parallel using supplied parameters")
#   }
#   if (parallel & is.null(parallelParam)) {
#     message("Step 1: Computation will run in parallel using BiocParallel::bpparam()")
# #     parallelParam = BiocParallel::bpparam()
#   }
#   if (!parallel | is.null(parallelParam)) {
#     message("Step 1: Computation will run in serial")
# #     parallelParam = BiocParallel::SerialParam()
#   }
#   t1 <- Sys.time()
#   repMat <- scMerge::scReplicate(sce_combine = sce_combine, batch = batch, 
#                                  kmeansK = kmeansK, 
#                                  exprs = exprs, 
#                                  hvg_exprs = hvg_exprs, 
#                                  marker = marker, 
#                                  marker_list = marker_list, 
#                                  replicate_prop = replicate_prop, 
#                                  cell_type = cell_type, 
#                                  cell_type_match = cell_type_match, 
#                                  cell_type_inc = cell_type_inc, 
#                                  dist = dist, 
#                                  WV = WV, 
#                                  WV_marker = WV_marker, 
# #                                  parallelParam = parallelParam, 
                                 
#                                  fast_svd = fast_svd, verbose = verbose)
#   t2 <- Sys.time()
#   timeReplicates <- t2 - t1
#   cat("Dimension of the replicates mapping matrix: \n")
#   print(dim(repMat))
  
#   message("Step 2: Performing RUV normalisation. This will take minutes to hours. \n")
#   ruv3res <- scRUVIII(Y = t(exprs_mat), M = repMat, ctl = ctl, 
#                       k = ruvK, batch = batch, fullalpha = NULL, cell_type = cell_type, 
#                       return_all_RUV = return_all_RUV, fast_svd = fast_svd, 
#                       rsvd_prop = rsvd_prop)
#   t3 <- Sys.time()
#   timeRuv <- t3 - t2
  
#   return(list(repMat,ruv3res))
# } 
