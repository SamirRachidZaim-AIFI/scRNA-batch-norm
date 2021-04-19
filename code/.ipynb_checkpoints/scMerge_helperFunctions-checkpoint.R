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
################################################################################################


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
