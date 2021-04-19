
source('scMerge_helperFunctions.R')
library(SingleCellExperiment)



scMerge.mod <- function (sce_combine, ctl = NULL, kmeansK = NULL, exprs = "logcounts", 
                        hvg_exprs = "counts", marker = NULL, marker_list = NULL, 
                        ruvK = 20, replicate_prop = 0.5, cell_type = NULL, cell_type_match = FALSE, 
                        cell_type_inc = NULL, 
#                       fast_svd = FALSE,
                        rsvd_prop = 0.1, 
                        dist = "cor", WV = NULL, WV_marker = NULL, 
                        parallel = FALSE, 
                        parallelParam = NULL, 
                         return_all_RUV = FALSE, assay_name = NULL, 
                        verbose = FALSE) 
{

 t1 <- Sys.time()
  repMat <- scMerge::scReplicate(sce_combine = sce_combine, batch = batch, 
                                 kmeansK = kmeansK, 
                                 exprs = exprs, 
                                 hvg_exprs = hvg_exprs, 
                                 marker = marker, 
                                 marker_list = marker_list, 
                                 replicate_prop = replicate_prop, 
                                 cell_type = cell_type, 
                                 cell_type_match = cell_type_match, 
                                 cell_type_inc = cell_type_inc, 
                                 dist = dist, 
                                 WV = WV, 
                                 WV_marker = WV_marker, 
#                                parallelParam = parallelParam, 
                                 
#                                  fast_svd = fast_svd, 
                                 verbose = verbose)
  t2 <- Sys.time()
  timeReplicates <- t2 - t1
  cat("Dimension of the replicates mapping matrix: \n")
  print(dim(repMat))
  
  message("Step 2: Performing RUV normalisation. This will take minutes to hours. \n")
  ruv3res <- scRUVIII(Y = t(exprs_mat), M = repMat, ctl = ctl, 
                      k = ruvK, batch = batch, fullalpha = NULL, cell_type = cell_type, 
                      return_all_RUV = return_all_RUV, 
#                       fast_svd = fast_svd, 
                      rsvd_prop = rsvd_prop)
  t3 <- Sys.time()
  timeRuv <- t3 - t2
  
  return(list(repMat,ruv3res))
}

break

#############

setwd('/home/')
sce <- readRDS('sceFltrd_kpGns.rds')

sceLP<-sce[,sce$subjid=="IMM19"] # dim(sceLP):  22758 x 9699; table(sceLP$batch, sceLP$cellType)
assay(sceLP,1)<-as.matrix(assay(sceLP,1))
assay(sceLP,2)<-as.matrix(assay(sceLP,2))
sceLP<-sceLP[which(rowSums(assay(sceLP,2)) != 0), which(colSums(assay(sceLP,2)) != 0)] #dim(sceLP): 19524  9699
table(sceLP$cellType,sceLP$batch); dim(sceLP)


#############
smpl.prop<-1

if(smpl.prop<1){
  tmp<-c()
  for(b in unique(colData(sceLP)$batch)){
    print(b)
    tmp.idx<-unlist(
      lapply(cellType,
             function(x){
               idx<-which(colData(sceLP)$batch==b & colData(sceLP)$cellType==x)
               i<-which(cellType==x)
               #print(i)
               set.seed(123*i)
               sample(idx,smpl.prop*length(idx),replace=FALSE)
             }
      )
    )
    tmp<-c(tmp,tmp.idx)
  }
  tmp<-sort(tmp)
  tmp <- sceLP[,tmp]
}else{tmp<-sceLP}

tmp<-tmp[which(rowSums(assay(tmp,2)) != 0), which(colSums(assay(tmp,2)) != 0)] 
n.ct<-length(unique(colData(tmp)$cellType))
n.batch<-length(unique(colData(tmp)$batch))
table(tmp$cellType,tmp$batch); dim(tmp)
which(duplicated(colnames(tmp))==TRUE)

#############
require(scMerge)
require(BiocParallel)

data("segList_ensemblGeneID", package = "scMerge")
seg_index <- segList_ensemblGeneID$human$human_scSEG
cmSEGs <- intersect(row.names(tmp), seg_index)
#tmp<-tmp[tmp$genes %in% intersect(tmp$genes,seg_index),]


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
#   t2 <- Sys.time()
#   timeReplicates <- t2 - t1
#   cat("Dimension of the replicates mapping matrix: \n")
#   print(dim(repMat))

  exprs_mat <- SummarizedExperiment::assay(tmp, 'logcpm')
  ctl <- which(row.names(tmp) %in% cmSEGs)
  ruvK <- 1:5
  
#   message("Step 2: Performing RUV normalisation. This will take minutes to hours. \n")

  ## if you subsample too small a number.. a few rows will be zeros
  ## which will lead to 0 denoms and NAs 
  ruv3res <- scMerge::scRUVIII(Y = exprs_mat, 
                               M = repMat, 
                               ctl = ctl, 
                               batch=tmp$batch,                          
                               k = ruvK,
                               cell_type = tmp$cellType)
# supLP <- scMerge.mod(
#   sce_combine = subsample.tmp, #[ , which(sceLP$batch %in% c("B001","B002"))],
#   ctl = as.character(cmSEGs), #segList_ensemblGeneID$human$human_scSEG,
#   kmeansK = rep(n.ct,n.batch),
#   exprs = "logcpm",
#   hvg_exprs = "counts", #"logcpm",
#   ruvK = 1:5, #rep(NA,20),
#   replicate_prop = 1,
#   cell_type = tmp$cellType,
# #   fast_svd = TRUE,
#   rsvd_prop = 0.2, #recommended 0.01
# #   BPPARAM = MulticoreParam(workers = 2),

# #   parallel = TRUE,
# #   parallelParam = MulticoreParam(), #workers = multicoreWorkers())
#   return_all_RUV = TRUE,
#   assay_name = paste0("ruv_sup",1:5), #"scMerge_supervised",      
#   verbose=TRUE
# )
# stop.time <- Sys.time()
# print(paste0("Supervised: ",round(difftime(stop.time, start.time, units = "min"), 3))) #"Supervised: 26.348"

k<-ruv3res$optimal_ruvK
assay(tmp,"normalized")<-t(ruv3res[[k]]$newY) #newY_mc
saveRDS(object=list(sce_object =tmp, scNormalization=ruv3res), '../data/tmp_scMerge.RDS')
