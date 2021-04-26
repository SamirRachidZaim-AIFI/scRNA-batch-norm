## Load libraries 
options(warn=-1)
source('../code/scMerge_helperFunctions.R')
source('../code/quant_scMerge.R')
library(SingleCellExperiment)
require(scMerge)
require(BiocParallel)
require(ggpubr)
data("segList_ensemblGeneID", package = "scMerge")


## load brigding data
bridge_controls <- readRDS('../../data/bridging_normalized.RDS')




## Load data
options(warn=-1)
setwd('/home/scMergeNormalization/scRNA-batch-norm/code/')

sce<-readRDS("../../data/sceFltrd_kpGns.rds")
sce<-sce[,sce$subjid!="IMM19"]
sce$subjid<-factor(sce$subjid,levels=sort(unique(sce$subjid)))






wrap_varPart <- function(sub.sce, bridging_controls, segList_ensemblGeneID,
                        frmla='~ cellType + batch', seed){
    set.seed(seed)
    ### Subsample the # of cells (i.e., by barcode)
    sub.sce <- sample_cellTypes(sce, 100, FALSE)
       
    sub.sce <- normalizeClinicalSamples(sub.sce, 
                                    bridging_controls,
                                    segList_ensemblGeneID$human$human_scSEG
                                   )

    ## plot before & after on clinical samples 

    ratio_df <- quant_varPart(sub.sce, frmla= frmla)
    return(ratio_df)
                              
}

ratio_list <- mclapply(1:12, function(x) wrap_varPart(sub.sce, bridging_controls, 
                                                      segList_ensemblGeneID, 
                                                      frmla='~ cellType + batch +Week', x),
                       mc.cores=detectCores()-1)
                       
                       


break



































































wks<-sort(unique(sce$wk))
names(fres)










#### Extract Elements from scMerge/RUV normalization
data("segList_ensemblGeneID", package = "scMerge")
seg_index <- segList_ensemblGeneID$human$human_scSEG
cmSEGs <- intersect(row.names(fres$sce_object), seg_index)

k<-fres[[2]]$optimal_ruvK ## k=1
falpha<-fres[[2]][[k]]$fullalpha # dim(falpha) 52 x 20240 (k x G)
colnames(falpha)<-rownames(fres$sce_object)
alpha <- falpha[seq_len(min(k, nrow(falpha))), , drop = FALSE] # dim(alpha) : 1 x 20240
ac <- alpha[,cmSEGs, drop = FALSE] # dim(ac): 1 x 52; sum(is.na(ac[1,]))
head(ac)
d=b=1
tmp <- fres$sce_object

#remove zeros
sce<-sce[which(rowSums(assay(sce,2)) != 0), which(colSums(assay(sce,2)) != 0)]

dim(sce)
dim(tmp)

idx <- subSmpl.lst[[d]][b, ]

### run on subsampled data
sub.sce <- sce[which(rownames(sce) %in% rownames(tmp)) , which(sce$wk==wks[d])[idx] ]
sub.sce <- normalizeClinicalSamples(sub.sce)

###
saveRDS(sub.sce, '../../data/sub.sce_norm.RDS')

    
