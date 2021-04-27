## Load libraries 
options(warn=-1)
source('../code/scMerge_helperFunctions.R')
source('../code/quant_scMerge.R')
library(SingleCellExperiment)
require(scMerge)
require(BiocParallel)
require(ggpubr)
require(variancePartition)
data("segList_ensemblGeneID", package = "scMerge")


## load brigding data
bridging_controls <- readRDS('../../data/bridging_normalized.RDS')




## Load data
options(warn=-1)
setwd('/home/scMergeNormalization/scRNA-batch-norm/code/')

sce<-readRDS("../../data/sceFltrd_kpGns.rds")
sce<-sce[,sce$subjid!="IMM19"]
sce$subjid<-factor(sce$subjid,levels=sort(unique(sce$subjid)))

wk2.sce <- sce[,sce$Week==2]


Nsim=5
subsampleN=500
sub.sceListEven <- lapply(1:Nsim, function(x) sample_cellTypes(wk2.sce, 
                                                            subsampleN, 
                                                            stratified=FALSE))

sub.sceListStrat <- lapply(1:Nsim, function(x) sample_cellTypes(wk2.sce, 
                                                            subsampleN, 
                                                            stratified=TRUE))

                           
seg_index <- segList_ensemblGeneID$human$human_scSEG 
cmSEGs <- intersect(seg_index, row.names(wk2.sce))                           
ctl.gn <- cmSEGs[1:10]
                           
                           
sce_tmp <- sub.sceListEven[[1]]@assays@.xData$data$logcpm                        
sce_tmp_norm <- normalizeClinicalSamples(sub.sceListEven[[1])                           
                

### Welch's T-test                         
welchs_ttest <- function(ctl.gn,sub.sceListEven){
    
    b1 <- sub.sceListEven[[1]][,sub.sceListEven[[1]]$batch=='X001']
    b2 <- sub.sceListEven[[1]][,sub.sceListEven[[1]]$batch=='X002']   
    
    v1<-b1@assays$data$logcpm[ctl.gn,]
    v2<-b2@assays$data$logcpm[ctl.gn,] 
    t.test(v1,v2, paired=FALSE, var.equal=FALSE)$p.value
    
}                           
                           
prenorm <- sapply(1:10, function(x) get_ttest(ctl.gn[x], sub.sceListEven))           
postnorm<- sapply(1:10)                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
###### ############ ###### ############ ###### ############ ######
###### ############ ###### ############ ###### ############ ######                           
                           
# Variance Partition Analysis 
# Compare ratio of pre- and post-normalization
# of batch effect contributions (absolute & percentage)
# using the variance partition R package

wrap_varPart <- function(sub.sce, bridging_controls, segList_ensemblGeneID,
                        frmla='~ cellType + batch'){
    # normalize
    sub.sce <- normalizeClinicalSamples(sub.sce, 
                                    bridging_controls,
                                    segList_ensemblGeneID$human$human_scSEG
                                   )
    rm(bridging_controls)
    print(pryr::mem_used())

    ## calculate pre- post varpart
    ratio_df <- quant_varPart(sub.sce, frmla= frmla)
    return(ratio_df)
                              
}

require(parallel)
## Look at how batch effects look in stratified vs. even sampling
even_smpl_ratio_list <- mclapply(1:Nsim, function(x) wrap_varPart(sub.sce=sub.sceListEven[[x]], 
                                                                bridging_controls=bridging_controls, 
                                                                segList_ensemblGeneID, 
                                                                frmla='~ cellType * batch'
                                                                ),
                       mc.cores=detectCores()-1)
                   

str_smpl_ratio_list <- mclapply(1:Nsim, function(x) wrap_varPart(sub.sce=sub.sceListStrat[[x]], 
                                                                bridging_controls=bridging_controls, 
                                                                segList_ensemblGeneID, 
                                                                frmla='~ cellType * batch'
                                                                ),
                       mc.cores=detectCores()-1)




even<- do.call(rbind, even_smpl_ratio_list)
strat<-do.call(rbind, str_smpl_ratio_list)

ratio_even <- even[grep('Ratio',row.names(even)),]
ratio_strt <- strat[grep('Ratio',row.names(strat)),]
                           
dfstrt <- data.table::melt(ratio_strt)
dfeven <- data.table::melt(ratio_even)                           

pdf('strat.pdf')
ggplot(dfstrt, aes(x=variable, y=value, col=variable))+geom_boxplot()+ylim(1,1.5)
dev.off()
                  
write.csv(even,'even.csv')
write.csv(strat,'strat.csv')                             
                  
pdf('even.pdf')                  
ggplot(dfeven, aes(x=variable, y=value, col=variable))+geom_boxplot()      +ylim(1,1.5)            
dev.off()
### plot the ratio of variance conservation 

### 

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

    
