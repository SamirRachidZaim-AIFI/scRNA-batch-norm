####################################################################
####################################################################

# The deep & wide analysis is a study to explore how to 
# anchor scRNA batch correction on the deep & wide
# replicates. The analysis consists of assessing which of 
# the following approaches leads to the best batch correction
# results across elements in different batches: 

### load libraries 
source('/home/jupyter/scRNA-batch-norm/code/scMerge_helperFunctions.R')
source('/home/jupyter/scRNA-batch-norm/code/ruviii_bridging_controls.R')
source('/home/jupyter/scRNA-batch-norm/code/quant_scMerge.R')
source('/home/jupyter/scRNA-batch-norm/code/normalizeClinicalSamples.R')
source('/home/jupyter/scRNA-batch-norm/code/quant_scMerge.R')
source('/home/jupyter/scRNA-batch-norm/code/pvca.R')
require(scMerge)
require(BiocParallel)
require(Seurat)
setwd('/home/jupyter/deepNwide/')
require(H5weaver)
require(SingleCellExperiment)
require(ggpubr)
require(dplyr)
require(magrittr)

####################################################################
####################################################################

####################################################################
####################################################################

# Now to show benefits of batch correction we need to identify samples
# with potential batch effects so that when we apply the 
# batch correction procedure we can detect "measurable" changes.

# Patients in batch 002, 004, and 005 are great starter candidates. 
# # and in batches 43 & 39 respectively for day 0/7. 
# So we can look at bridging controls for batches 43, 39, 50
# and see if the batch effects are present there, and then 
# use those to correct. 

### load br1 samples 
br1Lookup <- read.csv('/home/jupyter/scRNA-batch-norm/bri_analysis/lookupTable_br1.csv')
clinical_batch_idx <- which(br1Lookup$subject.subjectGuid!='' &
                           br1Lookup$file.batchID %in% c('B004','B005'))
                            
clin_meta <- br1Lookup[clinical_batch_idx,]
sub_clin_meta <- clin_meta[seq(1, 23, 2), ]

### select "latest" sample per sample.samplekitguid
### since some samples can have multiple 'entries'

sub_clin_meta <- sub_clin_meta %>%
 arrange(sample.sampleKitGuid,lastUpdated) %>% 
 group_by(sample.sampleKitGuid) %>% 
 summarise(across(everything(), last))

sub_clin_meta <- as.data.frame(sub_clin_meta)
setwd('/home/jupyter/scRNA-batch-norm/bri_analysis')
fnames <- sub_clin_meta$filePath
br1_clin_smps <- lapply(fnames, function(x) read_h5_sce(x)
                   )

####################################################################
####################################################################


####################################################################
####################################################################

# Now we have to load the corresponding bridging controls to 
# conduct the ruviii normalization and use that to 'anchor' our 
# batch correction procedure. 

# load bridging samples from batch 2, 4,5
brgd2_4_5.idx <- br1Lookup$sample.bridgingControl =="" & br1Lookup$file.batchID %in% c('B004','B005')
fnames <- br1Lookup$filePath[brgd2_4_5.idx]
br1_ctls <- lapply(fnames, function(x) read_h5_sce(x)
                   )

# Most bridging controls contain ~ 4,5k cells, whereas
# others contain ~20k. Here we test for cell
# abundance in the bridging controls 

br1_ctls <- br1_ctls[c(1,3)]
br1_ctls <- do.call(cbind, br1_ctls)
br1_ctls <- scater::logNormCounts(br1_ctls)

numComponents=50
cellType_name ='seurat_pbmc_type'
batchColName  ='batch_id'


## plot batch effects
numpComponents=50
plot_batchEffects(br1_ctls, 
                  exprs_values='logcounts', 
                  numComponents=numpComponents,
                  type='bridge_20k'
                 )

### After log-transforming the raw counts, 
### obtain the ruviii object with the bridge
### control normalization

br1_ctls <- ruviii_bridging_controls(br1_ctls, 
                                      exprsname='logcounts', 
                                      cellType='seurat_pbmc_type',
                                      batchColname='batch_id')

####################################################################
####################################################################


####################################################################
####################################################################

#### compare visuals before & after 
#### for bridge controls 

plot_batchEffects(br1_ctls$normed_sce, 
                  exprs_values='normalized', 
                  numComponents=numpComponents,
                  type='bridge'
                 )
####################################################################
####################################################################

####################################################################
####################################################################
####
####
#### batch correct and evaluate 
#### visuals on clinical samples

### batch correct clinical samples 
br1_clin_smps <- do.call(cbind, br1_clin_smps)
br1_clin_smps <- scater::logNormCounts(br1_clin_smps)

### run normalization on a 50k samples/cells (20% of data) runs in ~ 2mins

sub.sce <- br1_clin_smps[,br1_clin_smps$barcodes %in% sample(br1_clin_smps$barcodes, 50000)]
system.time(sub.sce <- transform_and_merge_sce(sub.sce, sub_clin_meta, br1_ctls))



### plot batch effects (tsne, umap, pca)
plot_batchEffects(sub.sce, 
                  exprs_values='logcounts', 
                  numComponents=numpComponents,
                  type='clin_sub'                  
                 )

plot_batchEffects(sub.sce, 
                  exprs_values='normalized', 
                  numComponents=numpComponents,
                  type='clin_sub'                  
                 )

####################################################################
####################################################################


####################################################################
####################################################################
require(variancePartition)
suppressPackageStartupMessages(library(SingleCellExperiment))
require(scMerge)
require(BiocParallel)


#### quantify variance partition bridging controls 

# gdata_logct <- data.frame(assay(br1_ctls$normed_sce, 'logcounts')) ; names(gdata_logct) <- gsub('X','', names(gdata_logct))
# gdata_norm <- data.frame(assay(br1_ctls$normed_sce, 'normalized')) ; names(gdata_norm) <- gsub('X','', names(gdata_norm))
# mdata <- info <- br1_ctls$normed_sce@colData
# sid <- 'barcodes'
# batch.factors <- c('seurat_pbmc_type','batch_id')
# threshold <- 0.8
# interaction=FALSE

estimatePVCA <- function(ncells = 200){
    idx <- sample(ncol(gdata_logct),ncells)
    pvca_logct<-runPVCA(gdata=gdata_logct[,idx], mdata=mdata, sid=sid, factors = batch.factors, 0.2, interaction)
    pvca_norm<-runPVCA(gdata=gdata_norm[,idx], mdata=mdata, sid=sid, factors = batch.factors, 0.2, interaction)

    df <- rbind(pvca_logct,
            pvca_norm,
            pvca_logct/pvca_norm
           )
    df <- round(df,2)
    row.names(df) <- c('Pre-norm','Normalized','Ratio')
    return(df)
}

# res <- mclapply(1:10, function(x) estimatePVCA(),
#                 mc.cores = 25)

# res.df <- do.call(rbind,res )
# ratios <- res.df[row.names(res.df)%in% 'Ratio',]
# prenorm <- res.df[row.names(res.df)%in% 'Pre-norm',]
# postnorm <- res.df[row.names(res.df)%in% 'Normalized',]

####################################################################
####################################################################






####################################################################
####################################################################

## Define function 

#### quantify variance partition clinical samples

sample_pvcaEstimates <- function(sub.sce =sub.sce,
                                 sid = 'barcodes',
                                 batch.factors = c('batch_id','subject.subjectGuid'),
                                 ncells = 500,
                                 nSims = 25,
                                 mc.cores=30
                                 
                                ){
    mdata <- sub.sce@colData
    br1_gdata_logct <- data.frame(assay(sub.sce, 'logcounts')); colnames(br1_gdata_logct) <- mdata$barcodes
    br1_gdata_norm <- data.frame(assay(sub.sce, 'normalized')); colnames(br1_gdata_norm) <- mdata$barcodes  
    


    interaction=FALSE

    estimatePVCA <- function(ncells ,
                             gdata_logct,
                             gdata_norm 
                            ){
        idx <- sample(ncol(gdata_logct),ncells)
        pvca_logct<-runPVCA(gdata=gdata_logct[,idx], mdata=mdata, 
                            sid=sid, factors = batch.factors, 
                            0.75, interaction)
        pvca_norm<-runPVCA(gdata=gdata_norm[,idx], mdata=mdata, 
                           sid=sid, factors = batch.factors,
                           0.75, interaction)

        df <- rbind(pvca_logct,
                pvca_norm,
                pvca_logct/pvca_norm
               )
        df <- round(df,2)
        row.names(df) <- c('Pre-norm','Normalized','Ratio')
        return(df)
    }

    res <- mclapply(1:nSims, function(x) estimatePVCA(ncells=ncells, 
                                                   gdata_logct =br1_gdata_logct, 
                                                   gdata_norm =br1_gdata_norm),
                    mc.cores=mc.cores
                    )


    res.df <- do.call(rbind,res )
    ratios <- res.df[row.names(res.df)%in% 'Ratio',]
    prenorm <- res.df[row.names(res.df)%in% 'Pre-norm',]
    postnorm <- res.df[row.names(res.df)%in% 'Normalized',]
    
    return(list(Results= res.df, 
                RatioDF= ratios,
                Prenorm= prenorm, 
                Postnorm=postnorm)
          )
    
}


####################################################################
####################################################################

####################################################################
####################################################################

require(variancePartition)
### Run variance Partition 
ncells = 1000
nsims=5

tmp <- br1_ctls$normed_sce[,br1_ctls$normed_sce$barcodes %in% sample(br1_ctls$normed_sce$barcodes, 5000)]
tmp<-tmp[which(rowSums(assay(tmp,2)) != 0), which(colSums(assay(tmp,2)) != 0)] 


param <- SnowParam(workers = 10, type = "SOCK")
quant_varPart(tmp,
             frmla= '~ (1|batch_id)',
             numGenesToAnalyze=1000,
             numCells=5000,
             param=param)



tmp2 <- sub.sce[,sub.sce$barcodes %in% sample(sub.sce$barcodes, 5000)]
tmp2<-tmp2[which(rowSums(assay(tmp2,2)) != 0), which(colSums(assay(tmp2,2)) != 0)] 


param <- SnowParam(workers = 10, type = "SOCK")
quant_varPart(tmp2,
             frmla= '~(1+seurat_cell_type + |batch_id)',
             numGenesToAnalyze=1000,
             param=param)



####################################################################
####################################################################


####################################################################
####################################################################

### Run PVCA estimates 
ncells = 1000
nsims=5

pvca0_br = sample_pvcaEstimates(sub.sce =br1_ctls$normed_sce,
                             sid = 'barcodes',
                             batch.factors = c('batch_id'),
                             ncells = ncells,
                             nSims = nsims,
                             mc.cores=5
                             )


### run 0th = 'batch id'
pvca0 = sample_pvcaEstimates(sub.sce =sub.sce,
                             sid = 'barcodes',
                             batch.factors = c('batch_id'),
                             ncells = ncells,
                             nSims = nsims,
                             mc.cores=5
                             )
pryr::mem_used()

df = data.frame(rbind( colMeans(pvca0$Prenorm),
            colMeans(pvca0$Postnorm),
            colMeans(pvca0$RatioDF)
           ))
row.names(df) = c('Prenorm','Postnorm','Ratio')
colnames(df) <- c('Batch_id','Residuals')
write.csv(df, file='/home/jupyter/scRNA-batch-norm/bri_analysis/results/pvca0.csv')



### run 1 = 'batch id' + subject ID 
pvca1 = sample_pvcaEstimates(sub.sce =sub.sce,
                             sid = 'barcodes',
                             batch.factors = c('batch_id','subject.subjectGuid'),
                             ncells = ncells,
                             nSims = nsims,
                             mc.cores=5
                             )

pryr::mem_used()
df1 = data.frame(rbind( colMeans(pvca1$Prenorm),
            colMeans(pvca1$Postnorm),
            colMeans(pvca1$RatioDF)
           ))
row.names(df1) = c('Prenorm','Postnorm','Ratio')
colnames(df1) <- c('Batch_id','Subject','Residuals')
write.csv(df1, file='/home/jupyter/scRNA-batch-norm/bri_analysis/results/pvca1.csv')




### run 2= 'batch id' + 'subject' + 'sex'

pvca2 = sample_pvcaEstimates(sub.sce =sub.sce,
                             sid = 'barcodes',
                             batch.factors = c('batch_id','subject.subjectGuid','subject.biologicalSex'),
                             ncells = ncells,
                             nSims = nsims,
                             mc.cores=5
                             )
pryr::mem_used()
df2 = data.frame(rbind( colMeans(pvca2$Prenorm),
            colMeans(pvca2$Postnorm),
            colMeans(pvca2$RatioDF)
           ))
row.names(df2) = c('Prenorm','Postnorm','Ratio')
colnames(df2) <- c('Batch_id','Subject','Sex', 'Residuals')
write.csv(df2, file='/home/jupyter/scRNA-batch-norm/bri_analysis/results/pvca2.csv')

df
df1
df2

####################################################################
####################################################################


####################################################################
####################################################################

### try normalizing on data from 10x genommics 

# require(TENxPBMCData)

# sce <- TENxPBMCData(dataset='pbmc3k')
# sce2 <- TENxPBMCData(dataset='pbmc4k')

# universe <- intersect(rownames(sce), rownames(sce2))
# length(universe)
# ## [1] 31232
# # Subsetting the SingleCellExperiment object.
# sce <- sce[universe,]
# sce2 <- sce2[universe,]
                     
# rowData(sce) <- rowData(sce2)
# sce$batch <- "3k"
# sce2$batch <- "4k"
# uncorrected <- cbind(sce, sce2)
# uncorrected <- scater::logNormCounts(uncorrected)

# library(scater)
# set.seed(0010101010)
# uncorrected <- runPCA(uncorrected, BSPARAM=BiocSingular::RandomParam())

# set.seed(1111001)
# uncorrected <- scater::runTSNE(uncorrected, dimred="PCA")
# png('/home/jupyter/scRNA-batch-norm/bri_analysis/plots/10x.png')
# scater::plotPCA(uncorrected, colour_by="batch")
# dev.off()

# png('/home/jupyter/scRNA-batch-norm/bri_analysis/plots/10x2.png')
# scater::plotTSNE(uncorrected, colour_by="batch")
# dev.off()

# ctl = segList_ensemblGeneID$human$human_scSEG

# colnames(uncorrected) <- colData(uncorrected)$Barcode
# tenX <- scMerge(uncorrected, 
#                 kmeansK=1,
#                 ctl = ctl,
#                 assay_name='normalized',
#                 ruvK=5)

# tenX <- runPCA(tenX, BSPARAM=BiocSingular::RandomParam())



# BC.post.smpl = scater::runPCA(tenX, exprs_values = "normalized",ncomponents=50); 
    
# PCs.post<-data.frame(cbind(reducedDim(BC.post.smpl),"batch"=colData(tenX)[,'batch']))   
# PCs.post[,1:50]<-apply(PCs.post[,1:50],2,function(x){as.numeric(as.character(x))})
# p2<-scater::plotPCA(BC.post.smpl, colour_by = "batch")


# png('/home/jupyter/scRNA-batch-norm/bri_analysis/plots/10x_norm.png')
# print(p2)
# dev.off()























































### run normalization on entire set and plot normalization
br1_clin_smps_norm <- normalizeClinicalSamples(br1_clin_smps,
                                               'logcounts',
                                               'batch_id',
                                               br1_ctls
                                               )

br1_clin_smps_norm



















#### compare variance partition before & after 
require(parallel)



cl <- makeCluster(25)
registerDoParallel(cl)

df = quant_varPart(br1_ctls$normed_sce,
              frmla='~ n_genes + 1|batch_id',
              numGenesToAnalyze=2000
              )

stopCluster(cl)





































## make sure each sample has the same "genes"
gene_list <- lapply(br1_clin_smps, function(x) row.names(x)
                    )

cmnGenes <- Reduce(intersect, gene_list)
br1_clin_smps2 <- lapply(br1_clin_smps, function(x) x[row.names(x) %in% cmnGenes]
                         )
## concatenate all the sce objects
sce_clin <- do.call(cbind, br1_clin_smps)

sub.sce <- sce_clin[,sce_clin$barcodes %in% sample(sce_clin$barcodes, 2000)]


sub.sce <- transform_and_merge_sce(sub.sce, clin_meta, ruviii_15kcells)


####################################################################
####################################################################


### join and obtain metadata from clinical h5s to 


####################################################################
####################################################################


quant_varPart(sub.sce,
             frmla='~ 1|batch_id + 1|sample.visitName + 1| subject.subjectGuid')
                         
                         
                         