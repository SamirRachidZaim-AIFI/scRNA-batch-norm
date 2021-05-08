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
                           br1Lookup$file.batchID %in% c('B002','B004','B005'))
                            
clin_meta <- br1Lookup[clinical_batch_idx,]
sub_clin_meta <- clin_meta[seq(1, 32, 2), ]

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
brgd2_4_5.idx <- br1Lookup$sample.bridgingControl =="" & br1Lookup$file.batchID %in% c('B002','B004','B005')
fnames <- br1Lookup$filePath[brgd1_4_5.idx]
br1_ctls <- lapply(fnames, function(x) read_h5_sce(x)
                   )

# Most bridging controls contain ~ 4,5k cells, whereas
# others contain ~20k. For now we will do correction on 
# the samples with 5k cells, and later see if the 20k cells 
# make a bigger difference. 

br1_ctls <- br1_ctls[c(1,4,5)]
br1_ctls <- do.call(cbind, br1_ctls)
br1_ctls <- scater::logNormCounts(br1_ctls)

numComponents=50
cellType_name ='seurat_pbmc_type'
batchColName  ='batch_id'


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

numpComponents=50
plot_batchEffects(br1_ctls$normed_sce, 
                  exprs_values='logcounts', 
                  numComponents=numpComponents,
                  type='bridge'
                 )

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

### run normalization on a 50k samples/cells (20% of data)
sub.sce <- br1_clin_smps[,br1_clin_smps$barcodes %in% sample(br1_clin_smps$barcodes,50000)]
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

gdata_logct <- data.frame(assay(br1_ctls$normed_sce, 'logcounts')) ; names(gdata_logct) <- gsub('X','', names(gdata_logct))
gdata_norm <- data.frame(assay(br1_ctls$normed_sce, 'normalized')) ; names(gdata_norm) <- gsub('X','', names(gdata_norm))
mdata <- info <- br1_ctls$normed_sce@colData
sid <- 'barcodes'
batch.factors <- c('seurat_pbmc_type','batch_id')
threshold <- 0.8
interaction=FALSE

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

res <- mclapply(1:10, function(x) estimatePVCA(),
                mc.cores = 25)

res.df <- do.call(rbind,res )
ratios <- res.df[row.names(res.df)%in% 'Ratio',]
prenorm <- res.df[row.names(res.df)%in% 'Pre-norm',]
postnorm <- res.df[row.names(res.df)%in% 'Normalized',]




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
    


    threshold <- 0.2
    interaction=FALSE

    estimatePVCA <- function(ncells ,
                             gdata_logct,
                             gdata_norm 
                            ){
        idx <- sample(ncol(gdata_logct),ncells)
        pvca_logct<-runPVCA(gdata=gdata_logct[,idx], mdata=mdata, sid=sid, factors = batch.factors, 0.5, interaction)
        pvca_norm<-runPVCA(gdata=gdata_norm[,idx], mdata=mdata, sid=sid, factors = batch.factors, 0.5, interaction)

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

ncells = 5000
nsims=25

pvca0 = sample_pvcaEstimates(sub.sce =sub.sce,
                             sid = 'barcodes',
                             batch.factors = c('batch_id'),
                             ncells = ncells,
                             nSims = nsims,
                             mc.cores=25
                             )
pryr::mem_used()

pvca1 = sample_pvcaEstimates(sub.sce =sub.sce,
                             sid = 'barcodes',
                             batch.factors = c('batch_id','subject.subjectGuid'),
                             ncells = ncells,
                             nSims = nsims,
                             mc.cores=25
                             )

pryr::mem_used()

pvca2 = sample_pvcaEstimates(sub.sce =sub.sce,
                             sid = 'barcodes',
                             batch.factors = c('batch_id','subject.subjectGuid','subject.biologicalSex'),
                             ncells = ncells,
                             nSims = nsims,
                             mc.cores=25
                             )
pryr::mem_used()


####################################################################
####################################################################



























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
                         
                         
                         