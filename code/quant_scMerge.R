################################################################################################
################################################################################################

# Title: quant_scMerge 
# Author: Samir Rachid Zaim
# Last Modified: Samir Rachid Zaim 

# Desc: 
#     This file contain a list of auxiliary helper functions to quantify scMerge batch normalization 

################################################################################################
################################################################################################

####################################################################
####################################################################
###
###

# Name quant_varPart

# Desc 
#     variancePartition is a bioconductor R package that 
#     runs a linear model to quantify the proportion of variance 
#     explained by each of the components. The quant_varPart is a 
#     wrapper that quantifies the ratio of "conserved" variance
#     pre & post normalization using variancePartition

# Input
#     sce_obc a single-cell experiment object containing both logcpm 
#             and normalized assays 
#     frmla a string containing the formula to decompose the variance 

# Output
#     varianceRatios_Df containing the pre & post nrmalization variance 
#     decompositions as well as their ratio, to quantify the % of variance
#     reduction after normalizing the dat.a

quant_varPart <- function(sce_obc, 
                          frmla='~ n_genes + n_reads + 1|batch_id',
                          numGenesToAnalyze=100,
                          numCells =5000,
                          param
                           ){

#     library(BiocParallel)
#     register(SnowParam(25))  
    idx <- sample(unique(row.names(sce_obc)),numGenesToAnalyze)
    idx2<- sample(unique(sce_obc$barcodes), numCells)
    sce_obc <- sce_obc[idx, sce_obc$barcodes %in% idx2]
    geneExpr_norm <- as.matrix(sce_obc@assays@data$normalized)
    geneExpr_logcpm <- as.matrix(sce_obc@assays@data$logcounts)
    info <- data.frame(sce_obc@colData)
    varPartfrmla <-as.formula(frmla)
    
    varPart_post <- fitExtractVarPartModel( geneExpr_norm, 
                                           varPartfrmla, 
                                           info,
                                           BPPARAM=param
                                          )
    varPart_pre <- fitExtractVarPartModel(geneExpr_logcpm, 
                                          varPartfrmla,
                                          info,
                                         BPPARAM=param)

    vp_pre <- sortCols( varPart_pre )
    vp_post <- sortCols( varPart_post )


    df= data.frame(rbind(colMeans(vp_pre),
                     colMeans(vp_post)))

    df <- df*100

    df <- rbind(df, df[1,]/df[2,])               
    rownames(df) <- c('varPart_Pre','varPart_Post','Ratio')  
    df <- round(df,2)
    return(df)
}

####################################################################
####################################################################

####################################################################
####################################################################
###
###

# Name sample_cellTypes

# Desc 
#     This function allows you to sample cells by celltype using the
#     colData (metadata) data frame in an sce object using the 
#     sampling R package. This function allows you to both sample
#     cell types proportionally (i.e., stratified) or with 
#     even-counts. 

# Input
#     sce_obc a single-cell experiment object containing both logcpm 
#             and normalized assays 
#     size integer, total number of cells
#     stratified boolean, do you want celltypes stratified prop. or evenly

# Output
#     subsetted metadata data.frame

sample_cellTypes <- function(sce,
                             size=100, 
                             stratified=FALSE,
                             cellType_name='seurat_pbmc_type'
                            ){
    require(sampling)
    
    probs <- as.data.frame(table(sce@colData[cellType_name]))
    colnames(probs) <- c(cellType_name,'Frequency')
    probs$prob <- probs$Frequency/sum(probs$Frequency)
    
    prob_vec <- round(probs$prob *size,0)+round(size/1000,0) #ensures .1% of rare cell types (prevents 0s) 
    
    if(stratified){
       
        
        stratified_df <- strata(sce@colData, cellType_name, 
                                size=prob_vec, method='srswor')

    } else {
        stratified_df <- strata(sce@colData, cellType_name, 
                                size=rep(round(size/13),length(unique(sce@colData[cellType_name]))), 
                                method='srswor')
    }
    new_metadf <- getdata(sce@colData, stratified_df)
    return(sce[,sce$barcodes %in% new_metadf$barcodes])
}

# break


# cat('Batch ef')
# ##################  ##################  ##################
# ##################  ##################  ##################


# # library(pvca)
# source('/home/scMergeNormalization/scRNA-batch-norm/code/pvca.R')

# gdata_logcpm <- data.frame(assay(sce_obc,'logcpm'))
# gdata_norm <- data.frame(assay(sce_object,'normalized'))
# head(gdata)
# mdata <- info <- data.frame(sce_object@colData)

# head(mdata)
# names(mdata)
# names(gdata_logcpm) <- gsub('X','', names(gdata_logcpm))
# names(gdata_norm) <- gsub('X','', names(gdata_norm))

# sid <- 'barcodes'
# batch.factors <- c('cellType','batch')
# threshold <- 0.8
# interaction=T

# pvcaObj_logcpm <-runPVCA(gdata=gdata_logcpm, mdata=mdata, sid=sid, factors = batch.factors, 0.8, T)

# pvcaObj_norm<-runPVCA(gdata=gdata_norm, mdata=mdata, sid=sid, factors = batch.factors, 0.8, T)


# pvca_ratios <- data.frame(rbind(pvcaObj_logcpm,
#                                pvcaObj_norm)
#                          )

# pvca_ratios <- rbind(pvca_ratios, 
#                     pvca_ratios[2,] / pvca_ratios[1,]
                        
#                     )





# ##################  ##################  ##################
# ##################  ##################  ##################
# library(scMerge)
# require(sparseMatrixStats)
# data('segList_ensemblGeneID', package = 'scMerge')
# ctl <- segList_ensemblGeneID$human$human_scSEG
# ctl <- intersect(row.names(sce), ctl)
# smp_ctl <- sample(ctl, 10)

# hvg <- row.names(sce)[!row.names(sce)%in% ctl]

# smp_hvg <- sample(hvg,10)

# ### quantify variance reduction on stable genes 
# cor_across_genes <- function(sce, smp_ctl, assay='normalized'){

#     batch1 <- sce[, sce$batch=='X001']
#     batch2 <- sce[, sce$batch=='X002']

#     b1_df <-batch1@assays@data[assay][[1]]
#     b2_df <-batch2@assays@data[assay][[1]]

   
#     sm_df1 <- b1_df[row.names(b1_df)%in% smp_ctl,]
#     sm_df2 <- b2_df[row.names(b2_df)%in% smp_ctl,]

#     return(cor(rowMeans2(sm_df1), rowMeans2(sm_df2)))
# }
