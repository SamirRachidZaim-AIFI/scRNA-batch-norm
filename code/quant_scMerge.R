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

quant_varPart <- function(sce_obc, frmla='~ cellType *batch'){
    require(variancePartition)
    # head(tmp)
    geneExpr_norm <- as.matrix(sce_obc@assays@data$normalized)
    geneExpr_logcpm <- as.matrix(sce_obc@assays@data$logcpm)
    info <- data.frame(sce_obc@colData)
    info$cellType <- as.factor(info$cellType)
    info$batch <- as.factor(info$batch)
    # dim(geneExpr)
    # dim(info)
    # str(info)
    # str(geneExpr)
    # head(info)
    # str(info)
    varPartfrmla <-as.formula(frmla)
    varPart_post <- fitExtractVarPartModel( geneExpr_norm, varPartfrmla, info )
    varPart_pre <- fitExtractVarPartModel( geneExpr_logcpm, varPartfrmla, info )

    vp1 <- sortCols( varPart_pre )
    vp2 <- sortCols( varPart_post )


    df= data.frame(rbind(colMeans(vp1),
                     colMeans(vp2)))

    df <- df*100

    df <- rbind(df, df[2,]/df[1,])               
    rownames(df) <- c('varPart_Post','varPart_Pre','Ratio')  
    df <- round(df,2)
    return(df)
}

####################################################################
####################################################################



cat('Batch ef')
##################  ##################  ##################
##################  ##################  ##################


# library(pvca)
source('/home/scMergeNormalization/scRNA-batch-norm/code/pvca.R')

gdata_logcpm <- data.frame(assay(sce_obc,'logcpm'))
gdata_norm <- data.frame(assay(sce_object,'normalized'))
head(gdata)
mdata <- info <- data.frame(sce_object@colData)

head(mdata)
names(mdata)
names(gdata_logcpm) <- gsub('X','', names(gdata_logcpm))
names(gdata_norm) <- gsub('X','', names(gdata_norm))

sid <- 'barcodes'
batch.factors <- c('cellType','batch')
threshold <- 0.8
interaction=T

pvcaObj_logcpm <-runPVCA(gdata=gdata_logcpm, mdata=mdata, sid=sid, factors = batch.factors, 0.8, T)

pvcaObj_norm<-runPVCA(gdata=gdata_norm, mdata=mdata, sid=sid, factors = batch.factors, 0.8, T)


pvca_ratios <- data.frame(rbind(pvcaObj_logcpm,
                               pvcaObj_norm)
                         )

pvca_ratios <- rbind(pvca_ratios, 
                    pvca_ratios[2,] / pvca_ratios[1,]
                        
                    )





##################  ##################  ##################
##################  ##################  ##################
library(scMerge)
require(sparseMatrixStats)
data('segList_ensemblGeneID', package = 'scMerge')
ctl <- segList_ensemblGeneID$human$human_scSEG
ctl <- intersect(row.names(sce), ctl)
smp_ctl <- sample(ctl, 10)

hvg <- row.names(sce)[!row.names(sce)%in% ctl]

smp_hvg <- sample(hvg,10)

### quantify variance reduction on stable genes 
cor_across_genes <- function(sce, smp_ctl, assay='normalized'){

    batch1 <- sce[, sce$batch=='X001']
    batch2 <- sce[, sce$batch=='X002']

    b1_df <-batch1@assays@data[assay][[1]]
    b2_df <-batch2@assays@data[assay][[1]]

   
    sm_df1 <- b1_df[row.names(b1_df)%in% smp_ctl,]
    sm_df2 <- b2_df[row.names(b2_df)%in% smp_ctl,]

    return(cor(rowMeans2(sm_df1), rowMeans2(sm_df2)))
}
