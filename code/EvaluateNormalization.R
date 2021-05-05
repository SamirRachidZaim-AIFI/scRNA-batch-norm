################################################################################################
################################################################################################

# Title: EvaluateNormalization 
# Creator: Monica Chaudhari
# Last Modified: Samir Rachid Zaim 
# Date Modified: 5/05/21

# Desc: 
#     This file provides the script to evaluate the effects of ruv-iii correction on the bridging controls. 

#     Technical Details of scMerge are found here: https://www.pnas.org/content/116/20/9775
#     The bioconductor R package is found here: https://bioconductor.org/packages/release/bioc/html/scMerge.html

# Inputs
#     tmp = sce object containing normalized data
#     assayName= string containing assayname to normalize
#     batchColname= string containing column name with batch info
#     comb_bridg_corr= 2-d list object containing output from ruviii_bridging_control function
#            where the 1st component are the normalized bridging controls, and the 2nd are the parametrs for RUVIII

# Output
#     graphical summaries of the normalization

################################################################################################
################################################################################################

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

evaluateNormalization <- function(tmp, 
                                  prenorm_assayName='logcounts', 
                                  postnorm_assayName='normalized', 
                                  cellType_name = 'seurate_pmbc_type',
                                  batchColName = 'batch_id',
                                  numComponents=50
                                 ){
    
    BC.prior = scater::runPCA(tmp, exprs_values = prenorm_assayName,ncomponents=5); 
    #reducedDimNames(sceLP); head(reducedDim(sceLP)); dim(reducedDim(sceLP))
    PCs.prior<-data.frame(cbind(reducedDim(BC.prior),
                                "cellType"=colData(tmp)[,cellType_name],
                                "batch"=colData(tmp)[,batchColName]))
    
    PCs.prior[,1:5]<-apply(PCs.prior[,1:5],2,function(x){as.numeric(as.character(x))})
    op <- options(contrasts = c("contr.helmert", "contr.poly"))
    maov.prior<-manova(cbind(PC1,PC2,PC3,PC4,PC5) ~ as.factor(cellType)*as.factor(batch),data=PCs.prior)
    summary(maov.prior)

    outfile = "BCPrior.pca1.pdf"
    p<-scater::plotPCA(BC.prior, colour_by = "cellType", shape_by = batchColName)
    pdf(file = outfile, width = 7, height = 5)
    print(p)
    dev.off()  
    outfile = "BCPrior.pca2.pdf"
    p<-scater::plotPCA(BC.prior, colour_by = "batch")
    pdf(file = outfile, width = 7, height = 5)
    print(p)
    dev.off() 


    ## POST 
    BC.post = scater::runPCA(tmp, exprs_values = postnorm_assayName,ncomponents=5); 
    PCs.post<-data.frame(cbind(reducedDim(BC.post),"cellType"=colData(tmp)$cellType,"batch"=colData(tmp)$batch))
    PCs.post[,1:5]<-apply(PCs.post[,1:5],2,function(x){as.numeric(as.character(x))})
    op <- options(contrasts = c("contr.helmert", "contr.poly"))
    maov.post<-manova(cbind(PC1,PC2,PC3,PC4,PC5) ~ as.factor(cellType)*as.factor(batch),data=PCs.post)
    summary(maov.post)

    outfile = "BCPost.pca1.pdf"
    p<-scater::plotPCA(BC.post, colour_by = "cellType", shape_by = "batch")
    pdf(file = outfile, width = 7, height = 5)
    print(p)
    dev.off()   
    outfile = "BCPost.pca2.pdf"
    p<-scater::plotPCA(BC.post, colour_by = "batch")
    pdf(file = outfile, width = 7, height = 5)
    print(p)
    dev.off() 

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

    ###### T-SNE



    set.seed(1000)
    BC.prior <- scater::runTSNE(BC.prior, perplexity=50, dimred="PCA", exprs_values = "logcpm",  n_dimred=5)

    outfile = "BCPrior.tsne1.pdf" 
    p<-scater::plotTSNE(BC.prior, colour_by = "cellType", shape_by = "batch")
    pdf(file = outfile, width = 7, height = 5)
    print(p)
    dev.off()
    outfile = "BCPrior.tsne2.pdf" 
    p<-scater::plotTSNE(BC.prior, colour_by = "batch")
    pdf(file = outfile, width = 7, height = 5)
    print(p)
    dev.off()

    BC.post <- scater::runTSNE(BC.post, perplexity=50, dimred="PCA", exprs_values = "normalized",  n_dimred=5)
    outfile = "BCPost.tsne1.pdf" 
    p<-scater::plotTSNE(BC.post, colour_by = "cellType", shape_by = "batch")
    pdf(file = outfile, width = 7, height = 5)
    print(p)
    dev.off()
    outfile = "BCPost.tsne2.pdf" 
    p<-scater::plotTSNE(BC.post, colour_by = "batch")
    pdf(file = outfile, width = 7, height = 5)
    print(p)
    dev.off()

    ##################### LP: Run UMAP - PRIOR/POST #####################################
    BC.prior <- scater::runUMAP(BC.prior,  exprs_values = "logcpm",  n_dimred=5)

    outfile = "BCPrior.umap1.pdf" 
    p<-scater::plotUMAP(BC.prior, colour_by = "cellType", shape_by = "batch")
    pdf(file = outfile, width = 7, height = 5)
    print(p)
    dev.off()
    outfile = "BCPrior.umap2.pdf" 
    p<-scater::plotUMAP(BC.prior, colour_by = "batch")
    pdf(file = outfile, width = 7, height = 5)
    print(p)
    dev.off()

    BC.post <- scater::runUMAP(BC.post, exprs_values = "normalized", use_dimred="PCA",  n_dimred=5) 

    outfile = "BCPost.umap1.pdf" 
    p<-scater::plotUMAP(BC.post, colour_by = "cellType", shape_by = "batch")
    pdf(file = outfile, width = 7, height = 5)
    print(p)
    dev.off()
    outfile = "BCPost.umap2.pdf" 
    p<-scater::plotUMAP(BC.post, colour_by = "batch")
    pdf(file = outfile, width = 7, height = 5)
    print(p)
    dev.off()
}





#######################################  VARIATION ACCOUNTED BY EACH CELL TYPE #####################################################
#######################################  Each cell of the "Terms:" divided by the sum of the corresponding row #####################
#######################################  gives an unbiased estimate of the proportion of variance explained by #####################
#######################################  that predictor in the model.###############################################################
####################################################################################################################################
cellType=unique(colData(tmp)$cellType)

PCs.prior.wide<-PCs.prior
PCs.post.wide<-PCs.post
for(i in 1:length(cellType)){
  k<-ncol(PCs.prior.wide)+1
  PCs.prior.wide[,paste(trimws(cellType[i]))]<-PCs.post.wide[,paste(trimws(cellType[i]))]<-0
  ct.idx<-which(PCs.prior.wide$cellType==cellType[i])
  PCs.prior.wide[ct.idx,paste0(cellType[i])]<-PCs.post.wide[ct.idx,paste0(cellType[i])]<-1
}

colnames(PCs.prior.wide) <- gsub(" ", "", colnames(PCs.prior.wide))
colnames(PCs.prior.wide) <- gsub("-", "", colnames(PCs.prior.wide))
colnames(PCs.prior.wide)[colnames(PCs.prior.wide) %in% c("CD14+Monocytes","CD16+Monocytes")] <-c("CD14Monocytes" , "CD16Monocytes")
colnames(PCs.post.wide) <- gsub(" ", "", colnames(PCs.post.wide))
colnames(PCs.post.wide) <- gsub("-", "", colnames(PCs.post.wide))
colnames(PCs.post.wide)[colnames(PCs.post.wide) %in% c("CD14+Monocytes","CD16+Monocytes")] <-c("CD14Monocytes" , "CD16Monocytes")

op <- options(contrasts = c("contr.helmert", "contr.poly"))
frmla<-formula(paste("cbind(PC1,PC2,PC3,PC4,PC5) ~ ", paste0(colnames(PCs.prior.wide)[9:ncol(PCs.prior.wide)], collapse = " + ")," + batch"))
maov.prior<-manova(frmla,data=PCs.prior.wide)
summary(maov.prior)
frmla<-formula(paste("cbind(PC1,PC2,PC3,PC4,PC5) ~ ", paste0(colnames(PCs.post.wide)[9:ncol(PCs.post.wide)], collapse = " + ")," + batch"))
maov.post<-manova(frmla,data=PCs.post.wide)
summary(maov.post, test='Pillai')

frmla2 <- 'cbind(PC1,PC2,PC3,PC4,PC5)~ Platelets'
maov.post<-manova(frmla2,data=PCs.post.wide)
summary(maov.post)



    
    
}

