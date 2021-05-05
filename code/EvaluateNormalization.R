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

evaluateNormalization <- function(tmp, 
                                  prenorm_assayName='logcounts', 
                                  postnorm_assayName='normalized', 
                                  cellType_name = 'seurate_pmbc_type',
                                  batchColName = 'batch_id',
                                  numComponents=50
                                 ){

    ################################################################################
    ################################################################################
    
    cat('Calculating principal components for visualization and downstream analyses')
    
    #### Calculate Principal components based on pre-batch-correction
    BC.prior = scater::runPCA(tmp, exprs_values = 
                              prenorm_assayName,
                              ncomponents=numComponents); 
    PCs.prior<-data.frame(cbind(reducedDim(BC.prior),
                                "cellType"=colData(tmp)[,cellType_name],
                                "batch"=colData(tmp)[,batchColName]))
    
    PCs.prior[,1:numComponents]<-apply(PCs.prior[,1:NumComponents],
                                       2,
                                       function(x){as.numeric(as.character(x))}
                                      )
    options(repr.plot.width=10, repr.plot.height=10)
    p1<-scater::plotPCA(BC.prior, colour_by = "cellType", shape_by = "batch")
    p2<-scater::plotPCA(BC.prior, colour_by = "batch")
    ggarrange(p1,p2, ncol=1)
    
    #### Calculate Principal components based on pre-batch-correction

    BC.post = scater::runPCA(tmp, 
                             exprs_values = postnorm_assayName,
                             ncomponents=NumComponents); 
    
    PCs.post<-data.frame(cbind(reducedDim(BC.post),
                               "cellType"=colData(tmp)[,cellType_name],
                                "batch"=colData(tmp)[,batchColName]))
    
    PCs.post[,1:NumComponents]<-apply(PCs.post[,1:NumComponents],2,function(x){as.numeric(as.character(x))})
    options(repr.plot.width=10, repr.plot.height=10)

    p1<-scater::plotPCA(BC.post, colour_by = "cellType", shape_by = "batch")
    p2<-scater::plotPCA(BC.post, colour_by = "batch")
    ggarrange(p1,p2, ncol=1)
    
    ################################################################################
    ################################################################################
   
    
    ################################################################################
    ################################################################################

    cat('\n\nCalculating T-SNEs')

    ###### T-SNE
    set.seed(1000)
    BC.prior <- scater::runTSNE(BC.prior, perplexity=50, 
                                dimred="PCA", 
                                exprs_values = prenorm_assayName,  
                                n_dimred=numComponents)

    p1<-scater::plotTSNE(BC.prior, colour_by = "cellType", shape_by = "batch")
    p2<-scater::plotTSNE(BC.prior, colour_by = "batch")

    BC.post <- scater::runTSNE(BC.post, perplexity=50, 
                               dimred="PCA", 
                               exprs_values = postnorm_assayName,  
                               n_dimred=numComponents)

    p3<-scater::plotTSNE(BC.post, colour_by = "cellType", shape_by = "batch")
    p4<-scater::plotTSNE(BC.post, colour_by = "batch")
    
    ## Pre normalization
    figure1 <- ggpubr::ggarrange(p1,p2, ncol=2) 
    annotate_figure(figure1,
                    top = text_grob("T-SNE: Before Normalization", 
                                    color = "black", face = "bold", 
                                    size = 24))

    figure2 <- ggpubr::ggarrange(p3,p4, ncol=2) 
    annotate_figure(figure2,
                    top = text_grob("T-SNE: Post Normalization", 
                                    color = "black", face = "bold",
                                    size = 24))
    
    cat('\n\nfinished calculating T-snes')
    ################################################################################
    ################################################################################

        
    ################################################################################
    ################################################################################
    cat('\n\ncalculating umap')
       
    BC.prior <- scater::runUMAP(BC.prior,  
                                exprs_values = prenorm_assayName,    
                                n_dimred=NumComponents)

    p1<-scater::plotUMAP(BC.prior, colour_by = "cellType", shape_by = "batch")+ggtitle('Celltype X Batch')
    p2<-scater::plotUMAP(BC.prior, colour_by = "batch") + ggtitle('Batch')

    BC.post <- scater::runUMAP(BC.post, 
                               exprs_values = postnorm_assayName, 
                               n_dimred=NumComponents)
    
    p3<-scater::plotUMAP(BC.post, colour_by = "cellType", shape_by = "batch")
    p4<-scater::plotUMAP(BC.post, colour_by = "batch")
    
    figure1 <- ggpubr::ggarrange(p1,p2, ncol=2) 
    annotate_figure(figure1,
                    top = text_grob("UMAP: Before Normalization", 
                                    color = "black", face = "bold", 
                                    size = 24))

    cat('\n\n')
    figure2 <- ggpubr::ggarrange(p3,p4, ncol=2) 
    annotate_figure(figure2,
                    top = text_grob("UMAP: Post Normalization", 
                                    color = "black", face = "bold", 
                                    size = 24))
    
    ################################################################################
    ################################################################################
 }
