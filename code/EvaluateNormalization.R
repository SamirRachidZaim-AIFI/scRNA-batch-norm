### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

BC.prior = scater::runPCA(tmp, exprs_values = "logcpm",ncomponents=5); 
#reducedDimNames(sceLP); head(reducedDim(sceLP)); dim(reducedDim(sceLP))
PCs.prior<-data.frame(cbind(reducedDim(BC.prior),"cellType"=colData(tmp)$cellType,"batch"=colData(tmp)$batch))
PCs.prior[,1:5]<-apply(PCs.prior[,1:5],2,function(x){as.numeric(as.character(x))})
op <- options(contrasts = c("contr.helmert", "contr.poly"))
maov.prior<-manova(cbind(PC1,PC2,PC3,PC4,PC5) ~ as.factor(cellType)*as.factor(batch),data=PCs.prior)
summary(maov.prior)

outfile = "BCPrior.pca1.pdf"
p<-scater::plotPCA(BC.prior, colour_by = "cellType", shape_by = "batch")
pdf(file = outfile, width = 7, height = 5)
print(p)
dev.off()  
outfile = "BCPrior.pca2.pdf"
p<-scater::plotPCA(BC.prior, colour_by = "batch")
pdf(file = outfile, width = 7, height = 5)
print(p)
dev.off() 


## POST 
BC.post = scater::runPCA(tmp, exprs_values = "normalized",ncomponents=5); 
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


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### quantify unwanted variation 

seg_index <- segList_ensemblGeneID$human$human_scSEG
cmSEGs <- intersect(row.names(tmp), seg_index)

# ctl.gn <- as.character(SEGs$geneid) ; # length(ctl.gn): 52
ctl.gn <- cmSEGs

k<-ruv3res$optimal_ruvK ## k=1
assay(tmp,"normalized")<-t(ruv3res[[k]]$newY) 
falpha<-ruv3res[[k]]$fullalpha # dim(falpha) 52 x 20240 (k x G)
colnames(falpha)<-rownames(tmp)
alpha <- falpha[seq_len(min(k, nrow(falpha))), , drop = FALSE] # dim(alpha) : 1 x 20240
ac <- alpha[,ctl.gn, drop = FALSE] # dim(ac): 1 x 52; sum(is.na(ac[1,]))

########### STANDARDIZE LEUKOPAK AND GET ESTIMATES  ########################
scale_res <- standardize2_mod(Y=assay(tmp,"logcpm"), batch=tmp$batch)
stand_tY <- DelayedArray::t(scale_res$stand_Y) # dim(stand_tY) : 9699 20240
stand_sd <- sqrt(scale_res$stand_var)
stand_mean <- scale_res$stand_mean

W <- stand_tY[, ctl.gn] %*% (DelayedArray::t(ac) %*% solve(ac %*% DelayedArray::t(ac)))
newY_LP <- stand_tY - W %*% alpha
## Add back the mean and sd to the normalised data
all.equal(t((t(newY_LP) * stand_sd + stand_mean)),ruv3res[[k]]$newY) ## TRUE
newY_LP <-(t(newY_LP) * stand_sd ) + stand_mean





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


