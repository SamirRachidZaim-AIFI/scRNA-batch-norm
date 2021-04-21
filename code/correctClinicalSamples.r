## Load data
options(warn=-1)

sce<-readRDS("../../data/sceFltrd_kpGns.rds")
sce<-sce[,sce$subjid!="IMM19"]
sce$subjid<-factor(sce$subjid,levels=sort(unique(sce$subjid)))
subSmpl.lst<-readRDS("../../data/subSmpl.lst.b15p10.rds")
print('data loaded')
fres <- readRDS('../../data/tmp_scMerge.RDS')



## Load libraries 
options(warn=-1)
source('../code/scMerge_helperFunctions.R')
library(SingleCellExperiment)
require(scMerge)
require(BiocParallel)
require(ggpubr)

length(subSmpl.lst) # length(subSmpl.lst): 6  ; dim(subSmpl.lst[[1]]) (bootstrap x cells: 15 x 7353)
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
sub.sce <- sce[which(rownames(sce) %in% rownames(tmp)) , which(sce$wk==wks[d])[idx] ]
ctl.gn <- intersect(cmSEGs, row.names(sub.sce))

ac = ac[,colnames(ac) %in% ctl.gn]
    
scale_res <- standardize2_mod(Y=assay(sub.sce,"logcpm"), batch=sub.sce$batch)
stand_tY <- DelayedArray::t(scale_res$stand_Y)
stand_sd <- sqrt(scale_res$stand_var)
stand_mean <- scale_res$stand_mean


ac_t = DelayedArray::as.matrix(ac)
ac = t(ac_t)

W <- stand_tY[, ctl.gn] %*% ac_t %*% solve(ac %*% ac_t)
noise = W %*% alpha
newY_mc <- stand_tY - noise[, which(colnames(noise) %in% colnames(stand_tY) )]

## Add back the mean and sd to the normalised data
newY_mc <-(t(newY_mc) * stand_sd) + stand_mean
head(newY_mc)
assay(sub.sce,"normalized")<-newY_mc

cat('\nClinical samples have been normalized for week 2\n\n')
sub.sce
cat('\n\n')
head(sub.sce@assays@data$normalized[1:10,1:10])

##################### LP: Run PCA - PRIOR/POST #####################################
BC.prior.smpl = scater::runPCA(sub.sce, exprs_values = "logcpm",ncomponents=5); 
#reducedDimNames(sceLP); head(reducedDim(sceLP)); dim(reducedDim(sceLP))
PCs.prior.smpl<data.frame(cbind(reducedDim(BC.prior.smpl),
                        "cellType"=colData(sub.sce)$cellType,"batch"=colData(sub.sce)$batch))
PCs.prior.smpl[,1:5]<-apply(PCs.prior.smpl[,1:5],2,function(x){as.numeric(as.character(x))})
op <- options(contrasts = c("contr.helmert", "contr.poly"))
maov.prior.smpl<-manova(cbind(PC1,PC2,PC3,PC4,PC5) ~
                        as.factor(cellType)*as.factor(batch),data=PCs.prior.smpl)   

summary(maov.prior.smpl)

fname <- paste('../../plots/clinical_PCAs.pdf')
pdf(fname) 
p1<-scater::plotPCA(BC.prior.smpl, colour_by = "cellType", shape_by = "batch")
p2<-scater::plotPCA(BC.prior.smpl, colour_by = "batch")

ggarrange(p1,p2, ncol=1)

BC.post.smpl = scater::runPCA(sub.sce, exprs_values = "normalized",ncomponents=5); 
PCs.post.smpl<-data.frame(cbind(reducedDim(BC.post.smpl),
                                    "cellType"=colData(sub.sce)$cellType,
                                    "batch"=colData(sub.sce)$batch))
PCs.post.smpl[,1:5]<-apply(PCs.post.smpl[,1:5],2,function(x){as.numeric(as.character(x))})
op <- options(contrasts = c("contr.helmert", "contr.poly"))
maov.post.smpl<-manova(cbind(PC1,PC2,PC3,PC4,PC5) ~ as.factor(cellType)*as.factor(batch),
                           data=PCs.post.smpl)
summary(maov.post.smpl)
    
p1<-scater::plotPCA(BC.post.smpl, colour_by = "cellType", shape_by = "batch")
p2<-scater::plotPCA(BC.post.smpl, colour_by = "batch")
     
ggarrange(p1,p2, ncol=1) 
dev.off()

    ##################### LP: Run TSNE - PRIOR/POST #####################################
fname <- paste('../../plots/clinical_tsne.pdf')
pdf(fname) 
set.seed(1000)
BC.prior.smpl <- scater::runTSNE(BC.prior.smpl, perplexity=50, 
                                 dimred="PCA", exprs_values ="logcpm",  n_dimred=5)
    
p1<-scater::plotTSNE(BC.prior.smpl, colour_by = "cellType", shape_by = "batch")
p2<-scater::plotTSNE(BC.prior.smpl, colour_by = "batch")
ggarrange(p1,p2, ncol=1) 


BC.post.smpl <- scater::runTSNE(BC.post.smpl, perplexity=50, dimred="PCA", exprs_values = "normalized",  n_dimred=5)
p3<-scater::plotTSNE(BC.post.smpl, colour_by = "cellType", shape_by = "batch")
p4<-scater::plotTSNE(BC.post.smpl, colour_by = "batch")

ggarrange(p3,p4, ncol=1) 
dev.off()
    
##################### LP: Run UMAP - PRIOR/POST #####################################
BC.prior.smpl <- scater::runUMAP(BC.prior.smpl,  exprs_values = "logcpm",  
                                     n_dimred=5)
BC.post.smpl <- scater::runUMAP(BC.post.smpl, exprs_values = "normalized", 
                                     n_dimred=5)    

outfile = "../../plots/clinical_UMAPs.pdf" 
pdf(file = outfile, width = 7, height = 5)
p1<-scater::plotUMAP(BC.prior.smpl, colour_by = "cellType", shape_by = "batch")
p2<-scater::plotUMAP(BC.prior.smpl, colour_by = "batch")
ggarrange(p1,p2, ncol=1) 
   
p3<-scater::plotUMAP(BC.post.smpl, colour_by = "cellType", shape_by = "batch")
p4<-scater::plotUMAP(BC.post.smpl, colour_by = "batch")
ggarrange(p3,p4, ncol=1) 

dev.off()






