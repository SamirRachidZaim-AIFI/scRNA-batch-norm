################################################################################################
################################################################################################

# Title: normalizeClinicalSamples 
# Creator: Monica Chaudhari
# Last Modified: Samir Rachid Zaim 
# Date Modified: 4/19/21

# Desc: 
#     This file provides the script to do an ruv-iii correction on the bridging controls. 

#     Technical Details of scMerge are found here: https://www.pnas.org/content/116/20/9775
#     The bioconductor R package is found here: https://bioconductor.org/packages/release/bioc/html/scMerge.html

# Inputs
#     sub.sce = sce object containing clinical samples to normalize
#     assayName= string containing assayname to normalize
#     batchColname= string containing column name with batch info
#     comb_bridg_corr= 2-d list object containing output from ruviii_bridging_control function
#            where the 1st component are the normalized bridging controls, and the 2nd are the parametrs for RUVIII

# Output
#     sub.sce = returns same input sce object with an added normalized assay

################################################################################################
################################################################################################

normalizeClinicalSamples <- function(sub.sce,
                                     assayName='logcounts',
                                     batchColName='batch_id',
                                     comb_bridg_corr){
    
        ### USE scMERGE common genes 
        data("segList_ensemblGeneID", package = "scMerge")
        seg_index <- segList_ensemblGeneID$human$human_scSEG
        symbols <- mapIds(org.Hs.eg.db, keys = seg_index, keytype = "ENSEMBL", column="SYMBOL")

        ctl.gn <- Reduce(intersect, list(row.names(sub.sce), 
                                         unique(symbols),
                                         row.names(comb_bridg_corr$normed_sce)
                                        )
                         )
        ## make sure control genes & SEGs from scMerge present in data
#         ctl.gn <- intersect(ctl.gn, row.names(sub.sce))
#         ctl.gn <- intersect(ctl.gn, row.names(comb_bridg_corr$normed_sce))
    
        ### extract objects for scMerge normalization
        k<-comb_bridg_corr[[2]]$optimal_ruvK ## k=1
        falpha<-comb_bridg_corr[[2]][[k]]$fullalpha # dim(falpha) 52 x 20240 (k x G)
        colnames(falpha)<-rownames(comb_bridg_corr[[1]])
        alpha <- falpha[seq_len(min(k, nrow(falpha))), , drop = FALSE] # dim(alpha) : 1 x 20240
        ac <- alpha[,ctl.gn, drop = FALSE] # dim(ac): 1 x 52; sum(is.na(ac[1,]))
    
        ## Remove any rows & columns that are all 0s
        sub.sce<-sub.sce[which(rowSums(assay(sub.sce,2)) != 0), 
                         which(colSums(assay(sub.sce,2)) != 0)]

       

        ## subset alpha_c (unwanted variation vectors)
        ## contain the control genes to correct against
        ac = ac[,colnames(ac) %in% ctl.gn]

        ## Standardize results by batch 
        scale_res <- standardize2_mod(Y=assay(sub.sce,assayName), 
                                      batch=colData(sub.sce)[,batchColName])
        stand_tY <- DelayedArray::t(scale_res$stand_Y)
        stand_sd <- sqrt(scale_res$stand_var)
        stand_mean <- scale_res$stand_mean
    
        ## Manually apply RUV-III by estimating W and then Y_hat
        ##  --> W = Y'a_c(a_c a_c')-1
        ##  --> Y_hat = Y-W_alpha
        Y_c = stand_tY[, ctl.gn]
        ac_t = t(ac)

        W <-  Y_c %*% ac_t %*% solve(ac %*% ac_t)
        noise = W %*% alpha
    
        cmnGenes <- intersect(colnames(noise),colnames(stand_tY))
    
        newY_mc <- stand_tY
        newY_mc[,cmnGenes] <- stand_tY[,cmnGenes] - noise[, cmnGenes]

        ## Add back the mean and sd to the normalised data
        newY_mc <-(t(newY_mc) * stand_sd) + stand_mean

        ## add normalized assay onto sub.sce object (summarized exp object)
        assay(sub.sce,"normalized")<-newY_mc

        cat('\nClinical samples have been normalized\n\n')
        return(sub.sce)
}
