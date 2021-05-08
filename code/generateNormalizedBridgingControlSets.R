 #### create various bridging control sets

####################################################################
####################################################################

# The deep & wide analysis is a study to explore how to 
# anchor scRNA batch correction on the deep & wide
# replicates. The analysis consists of assessing which of 
# the following approaches leads to the best batch correction
# results across elements in different batches: 
# - anchoring deep & wide 
# - tmp_

### load deep & wide bridging controls 
source('/home/jupyter/scRNA-batch-norm/code/scMerge_helperFunctions.R')
source('/home/jupyter/scRNA-batch-norm/code/ruviii_bridging_controls.R')
source('/home/jupyter/scRNA-batch-norm/code/quant_scMerge.R')
require(scMerge)
require(BiocParallel)
require(Seurat)
setwd('/home/jupyter/deepNwide/')
require(H5weaver)
require(SingleCellExperiment)

####################################################################
####################################################################

fnames <- dir('h5/')
fnames <- fnames[-grep('multiplet',fnames)]

deep_wide_ctls <- lapply(fnames,
                function(x) 
                    H5weaver::read_h5_sce(paste('h5/',x,sep='')
                                         )
                )


deep_ctls <- lapply(1:4, function(x) deep_wide_ctls[[x]]
               )


####################################################################
####################################################################
                

####################################################################
####################################################################

### load br1 bridging controls 
br1Lookup <- read.csv('/home/jupyter/scRNA-batch-norm/bri_analysis/lookupTable_br1.csv')

# Since the bridging controls used in the deep & wide study correspond 
# to batch 4, the below filters out for non-batch 4 bridging control
# in order to get a diversity of batch behavior. 

bridge_other_batch_idx <- which(br1Lookup$subject.subjectGuid=='' & br1Lookup$file.batchID != 'B004')
bridging_meta <- br1Lookup[bridge_other_batch_idx,]
setwd('/home/jupyter/scRNA-batch-norm/bri_analysis')
fnames <- bridging_meta$filePath
br1_ctls <- lapply(fnames, function(x) read_h5_sce(x)
                   )

####################################################################
####################################################################

# Combine data from different batches and run scMerge normalization
# will result in a superset of bridging controls cells to sample from

comb_bridg <- cbind(do.call(cbind, br1_ctls),
                    do.call(cbind, deep_ctls)
                    )


# To run scMerge you must log normalize the raw count data 
# obtained in the H5 files. One option is using scater's logNormCounts
# function to log normalize the raw counts. 

comb_bridg <- scater::logNormCounts(comb_bridg)

####################################################################
####################################################################



####################################################################
####################################################################

# To evaluate the quality of the normalization we will run two sets 
# of experiments: 
# 1. naive subsampling of all cells across the corpora of bridging controls
# 2. cell-type subsampling to ensure cell representation in the final bridging control set


## Experiment 1: naively subsample 5k, 10k, 15k, 20k, 30k cells from bridging controls
## since deep & wide = 60% of data, this effectively "anchors" normalization
## on the deep & wide via majority sampling

bridging_control_sizes <- c(5000, 10000,15000, 20000, 30000)

tmp_list <- lapply(bridging_control_sizes, 
                   function(x) comb_bridg[ ,comb_bridg$barcodes %in% sample(comb_bridg$barcodes, x)]
                   )

brdg_ruviii_list <- mclapply(tmp_list, 
                             function(x) ruviii_bridging_controls(x),
                             mc.cores = detectCores()-10
                             )

mclapply(1:length(bridging_control_sizes),
        function(x) saveRDS(brdg_ruviii_list[[x]], 
                            paste('/home/jupyter/scRNA-batch-norm/bridging_controls_ruviii/ruviii_',
                                  bridging_control_sizes[x],
                                  'cells.RDS', sep='')
                            )
         
        )


# Experiment 2: subsample 5k, 10k, 15k, 20k, 30k cells from bridging controls
#     but oversample rare cell types in order to ensure they are equally represented
#     across rare cell types. 

# Experiment 3: look at using KS test to batch correct individual cell types 
#     that are not well represented. (I.e., if PDC's dist'n can be approx via
#     B cells, then we can batch correct PDCs and B cells together). 
#     Need to look at UMAP 

#     Would require pairwise KS-test between rare cells and major cell types 
#     to see which "pair" is 'closest' distributionally, and then compare batch
#     correction on rare cells using their own, vs. an auxiliary set to help. 

         
         