# scRNA-batch-norm
This repository contains code and helper functions to develop &amp; deploy functionality for scRNA batch correction issues that may arise across our different studies and cohorts. The main components being studied are the roles of batch effects (by sample and by cell type) and to understand how cell composition may affect our batch effect across bridging controls with 5k vs. 20k cells. 

### Batch Effects 
- The first sets of analyses are conducted to evaluate the presence of batch effects and the role that the leukopak (bridging control) samples have in batch correction. Using the scMerge bioconductor library, the leukopak samples are used to estimate unwanted variance and batch effect, and are then removed from the clinical samples. 

### Batch Effects by Cell type
- 
### Cell Composition by Batch Effect
-TBD 

# Future
Adjusting for time correction and drifting in batch effects. 
