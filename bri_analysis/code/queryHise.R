fileDescToDataframe <- function(descriptors,
                                keep_labs = FALSE) {
  
  assertthat::assert_that(typeof(descriptors) == "list")
  assertthat::assert_that(typeof(keep_labs) == "logical")
  
  do.call(
    rbind,
    lapply(
      descriptors,
      function(desc) {
        desc <- unlist(desc)
        desc <- desc[!grepl("scheme", names(desc))]
        names(desc) <- sub("^descriptors.","",names(desc))
        if(!keep_labs) {
          desc <- desc[!grepl("^lab", names(desc))]
        }
        desc <- as.list(desc)
        df <- as.data.frame(desc)
        df
      }
    )
  )
}

br1_filter_list <- list(
    cohort.cohortGuid = "BR1"
)

br1_rna_desc <- getFileDescriptors(
    fileType = "scRNA-seq-labeled", 
    filter = br1_filter_list)
    
br1_rna_desc <- fileDescToDataframe(br1_rna_desc)    