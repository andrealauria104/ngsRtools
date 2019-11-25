# Quality Controls ====

GrepStats <- function(f, pipeline="Hisat2") {
  
  # Pipeline Hisat2
  if( pipeline=="Hisat2" ) {
    grp <- paste0("grep 'reads\\|alignment\\|aligned' ",f,"  | cut -d ' ' -f1,5")  
  } else {
    stop(message("[-] Please, provide a valid RNA-seq pipeline."))
  }
  
}

AlignStat <- function(STATDIR, ...) {
  
  fl <- list.files(STATDIR
                   , pattern = "txt"
                   , full.names = T)
  
  grps <- vapply(fl, GrepStats, FUN.VALUE = character(1), ...)
  
  tmp <- lapply(grps, system, intern=T)
  names(tmp) <- basename(fl)
  tmp <- do.call(rbind, tmp)
  
  nm <- rownames(tmp)
  
  tmp <- gsub("%| ","",tmp)
  tmp <- as.data.frame(apply(tmp, 2, as.numeric))
  colnames(tmp) <- c("library_size","unmapped","uniquely_map","multi_map","mapping_rate")
  tmp$sample <- gsub("\\.txt","",nm)
  
  return(tmp)
  
}
