# Quality check ----
get_chipseq_qc <- function(DATADIR)
{
  qc_path <- paste0(DATADIR, "/BAM/qc_stats")
  f1 <- list.files(qc_path
                   , pattern = ".fstat.qc"
                   , full.names = T)
  
  f1 <- f1[grep("filt", f1, invert = T)]
  
  tot_reads <- unlist(lapply(f1, function(x)
  {
    as.numeric(system(paste0("awk '{if(NR==1) print $1}' ", x), intern = T))
  }))
  names(tot_reads) <- sapply(strsplit(basename(f1), "\\."), "[[", 1)
  
  f2 <- list.files(qc_path
                   , pattern = ".pbc.qc"
                   , full.names = T)
  
  pbc <- lapply(f2, read.delim2, header = F)
  names(pbc) <- sapply(strsplit(basename(f2), "\\."), "[[", 1)
  pbc <- do.call(rbind, pbc)
  pbc <- pbc[,c(1,2,5)]
  colnames(pbc) <- c("Uniquely mapped reads", "Non-redundant reads", "NRF")
  pbc$tot <- tot_reads[rownames(pbc)]
  pbc$sample <- rownames(pbc)
  pbc <- cbind.data.frame(pbc[,c(5,4)],pbc[,1:3])
  colnames(pbc)[2] <- "Sequenced reads"
  rownames(pbc) <- NULL
  pbc$uniq_mrate <- pbc[,3]/pbc[,2]
  colnames(pbc)[6] <- "Uniquely-mapping ratio"
  
  pbc <- pbc[,c(1:3,6,4:5)]
  return(pbc)
}
