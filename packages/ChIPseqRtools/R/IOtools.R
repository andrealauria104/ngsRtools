# I/O tools ----
write_peaks_bed <- function(peaks, bedfile = NULL)
{
  options(scipen=999) # disable scientific notation
  
  message("[+] Writing peaks to file: ", bedfile)
  
  tobed <- cbind('chr' = peaks[,1], format(peaks[,2:3], scientific = FALSE), 'peak' = peaks[,4])
  tobed[,2:3] <- lapply(tobed[,2:3], as.numeric)
  
  if( is.null(bedfile) ) {
    stop(message("[!] Please, provide a path to output bed file."))
  }
  
  write.table(tobed,
              quote = F      ,
              row.names = F  ,
              col.names = F  ,
              sep = "\t",
              file = bedfile)
  
  options(scipen=0) # restore default
}
