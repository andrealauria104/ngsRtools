# General Utilities 
library(xlsx)

# Get OS info ====
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

# Source lines in file ====
source_lines <- function(file, start, end, ...) {
  
  file.lines <- scan(file
                     , what   = character()
                     , skip   = start-1
                     , nlines = end-start+1
                     , sep    = '\n')
  
  file.lines.collapsed <- paste(file.lines, collapse='\n')
  
  source(textConnection(file.lines.collapsed), ...)
  
}

# Read gmt ====
read.gmt.file = function(pathMsigDbFile) {
  inputFile <- pathMsigDbFile
  con  <- file(inputFile, open = "r")
  
  c = 1
  pathway.list <- vector(mode="list",length=0)
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0)
  {
    myVector <- do.call("rbind",strsplit(oneLine, "\t"))
    t = vector(mode="list",length=1)
    t[[1]] = myVector[3:length(myVector)]
    names(t) = myVector[1]
    pathway.list = c(pathway.list,t)
    c = c+1
  }
  
  close(con)
  return(pathway.list)
}

# Data summary ====
data_summary <- function(x) 
{
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# GeneId reference ====
prepare_idRef <- function(ANN){
  
  # Annotation generated from GENCODE gtf, column 9
  
  id_ref_mm9 <- read.delim2(ANN, stringsAsFactors = F,
                            header = F, sep = "\t", 
                            comment.char = "#")
  if(any(grepl("gene",id_ref_mm9$V3))) {
    id_ref_mm9 <- subset(id_ref_mm9, V3=="gene")
  } else if(any(grepl("transcript",id_ref_mm9$V3))) {
    id_ref_mm9 <- subset(id_ref_mm9, V3=="transcript")
  }
  if(ncol(id_ref_mm9)==9) id_ref_mm9 <- data.frame("V1" = id_ref_mm9$V9, stringsAsFactors = F)
  idx <- unlist(strsplit(id_ref_mm9[1,], "\\;"))
  i1 <- which(grepl("gene_id", idx))
  i2 <- which(grepl("gene_name", idx))
  i3 <- which(grepl("gene_type", idx))
  
  id_ref_mm9$gene_id   <- gsub(" ", "",gsub("gene_id|;"  ,"",sapply(strsplit(id_ref_mm9$V1, "\\;"), "[[", i1)))
  id_ref_mm9$gene_name <- gsub(" ", "",gsub("gene_name|;","",sapply(strsplit(id_ref_mm9$V1, "\\;"), "[[", i2)))
  id_ref_mm9$gene_type <- gsub(" ", "",gsub("gene_type|;","",sapply(strsplit(id_ref_mm9$V1, "\\;"), "[[", i3)))
  id_ref_mm9$V1 <- NULL
  return(id_ref_mm9)
}

# Create bedfile ====
write_bedfile <- function(tobed, bedfile) {
  
  write.table(tobed, 
              quote = F      , 
              row.names = F  , 
              col.names = F  ,
              sep = "\t",
              file = bedfile)
}

# Hypergeometric test ====
get_phyper <- function(x, y, z, w, type){
  # x = success_in_sample
  # y = success_in_bkgd
  # z = failure_in_bkgd
  # w = sample_size
  if ( type == 'enrichment'){
    phyper(q = x - 1, 
           m = y,
           n = z,
           k = w, lower.tail = F)
  } else if ( type == 'depletion') {
    phyper(q = x, 
           m = y,
           n = z,
           k = w, lower.tail = T)
  }
}
get_fisher <- function(x, y, z, w, type){
  # x = success_in_sample
  # y = success_in_bkgd - success_in_sample
  # z = sample_size - success_in_sample
  # w = failure_in_bkgd - sample_size + success_in_sample
  m <- matrix(data = c(x, y, 
                       z, w), 
              nrow = 2, ncol = 2)
  fisher.test(m, alternative = type)
  return(fisher.test(m, alternative = type)$p.value)
}

test_enrichment <- function(success_in_sample, success_in_bkgd, failure_in_bkgd, sample_size
                            , method = "hyper"
                            , retres = F) 
{
  if ( method == "hyper" ) {
    # Hypergeometric Test
    hyper <- lapply(c('enrichment', 'depletion'), 
                    function(x) get_phyper(success_in_sample,
                                           success_in_bkgd  ,
                                           failure_in_bkgd  ,
                                           sample_size      ,
                                           type = x)
    )
    names(hyper) <- c('enrichment', 'depletion')
    
    message("[*] Hypergeometric test [*]")
    message("    P-values")
    message("    Enrichment: ", formatC(hyper$enrichment, format = 'e', digits = 2))
    message("    Depletion : ", formatC(hyper$depletion, format = 'e', digits = 2))
    
    if(retres) return(hyper)
    
  } else if ( method == "fisher" ) {
    # Fisher's Exact Test
    fisher <- lapply(c('g', 'l'), 
                     function(x) get_fisher(success_in_sample,
                                            success_in_bkgd - success_in_sample,
                                            sample_size - success_in_sample,
                                            failure_in_bkgd - sample_size + success_in_sample,
                                            type = x)
    )
    
    names(fisher) <- c('enrichment', 'depletion')
    
    message("[*] Fisher's Exact test [*]")
    message("    P-values")
    message("    Enrichment: ", formatC(fisher$enrichment, format = 'e', digits = 2))
    message("    Depletion : ", formatC(fisher$depletion, format = 'e', digits = 2))
    
    if(retres) return(fisher)
    
  } else {
    
    stop(message("[!] Error: incorrect enrichment test method."))
  }
}

# Remove outiliers ====
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
