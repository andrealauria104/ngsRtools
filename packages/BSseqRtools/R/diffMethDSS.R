# Load data ====
prepareDSS <- function(path_files
                       , pipeline = "bismark"
                       , sample_names = NULL
                       , min.cov=5) 
{
  # Internal data import functions
  read_func <- function(x, pipeline) {
    
    read_bismark <- function(x)
    {
      methylc_data <- read.delim(x, header=F, stringsAsFactors=F)[,c(1:2,5:6)]
      methylc_data$N <- methylc_data[,3]+methylc_data[,4]
      methylc_data[,4] <- NULL
      methylc_data <- methylc_data[,c(1:2,4,3)]
      colnames(methylc_data) <- c("chr","pos","N","X")
      return(methylc_data)
    }
    
    read_bsmap <- function(x) 
    {
      methylc_data <- read.delim(x, header=F, stringsAsFactors=F)[,c(1:2,7:8)]
      methylc_data$X <- methylc_data[,4]-methylc_data[,3]
      methylc_data[,3] <- NULL
      colnames(methylc_data) <- c("chr","pos","N","X")
      return(methylc_data)
    }
    
    if(pipeline=="bismark") {
      return(read_bismark(x))
    } else if(pipeline=="bsmap"){
      return(read_bsmap(x))
    }
  }
  
  if(dir.exists(path_files)) {
    if(pipeline=="bismark") {
      file_pattern <- ".CpG_report.merged_CpG_evidence.cov.gz$"
    } else if(pipeline=="bsmap"){
      file_pattern <- ".txt$"
    }
    file_list <- list.files(path_files, pattern = file_pattern, full.names = T)
  } else {
    file_list <- path_files	
  }
  
  if(is.null(sample_names)) sample_names <- gsub("\\..*$","",basename(file_list))
  # Read data 
  methylc_data <- lapply(file_list, function(x) {
    message(" -- reading data from: ",x)
    mdata <- read_func(x, pipeline=pipeline)
    if(!is.null(min.cov)) mdata <- subset(mdata, N>=min.cov)
    return(mdata)
  })
  # Create BSobj
  message("\n -- creating BSobj")
  BSobj <- DSS::makeBSseqData(methylc_data, sample_names)
  
  return(BSobj)
}

# Run DSS ====
runTwoGroupDSS <- function(BSobj, group1, group2, smoothing
                           , equal.disp=FALSE
                           , smoothing.span=500
                           , dml.delta=0.1
                           , dml.p.threshold=0.001
                           , dmr.delta=0.1
                           , dmr.p.threshold=0.05
                           , dmr.minlen=50
                           , dmr.minCG=3
                           , dmr.dis.merge=100
                           , dmr.pct.sig=0.5
                           , n.cores=1) 
{	
  message(" [+] DSS for two-group comparisons")
  message(" [+] Running DMLtest, parameters:")
  message(" -- smoothing = ", smoothing)
  if(smoothing) message(" -- smoothing.span = ",smoothing.span)
  message(" -- comparison: ", group1, " vs ", group2)
  
  dmlTest <- DSS::DMLtest(BSobj
                          , group1=group1
                          , group2=group2
                          , equal.disp=equal.disp 
                          , smoothing=smoothing
                          , smoothing.span=smoothing.span
                          , BPPARAM=MulticoreParam(workers=n.cores, progressbar=TRUE))
  
  message(" [+] Running callDML, parameters:")
  message(" -- dml.delta = ",dml.delta)
  message(" -- dml.p.threshold = ", dml.p.threshold)
  message(" -- dmr.delta = ", dmr.delta)
  message(" -- dmr.p.threshold  = ", dmr.p.threshold)
  message(" -- dmr.minlen = ", dmr.minlen)
  message(" -- dmr.minCG = ", dmr.minCG)
  message(" -- dmr.dis.merge = ", dmr.dis.merge)
  message(" -- dmr.pct.sig = ", dmr.pct.sig)
  
  dmls <- DSS::callDML(dmlTest
                       , delta=dml.delta
                       , p.threshold=dml.p.threshold)
  
  message(" [+] Running callDMR, parameters:")
  message(" -- dmr.minlen = ", dmr.minlen)
  message(" -- dmr.minCG = ", dmr.minCG)
  message(" -- dmr.dis.merge = ", dmr.dis.merge)
  message(" -- dmr.pct.sig = ", dmr.pct.sig)
  
  dmrs <- callDMR(dmlTest
                  , delta=dmr.delta
                  , p.threshold=dmr.p.threshold
                  , minlen=dmr.minlen
                  , minCG=dmr.minCG
                  , dis.merge=dmr.dis.merge
                  , pct.sig=dmr.pct.sig)
  
  return(list("dmlTest"=dmlTest, "dmls"=dmls, "dmrs"=dmrs))
}