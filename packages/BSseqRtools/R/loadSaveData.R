# Load and prepare data ====
read_bsmap     <- function(bsmap_files
                           , context  = "CpG"
                           , assembly = "mm9"
                           , min  = NULL
                           , cv   = 10
                           , hiTh = 99.9
                           , normalize = F
                           , ...) {
  require(methylKit)
  
  message("[+] Read bsmap methratio \n")
  
  lapply(bsmap_files, function(i)
  {
    message(" [*] Processing file:", i)
  }
  )
  
  # n_pheno <- length(unique(gsub("_R+\\d", "", names(bsmap_files)))) - 1
  # n_rep   <- length(unique(gsub(".*_R"  , "", names(bsmap_files))))
  n_rep <- table(names(bsmap_files))
  pheno <- unique(names(bsmap_files))
  
  # treat <- rev(sapply(0:n_pheno, function(x) rep(x,n_rep)))
  i <- 0
  treat <- rev(sapply(pheno, function(x) {
    nr <- rep(i, n_rep[x])
    i <<- i + 1
    return(nr)
  }))
  
  obj <- methRead( location    = as.list(bsmap_files)
                   , sample.id  = as.list(names(bsmap_files))
                   , treatment  = treat
                   , assembly   = assembly
                   , header     = TRUE
                   , context    = context
                   , resolution = "base"
                   , mincov     = cv
                   , pipeline   = list(fraction       = TRUE
                                       , chr.col      = 1
                                       , start.col    = 2
                                       , end.col      = 2
                                       , coverage.col = 6
                                       , strand.col   = 3
                                       , freqC.col    = 5)
  )
  
  obj <- filterByCoverage(obj
                          , lo.count = cv
                          , lo.perc  = NULL
                          , hi.count = NULL
                          , hi.perc  = hiTh)
  
  if(normalize) {
    message('[+] Normalizing coverage ...')
    obj <- normalizeCoverage(obj)
  }
  
  meth <- unite(obj, destrand=F, min.per.group=min)
  
  return(meth)
}

write_to_bed <- function(mtr
                         , bedfile = NULL)
{
  options(scipen=999) # disable scientific notation
  
  message("[+] Writing methylkit object to file: ", bedfile)
  tobed <- as( mtr, "data.frame" )
  
  tobed <- cbind('chr' = tobed[,1], format(tobed[,2:3], scientific = FALSE))
  tobed[,2:3] <- lapply(tobed[,2:3], as.numeric)
  tobed$start <- tobed$start-1
  
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

create_begraph <- function(mtr.ratio, outfile, samp)
{
  message("[+] Creating bedgraph")
  message(" -- output: ", outfile)
  
  chr   <- sapply(strsplit(rownames(mtr.ratio), "[.]"), "[[", 1)
  start <- as.numeric(sapply(strsplit(rownames(mtr.ratio), "[.]"), "[[", 2))-1
  end   <- as.numeric(sapply(strsplit(rownames(mtr.ratio), "[.]"), "[[", 3))
  
  bg <- data.frame('chr' = chr, 'start' = start, 'end' = end)
  bg <- cbind.data.frame(bg, mtr.ratio)
  sidx <- grep(samp, colnames(bg))
  
  setstring <- paste0("track type=bedGraph name='",outfile,"' description='methratio' visibility=full color=204,0,0 altColor=0,0,153 maxHeightPixels=80:80:11")
  sink(file = "tmp_1.bg", append = F)
  cat(setstring, "\n")
  sink()
  
  options(scipen=10)
  
  write.table(bg[,c(1:3,sidx)],
              quote = F      ,
              row.names = F  ,
              col.names = F  ,
              sep = "\t",
              file = "tmp_2.bg")
  
  options(scipen=0)
  
  system(paste0('cat tmp_1.bg tmp_2.bg > ', outfile))
  unlink(c("tmp_1.bg", "tmp_2.bg"))
}

save_bedgraph <- function(mC, outbg)
{
  bedgraph(methylObj   = mC
           , col.name  = 'meth.diff'
           , file.name = outbg
           , unmeth    = T)
  
  cmd1 <- paste0("sed 's/color=255,0,255/color=204,0,0/g' ",outbg," > tmp.bedgraph")
  cmd2 <- "sed 's/altColor=102,205,170/altColor=0,0,153/g' tmp.bedgraph > tmp2.bedgraph"
  cmd3 <- paste0("rm tmp.bedgraph ",outbg)
  cmd4 <- paste0("mv tmp2.bedgraph ",outbg)
  
  cmds <- grep('cmd', ls(), value = T)
  cmds <- lapply(cmds, get, envir = environment())
  for(i in cmds) {
    system(command = i)
  }
}
