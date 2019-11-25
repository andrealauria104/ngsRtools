# Annotate CpG ====
associate_to_feature <- function(sites, feature)
{
  require(GenomicRanges)
  
  if(grepl('data.frame', class(sites))) {
    tranges <- makeGRangesFromDataFrame(sites)
  } else if(grepl('methyl', class(sites))) {
    tranges <- as(sites,"GRanges")
  } else {
    tranges <- sites
  }
  
  dist_feature <- distance2NearestFeature(tranges, tss = feature)
  sites$feature[dist_feature$target.row] <- dist_feature$feature.name
  sites$distance[dist_feature$target.row] <- dist_feature$dist.to.feature
  
  return(sites)
}
get_target_feature <- function(targets, mtr)
{
  require(plyr)
  mtr.ratio <- percMethylation(mtr, rowids = T)
  
  targets <- dlply(targets, ~feature, function(x) x[,c('names','feature')])
  n_sites <- unlist(lapply(targets, nrow))
  idx <- which(n_sites>=quantile(n_sites)['50%'])
  tmp <- targets[idx]
  j <- 1
  y <- lapply(tmp, function(x) {
    print(j)
    x$names <- gsub("\\.\\+","",gsub("_",".",x$names))
    x$names <- gsub(".random", "_random", x$names)
    idx     <- match(x$names, rownames(mtr.ratio))
    x  <- cbind(x, mtr.ratio[idx,])
    id <- unique(colnames(mtr.ratio))
    
    for(i in id) {
      idx <- which(colnames(x)==i)
      x$mean <- apply(x[,idx], 1, function(x) mean(na.omit(x)))
      colnames(x)[ncol(x)] <- paste0('mean_',i)
    }
    
    keep <- c('names', grep('mean', colnames(x), value = T))
    x <- reshape2::melt(x[,keep])
    x <- ddply(x, .(variable), summarize, mean = mean(value))
    colnames(x)[1] <- 'condition'
    x$condition <- gsub('mean_','', x$condition)
    j <<- j+1
    return(x)
  })
  
  for(i in names(y)) {
    y[[i]]$n_sites <- n_sites[i]
    y[[i]]$av_diff <- y[[i]][1,2] - y[[i]][2,2]
  }
  
  return(y)
}
get_mC_feature <- function(sites, feature_name)
{
  if(grepl('methyl', class(sites))) {
    sites <- getData(sites)
  }
  subset(sites, feature==feature_name)
}
build_genomic_annotation <- function(path_tss, path_exons, path_introns
                                     , path_names = NULL
                                     , up.flank = 1500
                                     , dw.flank = 1500)
{
  # Build genomic annotation from bed
  require(GenomicRanges)
  require(genomation)
  
  message("[+] Build genomic annotation from bed")
  
  message(" -- Reading TSS from: ", path_tss)
  tss       <- readBed(path_tss, zero.based = F)
  colnames(tss@elementMetadata)[2] <- 'gene_name'
  
  message(" -- Define promoters, boundaries: up = ", up.flank, " down = ", dw.flank)
  promoters <- extend(x = tss, upstream = up.flank, downstream = dw.flank)
  
  message(" -- Reading exons from: ", path_exons)
  exons   <- readBed(path_exons  , zero.based = T)
  
  message(" -- Reading introns from: ", path_introns)
  introns <- readBed(path_introns, zero.based = T)
  
  
  gen_ann <- list('tss'       = tss,
                  'promoters' = promoters,
                  'exons'     = exons,
                  'introns'   = introns)
  
  if( !is.null(path_names) ) {
    
    refseq_to_genename <- read.table(path_names, col.names = c('gene_name','refseq'))
    
    gen_ann[c('exons','introns')] <- lapply(gen_ann[c('exons','introns')], function(i)
    {
      i$refseq    <- unlist(regmatches(i$name, gregexpr('N+[M,R]_+\\d+',i$name)))
      
      rid <- intersect(i$refseq, refseq_to_genename$refseq)
      i <- subset(i, refseq%in%rid)
      
      idx <- match(i$refseq, refseq_to_genename$refseq)
      i$gene_name <- refseq_to_genename[idx,'gene_name']
      i$name <- NULL
      return(i)
    }
    )
  }
  
  return(gen_ann)
}
build_other_annotation <- function(paths)
{
  ann <- vector(mode = 'list', length = length(paths))
  get_path_names <- function(x) {
    gsub("\\.+.*$","",basename(x))
  }
  nm <- unlist(lapply(paths, get_path_names))
  names(ann) <- nm
  
  for(i in seq_along(paths)) {
    ann[[i]] <- genomation::readBed(paths[i])
  }
  
  return(ann)
}
prepare_maps <- function(mtr, gen_ann, ...)
{
  if(grepl('methyl', class(mtr))) {
    mtrranges <- as(mtr,"GRanges")
  } else if(grepl('data.frame', class(mtr))) {
    mtrranges <- makeGRangesFromDataFrame(mtr)
  } else if(grepl('GRanges', class(mtr))) {
    mtrranges <- mtr
  }
  
  if( is.character(gen_ann) ) {
    if(grepl('.bed', gen_ann)) {
      gen_ann <- readBed(gen_ann, ...)
    }
  }
  
  inputs <- list('mtrranges' = mtrranges
                 , 'gen_ann' = gen_ann)
  return(inputs)
}
map_cytosines <- function(mtr, gen_ann
                          , ordered = T
                          , precedence = c('promoters', 'exons', 'introns','intergenic')
                          , ...
)
{
  require(GenomicRanges)
  
  inputs <- prepare_maps(mtr = mtr, gen_ann = gen_ann, ...)
  
  mtrranges <- inputs$mtrranges
  gen_ann   <- inputs$gen_ann
  
  if(!is.vector(gen_ann)) ordered <- F
  
  if(!ordered) {
    
    message("[+] Mapping Cs")
    
    if(is.vector(gen_ann)) {
      nm <- names(gen_ann)
      maps        <- vector(mode = 'numeric', length = length(nm))
      names(maps) <- nm
      ovs         <- vector(mode = 'list'   , length = length(nm))
      names(ovs)  <- nm
      
      for(i in nm) {
        message(" --- region: ", i)
        
        ov <- findOverlaps( mtrranges
                            , gen_ann[[i]]
                            , ignore.strand=T
                            , type = 'any')
        
        hits    <- unique(queryHits(ov))
        maps[i] <- length(hits)
        ovs[[i]] <- mtrranges[hits]
        
      }
    } else {
      print('here')
      ov <- findOverlaps( mtrranges
                          , gen_ann
                          , ignore.strand=T
                          , type = 'any')
      
      hits <- unique(queryHits(ov))
      maps <- length(hits)
      ovs  <- mtrranges[hits]
    }
    
    # return(maps)
    
  } else if(ordered) {
    
    precedence <- ordered(factor(precedence))
    tomap <- mtrranges
    
    maps        <- vector(mode = 'numeric', length = length(precedence))
    names(maps) <- precedence
    ovs         <- vector(mode = 'list'   , length = length(precedence))
    names(ovs)  <- precedence
    
    message("[+] Mapping Cs with precedence: ", paste0(precedence, collapse = ' > '))
    
    for(i in precedence) {
      message(" --- region: ", i)
      
      if( i!='intergenic' ) {
        ov <- findOverlaps( tomap
                            , gen_ann[[i]]
                            , ignore.strand=T
                            , type = 'any')
        
        hits    <- unique(queryHits(ov))
        maps[i] <- length(hits)
        ovs[[i]] <- tomap[hits]
        tomap <- tomap[-hits,]
        
      } else {
        maps[i] <- length(tomap)
        ovs[[i]] <- tomap
      }
    }
  }
  return(list('maps' = maps, 'ovs' = ovs))
}
plot_map_pie <- function(maps, ptitle = NULL)
{
  #source("theme_setting.R")
  require(scales)
  
  mp <- cbind.data.frame('region' = names(maps), reshape2::melt(maps))
  if(any(grepl('promoters', mp$region))) {
    mp$region <- factor(mp$region, levels = c('promoters', 'exons','introns','intergenic'))
  }
  
  n <- length(mp$region)
  
  mp$ypos <- cumsum(mp$value)[4] - cumsum(mp$value) + mp$value/2
  mp$perc <- round(mp$value/sum(mp$value),4)
  
  pie <- ggplot(mp, aes(x="", y=value, fill=region)) +
    geom_bar(width = 1, stat = "identity", col = 'black') +
    coord_polar("y", start=0) +
    scale_fill_npg() +
    blank_theme + ggtitle(ptitle) +
    theme(axis.text.x=element_blank()
          , plot.title = element_text(size = 8, face = 'bold', hjust = 0.5, vjust = 0.7)
          , text = element_text(size = 8))+
    geom_label(aes(y = ypos,
                   label = percent(perc))
               , size=2.8
               , show.legend = F
               , fill = 'white')
  
  return(pie)
}
map_regions <- function(gen_ann, mtr, ...)
{
  inputs <- prepare_maps(mtr = mtr, gen_ann = gen_ann, ...)
  
  mtrranges <- inputs$mtrranges
  gen_ann   <- inputs$gen_ann
  
  cov <- countOverlaps( gen_ann
                        , mtrranges
                        , ignore.strand=T
                        , type = 'any')
  names(cov) <- as.character(gen_ann$gene_name)
  
  ov <- findOverlaps( gen_ann
                      , mtrranges
                      , ignore.strand=T
                      , type = 'any')
  
  out1 <- as.data.frame(gen_ann[queryHits(ov),])
  tmp <- mtrranges[subjectHits(ov),]
  names(tmp) <- NULL
  out2 <- as.data.frame(tmp)[,1:3]
  out <- cbind.data.frame(out1, out2)
  out$counts <- cov[as.character(out$gene_name)]
  
  if(any(grepl('refseq', colnames(out)))) {
    
    require(plyr)
    out$idx <- gsub(" ","", apply(out[,9:11], 1, paste0, collapse = '_'))
    out <- out[which(!duplicated(out[,'idx'])),]
    
    outc <- dlply(unique(out[,c('gene_name','refseq','counts')]), .(gene_name), function(x)
    {
      counts <- sum(x$counts)
      return(counts)
    })
    outc <- unlist(outc)
    
    out$counts <- outc[as.character(out$gene_name)]
    out$refseq <- NULL
    out$idx    <- NULL
  }
  
  return(out)
}
map_diff_mC <- function(diff_mC, gen_ann
                        , q_th  = 0.05
                        , m_th  = 10
                        , ...)
{
  if(grepl('methyl', class(diff_mC))) {
    diff_mC <- getData(diff_mC)
  }
  hyper <- subset(diff_mC, meth.diff>=m_th & qvalue<=q_th)
  hypo <- subset(diff_mC, meth.diff<=(-m_th) & qvalue<=q_th)
  
  hyper <- map_cytosines(mtr = hyper, gen_ann = gen_ann, ...)
  hypo <- map_cytosines(mtr = hypo, gen_ann = gen_ann, ...)
  
  return(list('hyper' = hyper, 'hypo' = hypo))
}
get_maps_mtr <- function(mtr, maps, average=T)
{
  # Calculate per sample methratio
  mtr.ratio <- percMethylation(mtr, rowids = T)
  rm(mtr)
  if(average) mtr.ratio <- get_average_methylation(mtr.ratio)
  
  message("[+] Get methylation levels by annotation")
  if(is.list(maps)) {
    maps <- lapply(maps, function(x)
    {
      paste0(seqnames(x),".",start(x),".",end(x))
    })
    
    maps.mtratio <- vector(mode = 'list', length(maps))
    names(maps.mtratio) <- names(maps)
    
    for(i in names(maps)) {
      message(" -- ",i)
      idx <- intersect(maps[[i]], rownames(mtr.ratio))
      maps.mtratio[[i]] <- mtr.ratio[idx,]
    }
  } else {
    maps <- paste0(seqnames(maps),".",start(maps),".",end(maps))
    idx <- intersect(maps, rownames(mtr.ratio))
    maps.mtratio <- mtr.ratio[idx,]
  }
  return(maps.mtratio)
}
prepare_enhancer <- function(path_enhancer)
{
  enhancer <- read.table(path_enhancer
                         , header = F
                         , nrows = 29134
                         , stringsAsFactors = F)
  
  enhancer <- makeGRangesFromDataFrame(enhancer
                                       , seqnames.field = 'V1'
                                       , start.field = 'V2'
                                       , end.field = 'V3')
  return(enhancer)
}
prepare_icr <- function(path_icr, mode = 'base')
{
  tmp_icr <- read.delim2(path_icr
                         , header = F
                         , col.names = c('chr', 'start', 'end', 'name')
                         , stringsAsFactors = F)
  
  tmp_icr <- makeGRangesFromDataFrame(tmp_icr
                                      , seqnames.field = 'chr'
                                      , start.field    = 'start'
                                      , end.field      = 'end'
                                      , keep.extra.columns = T)
  
  if(mode=='base') {
    nm  <- unique(tmp_icr$name)
    icr <- vector(mode = 'list', length = length(nm))
    names(icr) <- nm
    
    for(i in nm) {
      icr[[i]] <- tmp_icr[tmp_icr$name==i,]
    }
    
  } else if (mode=='region') {
    icr <- tmp_icr
  }
  
  return(icr)
}
prepare_repeats <- function(path_repeats, path_repinfo)
{
  message('[+] preparing repeats')
  ann <- genomation::readBed(path_repeats)
  info <- read.delim2(path_repinfo, header = T, stringsAsFactors = F)
  colnames(info) <- gsub('X.','',colnames(info))
  
  repClass    <- unique(info$repClass)
  repeats_ann <- vector(mode = 'list', length = length(repClass))
  names(repeats_ann) <- repClass
  
  for(i in repClass) {
    message(' -- class: ', i)
    cidx <- subset(info, repClass==i)[,'repName']
    midx <- which(ann$name%in%cidx)
    repeats_ann[[i]] <- ann[midx,]
  }
  
  return(repeats_ann)
}
get_rangestring <- function(grange)
{
  paste(seqnames(grange),start(grange),end(grange), sep = '.')
}
