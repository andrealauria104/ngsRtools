# Process DNA methylation levels ====
get_names <- function(x, y)
{
  y <- as.data.frame(y)
  x$idx <- gsub(" ","",apply(x[,c(1:4)],1, paste0, collapse="_"))
  y$idx <- gsub(" ","",apply(y[,c(1:3,5)],1, paste0, collapse="_"))
  name_id <- grep('gene.*name|^name', colnames(y))
  x$names <- y[match(x$idx,y$idx),name_id]
  x$idx <- NULL
  return(x)
}

get_region_methylation <- function(mtr, region, ...)
{
  cregion <- regionCounts(mtr, region, strand.aware = F, ...)
  tmp <- getData(cregion)
  tmp <- get_names(tmp, region)
  
  idx <- which(duplicated(tmp$names))
  if(length(idx)>0) {
    tmp <- tmp[-idx,]
    mratio <- percMethylation(cregion[-idx], rowids = T)
  } else {
    mratio <- percMethylation(cregion, rowids = T)
  }
  
  if(length(tmp$names)>0) rownames(mratio) <- tmp$names
  
  return(mratio)
}

get_average_methylation <- function(mratio)
{
  
  if(grepl('methyl', class(mratio))) {
    mratio <- percMethylation(mratio, rowids = T)
  }
  
  idx <- unique(colnames(mratio))
  idx <- unique(gsub("_rep_+\\d+","",idx))
  mratio <- lapply(idx, function(x)
  {
    i <- grep(x, colnames(mratio))
    rowMeans(mratio[,i], na.rm = T)
  }
  )
  names(mratio) <- idx
  mratio <- do.call(cbind, mratio)
  return(mratio)
}

get_average_coverage <- function(mtr)
{
  mtr <- getData(mtr)
  mtr$position <- gsub(" ","", apply(mtr[,1:3], 1, paste0, collapse = '.'))
  i <- grep('coverage', colnames(mtr))
  cv <- rowMeans(mtr[,i], na.rm = T)
  
  mtr$average_cv <- cv
  
  return(mtr[,c('position','average_cv')])
}

get_qmtr_ratio <- function(mtr.ratio)
{
  apply(mtr.ratio, 2, function(x) {
    n <- list()
    nn <- vector(mode = 'character')
    qn <- seq(0,100,20)
    nstep <- length(qn)-1
    for(i in 1:nstep) {
      if(i!=nstep) {
        n[[i]] <- sum( x>=qn[i] & x<qn[i+1] )
      } else {
        n[[i]] <- sum( x>=qn[i] & x<=qn[i+1] )
      }
      nn[i] <- paste0(qn[i],"-",qn[i+1], "%")
    }
    n <- do.call(rbind, n)
    names(n) <- nn
    return(n)
  } )
}

plot_qmtr_ratio <- function(q.mtr.ratio, fann = NULL)
{
  qm <- reshape2::melt(q.mtr.ratio, varnames = c('methylation','sample'))
  nm <- rev(levels(qm$methylation))
  qm$methylation <- factor(qm$methylation, levels = nm)
  qm$sample <- factor(qm$sample)
  
  if(!is.null(fann) & is.character(fann)) {
    qm$sample <- factor(qm$sample, levels =  names(fann))
    qm$fann <- fann[qm$sample]
  }
  
  nref <- max(apply(q.mtr.ratio, 2, function(i) max(cumsum(i))))
  om   <- floor(log10(nref))
  br <- c(0, floor(nref/10**om)*10**om/2, floor(nref/10**om)*10**om)
  
  pal <- rev(RColorBrewer::brewer.pal(5, 'Reds'))
  bars <- ggplot(qm, aes(x=sample, y=value, fill=methylation)) +
    geom_bar(width = 0.8, stat = "identity", col = 'black', size = 0.25) +
    scale_fill_manual(values = pal) + my_theme + theme_bw() +
    scale_y_continuous(labels = br/10**om, breaks = br) + 
    ylab(bquote(CG~count~"("~x~10^.(om)~")")) + xlab(NULL) + 
    theme(strip.background = element_blank()
          , strip.text = element_text(size = 8)
          , text = element_text(size = 8)
          , plot.title = element_text(size = 8, hjust = 0.5, face = "bold")
          , plot.background = element_rect(size = 0.25)
          , panel.background = element_rect(size = 0.25)
          , panel.border = element_rect(size=0.25)
          , panel.grid = element_line(size = 0.25)
          , axis.ticks = element_line(size = 0.25))
  
  if(!is.null(fann) & is.character(fann)) {
    bars <- bars + facet_grid(~fann, scales = 'free_x') 
  }
  
  return(bars)
}

plot_methratio_v2 <- function(mtr
                              , type='boxplot' # violin, bar, density, histogram
                              , pal = NULL
                              , precedence = NULL
                              , time.pos.id = NULL
                              , cond.pos.id = NULL
                              , sample.precedence = NULL)
{
  if(grepl('methyl', class(mtr))) {
    mtratio  <- methylKit::percMethylation(mtr, rowids = T)
    if(sum(table(colnames(mtr.ratio)))!=ncol(mtr.ratio)) mtratio <- get_average_methylation(mtratio)
  } else {
    mtratio  <- mtr
  }
  rm(mtr)
  
  mtratio <- reshape2::melt(mtratio, varnames = c('position','sample'), value.name = 'methratio')
  colnames(mtratio)[1] <- 'position'
  
  if(any(grepl("L1", colnames(mtratio))) & !is.null(precedence)) {
    mtratio$L1 <- factor(mtratio$L1, levels = precedence) 
  }
  if(!is.null(time.pos.id)) {
    mtratio$time <- sapply(strsplit(as.character(mtratio$sample), "\\_"), "[[", time.pos.id) 
    if(!is.null(precedence)) {
      tryCatch({mtratio$time <- factor(mtratio$time, levels = precedence)}
               , error = function(e) message("[!] Invalid precedence!"))
    }
  } else {
    mtratio$time <- mtratio$sample
  }
  
  if(!is.null(cond.pos.id)) {
    mtratio$condition <- sapply(strsplit(as.character(mtratio$sample), "\\_"), "[[", cond.pos.id) 
  } else {
    mtratio$condition <- mtratio$sample
  }
  
  
  if(any(grepl("_rep_+\\d+|_replicate_", unique(mtratio$sample))) & !(type%in%c("density","histogram"))){
    
    mtratio$tmp <- mtratio$condition
    mtratio$condition <- mtratio$sample
    mtratio$sample <- mtratio$tmp
    mtratio$sample <- gsub("_rep_+\\d+|_replicate_","",mtratio$sample)
    mtratio$tmp <- NULL
    mtratio$condition <- gsub("_","-",mtratio$condition)
    
  } else {
    stop(message("[!] Please provide average values for density/histogram plots with replicates."))
  }
  
  if(!is.null(sample.precedence)) {
    if(length(sample.precedence)==length(unique(mtratio$condition))) {
      mtratio$condition <- factor(mtratio$condition, levels = gsub("_","-",sample.precedence))
    } 
  }
  
  mtratio$sample <- gsub("_","-",mtratio$sample)
  
  if(type=='boxplot') {
    p <- ggplot(mtratio, aes(x=condition, y=methratio, fill=sample)) +  
      geom_boxplot(notch = T, outlier.shape = NA, width=0.6) + theme_bw() + my_theme + 
      theme(plot.title = element_text(face="bold", hjust = 0.5, size = 8)) + ylim(c(0,100)) +
      guides(col = guide_legend(nrow=2), fill = guide_legend(nrow = 2)) +
      scale_fill_manual(values = pal) + ylab("% CG methylation") 
    
  } else if(type=='violin') {
    if(any(grepl("L1", colnames(mtratio)))) {
      facet_formula <- paste0("L1~time")
      p0 <- ggplot(mtratio, aes(x=condition, y=methratio, fill=sample)) + facet_grid(as.formula(facet_formula))
    } else if(!is.null(time.pos.id)){
      facet_formula <- "~time"
      p0 <- ggplot(mtratio, aes(x=condition, y=methratio, fill=sample)) + facet_grid(as.formula(facet_formula))
    } else {
      p0 <- ggplot(mtratio, aes(x=condition, y=methratio, fill=sample))
    }
    dodge <- position_dodge(width = 0.8)
    p <- p0 +
      geom_violin(trim = T, position = dodge, lwd = 0.25, alpha = 0.7) + 
      geom_boxplot(width=0.08, outlier.color = NA, position = dodge, lwd = 0.25, show.legend = F, alpha = 0.7) + 
      stat_summary(fun.y=median, geom="point", size=0.5, color="black", position = dodge) +
      theme_classic() + my_theme +
      theme(plot.title = element_text(face="bold", hjust = 0.5, size = 8)
            , panel.grid = element_blank()
            , strip.background.y = element_blank()
            , strip.background.x = element_rect(size=0.25)) + guides(col = guide_legend(nrow=2), fill = guide_legend(nrow = 2)) +
      scale_fill_manual(values = pal) + ylab("% CG methylation")  + scale_y_continuous(breaks = c(0,50,100), limits = c(0,100))  
  }else if(type=='density'){
    if(any(grepl("L1", colnames(mtratio)))) {
      facet_formula <- paste0("L1~time")
    } else {
      facet_formula <- "~time"
    }
    dodge <- position_dodge(width = 0.8)
    p <- ggplot(mtratio, aes(x=methratio, fill=sample)) + facet_grid(as.formula(facet_formula), scales = 'free_y') +
      geom_density(lwd = 0.25, alpha = 0.6) +
      theme_classic() + my_theme +
      theme(plot.title = element_text(face="bold", hjust = 0.5, size = 8)
            , panel.grid = element_blank()
            , strip.background.y = element_blank()
            , strip.background.x = element_rect(size=0.25)
            # , axis.title.y = element_blank()
            , axis.text.y = element_blank()) + guides(fill = guide_legend(nrow=2)) +
      scale_fill_manual(values = pal) + xlab("% CG methylation") 
    # scale_x_continuous(breaks = c(0,50,100), limits = c(0,100))  
  }else if(type=='histogram') {
    if(any(grepl("L1", colnames(mtratio)))) {
      facet_formula <- paste0("L1~time")
    } else {
      facet_formula <- "~time"
    }
    dodge <- position_dodge(width = 0.8)
    p <- ggplot(mtratio, aes(x=methratio, fill=sample)) + facet_grid(as.formula(facet_formula), scales = 'free_y') +
      geom_histogram(lwd = 0.25, alpha = 0.6, binwidth = 10, col = 'black') +
      theme_classic() + my_theme +
      theme(plot.title = element_text(face="bold", hjust = 0.5, size = 8)
            , panel.grid = element_blank()
            , strip.background.y = element_blank()
            , strip.background.x = element_rect(size=0.25)
            # , axis.title.y = element_blank()
            # , axis.text.y = element_blank()
      ) + guides(fill = guide_legend(nrow=2)) +
      scale_fill_manual(values = pal) + xlab("% CG methylation") + ylab("CG count")
  }else if(type=="summary"){
    mtratio.sm <- ddply(mtratio, .(time, condition, sample)
                        , summarize
                        , mean = mean(methratio)
                        , median = median(methratio)
                        , min = min(methratio)
                        , max = max(methratio)
                        , sd  = sd(methratio)
                        , l = length(methratio)
                        , se = sd(methratio)/sqrt(length(methratio))
                        , se.md = 1.253*sd(methratio)/sqrt(length(methratio))
    )
    
    p <- ggplot(mtratio.sm, aes(x=time, y=median, fill=sample)) +  
      geom_col(position = 'dodge', size = 0.1, col = 'black', width = 0.6) + 
      geom_errorbar(
        # aes(ymin=mean-se, ymax=mean+se),
        aes(ymin=median-se.md, ymax=median+se.md)
        , size=.25
        , width=.2
        , position=position_dodge(.6)) + 
      theme_bw() + my_theme + 
      theme(plot.title = element_text(face="bold", hjust = 0.5, size = 8)) +
      scale_fill_manual(values = pal) + ylab("% of methylation")  + ylim(c(0,100))
    
    
  } else {
    stop(message("[!] Invalid plot type. Please, provide one of # violin, bar, density, histogram"))
  }
  
  return(p)
  
}

intersect_mtr <- function(mtr_x, mtr_y
                          , return.obj = "GRanges"
                          , return.merged = F)
{
  convert_input <- function(mtr)
  {
    if(grepl('methyl', class(mtr))) {
      mtrranges <- as(mtr,"GRanges")
    } else if(grepl('data.frame', class(mtr))) {
      mtrranges <- GenomicRanges::makeGRangesFromDataFrame(mtr)
    } else if(grepl('GRanges', class(mtr))) {
      mtrranges <- mtr
    }
    return(mtrranges)
  }
  mtr_x <- convert_input(mtr_x)
  mtr_y <- convert_input(mtr_y)
  
  ov <- IRanges::findOverlaps(mtr_x, mtr_y, ignore.strand=T, type = 'any')
  x <- mtr_x[S4Vectors::queryHits(ov),]
  y <- mtr_y[S4Vectors::subjectHits(ov),]
  
  if(return.merged) {
    colnames(x@elementMetadata) <- paste0(colnames(x@elementMetadata),"_x")
    colnames(y@elementMetadata) <- paste0(colnames(y@elementMetadata),"_y")
    x@elementMetadata <- cbind.DataFrame(x@elementMetadata, y@elementMetadata)
    m <- x
    rm(x,y)
    if(return.obj == "GRanges") {
      return(m)
    } else if(return.obj == "data.frame") {
      return(as.data.frame(m))
    } else {
      stop(message("[!] Invalid return object (GRanges/data.frame)."))
    }
  } else {
    if(return.obj == "GRanges") {
      return(list("mtr_x"=x,"mtr_y"=y))
    } else if(return.obj == "data.frame") {
      return(list("mtr_x"=as.data.frame(x),"mtr_y"=as.data.frame(y)))
    } else {
      stop(message("[!] Invalid return object (GRanges/data.frame)."))
    }
  }
}
