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
  idx <- unique(gsub("_rep.*\\d+","",idx, ignore.case = T))
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

get_average_methylation_2 <- function(mdata_ratio, metadata, metadata_features_summarize)
{
  metadata_i <- metadata[,c("sample.id",metadata_features_summarize)]
  idx <- list()
  for(i in unique(metadata_i[,metadata_features_summarize])) {
    idx[[i]] <- which(metadata_i[,metadata_features_summarize]==i)
  }
  m <- lapply(idx, function(i) {
   if(length(i)>1) {
     rowMeans(mdata_ratio[,i], na.rm = T)
   } else{
     mdata_ratio[,i]
   }
  })
  m <- do.call(cbind, m)
  return(m)
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

get_qmtr_ratio <- function(mtr.ratio
                           , qn = seq(0,100,20))
{
  apply(mtr.ratio, 2, function(x) {
    n <- list()
    nn <- vector(mode = 'character')
    
    nstep <- length(qn)-1
    for(i in 1:nstep) {
      if(i!=nstep) {
        n[[i]] <- sum( x>=qn[i] & x<qn[i+1], na.rm = T )
      } else {
        n[[i]] <- sum( x>=qn[i] & x<=qn[i+1], na.rm = T )
      }
      nn[i] <- paste0(qn[i],"-",qn[i+1], "%")
    }
    n <- do.call(rbind, n)
    names(n) <- nn
    return(n)
  } )
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

# Plotting =====
plot_qmtr_ratio <- function(q.mtr.ratio, fann = NULL)
{
  qm <- reshape2::melt(q.mtr.ratio, varnames = c('methylation','sample'))
  nm <- rev(levels(qm$methylation))
  qm$methylation <- factor(qm$methylation, levels = nm)
  qm$sample <- factor(qm$sample)
  
  if(!is.null(fann) & is.character(fann)) {
    qm$sample <- factor(qm$sample, levels =  names(fann))
    qm$fann <- fann[qm$sample]
    qm$fann <- factor(qm$fann, levels = unique(fann))
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
          , text = element_text(size = 6)
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

# Quantized distributions ---
plot_qmtr_ratio_2 <- function(mdata_ratio
                              , metadata
                              , metadata_features_summarize
                              , summary.precedence = NULL
                              , metadata_features_facet = NULL
                              , facet.precedence = NULL
                              , plot.title = NULL
                              , qn = seq(0,100,20)
                              , plot_type = "bar"
                              , ylab_stack_suffix = "CpG sites"
                              , bar_width = 0.7
                              , col_width = 0.9) 
{
  qmtr <- get_qmtr_ratio(mdata_ratio, qn = qn)
  
  qmtr_molten <- reshape2::melt(qmtr, varnames=c("Ranges","sample.id"))
  qmtr_molten <- merge(qmtr_molten,metadata, by = "sample.id")
  
  qmtr_molten <- plyr::ddply(qmtr_molten, as.formula(paste0("~sample.id"))
                             , mutate
                             , perc_value = round(100*value/sum(value),2))
  
  qmtr_molten <- plyr::ddply(qmtr_molten, as.formula(paste0("~",metadata_features_summarize,"+Ranges"))
                             , mutate
                             , average_value = mean(value))
  
  qmtr_molten <- plyr::ddply(qmtr_molten, as.formula(paste0("~",metadata_features_summarize))
                             , mutate
                             , perc_average_value = round(100*average_value/sum(unique(average_value)),2))
  
  if(!is.null(summary.precedence)) {
    qmtr_molten[,metadata_features_summarize] <- factor(qmtr_molten[,metadata_features_summarize], levels=summary.precedence) 
  }
  if(!is.null(metadata_features_facet) & !is.null(facet.precedence)) {
    qmtr_molten[,metadata_features_facet] <- factor(qmtr_molten[,metadata_features_facet], levels=facet.precedence)  
  }
  
  if(plot_type=="bar") {
    pal <- RColorBrewer::brewer.pal(length(qn)-1, 'Reds')
    
    nref <- max(qmtr_molten$value)
    om   <- floor(log10(nref))
    br <- c(0, floor(nref/10**om)*10**om/4
            , 2*floor(nref/10**om)*10**om/4
            , 3*floor(nref/10**om)*10**om/4
            , floor(nref/10**om)*10**om)
    
    pqmtr <- ggplot(qmtr_molten, aes_string(x=metadata_features_summarize,y="average_value", fill="Ranges")) + 
      geom_col(position = "dodge",width=col_width,col="black") +
      guides(fill=guide_legend(title = "Methylation range")) +
      geom_point(data=qmtr_molten, aes_string(x=metadata_features_summarize,y="value", fill="Ranges") 
                 , position = position_jitterdodge(dodge.width = .9, jitter.width = .2)
                 , show.legend = F
                 , size = .8) +  scale_fill_manual(values = pal)
    
    if(!is.null(metadata_features_facet)) pqmtr <- pqmtr + facet_wrap(as.formula(paste0("~",metadata_features_facet)), scales = "free_x", nrow=1)
    
    pqmtr <- pqmtr +scale_y_continuous(labels = br/10**om, breaks = br) + 
      theme_bw() + my_theme_2 + theme(plot.title = element_text(size = 8, face = "plain"), legend.key.size = unit(4,'mm')) +
      ylab(bquote(Count~"("~x~10^.(om)~")")) + ggtitle(plot.title)
    
  } else if(plot_type=="stacked_bar") {
    pal <- rev(RColorBrewer::brewer.pal(length(qn)-1, 'Reds'))
    
    qmtr_molten_stack <- unique(qmtr_molten[,c(metadata_features_summarize,"perc_average_value","Ranges",metadata_features_facet)])
    qmtr_molten_stack$Ranges <- factor(qmtr_molten_stack$Ranges
                                       , levels = rev(levels(qmtr_molten_stack$Ranges)))
    
    pqmtr <- ggplot(qmtr_molten_stack, aes_string(x=metadata_features_summarize,y="perc_average_value", fill="Ranges")) +
      geom_bar(width = bar_width, stat = "identity", col = 'black', size = 0.25) +
      scale_fill_manual(values = pal) + theme_bw() +  my_theme + ylab(paste0("% of ",ylab_stack_suffix)) + 
      guides(fill=guide_legend(title = "Methylation range")) +
      theme(legend.position = "right", legend.key.size = unit(4,'mm'), axis.text.x = element_text(angle=45, hjust = 1,vjust = 1, size=8))
    
    if(!is.null(metadata_features_facet)) pqmtr <- pqmtr + facet_wrap(as.formula(paste0("~",metadata_features_facet)), scales = "free_x", nrow=1)
  }
  
  return(pqmtr)
}

plot_coverage_stats <- function(m)
{
  if(grepl('methyl', class(m))) {
    cvstats <- getData(m)[,grep('coverage', colnames(m))]
    colnames(cvstats) <- m@sample.ids
  } else {
    cvstats <- m
  }
  cvstats <- reshape2::melt(cvstats, variable.name = "sample", value.name = "coverage") 
  pcvstats <- ggplot(cvstats, aes(x = coverage)) + facet_wrap(~sample, ncol = 2) +
    geom_histogram() + scale_x_log10() + xlab("Coverage Depth") + ylab("CG count") +
    theme_bw() + my_theme_2
  return(pcvstats)
}

plot_methratio_v2 <- function(mtr
                              , type='boxplot' # violin, bar, density, histogram
                              , pal = NULL
                              , precedence = NULL
                              , time.pos.id = NULL
                              , cond.pos.id = NULL
                              , sample.precedence = NULL
                              , violin.style = 1)
{
  if(grepl('methyl', class(mtr))) {
    mtratio  <- methylKit::percMethylation(mtr, rowids = T)
    if(sum(table(colnames(mtratio)))!=ncol(mtratio)) mtratio <- get_average_methylation(mtratio)
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
  
  
  if(any(grepl("_rep_+\\d+|_replicate_", unique(mtratio$sample))) && !(type%in%c("density","histogram"))){
    
    mtratio$tmp <- mtratio$condition
    mtratio$condition <- mtratio$sample
    mtratio$sample <- mtratio$tmp
    mtratio$sample <- gsub("_rep_+\\d+|_replicate_","",mtratio$sample)
    mtratio$tmp <- NULL
    mtratio$condition <- gsub("_","-",mtratio$condition)
    
  } else if(any(grepl("_rep_+\\d+|_replicate_", unique(mtratio$sample))) && type%in%c("density","histogram")) {
    stop(message("[!] Please provide average values for density/histogram plots with replicates."))
  }
  
  if(!is.null(sample.precedence)) {
    if(length(sample.precedence)==length(unique(mtratio$condition))) {
      mtratio$condition <- factor(mtratio$condition, levels = sample.precedence)
      levels(mtratio$condition) <- gsub("_","-",levels(mtratio$condition))
    }
  }
  
  mtratio$sample <- gsub("_","-",mtratio$sample)
  
  if(type=='boxplot') {
    p <- ggplot(mtratio, aes(x=condition, y=methratio, fill=sample)) +  
      geom_boxplot(notch = T, outlier.shape = NA, width=0.6, lwd=0.25) + theme_bw() + my_theme + 
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
    if(violin.style==1) {
      p <- p0 +
        geom_violin(trim = T, position = dodge, lwd = 0.25, alpha = 0.7) + 
        geom_boxplot(width=0.08, outlier.color = NA, position = dodge, lwd = 0.25, show.legend = F, alpha = 0.7) + 
        stat_summary(fun.y=median, geom="point", size=0.5, color="black", position = dodge, show.legend = F) +
        theme_classic() + my_theme +
        theme(plot.title = element_text(face="bold", hjust = 0.5, size = 8)
              , panel.grid = element_blank()
              , strip.background.y = element_blank()
              , strip.background.x = element_rect(size=0.25)) + guides(col = guide_legend(nrow=2), fill = guide_legend(nrow = 2)) +
        scale_fill_manual(values = pal) + ylab("% CG methylation")  + scale_y_continuous(breaks = c(0,50,100), limits = c(0,100))  
    } else if(violin.style==2) {
      p <- p0 +
        geom_violin(aes(col=sample), trim=T, fill="white") + 
        theme_bw() + my_theme_2 + scale_color_manual(values = pal) + 
        stat_summary(fun.y=median, geom="point", size=1, color="red", show.legend = F) +
        guides(col=guide_legend(nrow=2)) +
        ylab("% CG methylation")  + scale_y_continuous(breaks = c(0,50,100), limits = c(0,100))   
    }
    
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

plot_methratio_v3 <- function(mtratio
                              , type='boxplot' # violin, bar, density, histogram
                              , average = T
                              , pal = NULL
                              , experimental_info = NULL
                              , col_by = NULL
                              , facet_by = NULL
                              , compare_dist = F
                              , my_comparisons = NULL
                              , facet.precedence = NULL
                              , col.precedence = NULL
                              , violin.style = 1)
{
  if(grepl('methyl', class(mtratio))) {
    mtratio  <- methylKit::percMethylation(mtratio, rowids = T)
  } 
  
  if(average & !is.list(mtratio)) mtratio <- get_average_methylation(mtratio)
  
  mtratio <- reshape2::melt(mtratio, varnames = c('position','sample'), value.name = 'methratio')
  mtratio$sample <- factor(mtratio$sample, levels =unique(mtratio$sample))
  
  if(any(grepl("L1", colnames(mtratio))) & !is.null(facet.precedence)) {
    mtratio$L1 <- factor(mtratio$L1, levels = facet.precedence) 
  }
  if(!is.null(facet_by)) {
    mtratio$time <- sapply(strsplit(as.character(mtratio$sample), "\\_"), "[[", facet_by) 
    if(!is.null(facet.precedence)) {
      tryCatch({mtratio$time <- factor(mtratio$time, levels = facet.precedence)}
               , error = function(e) message("[!] Invalid facet.precedence!"))
    }
  } else {
    mtratio$time <- mtratio$sample
  }
  
  if(!is.null(col_by)) {
    if(length(col_by)>1) {
      condition <- sapply(1:length(col_by), function(x) sapply(strsplit(as.character(mtratio$sample), "\\_"), "[[", col_by[x]))
      mtratio$condition <- apply(condition,1,paste0,collapse="_")
    } else {
      mtratio$condition <- sapply(strsplit(as.character(mtratio$sample), "\\_"), "[[", col_by)   
    }
  } else {
    mtratio$condition <- mtratio$sample
  }
  
  if(!average & any(grepl("_rep(\\_)?\\d+|_replicate_", unique(mtratio$sample), ignore.case = T)) && !(type%in%c("density","histogram"))){
    mtratio$tmp <- mtratio$condition
    mtratio$condition <- mtratio$sample
    mtratio$sample <- mtratio$tmp
    mtratio$sample <- gsub("_rep(\\_)?\\d+|_replicate_","",mtratio$sample)
    mtratio$tmp <- NULL
    mtratio$condition <- gsub("_","-",mtratio$condition)
    mtratio$condition <- factor(mtratio$condition, levels = unique(mtratio$condition))
    dist_theme <- theme(axis.text.x = element_text(angle=45,hjust = 1,vjust=1), axis.title.x = element_blank(), legend.key.size = unit(4,'mm'))
    mtratio$sample <- factor(mtratio$sample, levels = unique(mtratio$sample))
    levels(mtratio$sample) <- gsub("_","-",levels(mtratio$sample))
  } else if(any(grepl("_rep(\\_)?\\d+|_replicate_", unique(mtratio$sample), ignore.case = T)) && type%in%c("density","histogram")) {
    stop(message("[!] Please provide average values for density/histogram plots with replicates."))
  } else {
    levels(mtratio$sample) <- gsub("_","-",levels(mtratio$sample))
  }
  
  if(!is.null(col.precedence)) {
    if(length(col.precedence)==length(unique(mtratio$condition))) {
      mtratio$condition <- factor(mtratio$condition, levels = col.precedence)
      levels(mtratio$condition) <- gsub("_","-",levels(mtratio$condition))
    }
  }
  
  if(type=='boxplot') {
    p <- ggplot(mtratio, aes(x=condition, y=methratio, fill=sample)) +  
      geom_boxplot(notch = T, outlier.shape = NA, width=0.6, lwd=0.25) + theme_bw() + my_theme + 
      theme(plot.title = element_text(face="plain", hjust = 0.5, size = 8)) + scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
      guides(col = guide_legend(nrow=2), fill = guide_legend(nrow = 2)) +
      scale_fill_manual(values = pal) + ylab("Methylation [%]") 
    if(any(grepl("(\\_|\\-)rep(\\_|\\-)?\\d+|(\\_|\\-)replicate(\\_|\\-)", levels(mtratio$condition), ignore.case = T))) p <- p + dist_theme
    if(compare_dist) p <- p + stat_compare_means(comparisons = my_comparisons, na.rm = T, size=2, label = "p.signif", method = "wilcox.test") 
  } else if(type=='violin') {
    if(any(grepl("L1", colnames(mtratio)))) {
      facet_formula <- paste0("L1~time")
      p0 <- ggplot(mtratio, aes(x=condition, y=methratio, fill=sample)) + facet_wrap(as.formula(facet_formula), scales = "free_x")
    } else if(!is.null(facet_by)){
      facet_formula <- "~time"
      p0 <- ggplot(mtratio, aes(x=condition, y=methratio, fill=sample)) + facet_wrap(as.formula(facet_formula), scales = "free_x")
    } else {
      p0 <- ggplot(mtratio, aes(x=condition, y=methratio, fill=sample))
    }
    dodge <- position_dodge(width = 0.8)
    if(violin.style==1) {
      p <- p0 +
        geom_violin(trim = T, position = dodge, lwd = 0.25, alpha = 0.7) + 
        geom_boxplot(width=0.08, outlier.color = NA, position = dodge, lwd = 0.25, show.legend = F, alpha = 0.7) + 
        stat_summary(fun.y=median, geom="point", size=0.5, color="black", position = dodge, show.legend = F) +
        theme_classic() + my_theme +
        theme(plot.title = element_text(face="plain", hjust = 0.5, size = 8)
              , panel.grid = element_blank()
              , strip.background.y = element_blank()
              , strip.background.x = element_rect(size=0.25)) + guides(col = guide_legend(nrow=2), fill = guide_legend(nrow = 2)) +
        scale_fill_manual(values = pal) + ylab("% CG methylation")  + scale_y_continuous(breaks = c(0,50,100), limits = c(0,100))  
    } else if(violin.style==2) {
      p <- p0 +
        geom_violin(aes(col=sample), trim=T, fill="white") + 
        theme_bw() + my_theme_2 + scale_color_manual(values = pal) + 
        stat_summary(fun.y=median, geom="point", size=1, color="red",  position = dodge, show.legend = F) +
        guides(col=guide_legend(nrow=2)) + theme(plot.title = element_text(face="plain", hjust = 0.5, size = 8), legend.key.size = unit(4,'mm')) +
        ylab("Methylation [%]")  + scale_y_continuous(breaks = c(0,50,100))   
      if(any(grepl("(\\_|\\-)rep(\\_|\\-)?\\d+|(\\_|\\-)replicate(\\_|\\-)", levels(mtratio$condition), ignore.case = T))) p <- p + dist_theme
      if(compare_dist) p <- p + stat_compare_means(comparisons = my_comparisons, na.rm = T, size=2, label = "p.signif", method = "wilcox.test")
    }
    
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
      theme(plot.title = element_text(face="plain", hjust = 0.5, size = 8)
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
      theme(plot.title = element_text(face="plain", hjust = 0.5, size = 8)
            , panel.grid = element_blank()
            , strip.background.y = element_blank()
            , strip.background.x = element_rect(size=0.25)
            # , axis.title.y = element_blank()
            # , axis.text.y = element_blank()
      ) + guides(fill = guide_legend(nrow=2)) +
      scale_fill_manual(values = pal) + xlab("% CG methylation") + ylab("CG count")
  }else if(type=="summary"){
    
    if(!average & any(grepl("-rep(\\-)?\\d+|-replicate-", levels(mtratio$condition), ignore.case = T))) {
      mtratio.sm <- ddply(mtratio, as.formula("~condition")
                          , summarize
                          , mean = mean(methratio, na.rm=T)
                          , median = median(methratio, na.rm=T)
                          , min = min(methratio, na.rm=T)
                          , max = max(methratio, na.rm=T)
                          , sd  = sd(methratio, na.rm=T)
                          , l = length(methratio)
                          , se = sd(methratio, na.rm=T)/sqrt(length(methratio))
                          , se.md = 1.253*sd(methratio, na.rm=T)/sqrt(length(methratio)))
      mtratio.sm$group <- gsub("-Rep.*","",mtratio.sm$condition,ignore.case = T)
      
      mtratio.sm <- ddply(mtratio.sm, .(group)
                          , mutate
                          , mean.mean = mean(mean, na.rm=T)
                          , se.mean = sd(mean, na.rm=T)/sqrt(length(mean))
                          , sd.mean = sd(mean, na.rm=T))
      mtratio.sm$group <- factor(mtratio.sm$group, levels =  gsub("_","-",col.precedence))
      
      p <- ggplot(mtratio.sm,aes(x=group, y=mean, fill=group)) +  
        geom_col(aes(x=group, y=mean.mean, fill=group),position = 'dodge', size = 0.1, col = 'black', width = 0.6) +
        geom_point(aes(x=group, y=mean, fill=group), show.legend=F, size = .5, position = position_jitterdodge(dodge.width = 0.3))+
        theme_bw() + my_theme + geom_text(aes(x=group, y=mean.mean,label=round(mean.mean,2)), size=2, nudge_y = 7)+ 
        theme(axis.text.x = element_text(angle=45,hjust = 1,vjust=1),plot.title = element_text(face="plain", hjust = 0.5, size = 8), axis.title.x = element_blank(), legend.key.size = unit(4,'mm'), legend.position = "none") +
        scale_fill_manual(values = pal) + ylab("Average methylation [%]")  + ylim(c(0,100))
      if(compare_dist) p <- p + stat_compare_means(data=mtratio.sm, mapping = aes(x=group, y=mean, fill=group),comparisons = my_comparisons,size=2, method = "t.test",inherit.aes = F, label.y = c(78, 85, 90,95))
    } else {
      # mtratio.sm <- ddply(mtratio, as.formula("~condition")
      #                     , summarize
      #                     , mean = mean(methratio, na.rm=T)
      #                     , median = median(methratio, na.rm=T)
      #                     , min = min(methratio, na.rm=T)
      #                     , max = max(methratio, na.rm=T)
      #                     , sd  = sd(methratio, na.rm=T)
      #                     , l = length(methratio)
      #                     , se = sd(methratio, na.rm=T)/sqrt(length(methratio))
      #                     , se.md = 1.253*sd(methratio, na.rm=T)/sqrt(length(methratio)))
      # mtratio.sm$group <- mtratio.sm$condition
      # p <- ggplot(mtratio.sm, aes(x=group, y=mean, fill=group)) +  
      #   geom_col(position = 'dodge', size = 0.1, col = 'black', width = 0.6) + 
      #   geom_errorbar(
      #     aes(ymin=mean-se, ymax=mean+se)
      #     # aes(ymin=median-se.md, ymax=median+se.md)
      #     , size=.25
      #     , width=.2
      #     , position=position_dodge(.6)) + 
      #   theme_bw() + my_theme + 
      #   theme(plot.title = element_text(face="bold", hjust = 0.5, size = 8)) +
      #   scale_fill_manual(values = pal) + ylab("Average methylation [%]")  
      stop(message("[!] Invalid summary plot type.")) 
    }
    
  } else {
    stop(message("[!] Invalid plot type. Please, provide one of # violin, boxplot, summary, density, histogram"))
  }
  
  return(p)
  
}


plot_methratio_v4 <- function(mratio
                              , plot_type='boxplot' # violin, bar, density, histogram
                              , average = T
                              , pal = NULL
                              , metadata = NULL
                              , col_by = NULL
                              , facet_by = NULL
                              , summarize_by = NULL
                              , facet_wrap_nrow = NULL
                              , facet_wrap_ncol = NULL
                              , x_lab = "Sample"
                              , y_lab = "Methylation [%]"
                              , y_breaks = c(0,20,40,60,80,100)
                              , y_limits = c(0,100)
                              , compare_dist = F
                              , my_comparisons = NULL
                              , plot.title = "BS-seq"
                              , facet.precedence = NULL
                              , sample.precedence = NULL
                              , col.precedence = NULL
                              , summary.precedence = NULL
                              , violin.style = 1
                              , violin.point_color = "black"
                              , violin.trim = T
                              , dist_theme = NULL)
{
  # Theme ---
  if(is.null(dist_theme)) {
    dist_theme <- theme_bw() + my_theme + 
      theme(plot.title = element_text(face="plain", hjust = 0.5, size = 8)
            , axis.text.x = element_text(angle=45,hjust = 1,vjust=1)
            , axis.title.x = element_blank()
            , legend.key.size = unit(4,'mm'))
  }
  
  # Reshaping ---
  mratio_molten <- reshape2::melt(mratio)
  colnames(mratio_molten) <- c("pos","sample","value")
  mratio_molten <- cbind.data.frame(mratio_molten, metadata[match(mratio_molten$sample,metadata$sample.id),])
  
  if(!is.null(sample.precedence)) mratio_molten$sample.id <- factor(mratio_molten$sample.id, levels = sample.precedence)
  if(!is.null(facet.precedence) & is.character(facet_by)) mratio_molten[,facet_by] <- factor(mratio_molten[,facet_by], levels = facet.precedence)
  if(!is.null(col.precedence)) mratio_molten[,col_by] <- factor(mratio_molten[,col_by], levels = col.precedence)
  
  # Colors ---
  if(!is.null(col_by)) {
    if(violin.style==1) {
      p0 <- ggplot(mratio_molten, aes_string(x="sample.id", y="value", fill=col_by))   
    } else {
      p0 <- ggplot(mratio_molten, aes_string(x="sample.id", y="value", col=col_by)) 
    }
    
  } else {
    p0 <- ggplot(mratio_molten, aes_string(x="sample.id", y="value")) 
  }
  # Faceting ---
  if(!is.null(facet_by)) {
    if(is.character(facet_by) & !any(grepl("~",facet_by)) & length(facet_by)==1) facet_by <- paste0("~",facet_by)
    p0 <- p0 + facet_wrap(as.formula(facet_by), scales = "free_x", nrow = facet_wrap_nrow,ncol = facet_wrap_ncol)
  }
  
  # Plots ---
  if(plot_type=='boxplot') {
    p <- p0 + geom_boxplot(notch = T, outlier.shape = NA, width=0.6, lwd=0.25) + 
      scale_y_continuous(breaks = y_breaks, limits = y_limits) + dist_theme +
      guides(col = guide_legend(nrow=2), fill = guide_legend(nrow = 2)) + ggtitle(plot.title) +
      scale_fill_manual(values = pal) + ylab(y_lab) + xlab(x_lab)
    
    if(compare_dist) p <- p + stat_compare_means(comparisons = my_comparisons, na.rm = T, size=2, label = "p.signif", method = "wilcox.test") 
    
  } else if(plot_type=='violin') {
    dodge <- position_dodge(width = 0.8)
    if(violin.style==1) {
      p <- p0 +
        geom_violin(trim = violin.trim, position = dodge, lwd = 0.25, alpha = 1) + 
        geom_boxplot(width=0.08, outlier.color = NA, position = dodge, lwd = 0.25, show.legend = F, alpha = 1) + 
        stat_summary(fun.y=median, geom="point", size=0.5, color=violin.point_color, position = dodge, show.legend = F) +
        dist_theme + guides(col = guide_legend(nrow=2), fill = guide_legend(nrow = 2)) +ggtitle(plot.title) +
        scale_fill_manual(values = pal) + xlab(x_lab) + ylab(y_lab) + scale_y_continuous(breaks = y_breaks, limits = y_limits) 
    } else if(violin.style==2) {
      p <- p0 +
        geom_violin(trim=violin.trim, fill="white") +
        dist_theme + scale_color_manual(values = pal) + ggtitle(plot.title) +
        stat_summary(fun.y=median, geom="point", size=1, color=violin.point_color, position = dodge, show.legend = F) +
        # stat_summary(fun.data ="mean_sdl",mult=1, geom="pointrange", color="red", position = dodge, show.legend = F) +
        guides(col=guide_legend(nrow=2)) + theme(plot.title = element_text(face="plain", hjust = 0.5, size = 8), legend.key.size = unit(4,'mm')) +
        ylab(y_lab) + xlab(x_lab) + scale_y_continuous(breaks = y_breaks, limits = y_limits)
      
      if(compare_dist) p <- p + stat_compare_means(comparisons = my_comparisons, na.rm = T, size=2, label = "p.signif", method = "wilcox.test")
    }
    
  } else if(plot_type=='density'){
    
  } else if(plot_type=='histogram') {
    
  } else if(plot_type=="summary" & !is.null(summarize_by)){
    
    mratio_sm <- plyr::ddply(mratio_molten, as.formula("~sample.id")
                             , summarize
                             , mean = mean(value, na.rm=T)
                             , median = median(value, na.rm=T)
                             , min = min(value, na.rm=T)
                             , max = max(value, na.rm=T)
                             , sd  = sd(value, na.rm=T)
                             , l = length(value)
                             , se = sd(value, na.rm=T)/sqrt(length(value))
                             , se.md = 1.253*sd(value, na.rm=T)/sqrt(length(value)))
    mratio_sm <- merge(mratio_sm, metadata, by="sample.id")
    mratio_sm <- ddply(mratio_sm, as.formula(paste0("~",summarize_by))
                       , mutate
                       , mean.mean = mean(mean, na.rm=T)
                       , se.mean = sd(mean, na.rm=T)/sqrt(length(mean))
                       , sd.mean = sd(mean, na.rm=T))
    mratio_sm$summarize_by <- mratio_sm[,summarize_by]
    
    if(!is.null(summary.precedence)) mratio_sm$summarize_by <- factor(mratio_sm$summarize_by, levels = summary.precedence)
    if(!is.null(col.precedence)) mratio_sm[,col_by] <- factor(mratio_sm[,col_by], levels = col.precedence)
    
    p <- ggplot(mratio_sm,aes_string(x="summarize_by", y="mean", fill=col_by)) +  
      geom_col(aes_string(x="summarize_by", y="mean.mean", fill=col_by),position = 'dodge', size = 0.1, col = 'black', width = 0.6) +
      geom_text(aes(x=summarize_by, y=mean.mean,label=round(mean.mean,2)), size=2, nudge_y = 8)+ 
      geom_point(aes_string(x="summarize_by", y="mean", fill=col_by), show.legend=F, size = .5, position = position_jitterdodge(dodge.width = 1.2, jitter.width = 1))+
      dist_theme + scale_fill_manual(values = pal) + ylab("Average methylation [%]") + 
      scale_y_continuous(breaks = y_breaks, limits = y_limits) + ggtitle(plot.title)
    
    if(compare_dist) p <- p + stat_compare_means(data=mratio_sm, mapping = aes_string(x=summarize_by, y="mean", fill=col_by)
                                                 , comparisons = my_comparisons,size=2
                                                 , method = "t.test"
                                                 , inherit.aes = F
                                                 , label.y = c(78, 85, 90,95))
    
  } else {
    stop(message("[!] Invalid plot type. Please, provide one of # violin, boxplot, summary, density, histogram"))
  }
  return(p)
}
