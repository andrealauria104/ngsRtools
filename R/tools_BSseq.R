# Tools for BS-seq data analysis

# Load external resources and libraries
#source("tools_sequence_analysis.R")
require(genomation)
require(xlsx)
require(methylKit)
# Load and save data ====
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
  #source("theme_setting.R")
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
  require(plyr)
  
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
  require(GenomicRanges)
  
  convert_input <- function(mtr)
  {
    if(grepl('methyl', class(mtr))) {
      mtrranges <- as(mtr,"GRanges")
    } else if(grepl('data.frame', class(mtr))) {
      mtrranges <- makeGRangesFromDataFrame(mtr)
    } else if(grepl('GRanges', class(mtr))) {
      mtrranges <- mtr
    }
    return(mtrranges)
  }
  mtr_x <- convert_input(mtr_x)
  mtr_y <- convert_input(mtr_y)
  
  ov <- findOverlaps(mtr_x, mtr_y, ignore.strand=T, type = 'any')
  x <- mtr_x[queryHits(ov),]
  y <- mtr_y[subjectHits(ov),]
  
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

# Principal Component Analysis ====
get_meth_pca   <- function(meth
                           , pal = NULL
                           , labels = T
                           , scree_plot = F
                           , max.pcs = NULL
                           , scree_plot_type = "explained_variance"
                           , point_size = 2){

  #source("theme_setting.R")

  require(ggplot2)
  require(ggrepel)
  require(reshape2)
  require(methylKit)

  pca  <- PCASamples(meth,
                     screeplot=F,
                     filterByQuantile = T,
                     sd.threshold=0.5,
                     obj.return = T)
  if((is.null(max.pcs) || max.pcs > ncol(RNA_pca$x))) {
    max.pcs <- ncol(pca$x)
  }
  
  eigs <- pca$sdev^2
  pca_var  <- rapply(as.list(eigs), function(x) sum(x/sum(eigs)))
  pca_plot <- as.data.frame(pca$x[,c("PC1", "PC2", "PC3", "PC4")])
  pca_plot$sample <- gsub("_R+\\d|\\.+\\d+$|_rep_+\\d+","",rownames(pca_plot))
  pca_plot$sample <- gsub("_"," ",pca_plot$sample)
  pca_summary <- summary(pca)$importance
  
  if(is.null(pal)) {
    pal <- ggsci::pal_npg()
    pal <- pal(length(unique(pca_plot$sample)))
  }
  
  p <- ggplot(pca_plot, aes(x=PC1, y=PC2, col=sample)) + geom_point(size=point_size) +
    xlab(paste0("PC1 (",round(pca_var[1]*100,1),"%)")) +
    ylab(paste0("PC2 (",round(pca_var[2]*100,1),"%)")) +
    theme_bw() + my_theme + ggtitle("CpG Methylation PCA") +
    theme(panel.grid.minor = element_blank()
          , plot.title = element_text(face="bold", hjust = 0.5, size=10)
          , aspect.ratio = 1) +
    scale_color_manual(values=pal) 
  
  if(labels) {
    rownames(pca_plot) <- gsub("_"," ",rownames(pca_plot))
    rownames(pca_plot) <- gsub(" rep "," - rep",rownames(pca_plot))
    p <- p + geom_label_repel(aes(label = rownames(pca_plot), col=sample),
                              fontface = 'bold'
                              , show.legend = F
                              # , color = 'black'
                              , size=2,
                              box.padding = 0.35, point.padding = 0.5,
                              segment.color = 'grey50') 
  }
  
  if(scree_plot) {
    
    if(scree_plot_type=='explained_variance') {
      spdata      <- reshape2::melt(pca_summary[c(2:3),])
      spdata$Var2 <- as.numeric(gsub("PC", "",spdata$Var2))
      spdata      <- spdata[spdata$Var2<=max.pcs,]
      sp <- ggplot(spdata, aes(x=Var2, y=value, group=Var1, col=Var1)) + geom_point(size=1) + geom_line() +
        theme_bw() + my_theme + scale_x_continuous(breaks = spdata$Var2) + 
        geom_hline(yintercept = 0.9, lwd=0.2, linetype="dashed", col = "darkred") +
        xlab("Principal Component") + ylab("Explained Variance") + ggtitle("PCA - Scree Plot") + 
        theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title = element_blank(), panel.grid = element_blank()) +
        scale_y_continuous(breaks = seq(0,1,0.1),labels = scales::percent, limits = c(0,1)) + scale_color_manual(values=c("black","grey"))
    } else if(scree_plot_type=='standard_deviation') {
      spdata      <- reshape2::melt(pca_summary[1,])
      
      spdata$Var2 <- as.numeric(gsub("PC", "",rownames(spdata)))
      spdata      <- spdata[spdata$Var2<=max.pcs,]  
      
      sp <- ggplot(spdata, aes(x=Var2, y=value), col='black') + geom_point(size=1) + geom_line() +
        theme_bw() + my_theme + scale_x_continuous(breaks = spdata$Var2) + 
        geom_hline(yintercept = 0.9, lwd=0.2, linetype="dashed", col = "darkred") +
        xlab("Principal Component") + ylab("Standard Deviation") + ggtitle("PCA - Scree Plot") + 
        theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title = element_blank(), panel.grid = element_blank())
      
    }
    
    
    
    return(list("pca" = p, "screeplot" = sp))
  } else {
    return(p) 
  }
}

# Correlation analysis ----
plot_pairwise_correlation <- function(x
                                      , method = "pearson"
                                      , pal    = NULL
                                      , ...) {

  require(ComplexHeatmap)
  require(RColorBrewer)
  mcor <- cor(x, method=method, ...)

  if(is.null(pal)) {
    myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
    pal <- myPalette(6)
  }
  hname <- paste0(method, " correlation")
  cHM <- Heatmap(mcor
                 , col  = pal
                 , cell_fun = function(j, i, x, y, w, h, col) {
                   grid.text(round(mcor[i, j], digits = 2), x, y,gp = gpar(col='black', fontsize=5))
                 }
                 , name = hname
                 , row_names_gp = gpar(fontsize = 6)
                 , column_names_gp = gpar(fontsize = 6)
                 , heatmap_legend_param = list(title_position = "topcenter",
                                               # legend_width  = unit(4, "cm"),
                                               # legend_height = unit(0.5, "mm"),
                                               values_gp     = gpar(fontsize=8),
                                               legend_direction = "horizontal"))
  return(draw(cHM, heatmap_legend_side = "bottom"))
}

analyze_pairwise_methyl_correlation <- function(mratio, xvar, yvar
                                                , ptitle         = NULL
                                                , pal            = "RdYlBu"
                                                , plot.type      = "density"
                                                , density.limits = c(0,0.001)
                                                , density.na.value = "red")
{
  #source("theme_setting.R")
  cor.test.res <- cor.test(x = mratio[,xvar], y = mratio[,yvar])
  if(plot.type=="density") {
    p0 <- ggplot(as.data.frame(mratio), aes_string(x=xvar, y=yvar) ) +
      stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_distiller(palette=pal, direction=-1, limits=density.limits, na.value = density.na.value) +
      my_theme
  } else if(plot.type=="scatter") {
    p0 <- ggplot(as.data.frame(mratio), aes_string(x=xvar, y=yvar) ) +
      geom_point(size = 1, alpha = 0.4) + theme_bw() + my_theme_2
  }
  p <- p0 + ggtitle(ptitle) +
    xlab(paste0("% CG methylation in ",gsub(".*_","#",xvar))) +
    ylab(paste0("% CG methylation in ",gsub(".*_","#",yvar))) +
    theme(plot.title = element_text(size = 8, hjust = 0.5)
          , aspect.ratio = 1
          , legend.position = 'none'
    ) +
    annotate("text", label = paste0("r = ",round(cor.test.res$estimate,2)), x = 30, y = 70, size = 2, col = "white")
  print(cor.test.res)
  return(list("p" = p,"cor.test.res"=cor.test.res))
}
# Differential Methylation Analysis ====
# Volcano plot
plot_diffmC_volcano <- function(diffMeth
                                , title = NULL
                                , q_th  = 0.05
                                , m_th  = 10
                                , pal = NULL
                                )
{

  #source("theme_setting.R")
  require(ggrepel)

  if(is.null(pal)) {
    pal    = c("orange","darkblue","grey")
  }

  diffMeth$status <- "none"
  diffMeth$status[which(diffMeth$qvalue<=q_th & diffMeth$meth.diff>=m_th)]   <- "Gain"
  diffMeth$status[which(diffMeth$qvalue<=q_th & diffMeth$meth.diff<=(-m_th))] <- "Loss"
  diffMeth$status <- as.factor(diffMeth$status)

  counts <- table(diffMeth$status)

  diffMeth$label <- as.character(counts['none'])
  diffMeth$label[which(diffMeth$status=='Gain')] <- as.character(counts['Gain'])
  diffMeth$label[which(diffMeth$status=='Loss')] <- as.character(counts['Loss'])

  if(grepl('methyl', class(diffMeth))) {
    diffMeth <- getData(x = diffMeth)
  }

  vp <- ggplot(data=subset(diffMeth,status=='none'), aes(x=meth.diff, y=-log10(qvalue), col=status)) + geom_point(size=0.5, alpha=0.1) +
    geom_point(data=subset(diffMeth,status!='none'), aes(x=meth.diff, y=-log10(qvalue), col=status),size=0.5, alpha=0.5) +
    geom_hline(yintercept = -log10(q_th), linetype = 'dashed', lwd = 0.25) +
    geom_vline(xintercept = m_th, linetype = 'dashed', lwd = 0.25) +
    geom_vline(xintercept = -m_th, linetype = 'dashed', lwd = 0.25) +
    theme_bw() + my_theme + theme(plot.title = element_text(hjust = 0.5, size = 8)) +
    scale_color_manual(values = pal) + ggtitle(title) + xlab("Methylation difference (%)") +
    # scale_y_log10(limits=c(1,min(log10(diffMeth$qvalue)))) +
    scale_y_log10() +
    scale_x_continuous(breaks = c(-50,-10, 10, 50), limits=c(-100,100)) +
    annotate("label"
             , label = as.character(counts['Loss'])
             , x = -30, y = 100
             , size = 3
             , colour = pal[2]) +
    annotate("label"
             , label = as.character(counts['Gain'])
             , x = 30, y = 100
             , size = 3, colour = pal[1])

  return(vp)
}

# analyze regions
get_diffmeth_region <- function(x, region
                                , m_th = 10
                                , q_th = 0.05
                                , diff_only = F
                                , ...)
{
  require(methylKit)

  # Count C sites per region
  cregion <- regionCounts(x, region, strand.aware = F)
  myDiff    <- calculateDiffMeth(cregion, ...)

  if(diff_only){
    myDiff10p <- getMethylDiff(myDiff           ,
                               difference = m_th,
                               qvalue     = q_th )
    diffregion <- as(myDiff10p, "data.frame")

  } else {
    diffregion <- as(myDiff, "data.frame")

  }

  diffregion <- get_names(diffregion, region)
  idx <- which(duplicated(diffregion$names))

  if(length(idx)>0) diffregion <- diffregion[-idx,]

  return(diffregion)
}

get_dmr <- function(dmr
                    , pth = 0.05
                    , mth = 10
                    , type = 'all' # 'hyper', 'hypo'
)
{

  get_dmt <- function(x, pth, mth, type)
  {
    if(type == 'all') {
      subset(getData(x), qvalue<=pth & abs(meth.diff)>mth)
    } else if (type == 'hyper') {
      subset(getData(x), qvalue<=pth & meth.diff>mth)
    } else if (type == 'hypo') {
      subset(getData(x), qvalue<=pth & meth.diff<(-mth))
    } else {
      stop(message('[-] Incorrect dmr setting (all/hyper/hypo)'))
    }
  }
  message('[+] Get DMR')
  message(' - meth difference = ', mth)
  message(' - qvalue = ', pth)

  if( is.list(dmr) & !is.data.frame(dmr) ) {
    lapply(dmr, function(i) {
      if( is.list(i) & !is.data.frame(i) ) {
        lapply(i, get_dmt, pth = pth, mth = mth, type = type)
      } else {
        get_dmt(i, pth = pth, mth = mth, type = type)
      }
    })
  } else {
    get_dmt(dmr, pth = pth, mth = mth, type = type)
  }
}
# Visalize Differentially methylated regions
get_heatmap_bs <- function(m
                         , annotDF  = NULL
                         , annotCol = NULL
                         , fig_out  = NULL
                         , retHm    = F
                         , myPalette = NULL
                         , ...)
{

  require(ComplexHeatmap)
  require(circlize)
  require(RColorBrewer)

  # base_mean <- rowMeans(m)
  # m_scaled <- t(apply(m, 1, scale))
  # colnames(m_scaled) <- colnames(m)
  # m <- m_scaled

  # bPalette <- colorRampPalette(brewer.pal(4, "Reds"))
  # bPalette <- colorRampPalette(brewer.pal(4, "YlOrBr"))
  # bPalette <- colorRampPalette(brewer.pal(5, "Blues"))(5)
  # bPalette <- colorRampPalette(brewer.pal(5, "Reds"))(5)
  bPalette <- c('black','red4','red3','red2','red1','red')
  # bPalette <- c('black','red3','red1')
  if(is.null(myPalette)) myPalette <- c("#7F00FF","black","#FFFF33")
  # ramp <- colorRamp2(c(0, 20, 40, 60, 100), c('white',bPalette(4)))
  # ramp <- colorRamp2(c(0, 50, 70, 90, 100), bPalette)
  ramp <- colorRamp2(c(0, 50, 60, 70, 90, 100), bPalette)
  # ramp <- colorRamp2(c(0,1,2), bPalette)

  if (!is.null(annotDF)) {
    if (!is.null(annotCol)) {
      ha_column <- HeatmapAnnotation(df  = annotDF,
                                     col = annotCol,
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8),
                                                                    legend_direction = "horizontal"))
    } else {
      ha_column <- HeatmapAnnotation(df  = annotDF,
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8),
                                                                    legend_direction = "horizontal"))
    }
  } else {
    ha_column <- new("HeatmapAnnotation")
  }

  hm <- Heatmap(m, col = ramp,
                # show_row_dend = T,
                # row_names_side = "left",
                row_names_gp = gpar(fontsize=8),
                column_names_gp = gpar(fontsize=8),
                column_title_gp = gpar(fontsize=10, fontface="bold"),
                heatmap_legend_param = list(title = "% methylation",
                                            title_gp = gpar(fontsize=8),
                                            title_position = "topcenter",
                                            # legend_width  = unit(4, "cm"),
                                            # legend_height = unit(0.5, "mm"),
                                            values_gp     = gpar(fontsize=8),
                                            legend_direction = "horizontal")
                , top_annotation = ha_column
                , top_annotation_height = unit(4, "mm")
                , width = unit(3,'cm')
                # , width = unit(2,'cm')
                , rect_gp = gpar(col = 'black', lwd=0.25 )
                , ...)

  hmOut <- hm

  if(!is.null(fig_out)){
    pdf(file = fig_out, useDingbats = F, h=8, w=3, paper = "a4")
    draw(hmOut, heatmap_legend_side = "bottom")
    dev.off()
  } else{
    draw(hmOut, heatmap_legend_side = "bottom")
  }

  if(retHm) return(hmOut)
}


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
