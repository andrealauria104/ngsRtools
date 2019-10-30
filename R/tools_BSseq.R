# Tools for BS-seq data analysis

# Load external re#sources and libraries
#source("tools_sequence_analysis.R")
require(genomation)
require(xlsx)
require(methylKit)
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
plot_qmtr_ratio <- function(q.mtr.ratio)
{
  #source("theme_setting.R")
  qm <- reshape2::melt(q.mtr.ratio, varnames = c('methylation','sample'))
  nm <- rev(levels(qm$methylation))
  qm$methylation <- factor(qm$methylation, levels = nm)
  qm$sample <- factor(qm$sample)
  # levels(qm$sample) <-  c('r_D12_2i-LIF','r_D7_2i-LIF'
  #                                 ,'S3 -/-_2i-LIF', 'S3 +/+_2i-LIF','S3 +/+_serum-LIF'
  #                                 ,'S3 +/+_2i')

  pal <- rev(RColorBrewer::brewer.pal(5, 'Reds'))
  bars <- ggplot(qm, aes(x=sample, y=value, fill=methylation)) +
    geom_bar(width = 0.8, stat = "identity", col = 'black') +
    scale_fill_manual(values = pal) + my_theme + theme_bw() +
    ylab('Number of CpG') + xlab(NULL)

  return(bars)
}

# Principal Component Analysis ====
get_meth_pca   <- function(meth){

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
  eigs <- pca$sdev^2
  pca_var  <- rapply(as.list(eigs), function(x) sum(x/sum(eigs)))
  pca_plot <- as.data.frame(pca$x[,c("PC1", "PC2", "PC3", "PC4")])
  pca_plot$sample <- gsub("_R+\\d|\\.+\\d+$","",rownames(pca_plot))
  pal <- function(x) (viridis::viridis(x))
  ggplot(pca_plot, aes(x=PC1, y=PC2, col=sample)) + geom_point(size=3) +
    xlab(paste0("PC1 (",round(pca_var[1]*100,1),"%)")) +
    ylab(paste0("PC2 (",round(pca_var[2]*100,1),"%)")) +
    theme_bw() + my_theme + ggtitle("CpG Methylation PCA") +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face="bold", hjust = 0.5),aspect.ratio = 1) +
    scale_color_manual(values=pal(length(unique(pca_plot$sample)))) +
    geom_label_repel(aes(label = pca_plot$sample),
                     fontface = 'bold', color = 'black', size=3,
                     box.padding = 0.35, point.padding = 0.5,
                     segment.color = 'grey50')
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

analyze_pairwise_methyl_correlation <- function(mratio, xvar, yvar, ptitle=NULL)
{
  #source("theme_setting.R")
  cor.test.res <- cor.test(x = mratio[,xvar], y = mratio[,yvar])
  p <- ggplot(as.data.frame(mratio), aes_string(x=xvar, y=yvar) ) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_distiller(palette="RdYlBu", direction=-1, limits=c(0,0.002), na.value = "red") +
    # scale_fill_brewer(palette = "RdBu")+
    # scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(6,"PuOr")[c(1,3,5)])
    #                       , guide = guide_colourbar(barheight = 0.6, title.vjust = 1 )) +
    my_theme + ggtitle(ptitle) +
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
