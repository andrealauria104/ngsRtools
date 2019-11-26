# Differential Methylation Analysis ====
# Volcano plot
plot_diffmC_volcano <- function(diffMeth
                                , title = NULL
                                , q_th  = 0.05
                                , m_th  = 10
                                , pal = NULL)
{
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
  # Count C sites per region
  cregion <- methylKit::regionCounts(x, region, strand.aware = F)
  myDiff    <- methylKit::calculateDiffMeth(cregion, ...)
  
  if(diff_only){
    myDiff10p <- methylKit::getMethylDiff(myDiff,
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
      subset(methylKit::getData(x), qvalue<=pth & abs(meth.diff)>mth)
    } else if (type == 'hyper') {
      subset(methylKit::getData(x), qvalue<=pth & meth.diff>mth)
    } else if (type == 'hypo') {
      subset(methylKit::getData(x), qvalue<=pth & meth.diff<(-mth))
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
  
  # base_mean <- rowMeans(m)
  # m_scaled <- t(apply(m, 1, scale))
  # colnames(m_scaled) <- colnames(m)
  # m <- m_scaled
  
  # bPalette <- colorRampPalette(RColorBrewer::brewer.pal(4, "Reds"))
  # bPalette <- colorRampPalette(RColorBrewer::brewer.pal(4, "YlOrBr"))
  # bPalette <- colorRampPalette(RColorBrewer::brewer.pal(5, "Blues"))(5)
  # bPalette <- colorRampPalette(RColorBrewer::brewer.pal(5, "Reds"))(5)
  bPalette <- c('black','red4','red3','red2','red1','red')
  # bPalette <- c('black','red3','red1')
  if(is.null(myPalette)) myPalette <- c("#7F00FF","black","#FFFF33")
  # ramp <- circlize::colorRamp2(c(0, 20, 40, 60, 100), c('white',bPalette(4)))
  # ramp <- circlize::colorRamp2(c(0, 50, 70, 90, 100), bPalette)
  ramp <- circlize::colorRamp2(c(0, 50, 60, 70, 90, 100), bPalette)
  # ramp <- circlize::colorRamp2(c(0,1,2), bPalette)
  
  if (!is.null(annotDF)) {
    if (!is.null(annotCol)) {
      ha_column <- ComplexHeatmap::HeatmapAnnotation(df  = annotDF,
                                     col = annotCol,
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8),
                                                                    legend_direction = "horizontal"))
    } else {
      ha_column <- ComplexHeatmap::HeatmapAnnotation(df  = annotDF,
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8),
                                                                    legend_direction = "horizontal"))
    }
  } else {
    ha_column <- new("HeatmapAnnotation")
  }
  
  hm <- ComplexHeatmap::Heatmap(m, col = ramp,
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
    ComplexHeatmap::draw(hmOut, heatmap_legend_side = "bottom")
    dev.off()
  } else{
    ComplexHeatmap::draw(hmOut, heatmap_legend_side = "bottom")
  }
  
  if(retHm) return(hmOut)
}

