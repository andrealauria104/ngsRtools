# Heatmap ====
get_heatmap3 <- function(m
                         , annotDF  = NULL
                         , annotCol = NULL
                         , fig_out  = NULL
                         , retHm    = F
                         , bm = T
                         , myPalette = NULL
                         , myZscale = NULL
                         , myLegend = NULL
                         , scale = T
                         , text_size = 8
                         , ...){
  
  base_mean <- rowMeans(m)
  if(scale) {
    m_scaled <- t(apply(m, 1, scale))
    colnames(m_scaled) <- colnames(m)
  } else {
    m_scaled <- m
  }
  
  bPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))
  
  if(is.null(myPalette)) {
    myPalette <- c("blue","black","red")
  }
  
  if(is.null(myZscale)) {
    myZscale <- c(-2, 0, 2)
  }
  ramp <- circlize::colorRamp2(myZscale, myPalette)
  
  
  if (!is.null(annotDF)) {
    if (!is.null(annotCol)) {
      ha_column <- ComplexHeatmap::HeatmapAnnotation(df  = annotDF, 
                                                     col = annotCol, 
                                                     annotation_legend_param = list(title_gp  = gpar(fontsize=text_size),
                                                                                    labels_gp = gpar(fontsize=text_size)))
    } else {
      ha_column <- ComplexHeatmap::HeatmapAnnotation(df  = annotDF, 
                                                     annotation_legend_param = list(title_gp  = gpar(fontsize=text_size),
                                                                                    labels_gp = gpar(fontsize=text_size)))
    }
  } else {
    ha_column <- new("HeatmapAnnotation")
  }
  
  if(is.null(myLegend)) myLegend <- "RPKM" 
  
  if(scale) {
    myLegend_title <- paste0("Z-score (",myLegend,")")
  } else {
    myLegend_title <- myLegend
  }
  
  hm <- ComplexHeatmap::Heatmap(m_scaled, col = ramp,
                                # show_row_dend = T,
                                # row_names_side = "left",
                                row_names_gp = gpar(fontsize=text_size),
                                column_names_gp = gpar(fontsize=text_size),
                                column_title_gp = gpar(fontsize=10, fontface="bold"),
                                heatmap_legend_param = list(title = myLegend_title,
                                                            title_gp = gpar(fontsize=text_size),
                                                            title_position = "topcenter",
                                                            legend_width  = unit(3, "cm"),
                                                            legend_height = unit(0.5, "mm"),
                                                            values_gp     = gpar(fontsize=text_size),
                                                            labels_gp     = gpar(fontsize=text_size),
                                                            # legend_direction = "vertical"
                                                            legend_direction = "horizontal"
                                )
                                , top_annotation = ha_column
                                , top_annotation_height = unit(4, "mm")
                                , ...)
  
  if(bm) {
    bmscale <- summary(base_mean)
    bmramp <- circlize::colorRamp2(c(bmscale[1],bmscale[3],bmscale[5]), bPalette(3))
    bmh <- ComplexHeatmap::Heatmap(base_mean 
                                   # , name = "Mean Expression"
                                   , column_names_gp = gpar(fontsize=text_size)
                                   , show_row_names = FALSE
                                   , width = unit(3, "mm") 
                                   , col = bmramp
                                   , heatmap_legend_param = list(title = paste0("Average ",myLegend)
                                                                 ,title_gp = gpar(fontsize=text_size),
                                                                 title_position = "topcenter",
                                                                 legend_width  = unit(3, "cm"),
                                                                 legend_height = unit(0.5, "mm"),
                                                                 values_gp     = gpar(fontsize=text_size),
                                                                 labels_gp     = gpar(fontsize=text_size),
                                                                 # legend_direction = "vertical"
                                                                 legend_direction = "horizontal"))
    
    hmOut <- hm + bmh
  } else {
    hmOut <- hm 
  }
  
  if(!is.null(fig_out)){
    pdf(file = fig_out, useDingbats = F, h=8, w=3, paper = "a4")
    ComplexHeatmap::draw(hmOut, heatmap_legend_side = "right")
    dev.off()
  } else{
    ComplexHeatmap::draw(hmOut, heatmap_legend_side = "right")
  }
  
  if(retHm) return(hmOut)
}

get_heatmap4 <- function(m
                         , annotDF  = NULL
                         , annotCol = NULL
                         , show_annotation_name = T
                         , annotation_width_heigth = 2
                         , annotation_nrow = 1
                         , annotation_ncol = NULL
                         , fig_out  = NULL
                         , retHm    = F
                         , bm = T
                         , myPalette = NULL
                         , myZscale = NULL
                         , myLegend = "CPM"
                         , scale = T
                         , text_size = 8
                         , use_raster = F
                         , raster_quality = 10
                         , ...){
  
  base_mean <- rowMeans(m)
  if(scale) {
    m_scaled <- t(apply(m, 1, scale))
    colnames(m_scaled) <- colnames(m)
  } else {
    m_scaled <- m
  }
  
  bPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))
  
  if(is.null(myPalette)) {
    myPalette <- c("blue","black","red")
  }
  
  if(is.null(myZscale)) {
    myZscale <- c(-2, 0, 2)
  }
  ramp <- circlize::colorRamp2(myZscale, myPalette)
  
  if(is.null(myLegend)) myLegend <- gsub("_"," ",assay.type)
  
  if(scale) {
    myLegend_title <- paste0("Z-score (",myLegend,")")
  } else {
    myLegend_title <- myLegend
  }
  if (!is.null(annotDF)) {
    if (!is.null(annotCol)) {
      ha_column <- ComplexHeatmap::HeatmapAnnotation(df  = annotDF 
                                                     , col = annotCol 
                                                     , annotation_legend_param = list(title_gp  = gpar(fontsize=text_size),
                                                                                      labels_gp = gpar(fontsize=text_size))
                                                     , height = unit(annotation_width_heigth, "mm")
                                                     , show_annotation_name = show_annotation_name
                                                     , which = "column"
                                                     , show_legend = T
                                                     , annotation_name_gp = gpar(fontsize=text_size))
    } else {
      ha_column <- ComplexHeatmap::HeatmapAnnotation(df  = annotDF 
                                                     , annotation_legend_param = list(title_gp  = gpar(fontsize=text_size),
                                                                                      labels_gp = gpar(fontsize=text_size))
                                                     , height = unit(annotation_width_heigth, "mm")
                                                     , show_annotation_name = show_annotation_name
                                                     , which = "column"
                                                     , show_legend = T
                                                     , annotation_name_gp = gpar(fontsize=text_size))
    }
  } else {
    ha_column <- NULL
  }
  hm <- ComplexHeatmap::Heatmap(m_scaled
                                , col = ramp
                                , heatmap_legend_param = list(title = myLegend_title,
                                                              title_gp = gpar(fontsize=text_size),
                                                              title_position = "topcenter",
                                                              legend_width  = unit(3, "cm"),
                                                              legend_height = unit(0.5, "mm"),
                                                              values_gp     = gpar(fontsize=text_size),
                                                              labels_gp     = gpar(fontsize=text_size),
                                                              legend_direction = "horizontal"
                                )
                                , top_annotation = ha_column
                                , use_raster = use_raster
                                , raster_device = "tiff"
                                , raster_quality = raster_quality
                                , row_names_gp = gpar(fontsize=text_size)
                                , column_names_gp = gpar(fontsize=text_size)
                                , column_title_gp = gpar(fontsize=text_size, fontface="plain")
                                , ...)
  
  if(bm) {
    bmscale <- summary(base_mean)
    bmramp <- circlize::colorRamp2(c(bmscale[1],bmscale[3],bmscale[5]), bPalette(3))
    bmh <- ComplexHeatmap::Heatmap(base_mean 
                                   , name = "Average Expression"
                                   , column_names_gp = gpar(fontsize=text_size)
                                   , show_row_names = FALSE
                                   , width = unit(3, "mm") 
                                   , col = bmramp
                                   , heatmap_legend_param = list(title = paste0("Average ",myLegend)
                                                                 ,title_gp = gpar(fontsize=text_size),
                                                                 title_position = "topcenter",
                                                                 legend_width  = unit(3, "cm"),
                                                                 legend_height = unit(0.5, "mm"),
                                                                 values_gp     = gpar(fontsize=text_size),
                                                                 labels_gp     = gpar(fontsize=text_size),
                                                                 # legend_direction = "vertical"
                                                                 legend_direction = "horizontal"))
    
    hmOut <- hm + bmh
  } else {
    hmOut <- hm 
  }
  
  if(!is.null(fig_out)){
    pdf(file = fig_out, useDingbats = F, h=8, w=3, paper = "a4")
    ComplexHeatmap::draw(hmOut, heatmap_legend_side = "right")
    dev.off()
  } else{
    ComplexHeatmap::draw(hmOut, heatmap_legend_side = "right")
  }
  
  if(retHm) return(hmOut)
}

get_clusters <- function(m, hm){
  # Retrieve clusters from K-means
  clusters <- lapply(ComplexHeatmap::row_order(hm), 
                     function(x){
                       rownames(m[x,])
                     }
  )
  names(clusters) <- paste0('cluster_', seq_along(clusters))
  return(clusters)
}


plot_cluster_expression <- function(m, cl, pal)
{
  toplot <- lapply(cl, function(x) y <- m[x,])
  toplot <- reshape2::melt(toplot)
  colnames(toplot) <- c("gene","sample","tpm","cluster")
  toplot <- plyr::ddply(toplot, .(gene, sample, cluster)
                        , summarize
                        , av_tpm = mean(tpm)
                        , log2_av_tpm = log2(mean(tpm)))
  
  toplot$state <- gsub("/.*","",as.character(toplot$sample))
  toplot$cond  <- gsub(".*/| #.*","",as.character(toplot$sample))
  
  toplot$cluster <- gsub("_"," ", toplot$cluster)
  p <- ggplot(toplot, aes(x=sample, y=log2_av_tpm, fill = cond)) + 
    geom_boxplot(notch = T, outlier.colour = "grey", outlier.size = 0.2, outlier.alpha = 0.5) +
    facet_wrap(~cluster, ncol = 1) + theme_bw() + my_theme + ylab("average log2[RPKM]") +
    # coord_cartesian(ylim = quantile(toplot$log2_av_tpm, c(0.1, 0.95))) +
    scale_fill_manual(values = pal) + 
    theme(panel.grid = element_blank()
          , axis.title.x = element_blank()
          , axis.text.x = element_text(angle = 90, hjust = 1)
          , strip.background = element_blank())
  
}

plot_hm_fc <- function(res_df
                       , myPalette = NULL
                       , myFCscale = NULL
                       , myTitle = NULL
                       , ...)
{
  
  
  if(any(grepl("FDR", colnames(res_df)))) {
    idx <- grep("logFC", colnames(res_df))
    fc <- as.matrix(res_df[,idx])
    colnames(fc) <- gsub("logFC.group|NT_|Act_","", colnames(fc))
  } else {
    fc <- res_df
  }
  if(is.null(myPalette)) {
    myPalette <- rev(colorRampPalette(RColorBrewer::brewer.pal(6, "RdBu"))(3))
  }
  
  if(is.null(myFCscale)) {
    myFCscale <- c(-1, 0, 1)
  }
  ramp <- circlize::colorRamp2(myFCscale, myPalette)
  
  if(is.null(myTitle)) {
    myTitle <- "log2[fold-change]"
  }
  hmfc <-  ComplexHeatmap::Heatmap(fc, col = ramp,
                                   # show_row_dend = T,
                                   # row_names_side = "left",
                                   # row_names_gp = gpar(fontsize=8),
                                   column_names_gp = gpar(fontsize=8),
                                   column_title_gp = gpar(fontsize=10, fontface="bold"),
                                   heatmap_legend_param = list(title = myTitle,
                                                               title_gp = gpar(fontsize=8),
                                                               title_position = "topcenter",
                                                               legend_width  = unit(3, "cm"),
                                                               legend_height = unit(0.5, "mm"),
                                                               values_gp     = gpar(fontsize=8),
                                                               # legend_direction = "vertical"
                                                               legend_direction = "horizontal")
                                   , ...)
  
  return(hmfc)
}