# Correlation Analysis ====
reorder_cormat <- function(mcor, reorder_cormat_method = "ward.D")
{
  # Use correlation between variables as distance
  message("   - hclust method: ", reorder_cormat_method)
  dd <- as.dist((1-mcor))
  hc <- hclust(dd, method = reorder_cormat_method)
  mcor <- mcor[hc$order, hc$order]
  return(list(mcor=mcor,hc=hc))
}

plotCorrelation <- function(x
                            , method = "pearson"
                            , reorder = TRUE
                            , reorder_cormat_by_cor_dist = FALSE
                            , reorder_cormat_by_cor_dist_method = "ward.D"
                            , myPalette = NULL
                            , myBrewerPal = "Reds"
                            , draw_cor_values = T
                            , cor_values_size = 6
                            , annotDF  = NULL
                            , annotCol = NULL
                            , show_annotation_name = T
                            , annotation_width_heigth = 2
                            , annotation_nrow = 1
                            , annotation_ncol = NULL
                            , show_column_annotation = T
                            , simple_anno_size_mm = 2
                            , text_size = 8
                            , use_raster = F
                            , split = NULL
                            , column_split = NULL
                            , column_title = character(0)
                            , ...) {
  
  mcor <- cor(x, method=method, ...)
  
  if(is.null(myPalette)) {
    myPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, myBrewerPal))
  }
  if(draw_cor_values) {
    cell_fun <- function(j, i, x, y, w, h, col) {
      grid.text(round(mcor[i, j], digits = 2), x, y,gp = gpar(col='black', fontsize=cor_values_size))
    }
  } else {
    cell_fun <- NULL
  }
  if (!is.null(annotDF)) {
    if (!is.null(annotCol)) { 
      if(show_column_annotation) {
        ha_column <- ComplexHeatmap::HeatmapAnnotation(df  = annotDF 
                                                       , col = annotCol 
                                                       , annotation_legend_param = list(title_gp  = gpar(fontsize=text_size),
                                                                                        labels_gp = gpar(fontsize=text_size))
                                                       , height = unit(annotation_width_heigth, "mm")
                                                       , show_annotation_name = show_annotation_name
                                                       , which = "column"
                                                       , show_legend = F
                                                       , annotation_name_gp = gpar(fontsize=text_size)
                                                       , simple_anno_size = unit(simple_anno_size_mm, "mm"))
      } else {
        ha_column <- NULL
      }
      
      ha_row <- ComplexHeatmap::rowAnnotation(df  = annotDF
                                              , col = annotCol
                                              , annotation_legend_param = list(title_gp  = gpar(fontsize=text_size),
                                                                             labels_gp = gpar(fontsize=text_size),
                                                                             nrow = annotation_nrow, ncol = annotation_ncol)
                                              , width = unit(annotation_width_heigth, "mm")
                                              , show_legend = T
                                              , show_annotation_name = show_annotation_name
                                              , annotation_name_gp = gpar(fontsize=text_size)
                                              , simple_anno_size = unit(simple_anno_size_mm, "mm")
                                              )
    } else {
      if(show_column_annotation) {
        ha_column <- ComplexHeatmap::HeatmapAnnotation(df  = annotDF 
                                                       , annotation_legend_param = list(title_gp  = gpar(fontsize=text_size),
                                                                                        labels_gp = gpar(fontsize=text_size))
                                                       , height = unit(annotation_width_heigth, "mm")
                                                       , show_annotation_name = show_annotation_name
                                                       , which = "column"
                                                       , show_legend = F
                                                       , annotation_name_gp = gpar(fontsize=text_size)
                                                       , simple_anno_size = unit(simple_anno_size_mm, "mm"))
      } else {
        ha_column <- NULL
      }
      
      ha_row <- ComplexHeatmap::rowAnnotation(df  = annotDF
                                              , annotation_legend_param = list(title_gp  = gpar(fontsize=text_size),
                                                                               labels_gp = gpar(fontsize=text_size),
                                                                                nrow = annotation_nrow, ncol = annotation_ncol)
                                              , width = unit(annotation_width_heigth, "mm")
                                              , show_legend = T
                                              , show_annotation_name = show_annotation_name
                                              , annotation_name_gp = gpar(fontsize=text_size)
                                              , simple_anno_size = unit(simple_anno_size_mm, "mm"))
      }
    } else {
      ha_column <- NULL
      ha_row <- NULL
  }
  
  if(reorder) {
    if(reorder_cormat_by_cor_dist) {
      message(" -- reordering correlation matrix by hclust on correlation distance (1-r).")
      reorder_cormat_out <- reorder_cormat(mcor, reorder_cormat_method = reorder_cormat_by_cor_dist_method)
      cluster_columns <- as.dendrogram(reorder_cormat_out$hc)
      cluster_rows <- as.dendrogram(reorder_cormat_out$hc)
    } else {
      message(" -- reordering correlation matrix by ComplexHeatmap default hclust on correlation values.")
      cluster_columns <- T
      cluster_rows <- T
    }
  } else {
    cluster_columns <- F
    cluster_rows <- F
  }
  
  hm_name <- paste0(method, " correlation")
  hm <- ComplexHeatmap::Heatmap(mcor
                                 , col  = myPalette(6)
                                 , cell_fun = cell_fun
                                 , name = hm_name
                                 , use_raster = use_raster
                                 , raster_device = "tiff"
                                 , raster_quality = 10
                                 , row_names_gp = gpar(fontsize = text_size)
                                 , column_names_gp = gpar(fontsize = text_size)
                                 , row_title_gp = gpar(fontsize = text_size, fontface = "plain")
                                 , column_title_gp = gpar(fontsize = text_size, fontface = "plain")
                                 , heatmap_legend_param = list(title_gp = gpar(fontsize = text_size),
                                                               title_position = "topcenter",
                                                               legend_width  = unit(2.5, "cm"),
                                                               legend_height = unit(0.5, "mm"),
                                                               values_gp = gpar(fontsize = text_size), 
                                                               labels_gp = gpar(fontsize = text_size),
                                                               legend_direction = "horizontal")
                                , split = split
                                , column_split = column_split
                                , top_annotation = ha_column
                                , left_annotation = ha_row
                                , column_title = column_title
                                , cluster_rows = cluster_rows
                                , cluster_columns = cluster_columns)
  
  return(ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom"))
}
