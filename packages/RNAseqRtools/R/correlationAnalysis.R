# Correlation Analysis ====
plotCorrelation <- function(x
                            , method = "pearson"
                            , myPalette=NULL
                            , ...) {
  
  mcor <- cor(x, method=method, ...)
  
  if(is.null(myPalette)) {
    myPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))
  }
  hname <- paste0(method, " correlation")
  cHM <- ComplexHeatmap::Heatmap(mcor
                                 , col  = myPalette(6)
                                 , cell_fun = function(j, i, x, y, w, h, col) {
                                   grid.text(round(mcor[i, j], digits = 2), x, y,gp = gpar(col='black', fontsize=6))
                                 }
                                 , name = hname
                                 , row_names_gp = gpar(fontsize = 8)
                                 , column_names_gp = gpar(fontsize = 8)
                                 , heatmap_legend_param = list(title_position = "topcenter",
                                                               # legend_width  = unit(4, "cm"),
                                                               # legend_height = unit(0.5, "mm"),
                                                               values_gp     = gpar(fontsize=8),
                                                               legend_direction = "horizontal"))
  return(ComplexHeatmap::draw(cHM, heatmap_legend_side = "bottom"))
}
