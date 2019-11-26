# Correlation analysis ----
plot_pairwise_correlation <- function(x
                                      , method = "pearson"
                                      , pal    = NULL
                                      , ...) {
  
  mcor <- cor(x, method=method, ...)
  
  if(is.null(pal)) {
    myPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))
    pal <- myPalette(6)
  }
  hname <- paste0(method, " correlation")
  cHM <- ComplexHeatmap::Heatmap(mcor
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
  return(ComplexHeatmap::draw(cHM, heatmap_legend_side = "bottom"))
}

analyze_pairwise_methyl_correlation <- function(mratio, xvar, yvar
                                                , ptitle         = NULL
                                                , pal            = "RdYlBu"
                                                , plot.type      = "density"
                                                , density.limits = c(0,0.001)
                                                , density.na.value = "red")
{
  
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