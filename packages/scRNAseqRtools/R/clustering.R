# Clustering ----
plot_cinfo <- function(cinfo, structure, feature, pal, info_type = "absolute", mode = "bar")
{
  if(length(feature)>1) {
    nm <- apply(cinfo[,feature], 1, paste, collapse = "_")
    cinfo$sample <- gsub("_","/",nm)
    feature   <- "sample"
  }
  
  if(info_type == "absolute") {
    p <- ggplot(cinfo, aes_string(x=structure, fill=feature)) + geom_bar() +
      theme_bw() + my_theme +
      # theme(legend.title = element_blank()) +
      scale_fill_manual(values = pal) + xlab(structure) + ylab("n. of cells")
  } else if(info_type == "relative") {
    print("here")
    cinfo     <- split(cinfo, cinfo[,feature])
    cinfo     <- lapply(cinfo, function(x)
    {
      x$tot_n_cells = nrow(x)
      return(x)
    } )
    cinfo <- do.call(rbind, cinfo)
    
    cinfo     <- split(cinfo, cinfo[,c(structure,feature)])
    cinfo     <- lapply(cinfo, function(x)
    {
      x$relative_n_cells = nrow(x)/x$tot_n_cells
      return(x)
    } )
    
    cinfo     <- do.call(rbind, cinfo)
    
    if(mode=="reverse_bar") {
      p0 <- ggplot(cinfo, aes_string(x=feature, y = "relative_n_cells", fill=structure)) +
        geom_bar(position = "fill",stat = "identity")
      x_lab <- feature
    } else if(mode=="pie") {
      require(grid)
      mp <- cinfo
      mp$value <- mp$relative_n_cells
      mp <- unique(mp[,c(structure,"sample","value")])
      mp$ypos <- cumsum(mp$value)[4] - cumsum(mp$value) + mp$value/2
      mp$perc <- round(mp$value/sum(mp$value),4)
      pie <- vector(mode = 'list')
      
      for(s in unique(mp[,structure])) {
        
        mps <- mp[mp[,structure]==s,]
        mps$ypos <- cumsum(mps$value)[nrow(mps)] - cumsum(mps$value) + mps$value/2
        mps$perc <- round(mps$value/sum(mps$value),4)
        
        pie[[s]] <- ggplot(mps, aes(x="", y=value, fill=sample)) +
          geom_bar(width = 1, stat = "identity", col = 'black') +
          coord_polar("y", start=0)  + scale_fill_manual(values = pal) +
          blank_theme + ggtitle(s) +
          theme(axis.text.x=element_blank()
                , plot.title = element_text(size = 8, face = "plain", hjust = 0.5, vjust = 0)
                , text = element_text(size = 8)
                , legend.position = 'none'
                , plot.margin=unit(c(0,0,0,0),"cm")
                , panel.spacing = unit(c(0, 0, 0, 0), "cm"),)
      }
      # dev.off()
      if(structure=="State") {
        pushViewport(viewport(layout = grid.layout(4, 3, widths = unit(2,"cm"), heights = unit(2,"cm"))))
        vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
        grid.text("Relative proportion of cells by state", vp = vplayout(1, 2), gp = gpar(fontsize = 8))
        print(pie$A, vp = vplayout(2, 1))
        print(pie$E, vp = vplayout(2, 3))
        print(pie$C, vp = vplayout(3, 2))
        print(pie$B, vp = vplayout(4, 1))
        print(pie$D, vp = vplayout(4, 3))
      } else if(structure=="Cluster") {
        pushViewport(viewport(layout = grid.layout(4, 3, widths = unit(2,"cm"), heights = unit(2,"cm"))))
        vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
        grid.text("Relative proportion of cells by cluster", vp = vplayout(1, 2), gp = gpar(fontsize = 8))
        print(pie[[1]], vp = vplayout(2, 1))
        print(pie[[2]], vp = vplayout(2, 2))
        print(pie[[3]], vp = vplayout(2, 3))
        print(pie[[6]], vp = vplayout(3, 1))
        print(pie[[4]], vp = vplayout(3, 2))
        print(pie[[5]], vp = vplayout(3, 3))
      }
      
      # gridExtra::grid.arrange(pie$A,pie$B,pie$C,pie$D,pie$E, nrow = 2)
      p0 <- NULL
    } else if(mode=="bar") {
      mp   <- cinfo
      mp   <- unique(mp[,c(structure,"sample","relative_n_cells")])
      refs <- unique(cinfo$sample)
      
      mp <- split(mp, mp[,structure])
      mp <- lapply(mp, function(x) {
        if(length(setdiff(refs, x$sample))>0) {
          sname <- setdiff(refs, x$sample)
          snum  <- length(sname)
          y <- data.frame("sample" = sname
                          , structure = rep(unique(x[,structure]), snum)
                          , "relative_n_cells" = rep(0,snum))
          colnames(y)[colnames(y)=="structure"] <- structure
          x <- rbind.data.frame(x, y)
          rownames(x) <- NULL
        }
        return(x)
      })
      
      mp <- do.call(rbind, mp)
      mp[,structure] <- factor(mp[,structure], levels = levels(mp[,structure])[order(levels(mp[,structure]), decreasing = F)])
      p0 <- ggplot(mp, aes_string(x=structure, y = "relative_n_cells", fill=feature)) +
        geom_col(position = 'dodge')
      # + coord_flip()
      x_lab <- structure
    }
    
    if(mode!="pie") {
      p <- p0 +
        theme_bw() + my_theme +
        scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
        scale_fill_manual(values = pal) + xlab(x_lab) + ylab("Percent of cells")
    } else {
      p <- p0
    }
  } else {
    stop(message("[!] Invalid info_type (absolute/relative)."))
  }
  
  return(p)
}

# Heatmap ----
plotsceExpressionHeatmap <- function(sce, assay.type
                                     , scale=T
                                     , ridx=NULL
                                     , cidx=NULL
                                     , myLegend=NULL
                                     , myPalette=NULL
                                     , myZscale=NULL
                                     , text_size = 8
                                     , annotDF = NULL
                                     , annotCol =  NULL) 
{
  m <- as.matrix(assay(sce,assay.type))
  
  if(!is.null(ridx))  m <- m[ridx,]
  if(!is.null(cidx))  m <- m[,cidx]
  
  if(scale) {
    m_scaled <- t(apply(m, 1, scale))
    colnames(m_scaled) <- colnames(m)
  } else {
    m_scaled <- m
  }
  
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
  hm <- ComplexHeatmap::Heatmap(m_scaled
                                , cluster_rows = F
                                , cluster_columns = F
                                , show_row_names = F
                                , show_column_names = F
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
                                , use_raster = T
                                , raster_device = "tiff"
                                , raster_quality = 10)
  ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom")
  return(hm)
}