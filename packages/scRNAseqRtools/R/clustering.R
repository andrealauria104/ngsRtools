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
          y <- data.frame("sample" = setdiff(refs, x$sample)
                          , structure = x[,structure]
                          , "relative_n_cells" = 0)
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
