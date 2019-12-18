# Principal Component Analysis ====
plotMethPCA <- function(meth
                    , filterByQuantile = T
                    , experimental_info = NULL
                    , labels   = T
                    , col_by   = NULL
                    , shape_by = NULL
                    , pal      = NULL
                    , scree_plot = F
                    , scree_plot_type = "explained_variance"
                    , point_size = 2
                    , max.pcs = NULL
                    , dim_1 = "PC1"
                    , dim_2 = "PC2"
                    , default.title = "CpG Methylation PCA") {
  
  pca  <- methylKit::PCASamples(meth,
                                screeplot = F,
                                filterByQuantile = filterByQuantile,
                                sd.threshold = 0.5,
                                sd.filter = T,
                                obj.return = T)
  
  pca_summary <- summary(pca)$importance
  
  if((is.null(max.pcs) || max.pcs > ncol(pca$x))) {
    max.pcs <- ncol(pca$x)
  }
  
  pca_plot <- as.data.frame(pca$x[,paste0("PC",1:max.pcs)])
  
  if( !is.null(experimental_info) ) {
    if(is.character(experimental_info)) {
      # for compatibility with previous version
      experimental_info <- reshape2::melt(experimental_info)
      if(!is.null(col_by)) {
        colnames(experimental_info)[1] <- col_by 
      }
    }
    if(is.null(col_by)) {
      col_by <- colnames(experimental_info)[1]
    }
    pca_plot[[col_by]] <- experimental_info[rownames(pca_plot), col_by]
    if(!is.null(shape_by)) {
      pca_plot[[shape_by]] <- experimental_info[rownames(pca_plot), shape_by]
    }
  } else {
    if(is.null(col_by)) col_by <- "samples"
    pca_plot[[col_by]] <- rownames(pca_plot) 
  }
  
  if(is.null(pal)) {
    pal <- ggsci::pal_d3()
    pal <- pal(length(unique(pca_plot[,col_by])))
  }
  pca_plot$repel_col_by <- pca_plot[,col_by]
  
  p <- ggplot(pca_plot, aes_string(x=dim_1, y=dim_2, col=col_by, shape=shape_by)) + geom_point(size=point_size) +
    xlab(paste0(dim_1," (",round(pca_summary[2,grep(paste0(dim_1,"$"), colnames(pca_summary))]*100,1),"%)")) +
    ylab(paste0(dim_2," (",round(pca_summary[2,grep(paste0(dim_2,"$"), colnames(pca_summary))]*100,1),"%)")) + 
    theme_bw() + ggtitle(default.title) + my_theme_2 +
    theme(panel.grid.minor = element_blank()
          , plot.title = element_text(face="bold", hjust = 0.5, size=10)
          , aspect.ratio = 1) +
    scale_color_manual(values=pal)
  
  
  if(labels) {
    p <- p + ggrepel::geom_label_repel(aes(label = rownames(pca_plot), col=repel_col_by),
                                       fontface = 'bold'
                                       # , color = 'black'
                                       , show.legend = F
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
        theme_bw() + my_theme_2 + scale_x_continuous(breaks = spdata$Var2) + 
        geom_hline(yintercept = 0.9, lwd=0.2, linetype="dashed", col = "darkred") +
        xlab("Principal Component") + ylab("Explained Variance") + ggtitle("PCA - Scree Plot") + 
        theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title = element_blank(), panel.grid = element_blank()) +
        scale_y_continuous(breaks = seq(0,1,0.1),labels = scales::percent, limits = c(0,1)) + scale_color_manual(values=c("black","grey"))
    } else if(scree_plot_type=='standard_deviation') {
      spdata      <- reshape2::melt(pca_summary[1,])
      
      spdata$Var2 <- as.numeric(gsub("PC", "",rownames(spdata)))
      spdata      <- spdata[spdata$Var2<=max.pcs,]  
      
      sp <- ggplot(spdata, aes(x=Var2, y=value), col='black') + geom_point(size=1) + geom_line() +
        theme_bw() + my_theme_2 + scale_x_continuous(breaks = spdata$Var2) + 
        geom_hline(yintercept = 0.9, lwd=0.2, linetype="dashed", col = "darkred") +
        xlab("Principal Component") + ylab("Standard Deviation") + ggtitle("PCA - Scree Plot") + 
        theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title = element_blank(), panel.grid = element_blank())
      
    }
    
    return(list("pca" = p, "screeplot" = sp))
  } else {
    return(p) 
  }
}