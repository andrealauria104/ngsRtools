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
