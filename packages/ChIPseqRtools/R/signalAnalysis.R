# Signal analysis ----
read_matrix <- function(path_matrix, pheno_1, pheno_2, ...)
{
  mat <- read.delim2(path_matrix
                     , header = F
                     , comment.char = "@"
                     , stringsAsFactors = F
                     , sep = "\t"
                     , ...)
  colnames(mat)[1:4] <- c("chr", "start", "end","peak")
  
  if(!missing(pheno_1) & !missing(pheno_2)) {
    idx <- 7:ncol(mat)
    nb  <- length(idx)/2
    mat[,idx] <- lapply(mat[,idx], as.numeric)
    colnames(mat)[idx] <- c(paste0(pheno_1,"_bin_",1:nb), paste0(pheno_2,"_bin_",1:nb))
  }
  
  return(mat)
}

read_matrix_v2 <- function(path_matrix, pheno, ...)
{
  mat <- read.delim2(path_matrix
                     , header = F
                     , comment.char = "@"
                     , stringsAsFactors = F
                     , sep = "\t"
                     , ...)
  colnames(mat)[1:4] <- c("chr", "start", "end","peak")
  
  if(!missing(pheno)) {
    npheno <- length(pheno)
    idx <- 7:ncol(mat)
    nb  <- length(idx)/npheno
    mat[,idx] <- lapply(mat[,idx], as.numeric)
    colnames(mat)[idx] <- unlist(sapply(pheno, function(x) paste0(x,"_bin_",1:nb), simplify = F))
  }
  
  return(mat)
}

analyze_peaks <- function(path_matrix, pheno_1, pheno_2, method = 'wilcoxon-signed-rank', fdr.th = 1e-5, diff.th = NULL, score = 'mean',...)
{
  # internal functions ---
  perform_wx_test <- function(peaks, mat, ipheno_1, ipheno_2, fdr.th = 1e-5, diff.th = NULL)
  {
    # Perform Wilcoxon Signed-Rank test for signal in peaks --
    get_wx_test <- function(i1, i2)
    {
      wx.g <- wilcox.test(x = as.numeric(i1), y=as.numeric(i2), alternative = 'g', paired = T)
      wx.l <- wilcox.test(x = as.numeric(i1), y=as.numeric(i2), alternative = 'l', paired = T)
      
      return(list("g"=wx.g$p.value,"l"=wx.l$p.value))
    }
    
    peaks$wx.g <- apply(mat, 1, function(x) get_wx_test(i1 = x[ipheno_1], i2 = x[ipheno_2])$g)
    peaks$wx.g <- p.adjust(peaks$wx.g, method = "BH")
    peaks$wx.l <- apply(mat, 1, function(x) get_wx_test(i1 = x[ipheno_1], i2 = x[ipheno_2])$l)
    peaks$wx.l <- p.adjust(peaks$wx.l, method = "BH")
    
    peaks$status <- "shared"
    
    if(is.null(diff.th)) {
      ge_idx <- which(peaks$wx.g<=fdr.th)
      le_idx <- which(peaks$wx.l<=fdr.th)
    } else {
      ge_idx <- which(peaks$wx.g<=fdr.th & peaks$score_1-peaks$score_2>=diff.th)
      le_idx <- which(peaks$wx.l<=fdr.th & peaks$score_2-peaks$score_1>=diff.th)
    }
    peaks$status[ge_idx] <- paste0(pheno_1,"-enriched")
    peaks$status[le_idx] <- paste0(pheno_2,"-enriched")
    
    return(peaks)
  }
  
  test_signal_distribution <- function(method, ...)
  {
    if(method=='wilcoxon-signed-rank') {
      return(perform_wx_test(...))
    }
  }
  
  # analysis ---
  message("[*] Analyze ChIP-seq signal in peaks")
  message(" -- Phenotype ",pheno_1, "-",pheno_2)
  
  if(is.character(path_matrix)) {
    message(" -- Read scaled matrix")
    mat <- read_matrix(path_matrix, pheno_1, pheno_2, ...)
  } else if(length(dim(path_matrix)) > 1) {
    message(" -- Loading provided scaled matrix")
    mat <- path_matrix
  }
  
  message(" -- Compute phenotype score")
  ipheno_1 <- grep(pheno_1, colnames(mat))
  ipheno_2 <- grep(pheno_2, colnames(mat))
  peaks <- mat[,1:4]
  
  if(score == 'mean') {
    peaks$score_1 <- rowMeans(mat[,ipheno_1])
    peaks$score_2 <- rowMeans(mat[,ipheno_2])
    
  } else if(score == 'median') {
    peaks$score_1 <- apply(mat[,ipheno_1], 1, median)
    peaks$score_2 <- apply(mat[,ipheno_2], 1, median)
    
  }
  
  message(" -- Test signal distribution, method = ", method)
  peaks <- test_signal_distribution(method = method, peaks, mat, ipheno_1, ipheno_2, fdr.th = fdr.th, diff.th = diff.th, ...)
  
  colnames(peaks)[colnames(peaks)=='score_1'] <- paste0("score_", pheno_1)
  colnames(peaks)[colnames(peaks)=='score_2'] <- paste0("score_", pheno_2)
  
  message("[*] all done.")
  return(peaks)
}

plot_peaks <- function(peaks, pheno_1, pheno_2, pal = NULL)
{
  require(ggplot2)
  #source("theme_setting.R")
  
  if(is.null(pal)) pal <- c("orange","black","grey")
  
  ax.max <- round(max(unlist(peaks[,grep('score', colnames(peaks))]))+0.5,1)
  ax.min <- round(min(unlist(peaks[,grep('score', colnames(peaks))]))-0.5,1)
  
  p0 <- ggplot(peaks, aes_string(x=paste0("score_",pheno_2), y=paste0("score_",pheno_1))) +
    geom_point(aes(col=status), size = 0.8, alpha = 0.1) +
    geom_point(data = subset(peaks, status!="shared"), aes(col=status), size = 0.8, alpha = 0.4)
  
  
  if(any(grepl("regulation", colnames(peaks)))) p0 <- p0 + facet_grid(~regulation)
  
  p <- p0 +
    # geom_hline(yintercept = c(1,0), lwd = 0.25, linetype = 'dashed', col = c('red','blue')) +
    # geom_vline(xintercept = c(1,0), lwd = 0.25, linetype = 'dashed', col = c('red','blue')) +
    xlab(paste0(pheno_2," score")) + ylab(paste0(pheno_1," score")) + theme_bw() + my_theme +
    theme(aspect.ratio = 1, panel.grid.major = element_blank(), legend.title = element_blank()) +
    xlim(c(ax.min,ax.max)) + ylim(c(ax.min,ax.max)) +
    scale_color_manual(values = pal)
  
  return(p)
}

analyze_profile <- function(mat, method = 'mean', distance = 3000)
{
  # internal functions ---
  get_profile <- function(mat, distance, method = 'mean')
  {
    require(reshape2)
    
    profile_func <- function(method, ...)
    {
      if(method=='mean') {
        return(mean(..., na.rm=T))
      } else if(method=='sum') {
        return(sum(..., na.rm = T))
      } else if(method=='median') {
        return(median(..., na.rm=T))
      } else if (is.null(method)) {
        stop(message("[!] Please, provide profile summary method."))
      }
    }
    
    cn  <- colnames(mat)
    bid <- grep("bin", cn)
    
    profile <- apply(mat[,bid],2, profile_func, method = method)
    profile <- melt(profile)
    profile$pheno <- sapply(strsplit(rownames(profile),"_"),"[[",1)
    profile$pheno[grep("KO", profile$pheno)] <- "3BKO"
    profile$bin   <- as.numeric(sapply(strsplit(rownames(profile),"_"),"[[",3))
    
    bsize <- 2*distance/length(unique(profile$bin))
    profile$distance <- NA
    profile$distance[grep("KO", profile$pheno)] <- seq(-distance,distance, bsize)
    profile$distance[grep("WT", profile$pheno)] <- seq(-distance,distance, bsize)
    
    return(profile)
  }
  plot_profile <- function(profile)
  {
    require(ggplot2)
    #source("theme_setting.R")
    
    p <- ggplot(profile, aes(x=distance, y=value, group=pheno)) + geom_line(aes(col=pheno)) + facet_grid(~status) +
      xlab("distance from midpoint \n (base-pairs)") + ylab("log2 [IP/input] ") + theme_bw() + my_theme +
      theme(aspect.ratio = 1, panel.grid.major = element_blank(), legend.title = element_blank()
            # , strip.text = element_text(size=8, face = 'bold')
      ) +
      scale_color_manual(values = c("orange","grey"))
    
    return(p)
  }
  # analyze profile ---
  if(any(grepl("status", colnames(mat)))) {
    mat <- split(mat, mat[,'status'])
    profile <- lapply(mat, function(x) {
      y <- get_profile(mat=x, distance = distance, method = method)
      y$status <- paste0(unique(x$status),"\n n = ", nrow(x))
      return(y)
    })
    profile <- do.call(rbind, profile)
  } else {
    profile <- get_profile(mat = mat, distance = distance, method = method)
    profile$status <- paste0('all peaks', nrow(mat))
  }
  
  rm(mat)
  gc(verbose = F)
  
  p <- plot_profile(profile)
  
  res <- list('profile' = profile, 'plot' = p)
  return(res)
}

analyze_profile_v2 <- function(mat, method = 'mean', distance = 3000, pal = NULL, show_numbers = T)
{
  # internal functions ---
  get_profile <- function(mat, distance, method = 'mean')
  {
    require(reshape2)
    
    profile_func <- function(method, ...)
    {
      if(method=='mean') {
        return(mean(..., na.rm=T))
      } else if(method=='sum') {
        return(sum(..., na.rm = T))
      } else if(method=='median') {
        return(median(..., na.rm=T))
      } else if (is.null(method)) {
        stop(message("[!] Please, provide profile summary method."))
      }
    }
    
    cn  <- colnames(mat)
    bid <- grep("bin", cn)
    
    profile <- apply(mat[,bid],2, profile_func, method = method)
    profile <- melt(profile)
    profile$pheno <- sapply(strsplit(rownames(profile),"_"),"[[",1)
    profile$bin   <- as.numeric(sapply(strsplit(rownames(profile),"_"),"[[",3))
    
    bsize <- 2*distance/length(unique(profile$bin))
    profile$distance <- NA
    for(i in unique(profile$pheno)) {
      profile$distance[grep(i, profile$pheno)] <- seq(-distance,distance, bsize)
    }
    
    return(profile)
  }
  plot_profile <- function(profile, pal)
  {
    require(ggplot2)
    #source("theme_setting.R")
    
    p <- ggplot(profile, aes(x=distance, y=value, group=pheno)) + geom_line(aes(col=pheno)) + facet_grid(~status) +
      xlab("distance from midpoint \n (base-pairs)") + ylab("log2 [IP/input] ") + theme_bw() + my_theme +
      theme(panel.grid.major = element_blank()
            , legend.title = element_blank()
            # , aspect.ratio = 1
            , strip.text = element_text(size=8, face = 'bold')
            , strip.background = element_blank()
      ) +
      scale_color_manual(values = pal)
    
    return(p)
  }
  # analyze profile ---
  if(any(grepl("status", colnames(mat)))) {
    mat <- split(mat, mat[,'status'])
    profile <- lapply(mat, function(x) {
      y <- get_profile(mat=x, distance = distance, method = method)
      if(show_numbers) {
        y$status <- paste0(unique(x$status),"\n n = ", nrow(x))
      } else {
        y$status <- unique(x$status)
      }
      return(y)
    })
    profile <- do.call(rbind, profile)
  } else {
    profile <- get_profile(mat = mat, distance = distance, method = method)
    if(show_numbers) {
      profile$status <- paste0('all peaks [', nrow(mat),']')
    } else {
      profile$status <- 'all peaks'
    }
  }
  
  rm(mat)
  gc(verbose = F)
  
  if(is.null(pal)) pal <- c("grey","darkred","purple")
  p <- plot_profile(profile, pal = pal)
  
  res <- list('profile' = profile, 'plot' = p)
  return(res)
}
