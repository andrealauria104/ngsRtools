# Modelling ====
find.variable.genes <- function(m
                                , n=1000
                                , x0 = c(-0.5, 0.5)
                                , ret.plot=T)
{
  
  log2_cv2 <- log2(apply(m, 1, function(r) (sd(r)/mean(r))**2))
  log2_mn  <- log2(apply(m, 1, function(r) mean(r)))
  
  idx      <- names(which(!is.na(log2_cv2)))
  
  log2_cv2 <- log2_cv2[idx]
  log2_mn  <- log2_mn[idx]
  
  noise.model <- function(log2_mn, log2_cv2)
  {
    function(x) sum(abs(log2((2 ** log2_mn) ** x[1] + x[2]) - log2_cv2))
  }
  
  xopt <- optim(x0, noise.model(log2_mn, log2_cv2))
  log2_cv2_fit  <- log2(( 2 ** log2_mn) ** xopt$par[1] + xopt$par[2])
  log2_cv2_diff <- log2_cv2 - log2_cv2_fit
  
  idx.ord <- order(log2_cv2_diff, decreasing = T)
  
  log2_cv2_diff[idx.ord][1:n] <- 'TRUE'
  log2_cv2_diff[idx.ord][(n+1):length(log2_cv2_diff)] <- 'FALSE'
  
  y <- cbind.data.frame(    'log2_cv2'      = log2_cv2[idx.ord]
                            , 'log2_mean'     = log2_mn[idx.ord]
                            , 'log2_cv2_fit'  = log2_cv2_fit[idx.ord]
                            , 'High'          = log2_cv2_diff[idx.ord] )
  
  genes <- rownames(y)[y$High=='TRUE']
  
  if(ret.plot) {
    
    p <- ggplot(y, aes(x=log2_mean, y=log2_cv2, col=High)) + 
      geom_point(size=1) +
      geom_line(aes(y=log2_cv2_fit), col='#CC0000') +
      theme_bw() + my_theme + scale_color_manual(values=c('black', '#0066CC')) +
      xlab("log2(mean)") + ylab("log2(cv^2)")
    
    return(list("genes" = genes,
                "plot" = p))
    
  } else {
    return(genes)
  }
}
