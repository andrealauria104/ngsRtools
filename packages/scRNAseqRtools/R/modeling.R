# Modeling ----
model_variance <- function(sce
                           , assay.type  = "logcounts"
                           , method      = "loess"
                           , ntop        = 2500
                           , fdr.cutoff  = 0.05
                           , bio.cutoff  = 0.5
                           , mean.cutoff = 0.1)
{
  
  plot_variance <- function(decomp)
  {
    p <- ggplot(as.data.frame(decomp), aes(x=mean, y=total, col=hv)) +
      geom_point(size=0.8, alpha=0.7) +
      geom_line(aes(y=tech), col='#CC0000') +
      theme_bw() + my_theme + scale_color_manual(values=c('black', '#0066CC')) +
      xlab("Mean - log2[CPM]") + ylab("Variance")
    
    return(p)
  }
  
  # Fit a mean-dependent trend to the gene-specific variances (technical variance)
  varfit <- scran::trendVar(sce
                     , parametric = T
                     , method     = method
                     , assay.type = assay.type
                     , use.spikes = F)
  
  # Decompose the gene-specific variance into biological and technical components
  decomp  <- scran::decomposeVar(sce, varfit)
  decomp  <- decomp[order(decomp$bio, decreasing=TRUE), ]
  
  hv      <- decomp
  
  if(!is.null(mean.cutoff)) {
    hv <- subset(hv, mean >= mean.cutoff)
  }
  if(!is.null(fdr.cutoff)) {
    hv <- subset(hv, FDR <= fdr.cutoff)
  }
  if(!is.null(bio.cutoff)) {
    hv <- subset(hv, bio >= bio.cutoff)
  }
  
  hv      <- hv[order(hv$bio, decreasing=TRUE), ]
  
  if(is.null(ntop) || ntop > nrow(hv)) {
    hvgenes <- rownames(hv)
  } else {
    hvgenes <- rownames(hv)[1:ntop]
  }
  
  decomp$hv <- F
  decomp$hv[match(hvgenes, rownames(decomp))] <- T
  
  p <- plot_variance(decomp)
  
  return(list("model" = decomp, "hvgenes" = hvgenes, "plot" = p))
}

regress_variation_out <- function(expression_matrix, covariates, model_formula)
{
  mod <- model.matrix(as.formula(model_formula)
                      , data = covariates
                      , drop.unused.levels = TRUE)
  fit <- limma::lmFit(expression_matrix, mod)
  beta <- fit$coefficients[, -1, drop = FALSE]
  beta[is.na(beta)] <- 0
  expression_matrix <- as.matrix(expression_matrix) - beta %*% t(mod[, -1])
  
  return(expression_matrix)
}