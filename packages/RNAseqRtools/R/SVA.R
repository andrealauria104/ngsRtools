# Surrogate Variable Analysis ====
cleaningP <- function(y, mod, svaobj, P=ncol(mod))
{
  # adapted from Jaffe et al. (BMC Bioinformatics 2015):
  # Practical impacts of genomic data “cleaning” on 
  # biological discovery using surrogate variable analysis
  
  X      <- cbind(mod,svaobj$sv)
  Hat    <- solve(t(X)%*%X)%*%t(X)
  beta   <- (Hat%*%t(y))
  
  cleany <- y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  
  return(cleany)
}
estimate_surrogate_vars <- function(y
                                    , mod = NULL # Custom model
                                    , clean = F
                                    , norm.var = "CPM") 
{
  if(class(y)!="DGEList") {
    stop(message("[!] Please, provide DGEList object."))
  }
  
  # build model --
  if(is.null(mod)) {
    group <- relevel(y$samples$group, ref = y$ref)
    mod   <- model.matrix(~group)
    colnames(mod)[-1] <- paste0(levels(group)[-1], "vs", levels(group)[1])
  }
  
  # estimate surrogate variables --
  svseq <- sva::svaseq(y[[norm.var]], mod, mod[,1])
  colnames(svseq$sv) <- paste0("sv", 1:ncol(svseq$sv)) 
  mod.sv <- cbind(mod, svseq$sv)
  
  res <- list("mod.sv" = mod.sv, "svseq" = svseq)
  
  if(clean) {
    # regress out surrogate variables from data --
    y$cidx <- cleaningP(y = y[[norm.var]], mod = mod, svaobj = svseq)
    names(y)[grep('cidx', names(y))] <- paste0('clean', norm.var) 
    
    res <- c(list("y" = y), res)
  }
  
  return(res)
}
