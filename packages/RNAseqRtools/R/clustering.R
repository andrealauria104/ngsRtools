# Clustering ====
# K-means ---
findClustSSE <- function(scaledata, nKM=20, ret=F)
{
  wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
  for (i in 2:nKM) wss[i] <- sum(kmeans(scaledata,
                                        centers=i)$withinss)
  wss.best <- as.numeric(which.max(wss))
  plot(1:nKM, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  abline(v = wss.best, lty = 2)
  
  if(ret) return(wss.best)
}

findClustASW <- function(scaledata, nKM=20, nstart = 100, ret=F)
{
  sil <- rep(0, nKM)
  #repeat k-means for 1:n and extract silhouette:
  for(i in 2:nKM){
    k1to20 <- kmeans(scaledata, centers = i, nstart = nstart, iter.max = 20)
    ss <- cluster::silhouette(k1to20$cluster, dist(scaledata))
    sil[i] <- mean(ss[, 3])
  }
  
  # Plot the  average silhouette width
  plot(1:nKM, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
  abline(v = which.max(sil), lty = 2)
  
  sil.best <- as.numeric(which.max(sil))
  cat("Average silhouette width optimal number of clusters:", sil.best, "\n")
  
  if(ret) return(sil.best)
}

findClustCal <- function(scaledata, nKM=20, ret=F)
{
  # Calinski-Harabasz index
  fit <- vegan::cascadeKM(scaledata, 1, nKM, iter = 100)
  plot(fit, sortg = TRUE, grpmts.plot = TRUE)
  calinski.best <- as.numeric(which.max(fit$results[2,]))
  cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")
  
  if(ret) return(calinski.best)
}

setKMClusters <- function(m, method="silhouette", ...)
{
  scaledata <- t(scale(t(m)))
  
  if( method=="sse" ) {
    return(findClustSSE(scaledata, ...))
  }
  
  if( method=="silhouette") {
    return(findClustASW(scaledata, ...))
  }
  
  if( method=="calinski") {
    return(findClustCal(scaledata, ...))
  } else {
    stop(message('[!] Invalid method, please provide one of sse/silhouette/calinksi'))
  }
}

getKMeans <- function(x, centers, scaledata = T, ...) {
  if(scaledata) x <- t(scale(t(x)))
  stats::kmeans(x = x, centers = centers, ...)
}

# partitioning aroung medoids (PAM) ---
setPAMCnumber <- function(m, metric = "euclidean", scaledata = T, nc = 12, ret = T)
{
  if(scaledata) m <- t(scale(t(m)))
  pams <- lapply(2:nc, function(k) cluster::pam(m, k = k, metric = metric))
  names(pams) <- 2:nc
  sw <- unlist(lapply(pams, function(i) i$silinfo$avg.width))
  oc <- as.numeric(names(pams)[which.max(sw)])
  message(" - PAM optimal number of clusters = ",  oc)
  
  if(ret) return(oc)
}

getPAM <- function(m, k, metric = "euclidean", scaledata = T, ...) {
  if(scaledata) m <- t(scale(t(m)))
  list('cluster' = cluster::pam(m,k, metric = metric, ...))
}

# Hierarchical clustering ---
get_hclust <- function(m
                       , distance = 'euclidean'
                       , method   = 'complete'
                       , cut_tree = 'static'
                       , cut_tree_h = NULL
                       , cut_tree_k = NULL
                       , return_d   = F
                       , return_hr  = T
                       , minClusterSize = 50
                       , ...)
{
  
  message("[*] Hierarchical clustering")
  message(" -- distance: ", distance)
  message(" -- method: ", method)
  
  m  <- t(scale(t(m), center = T, scale = T))
  
  res <- vector(mode = 'list')
  
  if( any(grepl("pearson|spearman|kendall", distance)) ) {
    c <- cor(t(y), method = distance) 
    d <- as.dist(1-c)
  } else {
    d <- dist(m, method = distance)
  }
  
  if(return_d) {
    res$d <- d
  }
  
  hr <- hclust(d, method = method, members=NULL)
  
  if(return_hr) {
    res$hr <- hr
  }
  
  if(cut_tree == 'static') {
    
    if(is.null(cut_tree_h) || is.null(cut_tree_k)) {
      cut_tree_h <- quantile(hr$height, probs = seq(0,1,0.01))['99%']  
    } 
    message(" -- Static tree cutting")
    message(" -- heigth = ", cut_tree_h)
    cl <- cutree(hr, h=cut_tree_h, cut_tree_k)
    
    res$cl <- cl
    rm(cl)
    
  } else if(cut_tree == 'dynamic') {
    
    message(" -- Dynamic tree cutting")
    
    cl <-  dynamicTreeCut::cutreeDynamic(dendro = hr
                         , distM = as.matrix(d)
                         , minClusterSize = minClusterSize
                         , deepSplit = 0)
    res$cl <- cl
    rm(cl)
    
  } else if(!return_d & !return_hr) {
    res$hr <- hr
  } 
  
  return(res)
}