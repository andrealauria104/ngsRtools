# Distance between strings ----
compute_hamming_distance <- function(s) 
{
  # input s must be a character vector
  # containg strings to compare
  slen <- nchar(s[1])
  n <- vector(mode = 'numeric', length = slen)
  for(i in 1:slen) {
    n[i] <- substr(s[1],i,i)!=substr(s[2],i,i)
  }
  return(sum(n))
}

get_all_substr <- function(x, n, sub.method) 
{
  # find all substring of length n
  # for string x
  if(sub.method == 'all') {
    ss <- unique(substring(x, 1:(nchar(x) - n + 1), n:nchar(x)))  
  } else if(sub.method == 'trim') {
    ss <- substring(x, 1,n)
  }
  return(ss)
}

calculate_string_distance <- function(s, sub.method = "all", split.by = "\\-", verbose = F)
{
  calc_distance <- function(s, sub.method, split.by, verbose) {
    
    if( is.character(s) & length(s)==1 ) {
      s <- unlist(strsplit(s, split.by))  
    }
    
    ln <- sapply(s, nchar)
    
    if(length(unique(ln))==1) {
      if(verbose) message("[+] Computing Hamming distance")
      d <- compute_hamming_distance(s)
    } else {
      if(verbose) message("[+] Unequal string lenghts, computing Hamming distance on substrings")
      gidx <- which.max(ln)
      lidx <- which.min(ln)
      
      gs <- get_all_substr(x = s[gidx], n = ln[lidx], sub.method = sub.method)
      d  <- c(gs, s[lidx])
      d  <- combn(x = d, 2, FUN=compute_hamming_distance)
      d  <- d[which.min(d)]
    }
    return(d)
  }
  
  if(is.list(s)) s <- unique(unlist(s))
  
  if(length(s) <= 2) {
    index_distance <- calc_distance(s = s, sub.method = sub.method, split.by = split.by, verbose = verbose)
  } else {
    index_comb <- combn(x = index_comb,2, FUN=paste, collapse = "-")
    index_distance <- lapply(index_comb, calculate_string_distance, sub.method = sub.method, split.by = split.by, verbose = verbose)
    names(index_distance) <- index_comb
  }
  return(index_distance)
}