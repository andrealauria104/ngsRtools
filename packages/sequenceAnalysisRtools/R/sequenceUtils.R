# sequence utils ====
sequence_to_list <- function(fasta)
{
  fa <- Biostrings::readDNAStringSet(fasta)
  
  sequence_list <- vector(mode = 'list', length = length(fa))
  names(sequence_list) <- names(fa)
  
  for( i in names(fa)) {
    sequence_list[[i]] <- toString(fa[[i]])
  }
  
  rm(fa)
  gc(verbose = F)
  
  return(sequence_list)
}
# GenomicRanges Tools ====
extend <- function(x, upstream=0, downstream=0)
{
  if (any(GenomicRanges::strand(x) == "*")) warning("'*' ranges were treated as '+'")
  
  on_plus <- GenomicRanges::strand(x) == "+" | GenomicRanges::strand(x) == "*"
  new_start <- GenomicRanges::start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- GenomicRanges::end(x) + ifelse(on_plus, downstream, upstream)
  GenomicRanges::ranges(x) <- IRanges::IRanges(new_start, new_end)
  GenomicRanges::trim(x)
}
