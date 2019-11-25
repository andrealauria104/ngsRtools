# Normalization ====
countToTpm <- function(counts, effLen)
{
  effLen <- effLen[rownames(counts)]
  # Transcripts Per Million (TPM)
  tpms <- apply(counts, 2,
                function(x) {
                  rate <- log(x) - log(effLen)
                  denom <- log(sum(exp(rate)))
                  c <- exp(rate - denom + log(1e6))
                  return(c)
                }
  )
  return(tpms)
}

countToFpkm <- function(counts, effLen, pcg)
{
  # Fragments Per Kilobases per Million (FPKM) 
  effLen <- effLen[rownames(counts)]
  pcg    <- intersect(pcg, rownames(counts))
  
  fpkms <- apply(counts, 2,
                 function(x) {
                   N <- sum(x[pcg])
                   c <- exp( log(x) + log(1e9) - log(effLen) - log(N) )
                   return(c)
                 }
  )
  return(fpkms)
}
