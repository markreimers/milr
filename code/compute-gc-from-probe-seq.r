##  script to compute gc content probe 50-mer probes sequences
##  to be passed to milr function

# read in annotation
setwd('data')
probe.annotation <- read.csv( "infinium 450 annotation.csv", 
                              as.is=T, header=T, quote='',row.names=1)

## compute gc content for probes (this should take a while)
library('Biostrings')
sequences <- BStringSet(probe.annotation$ALLELEA_PROBESEQ)
computeGC <- function(seq){
  sum( letterFrequency( seq, c('C','G') ) )/50 # sum G & C and divide by 
}                                             # 50 mer probe length
gc <- sapply(sequences, computeGC)
save(gc, file = '450k-probe-gc-content.rda')