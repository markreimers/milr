# code to create list data object containing info to be passed to milr:
# GC content
# probe IDs
# probe type
# X&Y indicator
# housekeeping indicator
# SNP indicator


# load data ---------------------------------------------------------------

setwd('data')
probe.annotation <- read.csv('infinium 450 annotation.csv', 
                             as.is=T, header=T, quote='',row.names=1)
load('450k-probe-gc-content.rda')
load('house-keeping-ILMNID.rda')

# put together the pieces -------------------------------------------------
milrData <- list()

# gc content
milrData$gc <- gc

# ILMNIDs
milrData$ILMNID <- probe.annotation$ILMNID

# probe type Ind
probeSetInd <- rep(NA,nrow(probe.annotation))
probeType <- c('II', 'Red', 'Grn')
probeSetInd[which(probe.annotation$INFINIUM_DESIGN_TYPE == 'II')]  <- 'II'
probeSetInd[which(probe.annotation$COLOR_CHANNEL == 'Red')] <- 'Red'
probeSetInd[which(probe.annotation$COLOR_CHANNEL == 'Grn')] <- 'Grn'
milrData$probeSetInd <- probeSetInd

# XYind
milrData$XYind <- which(probe.annotation$CHR %in% c('X','Y'))

# housekeeping indicator
milrData$HKind <- match(hkILMNIDs, probe.annotation$ILMNID)

# SNP rs probes
milrData$ILMNIDannotSNPs <- grep('rs', probe.annotation$ILMNID)

save(milrData, file='milr-data.rda')


