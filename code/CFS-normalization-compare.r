library(methylumi) # for milr
library(minfi)     # for SWAN
library(IlluminaHumanMethylation450kmanifest)
library(IMA)       # for IMA

# source milr code and read in various forms of data ----------------------

setwd('code')
source('milr-source.r') # source code
setwd('..')
setwd('data/')
load('milr-data.rda')   # load milrData object

# read in methylumi object for milr normalization
setwd('CFS-data/')
load('CFS-methylumi.rda')    # load methylumi object for milr
load('CFS-RGset.rda')        # load RGset object for minfi
load('CFS-exprmethy450.rda') # load expermethy450 object for IMA

###########################################################################
###############  PERFORM VARIOUS FORMS OF NORMALIZATION  ##################
###########################################################################


# perform normalization using IMA -----------------------------------------

# Normalize using peak correction and quantile normalization
IMAbetas <- IMA.methy450PP(CFSexprmethy450, na.omit = TRUE, normalization=TRUE, transfm =FALSE,
                           peakcorrection = TRUE, samplefilterdetectP = FALSE, sitefilterdetectP = FALSE, 
                           snpfilter = FALSE, locidiff = FALSE)


# perform normalization using minfi & SWAN --------------------------------

# perform SWAN on raw signals
rawMSet <- preprocessSWAN(CFSRGset)

# perform SWAN after using Illumina background correction
bgCorrMSet <- preprocessIllumina(CFSRGset, bg.correct = TRUE,
                                      normalize = "controls", reference = 2)

SWANbgcorrMSet <- preprocessSWAN(CFSRGset, CFSbgCorrMSet)


# perform normalization using milr ----------------------------------------

# milr without housekeeping genes
milrNoHK <- milr(CFSmethylumi, milrData, useHK = FALSE)
milrNoHKbetas <- computeBetas(milrNoHK)

# milr with housekeeping genes winsorizing probes out of prediciton range
milrHKwins <- milr(CFSmethylumi, milrData)
milrHKwinsBetas <- computeBetas(milrHKwins)

# milr with housekeeping genes extrapolating loess surface
# to probes out of prediction range
milrHKextrap <- milr(CFSmethylumi, milrData, edgeMethod = 'gridExtrap' )
milrHKextrapBetas <- computeBetas(milrHKextrap)


###########################################################################
##############  MAKE COMPARISONS BETWEEN DIFFERENT METHODS  ###############
###########################################################################


