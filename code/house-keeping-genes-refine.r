# refine list of housekeeping genes

library(methylumi)

# load annotation data --------------------------------------------------------
setwd('data')
probe.annotation <- read.csv('infinium 450 annotation.csv', 
							                as.is=T, header=T, quote='',row.names=1)

# load list of housekeeping gene refseq IDs
hkGenes <- read.csv('house-keeping-refseq-ids.csv')[,1]

# load CFS methylumi data object ----------------------------------------------
setwd('CFS-data')
load('CFS-methylumi.rda') # load lumiSet object for our method
CFSbetas <- betas(CFSmethylumi)
rm(CFSmethylumi)
setwd('..')
gc()

# load BS methylumi data -------------------------------------------
setwd('BrainSpan-data')
load('BrainSpan-methylumi.rda') # load BS methylumi object
badChips <- c(7,59,74)
BSbetas <- betas(BSmethylumi)[,-badChips]
rm(BSmethylumi)
setwd('..')
gc()


# map house keeping refseq ids to probes ----------------------------------

# get mutual set for both data sets with no NAs
BSna <- unique(which(is.na(BSbetas), arr.ind=TRUE)[,1])
CFSna <- unique(which(is.na(CFSbetas), arr.ind=TRUE)[,1])
NAind <- union(BSna,CFSna)

RefGeneGroup <- sapply(strsplit(probe.annotation$UCSC_REFGENE_GROUP,';'), function(x) x[1])[-NAind]
# indicate which are near TSS?
tssInd <- grep('TSS', RefGeneGroup)
# which are in gene body?
bodyInd <- grep('Body', RefGeneGroup)
# which are probes are housekeeping?
hkLoci <- c(unlist(lapply(hkGenes, function(x) grep(paste('^', x, '$', sep=''), probe.annotation$UCSC_REFGENE_ACCESSION[-NAind]))))
# get set of housekeeping controls
hkControls <- intersect( c(tssInd, bodyInd), hkLoci)
hkBody <- intersect(bodyInd, hkLoci)
hkTSS <- intersect(tssInd, hkLoci)

# look at housekeeping genes for BrainSpan data -----------------------------

BShkBodyRowMeans <- rowMeans(BSbetas[-NAind,][hkBody,])
BShkTSSRowMeans <-  rowMeans(BSbetas[-NAind,][hkTSS,])

par(mfcol=c(2,3), mar=c(3,3,3,3))
image(BSbetas[-NAind,][hkBody,][order(BShkBodyRowMeans),],main='BS gene body',
      axes=FALSE,zlim=c(0,1))
image(BSbetas[-NAind,][hkTSS,][order(BShkTSSRowMeans),],main='BS tss',
      axes=FALSE,zlim=c(0,1))

# which housekeeping sites are differentially methylated in BS data?
cutoff=.05
BSbodySD <- apply(BSbetas[-NAind,][hkBody,],1,sd)
BStssSD <- apply(BSbetas[-NAind,][hkTSS,],1,sd)
plot(density(BSbodySD),main='BS body sdevs')
abline(v=cutoff)
plot(density(BStssSD),main='BS tss sdevs')
abline(v=cutoff)

# remove and replot 
BSbodyKeep = which(BSbodySD[order(BShkBodyRowMeans)]<=cutoff)
BStssKeep = which(BStssSD[order(BShkTSSRowMeans)]<=cutoff)

BSBodyThrowOut = names(which(BSbodySD[order(BShkBodyRowMeans)]>cutoff))
BStssThrowOut =  names(which(BStssSD[order(BShkTSSRowMeans)]>cutoff))

image(BSbetas[-NAind,][hkBody,][order(BShkBodyRowMeans),][BSbodyKeep,],
      main='BS gene body less high var CpGs', axes=FALSE,zlim=c(0,1))

image(BSbetas[-NAind,][hkTSS,][order(BShkTSSRowMeans),][BStssKeep,],
      main='BS tss less high var CpGs', axes=FALSE,zlim=c(0,1))

# look at housekeeping genes for CFS data -----------------------------------------

par(mfcol=c(2,3), mar=c(3,3,3,3))
CFShkBodyRowMeans <- rowMeans(CFSbetas[-NAind,][hkBody,])
CFShkTSSRowMeans <-  rowMeans(CFSbetas[-NAind,][hkTSS,])

image(CFSbetas[-NAind,][hkBody,][order(CFShkBodyRowMeans),],main='CFS gene body',
      axes=FALSE,zlim=c(0,1))
image(CFSbetas[-NAind,][hkTSS,][order(CFShkTSSRowMeans),],main='CFS tss',
      axes=FALSE,zlim=c(0,1))

# which housekeeping sites are differentially methylated in BS data?
cutoff=.05
CFSbodySD <- apply(CFSbetas[-NAind,][hkBody,],1,sd)
CFStssSD <- apply(CFSbetas[-NAind,][hkTSS,],1,sd)
plot(density(CFSbodySD),main='CFS body sdevs')
abline(v=cutoff)
plot(density(CFStssSD),main='CFS tss sdevs')
abline(v=cutoff)

# remove and replot 
CFSbodyKeep = which(CFSbodySD[order(CFShkBodyRowMeans)]<=cutoff)
CFStssKeep = which(CFStssSD[order(CFShkTSSRowMeans)]<=cutoff)

CFSBodyThrowOut = names(which(CFSbodySD[order(CFShkBodyRowMeans)]>cutoff))
CFStssThrowOut =  names(which(CFStssSD[order(CFShkTSSRowMeans)]>cutoff))

image(CFSbetas[-NAind,][hkBody,][order(CFShkBodyRowMeans),][CFSbodyKeep,],main='CFS gene body',
      axes=FALSE,zlim=c(0,1))
image(CFSbetas[-NAind,][hkTSS,][order(CFShkTSSRowMeans),][CFStssKeep,],main='CFS tss',
      axes=FALSE,zlim=c(0,1))

save(CFSbodyKeep, CFStssKeep, CFSBodyThrowOut, CFStssThrowOut,
     BSbodyKeep, BStssKeep, BSBodyThrowOut, BStssThrowOut,
     file = 'edited-HK-lists.rda')

hkILMNIDs <- c(intersect(names(BSbodyKeep),names(CFSbodyKeep)), 
               intersect(names(BStssKeep), names(CFStssKeep)))

save(hkILMNIDs, file='house-keeping-ILMNID.rda')






