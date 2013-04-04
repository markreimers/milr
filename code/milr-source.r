milr <- function(lumiSet, milrData, useHK = TRUE,# data to be passed
                 sParm = .15, traceHat = 'approx', nn = 300, nslice = 5, #fit parameters
                 badChips = NULL, edgeMethod = 'winsorize',  # some options
                 autosomal = FALSE, SNPremove = TRUE){  ### other options 

  require('methylumi')
  # create object for methylated and unmethylated channels ------------------
  methTmp <- methylated(lumiSet)
  if(!is.null(badChips)){
    dimTmp <- dim(methTmp)
    dimTmp[2] <- dimTmp[2] - length(badChips)
    signals <- array(dim = c(dimTmp,2))
    signals[,,1] <- unmethylated(lumiSet)[,-badChips]
    signals[,,2] <- methylated(lumiSet)[,-badChips]    
    rm(methTmp, dimTmp)
  }else{
    signals = array(dim = c(dim(methTmp),2))
    signals[,,1] <- unmethylated(lumiSet)
    signals[,,2] <- methylated(lumiSet)
    rm(methTmp)
  }
  
  # quality control & probe selcetion ------------------------------------------
  
  probeType <- unique(milrData$probeSetInd)
  # which probes have zero signal intensity values?
  leaveOut <- unique(c(unlist(apply(signals, 3, function(x) which(rowSums(x==0)>0) ))))
  cat('Removing probes with zero values for signal intensity\n')
  
  # remove SNP probes?
  if(SNPremove){
    SNPind <- milrData$ILMNIDannotSNPs
    leaveOut <- union(leaveOut,SNPind)
    cat('Removing the 50 annotated SNP probes\n')
  }
  
  # only look at autosomes?
  if(autosomal){
    leaveOut <- union(leaveOut, milrData$XYind)
    cat('Removing X & Y chromosomes\n')
  }
  
  # remove these probes from array & annotation & gc content
  log2SignalsMLO <- log2(signals[-leaveOut,,]); rm(signals)
  GC <- milrData$gc[-leaveOut]
  probeIDs <- milrData$ILMNID[-leaveOut]
  
  #get bead type and color channel annotation so we can conduct normalizations separately
  probeSetInd <- milrData$probeSetInd[-leaveOut]
  
  # get set of hk controls --------------------------------------------------
  if(useHK){
    hkControls <- as.numeric(na.omit(match(milrData$HKind, (1:485577)[-leaveOut])))
    cat('Will use',length(hkControls), 'HK controls to fit loess curve \n')
  }
  
  # median center -----------------------------------------------------------
  cat('Median centering by color channel & probe type\n')
  
  log2Centered <- array(dim=dim(log2SignalsMLO))
  
  medCenter <- function(x, whichSet, whichControls, kk){
    sweep(x[whichSet,,kk], 2, tmp<-apply(x[whichControls,,kk],2,median)) + median(x[whichControls,,kk])
  }
  
  for(pp in 1:3){  # median center by probe type
    whichSet <- which(probeSetInd==probeType[pp])
    
    if(useHK){ 
      whichControls <- intersect(whichSet, hkControls)
    }else{
      whichControls <- whichSet
    }
    
    for(kk in 1:2){      
      log2Centered[whichSet,,kk] <- medCenter(log2SignalsMLO, whichSet, whichControls, kk)
    }
  }
  
  # compute robust expt average ---------------------------------------------
  cat('Computing robust experiment-wise average\n')
  log2Standard <- apply(log2Centered, c(1,3), mean, trim=.2)
  
  # compute deviations from average -----------------------------------------
  cat('Computing deviations from average\n')
  log2Deviations <- array(dim=dim(log2SignalsMLO))
  for(kk in 1:2) log2Deviations[,,kk] <- log2Centered[,,kk] - log2Standard[,kk]
  
  # winsorize by probe type if edgeMethod = winsorize -----------------------
  if(edgeMethod=='winsorize' & useHK){
    cat('Winsorizing probes out of prediction range\n')
    winsorizeBySubset <- function(x, whichSet, whichControls){
      x[whichSet][which(x[whichSet] > max(x[whichControls]))] <- max(x[whichControls])
      x[whichSet][which(x[whichSet] < min(x[whichControls]))] <- min(x[whichControls])
      x[whichSet]
    }
    
    for(pp in 1:3){
      whichSet <- which(probeSetInd==probeType[pp])
      whichControls <- intersect(whichSet, hkControls)
      # winsorize gc
      GC[whichSet] <- winsorizeBySubset(GC, whichSet, whichControls)
      # winsorize mean intensities
      for(kk in 1:2)
        log2Standard[whichSet,kk] <- winsorizeBySubset(log2Standard[,kk], whichSet, whichControls)  
    }
  }
  
  # winsorize GC for extrap method ------------------------------------------
  if(edgeMethod == 'gridExtrap' & useHK){
    cat('Winsorizing GC content\n')
    # winsorize GC
    winsorizeBySubset <- function(x, whichSet, whichControls){
      x[whichSet][which(x[whichSet] > max(x[whichControls]))] <- max(x[whichControls])
      x[whichSet][which(x[whichSet] < min(x[whichControls]))] <- min(x[whichControls])
      x[whichSet]
    }
    
    for(pp in 1:3){
      whichSet <- which(probeSetInd==probeType[pp])
      whichControls <- intersect(whichSet, hkControls)
      GC[whichSet] <- winsorizeBySubset(GC, whichSet, whichControls)
    }
  }
  
  # fit & subtract out loess ------------------------------------------------
  log2NormedDevs <- array(dim=dim(log2SignalsMLO))
  indepVars <- data.frame(GC=GC, UMavg=log2Standard[,1], Mavg=log2Standard[,2])
  
  cat('Fitting & subtracting out loess\n')
  for(pp in 1:3){
    whichSet <- which(probeSetInd==probeType[pp])
    if(useHK){ 
      whichControls <- intersect(whichSet, hkControls)
    }else{
      whichControls <- whichSet
    }
    
    for(jj in 1:dim(log2SignalsMLO)[2]){
      for(kk in 1:2){
        y <- log2Deviations[,jj,kk]
        modelDat <- as.data.frame(cbind(y, indepVars))
        tempFit <- loess(y ~ GC*Mavg*UMavg, trace.hat=traceHat,
                         span=sParm, modelDat, subset=whichControls)
        log2NormedDevs[whichSet,jj,kk] <- log2Deviations[whichSet,jj,kk] - predict(tempFit, modelDat[whichSet,])
      }
    }
  }
  
  # extrapolate by grid if method=gridExtrap --------------------------------
  NAind <- unique(which(is.na(log2NormedDevs[,,1]), arr.ind=T)[,1])
  if(length(NAind)==0) NAind <- NULL
  if(edgeMethod == 'gridExtrap' & useHK){
    cat('Extrapolating surface to missing values across ',nslice,'x',nslice,' grid \n' , sep='')
    
    modelLM <- y~GC+Mavg+UMavg+ GC*Mavg + GC*UMavg + Mavg*UMavg + GC*Mavg*UMavg
    indepVarsLM <- as.matrix(data.frame(GC=GC, UMavg=log2Standard[,1], Mavg=log2Standard[,2],
                                      GCxMavg=GC*log2Standard[,2], GCxUMavg=GC*log2Standard[,1],
                                      MavgxUMavg=log2Standard[,1]*log2Standard[,2],
                                      GCxMavgxUMavg = GC*log2Standard[,1]*log2Standard[,2],
                                      intercept=rep(1,nrow(log2Standard))))
    
    xRanges <- apply(indepVars, 2, range)[,2:3]
    gridMarks <- apply(xRanges, 2, function(x) seq(x[1],x[2],length.out=(nslice+1)))
    gridMarks[1,] <- rep(0, 2)
    for(ii in 1:nslice){
      for(jj in 1:nslice){
        inSlice <- Reduce(intersect, list(
          v1 = which((indepVars[,2]>gridMarks[ii,1])&(indepVars[,2]<=gridMarks[(ii+1),1])),
          v2 = which((indepVars[,3]>gridMarks[jj,2])&(indepVars[,3]<=gridMarks[(jj+1),2])))
        )
        inSliceNA <- intersect(inSlice, NAind)
        cat(ii,jj,'\n')
                  
        if(length(inSliceNA>0)){
          for(pp in 1:3){
            
            whichSet <- which(probeSetInd==probeType[pp])
            whichSetInSliceNA <- intersect(whichSet, inSliceNA)
            
            sliceCenter <- apply(indepVars[whichSetInSliceNA,2:3],2,mean)
            sliceCenterDist <- rowSums(sweep(sweep(indepVars[,2:3],2,sliceCenter), 2, (xRanges[2,]-xRanges[1,]), '/')^2)
            
            
            if(length(whichSetInSliceNA)>0){
              radius <- sqrt(2)/(4*nslice) + .00001
              repeat{
                insideCircleInd <- setdiff(intersect(whichSet, which(sliceCenterDist < radius^2)),whichSetInSliceNA)
                if(length(insideCircleInd)>=nn & 
                   sum(sliceCenterDist[whichSetInSliceNA]>radius)==0) break
                radius <- radius + .0005
              }
                  
              for(yy in 1:dim(log2SignalsMLO)[2]){
                for(zz in 1:2){
                  
                  y <- log2Deviations[,jj,kk]
                  betaHat <- lm.fit(indepVarsLM[insideCircleInd,], y[insideCircleInd])$coef
                  if(length(whichSetInSliceNA)>1) preds <- betaHat %*% t(indepVarsLM[whichSetInSliceNA,])
                  if(length(whichSetInSliceNA)==1) preds <-  betaHat %*% indepVarsLM[whichSetInSliceNA,]
                  log2NormedDevs[whichSetInSliceNA,yy,zz] <- log2Deviations[whichSetInSliceNA,yy,zz] - preds
                  
                }# kk in 1:2
              }# jj in dim
            }# length(whichSetNA) > 0
          }# pp in 1:3
        }# length(inSliceNA > 0)
      }
    }
  }
    
  # compute normalized log2 signals
  log2NormedSignals <- array(dim=dim(log2SignalsMLO))  
  for(kk in 1:2) log2NormedSignals[,,kk] <- log2NormedDevs[,,kk] + log2Standard[,kk]

  # put important variables into a list to output
  milrOutput = list()
  milrOutput$log2NormedSignals <- log2NormedSignals
  milrOutput$log2Deviations <- log2Deviations
  milrOutput$log2Standard <- log2Standard
  milrOutput$probeSetInd <- probeSetInd
  milrOutput$leaveOut <- leaveOut
  milrOutput$ILMNID <- milrData$ILMNID
  milrOutput$chipNames <- colnames(methylated(lumiSet))
  if(useHK){
    milrOutput$hkControls <- hkControls
  }else{
    milrOutput$hkControls <- 'Normalized with No Controls'
  }
  
  if(!is.null(badChips)){
    milrOutput$badChips <- badChips
  }else{
    milrOutput$badChips <- NULL
  }
  
  if(!is.null(NAind)){
    milrOutput$NAind <- NAind
  }else{
    milrOutput$NAind <- NULL
  }
  
  milrOutput
  
}


# small function to compute betas from normedSignals object ---------------

computeBetas = function(milrOutput, normed=TRUE, betaOffset=100){
  
  if(normed){
    x <- milrOutput$log2NormedSignals
  }else{
    x <- array(dim=dim(milrOutput$log2NormedSignals))
    x[,,1] <- milrOutput$log2Deviations[,,1] + milrOutput$log2Standard[,1]
    x[,,2] <- milrOutput$log2Deviations[,,2] + milrOutput$log2Standard[,2]
  }
  
  betas <- (2^x[,,2]) / (2^x[,,1] + 2^x[,,2] + betaOffset)
  rownames(betas) <- milrOutput$ILMNID[-milrOutput$leaveOut]
  
  if(!is.null(milrOutput$badChips)){
    colnames(betas) <- milrOutput$chipNames[-badChips]
  }else{
    colnames(betas) <- milrOutput$chipNames
  }
  
  betas
  
}

