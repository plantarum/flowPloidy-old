setOldClass("nls.lm")

setClass(
  Class = "flowHist4",
  representation = representation(
    raw = "flowFrame", ## raw data, object defined in flowCore
    channel = "character", ## data channel to use for histogram
    bins = "numeric", ## the number of bins to use
    histData = "data.frame", ## binned histogram data
    peaks = "matrix", ## peak coordinates for initial values
    comps = "list", ## TODO complete this
    model = "function", ## model to fit
    nls = "nls.lm" ## nls output
  ),
  prototype = prototype(
    ## TODO complete this
  )
)

setMethod(
  f = "initialize",
  signature = "flowHist4",
  definition = function(.Object, file, channel, bins = 256,
                        window = 20, smooth = 20, pick = FALSE,
                        ... ){
    .Object@raw <- read.FCS(file, dataset = 1, alter.names = TRUE)
    .Object@channel <- channel
    .Object <- setBins(.Object, bins)
    if(pick){
      .Object <- pickPeaks4(.Object)
    } else {
      .Object <- findPeaks4(.Object, window = window,
                                  smooth = smooth)
      .Object <- cleanPeaks4(.Object, window = window)
    }
    .Object <- addComponents4(.Object)
    .Object <- makeModel4(.Object)
    callNextMethod(.Object, 
                   ## window = window, smooth = smooth, pick = pick,
                   ...)
  })

setBins <- function(fh, bins){
  fh@bins = bins

  ## Extract the data channel
  chanDat <- exprs(fh@raw)[, fh@channel]

  ## remove the top bin - this contains clipped values representing all
  ## out-of-range data, not true values
  chanTrim <- chanDat[chanDat < max(chanDat)]

  metaData <- pData(parameters(fh@raw))
  maxBins <- metaData[which(metaData$name == fh@channel), "range"]
  
  ## aggregate bins: combine maxBins into bins via hist
  binAg <- floor(maxBins / bins)

  histBins <- hist(chanTrim, breaks = seq(from = 0, to = 1024, by = binAg),
                   plot = FALSE)

  intensity <- histBins$counts
  x <- 1:length(intensity)
  SCvals <- getSingleCutVals(intensity, x)
  fh@histData <- data.frame(x = x , intensity = intensity, SCvals = SCvals)
  
  ## NOTE!! add code to clear out out-dated model data when the hist
  ## changes.
  fh@peaks = matrix()
  fh@comps = list()
  fh@nls = structure(list(), class = "nls.lm")
  return(fh)
}

getSingleCutValsBase <- function(intensity, xx){
  ## compute the single cut debris model values
  
  ## Do not extend the model below/beyond the data
  ## Modfit appears to cut off the debris slightly above the lowest data,
  ## which gives a better fit. Perhaps set first.channel to 2-4? Need to
  ## test this and determine best fit. Possibly use an extra parameter to
  ## tune this for each data set individually.
  first.channel <- which(intensity > 0)[2]

  res <- 0
  if(xx >= first.channel & xx < length(intensity)){
    channels = (xx + 1):length(intensity)
    for(j in channels){
      res <- res + j^(1/3) * intensity[j] * 2 /
        (pi * j * sqrt(xx/j * (1 - xx/j)))
    }
  }
  res
}

getSingleCutVals <- Vectorize(getSingleCutValsBase, "xx")

findPeaks4 <- function(fh, window, smooth = window / 2){
  ## extract all peaks from data
  ## smoothing removes most of the noisy peaks
  dat <- fh@histData[, "intensity"]

  smDat <- runmean(dat, k = floor(smooth), endrule = "mean")
  localMax <- runmax(smDat, k = window)
  isMax <- localMax == smDat
  maxVals <- dat[isMax]                 # use the raw data for heights 
  res <- cbind(mean = (1:length(dat))[isMax], height = maxVals)
  fh@peaks <- res
  fh
}

cleanPeaks4 <- function(fh, window){
  ## Remove ties and multiple peaks for histogram analysis

  ## Screen out any ties - if two peaks have the same height, and are
  ## within the same 'window', we need to drop one.
  
  ## If a peak has a 'match' at half the size, use the smaller peak (ie.,
  ## take the G1 peak in cases where the G2 peak is higher) 

  ## After the first peak is selected, only consider peaks that are not a
  ## multiple of the size of this peak when selecting the next one.

  peaks <- fh@peaks
  peaks <- peaks[order(peaks[,2], decreasing = TRUE), ]

  ## eliminate the debris field?
  peaks <- peaks[which(peaks[, "mean"] > 40), ]

  drop <- numeric()
  for(i in 2: nrow(peaks)){
    if((peaks[i-1, "height"] == peaks[i, "height"]) &
       (abs(peaks[i-1, "mean"] - peaks[i, "mean"]) <= window)){ 
      ## It's a tie!
      drop <- c(drop, i)
    }
  }

  if(length(drop) > 0){                  # there was at least one tie 
    peaks <- peaks[-drop, ]
  }
  
  out <- matrix(NA, nrow = 0, ncol = 2)

  while(nrow(peaks) > 0){
    ## which peaks are half or double the size of the first peak:
    paircheck <-
      which(((peaks[, "mean"] < 0.53 * peaks[1, "mean"]) &
             (peaks[, "mean"] > 0.47 * peaks[1, "mean"])) |
            ((peaks[, "mean"] < 2.13 * peaks[1, "mean"]) &
             (peaks[, "mean"] > 1.89 * peaks[1, "mean"])))
    ## Add the first peak to that list:
    paircheck <- c(1, paircheck)
    if(length(paircheck) == 1){            # no pairs
      out <- rbind(out, peaks[1, ])
      peaks <- peaks[-1, , drop = FALSE]              # remove peak
    } else if (length(paircheck == 2)) {              # pick the smallest
                                        # of the pair 
      out <- rbind(out,
                   peaks[paircheck[which.min(peaks[paircheck, "mean"])], ])
      peaks <- peaks[-paircheck, , drop = FALSE]      # remove pair
    } else {
      warning("paircheck found more than 2 peaks")
    }

  }

  if(is.vector(peaks))
    out <- rbind(out, peaks)

  rownames(out) <- NULL

  out <- out[1:min(2, nrow(out)), , drop = FALSE]
  if(nrow(out) > 1){
    out <- out[order(out[, "mean"]), ]
  }
  fh@peaks <- out
  fh
}

pickPeaks4 <- function(fh){
  if(class(fh) != "flowHist4")
    stop("fh must be a flowHist object")
  message("plotting data...")
  plotFH4(fh)
  message("select peak A:")
  peakA <- unlist(locator(1))
  points(peakA[1], peakA[2], col = 2, cex = 3)
  message("select peak B:")
  peakB <- unlist(locator(1))
  points(peakB[1], peakB[2], col = 3, cex = 3)
  res <- rbind(peakA, peakB)
  colnames(res) <- c("mean", "height")
  rownames(res) <- NULL
  fh@peaks <- res
  fh
}

plotFH4 <- function(fh, ...){
  ## plots the raw data for a flowHist object
  plot(fh@histData$intensity, type = 'n', main = fh@raw@description$GUID,
       ylab = "Intensity", xlab = fh@channel, ...)
  polygon(x = c(fh@histData$x, max(fh@histData$x) + 1),
          y = c(fh@histData$intensity, 0),
          col = "lightgray", border = NA)
}


