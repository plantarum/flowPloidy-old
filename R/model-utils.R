## Utilities used in constructing and fitting models

#' findPeaks
#'
#' Locate potential peaks in histogram data
#'
#' Peaks are defined as local maxima in the vector of values, using a
#' moving window. Note that these are used in the context of finding
#' starting values - accuracy isn't important, we just need something
#' `close-enough' that the nls algorithm will be able to find the correct
#' value.
#'
#' Utility functions for use internally by flowPloidy; not exported and
#' won't be visible to users. Usually invoked from within \code{flowHist}.
#'
#' Note that there is a trade-off between accuracy in detected peaks, and
#' avoiding noise. Increasing the value of \code{smooth} will reduce the
#' amount of 'noise' that is included in the peak list. However, increasing
#' smoothing shifts the location of the actual peaks. Most of the time the
#' default values provide an acceptable compromise, given we only need to
#' get 'close enough' for the NLS optimization to find the true parameter
#' values. If you'd like to explore this, the internal (unexported)
#' function \code{fhPeakPlot} may be useful.
#' 
#' @param fh a \code{flowHist} object
#' @param window an integer, the width of the moving window to use in
#'   identifying local maxima via \code{caTools::runmax}
#' @param smooth an integer, the width of the moving window to use in
#'   removing noise via \code{caTools::runmean}
#' 
#' @return Returns a matrix with two columns:
#' \describe{
#' \item{mean}{the index position of each potential peak}
#' \item{height}{the height (intensity) of the peak at that index position}
#' }
#' 
#' @author Tyler Smith
#'
#' @seealso \code{\link{fhPeakPlot}}
#' 
#' @examples
#' \dontrun{
#' set.seed(123)
#' test.dat <-(cumsum(runif(1000, min = -1)))
#' plot(test.dat, type = 'l')
#' test.peaks <- flowPloidy::findPeaks(test.dat, window = 20)
#' points(test.peaks, col = 'red', cex = 2)
#' }
#'
#' @name findPeaks
findPeaks <- function(fh, window, smooth = window / 2){
  ## extract all peaks from data
  ## smoothing removes most of the noisy peaks
  dat <- fh$data[, "intensity"]

  smDat <- runmean(dat, k = floor(smooth), endrule = "mean")
  localMax <- runmax(smDat, k = window)
  isMax <- localMax == smDat
  maxVals <- dat[isMax]                 # use the raw data for heights 
  res <- cbind(mean = (1:length(dat))[isMax], height = maxVals)
  res
}

#' @rdname findPeaks
#'
#' @param peaks a matrix of peaks, as returned by \code{findPeaks}
#'
#' @details \code{cleanPeaks} filters the output of \code{findPeaks} to:  
#' \itemize{
#'
#' \item remove duplicates, ie., peaks with the same intensity that occur
#' within \code{window} positions of each other. Otherwise,
#' \code{findPeaks} will consider noisy peaks without a single highest
#' point to be multiple distinct peaks.
#'
#' \item drop G2 peaks. In some cases the G2 peak for one sample will have
#' greater intensity than the G1 peak for another sample. We correct for
#' this by removing detected peaks with means close to twice that of other
#' peaks.
#'
#' \item ignore noise, by removing peaks with \code{intensity} < 40. A
#' somewhat arbitrary value. It's tricky to deal with this issue when the
#' debris field is large.
#' }
#' 
cleanPeaks <- function(peaks, window){
  ## Remove ties and multiple peaks for histogram analysis

  ## Screen out any ties - if two peaks have the same height, and are
  ## within the same 'window', we need to drop one.
  
  ## If a peak has a 'match' at half the size, use the smaller peak (ie.,
  ## take the G1 peak in cases where the G2 peak is higher) 

  ## After the first peak is selected, only consider peaks that are not a
  ## multiple of the size of this peak when selecting the next one.

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
    } else if (length(paircheck == 2)) {              # pick the smallest of the pair
      out <- rbind(out, peaks[paircheck[which.min(peaks[paircheck, "mean"])], ])
      peaks <- peaks[-paircheck, , drop = FALSE]      # remove pair
    } else {
      warning("paircheck found more than 2 peaks")
    }

  }

  if(is.vector(peaks))
    out <- rbind(out, peaks)

  rownames(out) <- NULL

  out <- out[1:2, ]
  out <- out[order(out[, "mean"]), ]
  out
}

#' @title Interactively select model starting values
#'
#' @description Prompts the user to select the peaks to use as initial
#'   values for non-linear regression on a plot of the histogram data. 
#'
#' @details The raw histogram data are plotted, and the user is prompted to
#'   select the peak positions to use as starting values in the NLS
#'   procedure. This is useful when the automated peak-finding algorithm
#'   fails to discriminate between overlapping peaks, or is confused by
#'   noise.
#'
#' The normal use, \code{pickPeaks} is called from \code{pickInit}, rather
#'   than directly by the user.
#'
#' @param fh A \code{flowHist} object
#' 
#' @return \code{pickInit} returns the \code{flowHist} object with its
#'   initial value slot updated.
#'
#' \code{pickPeaks} returns a matrix with each peak as a row, with the mean
#' (position) in the first column, and the height (intensity) in the second
#' column.
#'
#' @author Tyler Smith
#'
#' @seealso \code{\link{flowInit}}
#'
#' @export
pickInit <- function(fh){
  fh$peaks <- pickPeaks(fh)
  fh$init <- flowInit(fh)
  fh
}

#' @describeIn pickInit Does the work of acutally plotting and selecting peaks for
#'   \code{pickInit}
pickPeaks <- function(fh){
  if(class(fh) != "flowHist")
    stop("fh must be a flowHist object")
  message("plotting data...")
  plotFH(fh)
  message("select peak A:")
  peakA <- unlist(locator(1))
  points(peakA[1], peakA[2], col = 2, cex = 3)
  message("select peak B:")
  peakB <- unlist(locator(1))
  points(peakB[1], peakB[2], col = 3, cex = 3)
  res <- rbind(peakA, peakB)
  colnames(res) <- c("mean", "height")
  rownames(res) <- NULL
  res
}

##' Plot detected peaks along with the raw and smoothed data.
##'
##' For examining/debugging the peak-finding algorithm. Allows users to see
##' the raw data, the smoothed data, and the peaks that are detected based
##' on the smoothed parameters.
##'
##' A convenience for examining how different values of \code{window} and
##' \code{smooth} work in peak detection. Not exported for regular use by
##' users.
##' 
##' @title fhPeakPlot
##' @inheritParams findPeaks
##' @return None, used for plotting
##' @author Tyler Smith
##' @seealso \code{findPeaks}
##' 
fhPeakPlot <- function(fh, window, smooth = window/2){
  dat <- fh$data[, "intensity"]
  smDat <- runmean(dat, k = floor(smooth), endrule = "mean")
  localMax <- runmax(smDat, k = window)
  isMax <- localMax == smDat
  maxVals <- dat[isMax]                 # use the raw data for heights 
  res <- cbind(mean = (1:length(dat))[isMax], height = maxVals)

  clean <- cleanPeaks(res, window)
  plot(fh$data$intensity, type = 'l')
  points(smDat, type = 'l', col = 'grey')
  points(res, cex = 2, col = 2)
  points(clean, cex = 2, col = 3)
}

