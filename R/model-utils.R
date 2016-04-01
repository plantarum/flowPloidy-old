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
#' This is a utility function for use internally by flowPloidy; it is not
#' exported and won't be visible to users.
#' 
#' @param dat a numeric vector containing the data to be searched
#' @param window an integer, the width of the moving window to use
#' @return Returns a matrix with two columns:
#' \describe{
#' \item{mean}{the index of each potential peak}
#' \item{height}{the height of the peak at that index position}
#' }
#' 
#' @author Tyler Smith
#' @examples
#' \dontrun{
#' set.seed(123)
#' test.dat <-(cumsum(runif(1000, min = -1)))
#' plot(test.dat, type = 'l')
#' test.peaks <- findPeaks(test.dat, window = 20)
#' points(test.peaks, col = 'red', cex = 2)
#' }
#' 
findPeaks <- function(dat, window){
  localMax <- runmax(dat, k = window)
  isMax <- localMax == dat
  maxVals <- dat[isMax]
  cbind(mean = (1:length(dat))[isMax], height = maxVals)
}
