## Contains the various model components used in the flowPloidy analysis.

#' Single-cut debris model
#'
#' Models debris using the single-cut model described by Bagwell et al.
#'   (1991).
#'
#' The model is:
#' S(x) = SCa ⱼ₌ₓ₊₁∑ⁿ ³√j YⱼPₛ(j, x)
#'
#' x is the histogram channel that we're estimating the debris value for
#' SCa is the amplitude parameter
#' Yⱼ is the histogram intensity for channel j.
#' 
#' where Pₛ(j, x) is the probability of a nuclei from channel j falling
#'   into channel x when cut. That is, for j > x, the probability that
#'   fragmenting a nuclei from channel j with a single cut will produce a
#'   fragment of size x. This probability is calculated as:
#'
#' Pₛ(j, x) = 2 / (πj √(x/j)(1 - x/j)
#'
#' This model involves a recursive calculation, since the fitted value for
#'   channel x depends not just on the intensity for channel x, but also
#'   the intensities at all channels > x. Consequently, this is coded with
#'   an internal loop, and then vectorized to produce a well-behaved
#'   function that we can use with the standard nls toolchain.
#'
#' Only the final vectorized function \code{singleCut} is exported for the
#'   use of users. The underlying code is in the internal function
#'   \code{singleCutBase}
#'
#' @name singleCut
#' 
#' @param SCa a numeric value, the single-cut amplitude parameter
#' @param intensity a numeric vector, the histogram intensity in each channel
#' @param xx an integer vector, the ordered channels corresponding to the
#'   values in `intensity'.
#' @return
#' 
#' @references Bagwell, C. B., Mayo, S. W., Whetstone, S. D., Hitchcox, S.
#'   A., Baker, D. R., Herbert, D. J., Weaver, D. L., Jones, M. A. and
#'   Lovett, E. J. (1991), DNA histogram debris theory and compensation.
#'   Cytometry, 12: 107–118. doi: 10.1002/cyto.990120203 
#' 
#' @author Tyler Smith
singleCutBase <- function(SCa, intensity, xx){
  ## unvectorized version, vectorized below.
  
  ## Modfit appears to cut off the debris slightly above the lowest data,
  ## which gives a better fit. The index here is the offset from the first
  ## channel with data in it to start fitting the singleCut model.
  ## Hard-coded to 2 it will always start at the *second* channel with data
  ## in it. Need to test this and determine best fit. Could be a better
  ## hard-coded value, or use an extra model parameter to dynamically fit
  ## each time.
  first.channel <- which(intensity > 0)[2]

  res <- 0
  if(xx >= first.channel & xx < length(intensity)){
    channels = (xx + 1):length(intensity)
    for(j in channels){
      res <- res + j^(1/3) * intensity[j] * 2 / (pi * j * sqrt(xx/j * (1 - xx/j)))
    }
  }
  return(SCa * res)
}

#' @rdname singleCut
#' @export
singleCut <- Vectorize(singleCutBase, "xx")

#' @rdname singleCut
#' @export
oneSampleSC <-
  function (xx, a1, Ma, Sa, a2, b1, Mb, Sb, b2, SCa, intensity) { 
    ## a1 == highest G1 peak
    a1/(sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 * Sa^2)) +
      ## a2 == G2 peak for a1
      a2/(sqrt(2 * pi) * Sa * 2) * exp(-((xx - Ma * 2)^2)/(2 * (Sa * 2)^2)) +
        ## b1 == second highest peak, assumed to be the other G1 peak
        b1/(sqrt(2 * pi) * Sb) * exp(-((xx - Mb)^2)/(2 * Sb^2)) +
        ## b2 == G2 peak for b1
        b2/(sqrt(2 * pi) * Sb * 2) * exp(-((xx - Mb * 2)^2)/(2 * (Sb * 2)^2)) +
        ## single cut debris curve
        singleCut(SCa, intensity, xx) 
}

#' @rdname singleCut
#' @export
oneSampleSCInit <- function(mCall, LHS, data) {
  ## Not sure we need this fancy stuff, given we have the data already in
  ## hand, and in order:
  
  xy <- sortedXyData(mCall[["xx"]], LHS, data)
  ##xy <- data[, "x", "intensity"]
  peaks <- findPeaks(xy[, "y"], 20)     # WARNING: window size hardcoded 
                                        # here!

  peaks <- peaks[order(peaks[, "height"], decreasing = TRUE),]
  Ma <- peaks[1, "mean"]
  Sa <- Ma / 20                         # assume CV = 0.05
  a1 <- peaks[1, "height"] * Sa / 0.4
  if((peaks[1, "mean"] * 2) > max(xy[1, ]))
    a2 <- 0
  else
    a2 <- xy[peaks[1, "mean"] * 2, "intensity"] * Sa * 2 / 0.4
  Mb <- peaks[2, "mean"]
  Sb <- Mb / 20
  b1 <- peaks[2, "height"] * Sb / 0.4
  if((peaks[2, "mean"] * 2) > max(xy[1, ]))
    b2 <- 0
  else
    b2 <- as.vector(xy[peaks[2, "mean"] * 2, "intensity"] * Sb * 2 / 0.4)

  SCa <- 0.1                            # just a wild guess for now
  value <- c(a1, Ma, Sa, a2, b1, Mb, Sb, b2, SCa)
  names(value) <-
    mCall[c("a1", "Ma", "Sa", "a2", "b1", "Mb", "Sb", "b2", "SCa")] 
  value
}

#' @rdname singleCut
#' @export
SSoneSSC <- selfStart(oneSampleSC, oneSampleSCInit,
                      c("a1", "Ma", "Sa", "a2", "b1", "Mb", "Sb", "b2",
                        "SCa"))  

#########################
## Gaussian components ##
#########################

#' @name gausModel
#' @title Gaussian model components
#'
#' Each component contains a two peaks, representing the G1 and G2
#' peaks. The G2 peak is constrained to have a mean and standard deviation 
#' twice that of the G1 peak values. 
#' 
#' @param a1,a2,b1,b2,c1,c2 numeric values, peak area
#' @param Ma,Mb,Mc numeric values, peak mean
#' @param Sa,Sb,Sc numeric values, peak standard deviation
#' @param xx numeric vector, intensity values for each channel
NULL

#' @rdname gausModel
#' @export
gaussA  <- function(a1, Ma, Sa, a2, xx){
    a1 / (sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 * Sa^2)) +
      a2 / (sqrt(2 * pi) * Sa * 2) * exp(-((xx - Ma * 2)^2)/(2 * (Sa * 2)^2))
}

#' @rdname gausModel
#' @export
gaussB  <- function(b1, Mb, Sb, b2, xx){
  b1 / (sqrt(2 * pi) * Sb) * exp(-((xx - Mb)^2)/(2 * Sb^2)) +
    b2 / (sqrt(2 * pi) * Sb * 2) * exp(-((xx - Mb * 2)^2)/(2 * (Sb * 2)^2))
}

#' @rdname gausModel
#' @export
gaussC  <- function(c1, Mc, Sc, c2, xx){
  c1 / (sqrt(2 * pi) * Sc) * exp(-((xx - Mc)^2)/(2 * Sc^2)) +
    c2 / (sqrt(2 * pi) * Sc * 2) * exp(-((xx - Mc * 2)^2)/(2 * (Sc * 2)^2))
}
