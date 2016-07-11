## Contains the various model components used in the flowPloidy analysis.

## Functions for building non-linear models for application to flowHist
## objects.

#' Gaussian model components
#'
#' Components for modeling Gaussian features in flow histograms
#'
#' Typically the complete models will contain fA1 and fB2, which model the
#' G1 peaks of the sample and the standard. In most cases, they will also
#' contain fA2 and fB2, which model the G2 peaks. The G2 peaks are linked
#' to the G1 peaks, in that they require some of the parameters from the
#' G1 peaks as well (mean and standard deviation).
#'
#' @param a1,a2,b1,b2 area parameters
#' @param Ma,Mb curve mean parameter
#' @param Sa,Sb curve standard deviation parameter
#' @param xx vector of histogram intensities
#' @return NA
#' @author Tyler Smith
#' @name gauss
fA1 <- function(a1, Ma, Sa, xx){
  (a1 / (sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 * Sa^2)))
}

#' @rdname gauss
fA2 <- function(a2, Sa, Ma, xx){
  (a2 / (sqrt(2 * pi) * Sa * 2) * exp(-((xx - Ma * 2)^2)/(2 * (Sa * 2)^2)))
}

#' @rdname gauss
fB1 <- function(b1, Mb, Sb, xx){
  (b1 / (sqrt(2 * pi) * Sb) * exp(-((xx - Mb)^2)/(2 * Sb^2)))
}

#' @rdname gauss
fB2 <- function(b2, Sb, Mb, xx){
  (b2 / (sqrt(2 * pi) * Sb * 2) * exp(-((xx - Mb * 2)^2)/(2 * (Sb * 2)^2)))
}

## Single-cut debris model
##
## S(x) = a \sum{j = x + 1}^{n} \sqrt[3]{j} Y_j P_s(j, x)
## P_s (j, x) = \frac{2}{(\pi j \sqrt{(x/j) (1 - x/j)}
##
## a = amplitude parameter
## Y_j = intensity in channel j
## P_s(j, x) = probability of a nuclei from channel j falling into channel x
## when cut.
##
## By this formula, the debris intensity in a channel/bin is a function of
## the intensity in all the subsequent bins. This recursive relationship is
## tricky to code; in order to take full advantage of all of the R tools
## that support nls, the mean function needs to return one fitted value for
## one predictor value. The following implementation of singleCut therefore
## takes the entire vector of the response vector (intensity), necessary to
## calculate the debris curve, and returns only the value for a single
## predictor value.


#' Single-cut debris model
#'
#' Models debris using the single-cut model described by Bagwell et al.
#' (1991).
#'
#' The model is:
#' \deqn{S(x) = a \sum{j = x + 1}^{n} \sqrt[3]{j} Y_j P_s(j, x)}
#'
#' x is the histogram channel that we're estimating the debris value for
#' SCa is the amplitude parameter
#' Y_j is the histogram intensity for channel j.
#'
#' where P_s(j, x) is the probability of a nuclei from channel j falling
#' into channel x when cut. That is, for j > x, the probability that
#' fragmenting a nuclei from channel j with a single cut will produce a
#' fragment of size x. This probability is calculated as:
#'
#' \deqn{P_s (j, x) = \frac{2}{(\pi j \sqrt{(x/j) (1 - x/j)}}}
#'
#' This model involves a recursive calculation, since the fitted value for
#' channel x depends not just on the intensity for channel x, but also the
#' intensities at all channels > x. Consequently, this is coded with an
#' internal loop, and then vectorized to produce a well-behaved function
#' that we can use with the standard nls toolchain.
#'
#' @name singleCut
#'
#' @param SCa a numeric value, the single-cut amplitude parameter
#' @param intensity a numeric vector, the histogram intensity in each channel
#' @param xx an integer vector, the ordered channels corresponding to the
#'   values in `intensity'.
#' @return NA
#'
#' @references Bagwell, C. B., Mayo, S. W., Whetstone, S. D., Hitchcox, S.
#'   A., Baker, D. R., Herbert, D. J., Weaver, D. L., Jones, M. A. and
#'   Lovett, E. J. (1991), DNA histogram debris theory and compensation.
#'   Cytometry, 12: 107-118. doi: 10.1002/cyto.990120203
#'
#' @author Tyler Smith

#' @rdname singleCut
singleCutBase <- function(SCa, intensity, xx){
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
      res <- res + j^(1/3) * intensity[j] * 2 / (pi * j * sqrt(xx/j * (1 - xx/j)))
    }
  }
  return(SCa * res)
}

## Need to Vectorize this so our comparisons in line 2 make sense. Without
## this, the comparisons try to test all values of the xx vector against
## first.channel and length(intensity). This raises a warning, and the
## results are not what we want. We need them tested one at a time, hence
## the vectorize here:
singleCutVect <- Vectorize(singleCutBase, "xx")

#' @rdname singleCut
singleCut <- function(SCa, SCvals){
  ## I have no idea why the 'flowPloidy:::' prefix is needed here. As far
  ## as I can tell, the package will load and function just fine without
  ## it. However, the function install_bitbucket(..., build_vignettes =
  ## true) fails without it. It's not doing any harm otherwise, except to
  ## my sense of decency, so I leave it here.
  SCa * SCvals
}

## singleCut <- function(SCa, intensity, xx){
##   ## I have no idea why the 'flowPloidy:::' prefix is needed here. As far
##   ## as I can tell, the package will load and function just fine without
##   ## it. However, the function install_bitbucket(..., build_vignettes =
##   ## true) fails without it. It's not doing any harm otherwise, except to
##   ## my sense of decency, so I leave it here.
##   flowPloidy:::singleCutVect(SCa, intensity, xx)
## }


#' Provide starting values for flowHist NLS models
#'
#' Given a flowHist object with peaks identified (as should be the case
#' after calls to \code{flowHist} or \code{pickInit} (and also the
#' internal functions \code{pickPeaks} or \code{findPeaks},
#' \code{flowInit} will provide rough guesses for starting parameter
#' values.
#'
#' In most cases, at least with reasonably clean histograms, these guesses
#' are sufficient for the NLS optimization routine to find appropriate
#' values.
#'
#' @param fh a \code{flowHist} object
#' @return Returns the \code{flowHist} object with the initial parameter
#'   estimates in the init slot.
#' @author Tyler Smith
#' @export
flowInit <- function(fh) {
  modChange <- FALSE
  xy <- fh$data
  if(class(fh$peaks) != "matrix")
    stop("no peaks identified for flowHist object -- see pickPeaks or findPeaks")

  peaks <- fh$peaks

  params <- names(formals(fh$model))
  params <- params[-which(params %in% c("", "xx", "intensity"))]
  value <- c()

  if("Ma" %in% params) {
    ## Any model with Ma will require all three of these parameters:
    Ma <- peaks[1, "mean"]
    Sa <- Ma / 20                         # assume CV = 0.05
    a1 <- peaks[1, "height"] * Sa / 0.4
    tmpval <- c(Ma, Sa, a1)
    names(tmpval) <- c("Ma", "Sa", "a1")
    value <- c(value, tmpval)
  }

  if("a2" %in% params) {
    ## Is a2 off the chart? It shouldn't be! Models with an a2 peak can
    ## break if the a2 peak is beyond the data range.
    if((peaks[1, "mean"] * 2) > max(xy[ ,"x"])){
      warning("a2 peak appears to be out of range, it has been removed from the model")
      params <- params[params != "a2"]
      fh$comps$fA2 <- NULL
      modChange <- TRUE
    } else {
      a2 <- xy[peaks[1, "mean"] * 2, "intensity"] * Sa * 2 / 0.4
      tmpval <- c(a2)
      names(tmpval) <- c("a2")
      value <- c(value, tmpval)
    }
  } else if(! "a2" %in% params) {
    if((peaks[1, "mean"] * 2) < max(xy[ ,"x"])){
      warning("a2 peak appears to be IN range, adding it to the model")
      fh$comps$fA2 <- fA2
      modChange <- TRUE
      a2 <- xy[peaks[1, "mean"] * 2, "intensity"] * Sa * 2 / 0.4
      tmpval <- c(a2)
      names(tmpval) <- c("a2")
      value <- c(value, tmpval)
    }
  }

  if("Mb" %in% params){
    Mb <- peaks[2, "mean"]
    Sb <- Mb / 20
    b1 <- peaks[2, "height"] * Sb / 0.4
    tmpval <- c(Mb, Sb, b1)
    names(tmpval) <- c("Mb", "Sb", "b1")
    value <- c(value, tmpval)
  }

  if("b2" %in% params) {
    if((peaks[2, "mean"] * 2) > max(xy[,"x"])){
      warning("b2 peak appears to be out of range, it has been removed from the model")
      params <- params[params != "b2"]
      fh$comps$fB2 <- NULL
      modChange <- TRUE
    } else {
      b2 <- as.vector(xy[peaks[2, "mean"] * 2, "intensity"] * Sb * 2 / 0.4)
      tmpval <- c(b2)
      names(tmpval) <- c("b2")
      value <- c(value, tmpval)
    }
  } else if(! "b2" %in% params) {
    if(nrow(peaks) > 1 && (peaks[2, "mean"] * 2) < max(xy[ ,"x"])){
      warning("b2 peak appears to be IN range, adding it to the model")
      fh$comps$fB2 <- fB2
      modChange <- TRUE
      b2 <- xy[peaks[1, "mean"] * 2, "intensity"] * Sb * 2 / 0.4
      tmpval <- c(b2)
      names(tmpval) <- c("b2")
      value <- c(value, tmpval)
    }
  }

  if("SCa" %in% params){
    SCa <- 0.1                            # just a wild guess for now
    tmpval <- c(SCa)
    names(tmpval) <- c("SCa")
    value <- c(value, tmpval)
  }
  fh$init <- as.list(value)
  if(modChange){
    fh$model <- makeModel(fh$comps)
  }

  return(fh)
}

#' Build an NLS model from a list of components
#'
#' Build an NLS model from a list of components
#'
#' @param components a list of model components to combine
#' @param env the parent frame. Not intended for use by users.
#' @return a function for use in the R nonlinear regression routine.
#' @author Tyler Smith
makeModel <- function(components, env = parent.frame()){

  names(components) <- NULL
  args <- unlist(lapply(components, formals))
  args <- args[unique(names(args))]

  bodList <- lapply(components, FUN = body)
  bod <- bodList[[1]]
  bodList <- bodList[-1]

  while(length(bodList) > 0){
    bod <- call("+", bod, bodList[[1]])
    bodList <- bodList[-1]
  }

  eval(call("function", as.pairlist(args), bod), env)

}
