## Functions for creating and viewing FlowHist objects.

#' @importFrom flowCore read.FCS exprs pData parameters
NULL

#' @importFrom knitr kable
NULL

#' @importFrom rmarkdown render
NULL

#' @importFrom graphics hist lines locator plot points polygon grconvertX
#'   grconvertY text abline
NULL

#' @importFrom stats as.formula coef integrate predict pnorm
NULL

#' @importFrom utils write.table
NULL

#' @importFrom methods setClass setMethod new callNextMethod
NULL

setOldClass("nls")

#' FlowHist
#'
#' Creates a \code{FlowHist} object from an FCS file, setting up the
#' histogram data for analysis.
#'
#' Starting with a \code{flowFrame} object, read from a FCS file,
#' \code{FlowHist} will:
#'
#' \enumerate{
#' \item Extract the intensity data from \code{channel}.
#'
#' \item Remove the top bin, which contains off-scale readings we ignore
#' in the analysis.
#'
#' \item aggregates the raw data into the desired number of bins, as
#' specified with the \code{bins} argument. The default is 256, but you may
#' also try 128 or 512. Any integer is technically acceptable, but I
#' wouldn't stray from the default without a good reason.
#'
#' \item identify model components to include. All \code{FlowHist} objects
#' will have the single-cut debris model and the G1 peak for sample A, and
#' the broadened rectangle for the S-phase of sample A. Depending on the
#' data, additional components for the G2 peak and sample B (G1, G2,
#' s-phase) may also be added.
#' 
#' \item Build the NLS model. All the components are combined into a single
#' model. 
#'
#' \item Identify starting values for Gaussian (G1 and G2 peaks) model
#' components. For reasonably clean data, the built-in peak detection is
#' fine. You can evaluate this by plotting the \code{FlowHist} object with
#' the argument \code{init = TRUE}. If it doesn't look good, you can play
#' with the \code{window} and \code{smooth} arguments (which is tedious!),
#' or pick the peaks visually yourself with \code{pick = TRUE}.
#' }
#' 
#' @name FlowHist
#'
#' @param file character, the name of the file to load
#' @param files character, a vector of file names to load
#' @param channel character, the name of the data column to use
#' @param bins integer, the number of bins to use to aggregate events into
#'   a histogram 
#' @param window the width of the moving window used to identify local
#'   maxima for peak detection via \code{caTools:runmax}
#' @param smooth the width of the moving window used to reduce noise in
#'   the histogram via \code{caTools::runmean}
#' @param pick boolean; if TRUE, the user will be prompted to select peaks
#'   to use for starting values. Otherwise (the default), starting values
#'   will be detected automatically.
#' @param verbose boolean; if TRUE, \code{histBatch} will list files as it
#' processes them. 
#' 
#' @slot raw a flowFrame object containing the raw data from the FCS file
#' @slot channel character, the name of the data column to use
#' @slot bins integer, the number of bins to use to aggregate events into a
#'   histogram
#' @slot histdata data.frame, the columns are the histogram bin number
#'   (xx), florescence intensity (intensity), and the raw single-cut debris
#'   model values (SCVals, used in model fitting). Additional columns may
#'   be added if/when I add gating, so refer to columns by name, not
#'   position. 
#' @slot peaks matrix, containing the coordinates used for peaks when
#'   calculcating initial parameter values.
#' @slot comps a list of \code{modelComponent} objects included for these
#'   data.
#' @slot model the function (built from \code{comps}) to fit to these
#' data.
#' @slot init a list of initial parameter estimates to use in fitting the
#'   model.
#' @slot nls the nls object produced by the model fitting
#' @slot counts a list of cells counted in each peak of the fitted model
#' @slot CV a list of the coefficients of variation for each peak in the
#'   fitted model.
#' @slot RCS numeric, the residual chi-square for the fitted model.
#'
#' @return \code{FlowHist} returns a \code{FlowHist} object.
#' @author Tyler Smith
setClass(
  Class = "FlowHist",
  representation = representation(
    raw = "flowFrame", ## raw data, object defined in flowCore
    channel = "character", ## data channel to use for histogram
    bins = "integer", ## the number of bins to use
    linearity = "character", ## "fixed" or "variable", to determine whether
    debris = "character", ## "SC" or "MC", to set the debris model. 
    ## or not linearity is fixed at 2, or allowed to vary as a model
    ## parameter 
    histData = "data.frame", ## binned histogram data
    peaks = "matrix", ## peak coordinates for initial values
    opts = "list",    ## flags for selecting model components
    comps = "list", ## model components
    model = "function", ## model to fit
    init = "list", ## inital parameter estimates
    nls = "nls", ## nls output
    counts = "list", ## cell counts in each peak
    CV = "list", ## CVs
    RCS = "numeric" ## residual chi-square
  ),
  prototype = prototype(
    ## TODO complete this?
  )
)

setMethod(
  f = "initialize",
  signature = "FlowHist",
  definition = function(.Object, file, channel, bins = 256,
                        window = 20, smooth = 20, pick = FALSE,
                        linearity = "variable", debris = "SC",
                        opts = list(), ... ){
    .Object@raw <- read.FCS(file, dataset = 1, alter.names = TRUE)
    .Object@channel <- channel
    .Object <- setBins(.Object, bins)
    if(pick){
      .Object <- pickPeaks(.Object)
    } else {
      .Object <- findPeaks(.Object, window = window,
                           smooth = smooth)
      .Object <- cleanPeaks(.Object, window = window)
    }
    .Object@linearity <- linearity
    .Object@debris = debris
    .Object@opts <- opts
    .Object <- addComponents(.Object)
    .Object <- makeModel(.Object)
    .Object <- getInit(.Object)
    callNextMethod(.Object, ...)
  })

resetFlowHist <- function(fh, from = "peaks"){
  ## Clear analysis slots
  ## Default is to clear everything from peaks onwards
  removeFrom <- c("peaks", "comps")
  ## Dependencies
  ## - changing peaks changes everything
  ## - changing comps changes model, init, nls

  ## coded to allow for further refinement, if/when additions to the
  ## FlowHist class makes it sensible to change the granularity of slot
  ## resetting. 
  
  removeNum <- which(removeFrom == from)

  rmF <- function(x)
    removeNum <= which(removeFrom == x)
  
  if(rmF("peaks"))
    fh@peaks = matrix()
  if(rmF("comps")){
    fh@comps = list()
    fh@model = function(){}
    fh@init = list()
    fh@nls = structure(list(), class = "nls")
    fh@counts = list()
    fh@CV = list()
    fh@RCS = NA_real_
  }
  fh
}


#' @rdname FlowHist
#' @examples
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' fh1
#' @export
FlowHist <- function(file, channel, bins = 256, window = 20, smooth = 20,
                     pick = FALSE, linearity = "variable", debris = "SC",
                     opts = list(), analyze = TRUE){
  fh <-  new("FlowHist", file = file, channel = channel,
             bins = as.integer(bins), window = window, smooth = smooth,
             pick = pick, linearity = linearity, debris = debris,
             opts = opts)
  if(analyze)
    fh <- fhAnalyze(fh)
  return(fh)
}

#' Displays the column names present in an FCS file
#'
#' A convenience function for viewing column names in an FCS, in order to
#'   select one for the \code{channel} argument in \code{FlowHist}.
#' 
#' @title viewFlowChannels
#' @param file character, the name of a FCS data file
#' @return A vector of column names from the FCS file.
#' @seealso \code{FlowHist}
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData) 
#' viewFlowChannels(flowPloidyFiles[1])
#' @export
viewFlowChannels <- function(file){
  tmp <- read.FCS(file, alter.names = TRUE, dataset = 1)
  cnames <- colnames(exprs(tmp))
  names(cnames) <- NULL
  cnames
}


#' @rdname FlowHist
#' @examples
#' library(flowPloidyData) 
#' batch1 <- batchFlowHist(flowPloidyFiles, channel = "FL3.INT.LIN")
#' batch1
#' @return \code{batchFlowHist} returns a list of \code{FlowHist} objects.
#' @export
batchFlowHist <- function(files, channel, bins = 256, verbose = TRUE,
                      window = 20, smooth = 20, linearity = "variable",
                      debris = "SC"){ 
  res <- list()
  for(i in seq_along(files)){
    if(verbose) message("processing ", files[i])
    tmpRes <- FlowHist(file = files[i], channel = channel, bins = bins,
                       window = window, smooth = smooth, pick = FALSE,
                       linearity = linearity, debris = debris)
    res[[getFHFile(tmpRes)]] <- tmpRes
  }              
  return(res)
}

getFHFile <- function(fh){
  fh@raw@description$GUID
}

setMethod(
  f = "show",
  signature = "FlowHist",
  def = function(object){
    cat("FlowHist object '")
    cat(getFHFile(object)); cat("'\n")
    cat("channel: "); cat(object@channel); cat("\n")
    cat("bins: "); cat(object@bins); cat("\n")
    cat("linearity: "); cat(object@linearity); cat("\n")
    cat("debris: "); cat(object@debris); cat("\n")
    cat(length(object@comps)); cat(" model components: ")
    cat(paste(names(object@comps), collapse = ", ")); cat("\n")
    pnames <- names(formals(object@model))
    specialP <- names(getSpecialParams(object))
    pnames <- pnames[which(! pnames %in% specialP)]
    cat(length(pnames)); cat(" parameters: ");
    cat(paste(pnames, collapse = ", ")); cat("\n")
    cat(length(specialP)); cat(" special parameters: ");
    cat(paste(specialP, collapse = ", ")); cat("\n")    
    if(length(object@nls) == 0)
      cat("Model fitting not complete\n")
    else
      cat("Model fit\n")

    if(length(object@counts) > 0){
      cat(paste("\nAnalysis\n========\n"))
      ## cat(paste("Modelled events: ",
      ##           round(object$counts$total$value, 1)))
      counts <- c(object@counts$firstPeak$value,
                  object@counts$secondPeak$value)
      size <- c(coef(object@nls)["Ma"],  coef(object@nls)["Mb"])
      if(is.na(size[2])) size <- size[1]
    }
  
    if(length(object@CV) > 0){
    cvs <- c(object@CV$CVa, object@CV$CVb)
    if(!is.null(object@CV$CVb)){
      cat(paste("\nRatio Peak A / Peak B: ", round(object@CV$CI[1], 3),
                ", SE: ", round(object@CV$CI[2], 5), sep = ""))
    }
  }

  if(length(object@counts) > 0 & length(object@CV) > 0){
    if(length(counts) == 2)
      rnames <- c("Peak A", "Peak B")
    else if (length(counts) == 1)
      rnames <- "Peak A"
    print(kable(data.frame(counts = counts, size = size, cvs = cvs,
                           row.names = rnames), format = "markdown",
                digits = 3))
  }
  
  if(!is.null(object@RCS)){
    cat(paste("\nRCS:", round(object@RCS, 3), "\n"))
  }

  }
)

########################
## Plotting functions ##
########################
#' Plot the raw data for a FlowHist object
#'
#' Creates a simple plot of the raw histogram data. Used as a utility for
#' other plotting functions, and perhaps useful for users who wish to
#' create their own plotting routines.
#' 
#' @param fh a \code{FlowHist} object
#' @param ... additional parameters passed to \code{plot}
#' @return Not applicable, used for plotting
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' plotFH(fh1)
#' @export
plotFH <- function(fh, ...){
  ## plots the raw data for a FlowHist object
  plot(fh@histData$intensity, type = 'n', main = getFHFile(fh),
       ylab = "Intensity", xlab = fh@channel, ...)
  polygon(x = c(fh@histData$xx, max(fh@histData$xx) + 1),
          y = c(fh@histData$intensity, 0),
          col = "lightgray", border = NA)
}

#' Plot histograms for FlowHist objects
#'
#' Plot histograms for FlowHist objects
#'
#' @param x a \code{FlowHist} object
#' @param init boolean; if TRUE, plot the regression model using the
#'   initial parameter estimates over the raw data. 
#' @param nls boolean; if TRUE, plot the fitted regression model over the
#'   raw data (i.e., using the final parameter values)
#' @param comps boolean; if TRUE, plot the individual model components
#'   over the raw data.
#' @param ... additional arguments passed on to plot()
#' @return Not applicable
#' @author Tyler Smith
#' @export
plot.FlowHist <- function(x, init = FALSE, nls = TRUE, comps = TRUE, ...){
  plotFH(x, ...)

  if(init){
    yy <- with(x@histData,
               do.call(x@model, args = c(getSpecialParams(x), x@init)))
    lines(x = x@histData$xx,
          y = yy, 
          col = "grey", lwd = 1, lty = 5)
    abline(v = x@init$Ma, col = "blue", lwd = 2)
    points(x = x@init$Ma, y  = x@histData$intensity[round(x@init$Ma, 0)],
           cex = 2, pch = 16, col = "blue")
    text(paste("Peak A: ", round(x@init$Ma, 0)), cex = 1,
         x = x@init$Ma, col = "blue", pos = 2,
         y = grconvertY(0.9, from = "npc", to = "user"))
    abline(v = 2 * x@init$Ma, col = "blue", lwd = 0.5)
    if(! is.null(x@init$Mb)){
      abline(v = x@init$Mb, col = "orange", lwd = 2)
      points(x = x@init$Mb, y  = x@histData$intensity[round(x@init$Mb, 0)],
             cex = 2, pch = 16, col = "orange")
      text(paste("Peak B: ", round(x@init$Mb, 0)), cex = 1,
           x = x@init$Mb, col = "orange", pos = 2,
           y = grconvertY(0.7, from = "npc", to = "user"))
      abline(v = 2 * x@init$Mb, col = "orange", lwd = 0.5)
    }
  }

  if(nls & (length(x@nls) > 0)){
    dat <- tabulateFlowHist(x)
    lines(x = x@histData$xx, y = predict(x@nls), col = 2)
    text(paste("RCS: ", round(dat$rcs, 3)), cex = 1, pos = 4,
         x = grconvertX(0.8, from = "npc", to = "user"),
         y = grconvertY(0.95, from = "npc", to = "user"))
    text(paste("Peak A: ", round(dat$sizeA, 1), "/",
               round(dat$countsA, 1), "/",
               round(dat$cvA, 3)),
         cex = 1, pos = 4, col = "blue",
         x = grconvertX(0.8, from = "npc", to = "user"),
         y = grconvertY(0.9, from = "npc", to = "user"))
    text(paste("Peak B: ", round(dat$sizeB, 1), "/",
               round(dat$countsB, 1), "/",
               round(dat$cvB, 3)),
         cex = 1, pos = 4, col = "orange",
         x = grconvertX(0.8, from = "npc", to = "user"),
         y = grconvertY(0.85, from = "npc", to = "user"))
    text(paste("A/B:  ", round(dat$ratioAB, 3)), cex = 1, pos = 4,
         x = grconvertX(0.8, from = "npc", to = "user"),
         y = grconvertY(0.8, from = "npc", to = "user"))
    text(paste("Linearity:  ", round(dat$linearity, 3)), cex = 1, pos = 4,
         x = grconvertX(0.8, from = "npc", to = "user"),
         y = grconvertY(0.75, from = "npc", to = "user"))

  }

  if(comps & (length(x@nls) > 0)){
    yy <- with(x@histData,
               do.call(x@model, args = c(getSpecialParams(x), x@init)))
  
  
    for(i in seq_along(x@comps)){
      params <-
        as.list(coef(x@nls)[names(formals(x@comps[[i]]@func))])
      params <- params[! is.na(names(params))]
      yy <-
        with(x@histData,
             do.call(x@comps[[i]]@func,
                      args = c(getSpecialParamsComp(x@comps[[i]]),
                               params)))
        lines(x = x@histData$xx, y = yy, col = x@comps[[i]]@color) 
    }
  }
}

####################
## Exporting Data ##
####################

#' Extract analysis results from a flowHist object
#'
#' A convenience function for extracting the results of the NLS
#'   curve-fitting analysis on a flowHist object.
#'
#' If \code{fh} is a single FlowHist object, a data.frame with a single
#' row is returned. If \code{fh} is a list of \code{FlowHist} objects, a
#' row for each object will be added to the data.frame.
#'
#' If a file name is provided, the data will be saved to that file.
#' 
#' @title exportFlowHist
#' @param fh a flowHist object, or a list of flowHist objects.
#' @param file character, the name of the file to save data to
#' @return a data frame 
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' fh1 <- fhAnalyze(fh1)
#' tabulateFlowHist(fh1)
#' @export
tabulateFlowHist <- function(fh, file = NULL){
  if(class(fh) == "FlowHist")
    res <- exFlowHist(fh)
  else if (class(fh) == "list" && all(sapply(fh, class) == "FlowHist")){
    res <- do.call(rbind, lapply(fh, exFlowHist))
  }
  row.names(res) <- res$file
  res$file <- NULL
  
  if(! is.null(file))
    write.table(x = res, file = file)

  res
}

exFlowHist <- function(fh){
  if(!is.null(fh@nls)){
    if(fh@linearity == "variable")
      linearity = coef(fh@nls)["d"]
    else
      linearity = NA
             
    data.frame(file = getFHFile(fh), channel = fh@channel,
             components = paste(names(fh@comps), collapse = ";"),
             totalEvents = sum(fh@histData$intensity),
             countsA = fh@counts$firstPeak$value,
             countsB = ifelse(is.null(fh@counts$secondPeak$value), NA,
                              fh@counts$secondPeak$value),
             sizeA = coef(fh@nls)["Ma"],
             sizeB = coef(fh@nls)["Mb"],
             cvA = fh@CV$CVa,
             cvB = ifelse(is.null(fh@CV$CVb), NA, fh@CV$CVb),
             ratioAB = unlist(ifelse(is.null(fh@CV$CI[1]), NA,
                                     fh@CV$CI[1])),
             ratioSE = unlist(ifelse(is.null(fh@CV$CI[2]), NA,
                                     fh@CV$CI[2])),
             rcs = fh@RCS,
             linearity = linearity,
             row.names = NULL)
  } else {
    data.frame(file = getFHFile(fh), channel = fh@channel,
               components = paste(names(fh@comps), collapse = ";"),
               totalEvents = sum(fh@data$intensity),
               countsA = NA,
               countsB = NA,
               sizeA = NA,
               sizeB = NA,
               cvA = NA,
               cvB = NA,
               ratioAB = NA,
               ratioSE = NA,
               rcs = NA,
               linearity = NA,
               row.names = NULL)
  }
}

#################################################
## Functions for initializing FlowHist objects ##
#################################################

#' (Re-) set the bins for a FlowHist object
#'
#' This function sets (or resets) the number of bins to use in aggregating
#' FCS data into a histogram, and generates the corresponding data matrix.
#'
#' The \code{histData} matrix also contains the \code{SCvals} column. This
#' is used to calculate the single-cut debris component in the NLS model.
#' 
#' @title setBins
#' @param fh a \code{FlowHist} object
#' @param bins integer, the number of bins to use in aggregating FCS data
#' @return a \code{FlowHist} object, with the \code{bins} slot set to
#'   \code{bins}, and the corresonding binned data stored in a matrix in
#'   the \code{histData} slot. Any previous analysis slots are removed:
#'   \code{peaks, comps, model, init, nls, counts, CV, RCS}.
#' @author Tyler Smith
#' @examples
#' ## defaults to 256 bins:
#' library(flowPloidyData) 
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' plot(fh1)
#' ## reset them to 512 bins:
#' fh1 <- setBins(fh1, 512)
#' plot(fh1)
#' @export
setBins <- function(fh, bins = 256){
  fh@bins = as.integer(bins)

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
  xx <- 1:length(intensity)
  startMax <- max(intensity[which(intensity != 0)][1:20])
  startBin <- which(intensity == startMax)[1]
  SCvals <- getSingleCutVals(intensity, xx, startBin)
  MCvals <- getMultipleCutVals(intensity, startBin)
  fh@histData <- data.frame(xx = xx, intensity = intensity,
                            SCvals = SCvals, MCvals = MCvals)
  fh <- resetFlowHist(fh)
  fh
}

#' @importFrom caTools runmean runmax
NULL

##############################
## Peak Detection/Selection ##
##############################

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
#' won't be visible to users. Usually invoked from within \code{FlowHist}.
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
#' @param fh a \code{FlowHist} object
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
  dat <- fh@histData[, "intensity"]

  smDat <- runmean(dat, k = floor(smooth), endrule = "mean")
  localMax <- runmax(smDat, k = window)
  isMax <- localMax == smDat
  maxVals <- dat[isMax]                 # use the raw data for heights 
  res <- cbind(mean = (1:length(dat))[isMax], height = maxVals)
  fh@peaks <- res
  fh
}

#' @rdname findPeaks
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
cleanPeaks <- function(fh, window){
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
#' @examples
#' library(flowPloidyData) 
#' fh2 <- FlowHist(file = flowPloidyFiles[12], channel = "FL3.INT.LIN")
#' plot(fh2, init = TRUE) ## automatic peak estimates
#' \dontrun{
#' fh2 <- pickInit(fh2)   ## hand-pick peak estimates
#' }
#' plot(fh2, init = TRUE) ## revised starting values
#' @export
pickInit <- function(fh){
  fh@peaks = matrix()
  fh@comps = list()
  fh@model = function(){}
  fh@init = list()
  fh@nls = structure(list(), class = "nls")
  fh@counts = list()
  fh@CV = list()
  fh@RCS = NA_real_

  fh <- pickPeaks(fh)
  fh <- addComponents(fh)
  fh <- makeModel(fh)
  fh <- getInit(fh)
  fh
}

pickPeaks <- function(fh){
  ## Does the work of actually plotting and selecting peaks for
  ##   \code{pickInit}
  if(class(fh) != "FlowHist")
    stop("fh must be a FlowHist object")
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
  fh@peaks <- res
  fh
}

##########################
## Change Model Options ##
##########################
updateFlowHist <- function(fh, linearity = NULL, debris = NULL,
                           analyze = TRUE){
  ## keep the existing peaks, as they may have already been tweaked by the
  ## user
  message("updating FlowHist")
  if(!is.null(linearity))
    if(linearity %in% c("fixed", "variable"))
      fh@linearity <- linearity
    else
      stop("Invalid linearity value")
  if(!is.null(debris))
    if(debris %in% c("SC", "MC"))
      fh@debris <- debris
    else
      stop("Invalid debris value")
  
  fh <- resetFlowHist(fh, from = "comps")
  fh <- addComponents(fh)
  fh <- makeModel(fh)
  fh <- getInit(fh)
  if(analyze)
    fh <- fhAnalyze(fh)
  fh
}
