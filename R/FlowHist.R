## Functions for creating and viewing FlowHist objects.

#' @importFrom flowCore read.FCS exprs pData parameters
NULL

#' @importFrom knitr kable
NULL

#' @importFrom rmarkdown render
NULL

#' @importFrom graphics hist lines locator plot points polygon grconvertX grconvertY text abline
NULL

#' @importFrom stats as.formula coef integrate predict pnorm
NULL

#' @importFrom utils write.table
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
#'
#' @slot raw a flowFrame object containing the raw data from the FCS file
#' @slot channel character, the name of the data column to use
#' @slot bins integer, the number of bins to use to aggregate events into a
#'   histogram 
#' @slot histdata data.frame, the columns are the histogram bin number (x),
#'   florescence intensity (intensity), and the raw single-cut debris model
#'   values (SCVals, used in model fitting). Additional columns may be
#'   added if/when I add gating, so refer to columns by name, not position.
#' @slot peaks matrix, containing the coordinates used for peaks when
#'   calculcating initial parameter values.
#' @slot comps a list of \code{modelComponent} objects included for these
#'   data.
#' @slot model the function (built from \code{comps}) to fit to these data.
#' @slot init a list of initial parameter estimates to use in fitting the
#'   model.
#' @slot nls the nls object produced by the model fitting
#' @slot counts a list of cells counted in each peak of the fitted model
#' @slot CV a list of the coefficients of variation for each peak in the
#'   fitted model.
#' @slot RCS numeric, the residual chi-square for the fitted model.
#'
#' @author Tyler Smith
setClass(
  Class = "FlowHist",
  representation = representation(
    raw = "flowFrame", ## raw data, object defined in flowCore
    channel = "character", ## data channel to use for histogram
    bins = "integer", ## the number of bins to use
    histData = "data.frame", ## binned histogram data
    peaks = "matrix", ## peak coordinates for initial values
    comps = "list", ## model components
    model = "function", ## model to fit
    init = "list", ## inital parameter estimates
    nls = "nls", ## nls output
    counts = "list", ## cell counts in each peak
    CV = "list", ## CVs
    RCS = "numeric" ## residual chi-square
  ),
  prototype = prototype(
    ## TODO complete this
  )
)

#' @rdname FlowHist
#' @export
FlowHist <- function(file, channel, bins = 256, window = 20, smooth = 20,
                     pick = FALSE){
  new("FlowHist", file = file, channel = channel, bins = as.integer(bins),
      window = window, smooth = smooth, pick = pick)
}

#' @rdname FlowHist 
#' @export
batchFlowHist <- function(files, channel, bins = 256, verbose = TRUE,
                      window = 20, smooth = 20){ 
  res <- list()
  for(i in seq_along(files)){
    if(verbose) message("processing ", files[i])
    tmpRes <- FlowHist(file = files[i], channel = channel, bins = bins,
                       window = window, smooth = smooth, pick = FALSE)
    res[[getFHFile(tmpRes)]] <- tmpRes
    res[[getFHFile(tmpRes)]] <- fhAnalyze(res[[getFHFile(tmpRes)]])
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
    cat(length(object@comps)); cat(" model components: ")
    cat(paste(names(object@comps), collapse = ", ")); cat("\n")
    pnames <- names(formals(object@model))
    pnames <- pnames[which(! pnames %in% c("xx", "SCvals"))]
    cat(length(pnames)); cat(" parameters: ");
    cat(paste(pnames, collapse = ", "))
    cat("\n")
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

setMethod(
  f = "initialize",
  signature = "FlowHist",
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
    .Object <- getInit4(.Object)
    callNextMethod(.Object, ...)
  })


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
#' @export
plotFH <- function(fh, ...){
  ## plots the raw data for a FlowHist object
  plot(fh@histData$intensity, type = 'n', main = getFHFile(fh),
       ylab = "Intensity", xlab = fh@channel, ...)
  polygon(x = c(fh@histData$x, max(fh@histData$x) + 1),
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
    yy <- do.call(x@model,
                  args = c(list(SCvals = x@histData$SCvals,
                                xx = x@histData$x),
                           x@init))

    lines(x = x@histData$x,
          y = yy, 
          col = 1, lwd = 3, lty = 5)

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
           y = grconvertY(0.9, from = "npc", to = "user"))
      abline(v = 2 * x@init$Mb, col = "orange", lwd = 0.5)
    }
  }

  if(nls & (length(x@nls) > 0)){
    lines(x = x@histData$x, y = predict(x@nls), col = 2)
    text(paste("RCS: ", round(x@RCS, 3)), cex = 1, pos = 2,
         x = grconvertX(0.9, from = "npc", to = "user"),
         y = grconvertY(0.9, from = "npc", to = "user"))
  }

  if(comps & (length(x@nls) > 0)){
    for(i in seq_along(x@comps)){
      if("SCvals" %in% names(formals(x@comps[[i]]@func))){
        params <-
          as.list(coef(x@nls)[names(formals(x@comps[[i]]@func))])
        params <- params[! is.na(names(params))]
        yy <- do.call(x@comps[[i]]@func,
                      args = c(list(SCvals = x@histData$SCvals),
                               params))
        lines(x = x@histData$x, y = yy, col = x@comps[[i]]@color) 
      } else {
        params <-
          as.list(coef(x@nls)[names(formals(x@comps[[i]]@func))])
        params <- params[! is.na(names(params))]
        yy <- do.call(x@comps[[i]]@func,
                      args = c(list(xx = x@histData$x), params))
        lines(x = x@histData$x, y = yy, col = x@comps[[i]]@color)
      }
    }
  }
}

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
             rcs = fh@RCS, row.names = NULL)
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
               rcs = NA, row.names = NULL)
  }
}

#################################################
## Functions for initializing FlowHist objects ##
#################################################

#' (Re-) set the bins for a FlowHist object
#'
#' This function sets the number of bins to use in aggregating FCS data
#'   into a histogram. The 
#' @title setBins
#' @param fh a \code{FlowHist} object
#' @param bins integer, the number of bins to use in aggregating FCS data
#' @return an \code{FlowHist} object, with the \code{bins} slot set to
#'   \code{bins}, and the corresonding binned data stored in a matrix in
#'   the \code{histData} slot. Any previous analysis slots are removed:
#'   \code{peaks, comps, model, init, nls, counts, CV, RCS}.
#' @author Tyler Smith
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
  x <- 1:length(intensity)
  SCvals <- getSingleCutVals(intensity, x)
  fh@histData <- data.frame(x = x , intensity = intensity, SCvals = SCvals)
  
  ## NOTE!! add code to clear out out-dated model data when the hist
  ## changes.
  fh@peaks = matrix()
  fh@comps = list()
  fh@model = function(){}
  fh@init = list()
  fh@nls = structure(list(), class = "nls")
  fh@counts = list()
  fh@CV = list()
  fh@RCS = NA_real_
  fh
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
  if(class(fh) != "FlowHist")
    stop("fh must be a FlowHist object")
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

