## Functions for creating and viewing flowHist objects.

#' @importFrom flowCore read.FCS exprs pData parameters
NULL

#' @importFrom knitr kable
NULL

#' @importFrom graphics hist lines locator plot points polygon
NULL

#' @importFrom stats as.formula coef integrate predict
NULL

#' @importFrom utils write.table
NULL

#' Create flowHist objects from an FCS file or a flowFrame object
#'
#' Creates a \code{flowHist} object from an FCS file, or a
#' \code{flowFrame} object already in memory
#'
#' Starting with a \code{flowFrame} object, either read from a FCS data
#' file or already in memory, \code{flowHist} will:
#'
#' \enumerate{
#' \item Extract the intensity data from CHANNEL. The actual channel to
#' use will depend on the original FCS file. In our case, we use
#' "FL3.INT.LIN". See the examples for hints on how to find the right
#' channel.  
#'
#' \item Remove the top bin, which contains off-scale readings we ignore
#' in the analysis.
#'
#' \item aggregates the raw data into the desired number of bins, as
#' specified with the 'bins' argument. The default is 256, but you may
#' also try 128 or 512. Any integer is technically acceptable, but I
#' wouldn't stray from the default without a good reason.
#'
#' \item Identify starting values for Gaussian model components. For
#' reasonably clean data, the built-in peak detection is fine. You can
#' evaluate this by plotting the \code{flowHist} object with the argument
#' \code{init = TRUE}. If it doesn't look good, you can play with the
#' \code{window} and \code{smooth} arguments, or pick the peaks yourself
#' with \code{pick = TRUE}.
#'
#' \item Build the NLS model. All \code{flowHist} objects will have
#' components for the single cut debris model and the G1 peaks for the
#' sample and the internal standard. If the G1 peaks are below half (i.e.,
#' the are at < 128 on a 256 bin histogram), their G2 peaks will be
#' modelled as well. All the components are combined into a single model.
#'
#' \item Initial values for the model parameters are determined, using the
#' function \code{flowInit}.
#' }
#' 
#' @title flowHist
#' @name flowHist
#' @param FCS a \code{flowFrame} object, as created by the \code{flowCore::read.FCS} function
#' @param FILE the path to an FCS file
#' @param CHANNEL the name of the channel containing the intensity data
#'   for the histogram
#' @param bins an integer specifying the number of bins to use in building
#'   the histogram
#' @param window the width of the moving window used to identify local
#'   maxima for peak detection via \code{caTools:runmax}
#' @param smooth the width of the moving window used to reduce noise in
#'   the histogram via \code{caTools::runmean}
#' @param pick boolean; if TRUE, the user will be prompted to select peaks
#'   to use for starting values. Otherwise (the default), starting values
#'   will be detected automatically.
#' @return A \code{flowHist} object. Initially, it has the following
#'   slots:
#'
#' \itemize{
#'
#' \item data: a dataframe containing the histogram data. Column 'x' is
#' the index position, and column 'intensity' contains the fluorescence
#' intensity at that position.
#'
#' \item file: the filename of the source data
#'
#' \item peaks: a matrix containing the location of the G1 peaks, as
#' determined by the user or the automated peak-finding routine.
#'
#' \item comps: a list of the individual model components included in the
#' model.
#'
#' \item model: the complete model used in the NLS regression
#'
#' \item init: the inital parameter estimates used in the NLS regression
#'
#' }
#'
#' Additionally, after the data is analyzed by supporting functions, the
#' following slots will be added:
#'
#' \itemize{
#'
#' \item nls: the fitted nls model for the data
#'
#' \item counts: the number of modelled events in each model component
#'
#' \item cv: the coefficients of variation for the modelled peaks, and the
#' confidence interval for the ratio of the sample and co-chopped standard
#' G1 peaks
#'
#' \item RCS: the regularized Chi-Square value for the regression. There
#' are some issues with the interpretation of this value. It can be quite
#' sensitive to values at the end of the histogram, as one example.
#'
#' }
#' 
#' @author Tyler Smith
#'
#' @examples
#' \dontrun{
#' ## To find possible values for CHANNEL:
#' fcs <- read.FCS("path-to-file.lmd", dataset = 1, alter.names = TRUE)
#' fcs
#' ## This will present a summary like this:
#' ##
#' ## flowFrame object 'filename'
#' ## with 11369 cells and 9 observables:
#' ##             name         desc range minRange maxRange
#' ## $P1   FS.INT.LIN   FS INT LIN  1024        0     1023
#' ## $P2   SS.INT.LIN   SS INT LIN  1024        0     1023
#' ## $P3         TIME         TIME  1024        0     1023
#' ## $P4  FL3.INT.LIN  FL3 INT LIN  1024        0     1023
#' ## $P5 FL3.PEAK.LIN FL3 PEAK LIN  1024        0     1023
#' ## $P6  FS.PEAK.LIN  FS PEAK LIN  1024        0     1023
#' ## $P7  SS.PEAK.LIN  SS PEAK LIN  1024        0     1023
#' ## $P8   FS.TOF.LIN   FS TOF LIN  1024        0     1023
#' ## $P9   SS.TOF.LIN   SS TOF LIN  1024        0     1023
#' ## 241 keywords are stored in the 'description' slot
#'
#' ## possible channels to use are listed in the 'name' column.
#'
#' }
#' @export
flowHist <- function(FCS = NULL, FILE = NULL, CHANNEL,
                     bins = 256, window = 20, smooth = 20, pick = FALSE){
  ## You probably want to subdivide the bins evenly. i.e., if there are
  ## 1024 bins in the data, use 128, 256, 512 bins
  if((1024 %% bins) != 0)
    warning("maxBins is not a multiple of bins!")

  if(sum(c(is.null(FCS), is.null(FILE))) != 1){
    stop("\nOne (and only one) of FCS or FILE must be set.")
  }

  if(!is.null(FILE))
    FCS <- read.FCS(FILE, dataset = 1, alter.names = TRUE)

  res <- list(channel = CHANNEL, nls = NULL, standard = NULL,
              file = FCS@description$GUID) 
  ## Build histogram
  res$data = buildHist(FCS, CHANNEL, bins)
  
  if(pick)
    res$peaks <- pickPeaks(res)
  else{
    res$peaks <- findPeaks(res, window = window, smooth = smooth)
    res$peaks <- cleanPeaks(res$peaks, window = window)  
  }

  ## Add model components
  res$comps <- list(singleCut = singleCut, fA1 = fA1)

  if(res$peaks[1, "mean"] * 2 <= nrow(res$data))
    res$comps <- c(res$comps, fA2 = fA2)

  if(nrow(res$peaks) > 1){
    res$comps <- c(res$comps, fB1 = fB1)
    if(res$peaks[2, "mean"] * 2 <= nrow(res$data))
      res$comps <- c(res$comps, fB2 = fB2)
  }
  
  res$model <- makeModel(res$comps)

  res <- flowInit(res)
  class(res) <- "flowHist"
  
  return(res)
}

#' @export
#' @rdname flowHist
buildHist <- function(FCS, CHANNEL, bins){
  ## Extract the data channel
  chanDat <- exprs(FCS)[, CHANNEL]

  ## remove the top bin - this contains clipped values representing all
  ## out-of-range data, not true values
  chanTrim <- chanDat[chanDat < max(chanDat)]

  metaData <- pData(parameters(FCS))
  maxBins <- metaData[which(metaData$name == CHANNEL), "range"]
  
  ## aggregate bins: combine maxBins into bins via hist
  binAg <- floor(maxBins / bins)

  histBins <- hist(chanTrim, breaks = seq(from = 0, to = 1024, by = binAg),
                   plot = FALSE)

  intensity <- histBins$counts

  return(data.frame(x = 1:length(intensity),
                    intensity = intensity))
}

#' @export
print.flowHist <- function(x, ...){
  cat(paste("flowHist object",
            "\n===============",
            "\nSource file: ", x$file,
            "\nChannel: ", x$channel,
            "\nValues: ", dim(x$data)[1],
            "\nTotal events: ", sum(x$data$intensity),
            "\nModel components: ",
            paste(names(x$comps), collapse = ", "), sep = "")) 
  
  if(is.null(x$standard)){
    cat("\nStandard not specified")
  } else {
    cat(paste("\nStandard GC value: ", substitute(x$standard)))
  }

  if(is.null(x$nls)){
    cat("\nNot fit")
  } else {
    ## replace this with some measure of goodness-of-fit
    cat(paste("\nFit!"))
  }

  if(!is.null(x$counts)){
    cat(paste("\n\nAnalysis\n========\n",
              paste("Modelled events: ", round(x$counts$total$value, 1)), sep = ""))
    counts <- c(x$counts$firstPeak$value,
                   x$counts$secondPeak$value)
    size <- c(coef(x$nls)["Ma"],  coef(x$nls)["Mb"])
    if(is.na(size[2])) size <- size[1]
  }

  if(!is.null(x$cv)){
    cvs <- c(x$cv$CVa, x$cv$CVb)
    if(!is.null(x$cv$CVb)){
      cat(paste("\nRatio Peak A / Peak B: ", round(x$cv$CI[1], 3), ", SE: ",
                    round(x$cv$CI[2], 5), sep = ""))
    }
  }

  if(!is.null(x$counts) & !is.null(x$cv)){
    if(length(counts) == 2)
      rnames <- c("Peak A", "Peak B")
    else if (length(counts) == 1)
      rnames <- "Peak A"
    print(kable(data.frame(counts = counts, size = size, cvs = cvs,
                           row.names = rnames), format = "markdown",
                digits = 3))
  }
  
  if(!is.null(x$RCS)){
    cat(paste("\nRCS:", round(x$RCS, 3), "\n"))
  }

}

#' Plot the raw data for a flowHist object
#'
#' Creates a simple plot of the raw histogram data. Used as a utility for
#' other plotting functions, and perhaps useful for users who wish to
#' create their own plotting routines.
#' 
#' @param self a flowHist object
#' @param ... additional parameters passed to \code{plot}
#' @return Not applicable, used for plotting
#' @author Tyler Smith
#' @export
plotFH <- function(self, ...){
  ## plots the raw data for a flowHist object
  plot(self$data$intensity, type = 'n', main = self$file,
       ylab = "Intensity", xlab = self$channel, ...)
  polygon(x = c(self$data$x, max(self$data$x) + 1), y = c(self$data$intensity, 0),
          col = "lightgray", border = NA)
}
  
#' Plot histograms for flowHist objects
#'
#' .. content for details ..
#' @title plot.flowHist
#' @param x a flowHist object
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
plot.flowHist <- function(x, init = FALSE, nls = TRUE, comps = TRUE, ...){
  plotFH(x, ...)
  
  if(init){
    yy <- do.call(x$model,
                  args = c(list(intensity = x$data$intensity,
                                xx = x$data$x),
                           x$init))

    lines(x = x$data$x,
          y = yy, 
          col = 1, lwd = 3, lty = 5)
  }
  
  if(nls & (! is.null(x$nls))){
    lines(x = x$data$x, y = predict(x$nls), col = 2)
  }

  coltab <- c(fA1 = "blue", fA2 = "blue", fB1 = "orange", fB2 = "orange",
                 singleCut = "green") 
  
  if(comps & (! is.null(x$nls))){
    for(i in seq_along(x$comps)){
      if("intensity" %in% names(formals(x$comps[[i]]))){
        params <-
          as.list(coef(x$nls)[names(formals(x$comps[[i]]))])
        params <- params[! is.na(names(params))]
        yy <- do.call(x$comps[[i]],
                      args = c(list(intensity = x$data$intensity,
                                    xx = x$data$x),
                               params))
        lines(x = x$data$x, y = yy, col = coltab[names(x$comps)[[i]]])
      } else {
        params <-
          as.list(coef(x$nls)[names(formals(x$comps[[i]]))])
        params <- params[! is.na(names(params))]
        yy <- do.call(x$comps[[i]],
                      args = c(list(xx = x$data$x), params))
        lines(x = x$data$x, y = yy, col = coltab[names(x$comps)[[i]]])
      }
    }
  }
}

#' Extract analysis results from a flowHist object
#'
#' A convenience function for extracting the results of the NLS
#'   curve-fitting analysis on a flowHist object.
#'
#' If \code{fh} is a single flowHist object, a data.frame with a single
#' row is returned. If \code{fh} is a list of \code{flowHist} objects, a
#' row for each object will be added to the data.frame.
#'
#' If a file name is provided, the data will be saved to that file.
#' 
#' @title exportFlowHist
#' @param fh a flowHist object, or a list of flowHist objects.
#' @param file character, the name of the file to save dat to
#' @return a data frame 
#' @author Tyler Smith
#' @export
exportFlowHist <- function(fh, file = NULL){
  if(class(fh) == "flowHist")
    res <- exFlowHist(fh)
  else if (class(fh) == "list" && all(sapply(fh, class) == "flowHist")){
    res <- do.call(rbind, lapply(fh, exFlowHist))
  }
  if(! is.null(file))
    write.table(x = res, file = file)

  res
}

exFlowHist <- function(fh){
  data.frame(file = fh$file, channel = fh$channel,
             components = paste(names(fh$comps), collapse = ";"),
             totalEvents = sum(fh$data$intensity),
             modelledEvents = fh$counts$total$value,
             countsA = fh$counts$firstPeak$value,
             countsB = ifelse(is.null(fh$counts$secondPeak$value), NA, fh$counts$secondPeak$value),
             sizeA = coef(fh$nls)["Ma"],
             sizeB = coef(fh$nls)["Mb"],
             cvA = fh$cv$CVa,
             cvB = ifelse(is.null(fh$cv$CVb), NA, fh$cv$CVb),
             ratioAB = unlist(ifelse(is.null(fh$cv$CI[1]), NA, fh$cv$CI[1])),
             ratioSE = unlist(ifelse(is.null(fh$cv$CI[2]), NA, fh$cv$CI[2])),
             rcs = fh$RCS, row.names = NULL)
}
