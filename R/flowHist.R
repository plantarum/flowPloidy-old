## Functions for creating and viewing flowHist objects.

#' @importFrom flowCore read.FCS exprs pData parameters
NULL

##' Create flowHist objects from an FCS file or a flowFrame object
##'
##' Creates a \code{flowHist} object from an FCS file, or a
##' \code{flowFrame} object already in memory
##'
##' Starting with a \code{flowFrame} object, either read from a FCS data
##' file or already in memory, \code{flowHist} will:
##'
##' \enumerate{
##' \item Extract the intensity data from CHANNEL. The actual channel to
##' use will depend on the original FCS file. In our case, we use
##' "FL3.INT.LIN". See the examples for hints on how to find the right
##' channel.  
##'
##' \item Remove the top bin, which contains off-scale readings we ignore
##' in the analysis.
##'
##' \item aggregates the raw data into the desired number of bins, as
##' specified with the 'bins' argument. The default is 256, but you may
##' also try 128 or 512. Any integer is technically acceptable, but I
##' wouldn't stray from the default without a good reason.
##'
##' \item Identify starting values for Gaussian model components. For
##' reasonably clean data, the built-in peak detection is fine. You can
##' evaluate this by plotting the \code{flowHist} object with the argument
##' \code{init = TRUE}. If it doesn't look good, you can play with the
##' \code{window} and \code{smooth} arguments, or pick the peaks yourself
##' with \code{pick = TRUE}.
##'
##' \item Build the NLS model. All \code{flowHist} objects will have
##' components for the single cut debris model and the G1 peaks for the
##' sample and the internal standard. If the G1 peaks are below half (i.e.,
##' the are at < 128 on a 256 bin histogram), their G2 peaks will be
##' modelled as well. All the components are combined into a single model.
##'
##' \item Initial values for the model parameters are determined, using the
##' function \code{flowInit}.
##' }
##' 
##' @title flowHist
##' @param FCS a \code{flowFrame} object, as created by the \code{flowCore::read.FCS} function
##' @param FILE the path to an FCS file
##' @param CHANNEL the name of the channel containing the intensity data
##'   for the histogram
##' @param bins an integer specifying the number of bins to use in building
##'   the histogram
##' @param window the width of the moving window used to identify local
##'   maxima for peak detection via \code{caTools:runmax}
##' @param smooth the width of the moving window used to reduce noise in
##'   the histogram via \code{caTools::runmean}
##' @param pick boolean; if TRUE, the user will be prompted to select peaks
##'   to use for starting values. Otherwise (the default), starting values
##'   will be detected automatically.
##' @return A \code{flowHist} object. Initially, it has the following
##'   slots:
##'
##' \itemize{
##'
##' \item data: a dataframe containing the histogram data. Column 'x' is
##' the index position, and column 'intensity' contains the fluorescence
##' intensity at that position.
##'
##' \item file: the filename of the source data
##'
##' \item peaks: a matrix containing the location of the G1 peaks, as
##' determined by the user or the automated peak-finding routine.
##'
##' \item comps: a list of the individual model components included in the
##' model.
##'
##' \item model: the complete model used in the NLS regression
##'
##' \item init: the inital parameter estimates used in the NLS regression
##'
##' }
##'
##' Additionally, after the data is analyzed by supporting functions, the
##' following slots will be added:
##'
##' \itemize{
##'
##' \item nls: the fitted nls model for the data
##'
##' \item counts: the number of modelled events in each model component
##'
##' \item cv: the coefficients of variation for the modelled peaks, and the
##' confidence interval for the ratio of the sample and co-chopped standard
##' G1 peaks
##'
##' \item RCS: the regularized Chi-Square value for the regression. There
##' are some issues with the interpretation of this value. It can be quite
##' sensitive to values at the end of the histogram, as one example.
##'
##' }
##' 
##' @author Tyler Smith
##'
##' @examples
##' \dontrun{
##' ## To find possible values for CHANNEL:
##' fcs <- read.FCS("path-to-file.lmd", dataset = 1, alter.names = TRUE)
##' fcs
##' ## This will present a summary like this:
##' ##
##' ## flowFrame object 'filename'
##' ## with 11369 cells and 9 observables:
##' ##             name         desc range minRange maxRange
##' ## $P1   FS.INT.LIN   FS INT LIN  1024        0     1023
##' ## $P2   SS.INT.LIN   SS INT LIN  1024        0     1023
##' ## $P3         TIME         TIME  1024        0     1023
##' ## $P4  FL3.INT.LIN  FL3 INT LIN  1024        0     1023
##' ## $P5 FL3.PEAK.LIN FL3 PEAK LIN  1024        0     1023
##' ## $P6  FS.PEAK.LIN  FS PEAK LIN  1024        0     1023
##' ## $P7  SS.PEAK.LIN  SS PEAK LIN  1024        0     1023
##' ## $P8   FS.TOF.LIN   FS TOF LIN  1024        0     1023
##' ## $P9   SS.TOF.LIN   SS TOF LIN  1024        0     1023
##' ## 241 keywords are stored in the 'description' slot
##'
##' ## possible channels to use are listed in the 'name' column.
##'
##' }
##' @export
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

  res <- list(channel = CHANNEL,
              data = data.frame(x = 1:length(intensity),
                                intensity = intensity),
              nls = NULL, standard = NULL,
              file = FCS@description$GUID)

  if(pick)
    res$peaks <- pickPeaks(res)
  else{
    res$peaks <- findPeaks(res, window = window, smooth = smooth)
    res$peaks <- cleanPeaks(res$peaks, window = window)  
  }
    ## res$peaks <- cleanPeaks(findPeaks(res, window = window, smooth = smooth),
    ##                         window = window)  

  res$comps <- list(singleCut, fA1)

  if(res$peaks[1, "mean"] * 2 <= nrow(res$data))
    res$comps <- c(res$comps, fA2)

  if(nrow(res$peaks) > 1){
    res$comps <- c(res$comps, fB1)
    if(res$peaks[2, "mean"] * 2 <= nrow(res$data))
      res$comps <- c(res$comps, fB2)
  }
  
  res$model <- makeModel(res$comps)
  ##res$model = makeModel(res$comps, env = globalenv())

  res$init <- flowInit(res)
  class(res) <- "flowHist"
  
  return(res)
}

##' @export
print.flowHist <- function(self){
  message("flowHist object")
  message("===============")
  message("Source file: ", self$file)
  message("Channel: ", self$channel)
  message("Values: ", dim(self$data)[1])
  message("Total events: ", sum(self$data$intensity))
  message("Model components: ",
          paste(unlist(lapply(fh1$comps,
                              FUN = function(x) attr(x, "compName"))),
                collapse = ", "))
  
  if(is.null(self$standard)){
    message("Standard not specified")
  } else {
    message(paste("Standard GC value: ", substitute(self$standard)))
  }

  if(is.null(self$nls)){
    message("Not fit")
  } else {
    ## replace this with some measure of goodness-of-fit
    message(paste("Fit!"))
  }

  if(!is.null(self$counts)){
    message("\nAnalysis\n========\n", paste("Modelled events:", round(self$counts$total$value, 1)))
    ## message(paste("Peak A:", round(self$counts$firstPeak$value, 1), " at ",
    ##               round(coef(self$nls)["Ma"], 1)))
    ## message(paste("Peak B:", round(self$counts$secondPeak$value, 1), " at ",
    ##               round(coef(self$nls)["Mb"], 1)))
    counts <- c(self$counts$firstPeak$value,
                   self$counts$secondPeak$value)
    ##if(length(counts) == 1) counts <- c(counts, NA)
    size <- c(coef(self$nls)["Ma"],  coef(self$nls)["Mb"])
    if(is.na(size[2])) size <- size[1]
  }

  if(!is.null(self$cv)){
    cvs <- c(self$cv$CVa, self$cv$CVb)
    ##if(length(cvs) == 1) cvs <- c(cvs, NA)
    if(!is.null(self$cv$CVb)){
      message(paste("Ratio Peak A / Peak B: ", round(self$cv$CI[1], 3), ", SE: ",
                    round(self$cv$CI[2], 5), sep = ""))
    }
  }

  if(!is.null(self$counts) & !is.null(self$cv)){
    message("\nPeak Data")
    cat("=========\n")
    if(length(counts) == 2)
      rnames <- c("Peak A", "Peak B")
    else if (length(counts) == 1)
      rnames <- "Peak A"
    peaktable <- kable(data.frame(counts = counts, size = size, cvs = cvs,
                                  row.names = rnames), format = "markdown",
                       digits = 3)
    
    for(i in 1:length(peaktable))
      message(peaktable[i])
  }
  
  if(!is.null(self$RCS)){
    message()
    message(paste("RCS:", round(self$RCS, 3)))
  }

}

##' Plot the raw data for a flowHist object
##'
##' Creates a simple plot of the raw histogram data. Used as a utility for
##' other plotting functions, and perhaps useful for users who wish to
##' create their own plotting routines.
##' 
##' @param self a flowHist object
##' @param ... additional parameters passed to \code{plot}
##' @return Not applicable, used for plotting
##' @author Tyler Smith
##' @export
plotFH <- function(self, ...){
  ## plots the raw data for a flowHist object
  plot(self$data$intensity, type = 'n', main = self$file,
       ylab = "Intensity", xlab = self$channel, ...)
  polygon(x = c(self$data$x, max(self$data$x) + 1), y = c(self$data$intensity, 0),
          col = "lightgray", border = NA)
}
  
##' Plot histograms for flowHist objects
##'
##' .. content for \details{} ..
##' @title plot.flowHist
##' @param self a flowHist object
##' @param init boolean; if TRUE, plot the regression model using the
##'   initial parameter estimates over the raw data. 
##' @param nls boolean; if TRUE, plot the fitted regression model over the
##'   raw data (i.e., using the final parameter values)
##' @param comps boolean; if TRUE, plot the individual model components
##'   over the raw data.
##' @return Not applicable
##' @author Tyler Smith
##' @export
plot.flowHist <- function(self, init = FALSE, nls = TRUE, comps = TRUE){
  plotFH(self)
  
  if(init){
    ##iv <- fHcall(self, "getInitial")
    yy <- do.call(self$model,
                  args = c(list(intensity = self$data$intensity,
                                xx = self$data$x),
                           self$init))

    lines(x = self$data$x,
          y = yy, 
          col = 1, lwd = 3, lty = 5)
  }
  
  if(nls & (! is.null(self$nls))){
    lines(x = self$data$x, y = predict(self$nls), col = 2)
  }
  
  if(comps & (! is.null(self$nls))){
    for(i in seq_along(self$comps)){
      if("intensity" %in% names(formals(self$comps[[i]]))){
        params <-
          as.list(coef(self$nls)[names(formals(self$comps[[i]]))])
        params <- params[! is.na(names(params))]
        yy <- do.call(self$comps[[i]],
                      args = c(list(intensity = self$data$intensity,
                                    xx = self$data$x),
                               params))
        lines(x = self$data$x, y = yy, col = i + 2)
      } else {
        params <-
          as.list(coef(self$nls)[names(formals(self$comps[[i]]))])
        params <- params[! is.na(names(params))]
        yy <- do.call(self$comps[[i]],
                      args = c(list(xx = self$data$x), params))
        lines(x = self$data$x, y = yy, col = i + 2)
      }
    }
  }
}

