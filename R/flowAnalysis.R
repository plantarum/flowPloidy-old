## Functions for analyzing FlowHist datasets

#' @importFrom car deltaMethod
NULL

#' @importFrom minpack.lm nlsLM
NULL

#' Complete non-linear regression analysis of FlowHist histogram data
#'
#' Completes the NLS analysis, and calculates the modelled events and CVs
#' for the result.
#' 
#' @title fhAnalyze
#' @param fh a \code{FlowHist} object
#' @return a \code{FlowHist} object with the analysis (nls, counts, cv,
#'   RCS) slots filled.
#' @seealso \code{\link{FlowHist}}
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData)
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' fh1 <- fhAnalyze(fh1)
#' @export
fhAnalyze <- function(fh){
  message("analyzing ", getFHFile(fh))
  tryVal <- try(fh <- fhNLS(fh))
  if(inherits(tryVal, "try-error")){
    message("\n*** Analysis Failed: ", getFHFile(fh), " ***\n")
  } else {
    fh <- fhCount(fh)
    fh <- fhCV(fh)
    fh <- fhRCS(fh)
  }
  return(fh)
}

fhNLS <- function(fh){
  model <- fh@model
  form1 <- paste("intensity ~ model(")
  args <- as.character(names(formals(fh@model)))
  args <- args[!args %in% c("", names(getSpecialParams(fh)))]
  args <- paste(args, collapse = ", ")
  form3 <- paste(", ", getSpecialParamArgs(fh), ")")
  form <- as.formula(paste(form1, args, form3))

  fh@nls <- eval(call("nlsLM", form, start = fh@init, data = fh@histData,
                      lower = rep(0, length = length(fh@init)),
                      control = list(ftol = .Machine$double.xmin,
                                     ptol = .Machine$double.xmin,
                                     maxiter = 1024)))
  return(fh)
}

fhCount <- function(fh){
  ## lower was originally an argument to fhCount, but I don't think it will
  ## ever be anything other than 0?
  lower = 0
  ## similarly, upper was an argument, but it should always be the number
  ## of bins 
  upper = nrow(fh@histData)
  ## I think anything >= the number of bins should be fine for
  ## subdivisions: 
  subdivisions = upper * 2

  ## This integration is the single slowest step in the analysis, requiring
  ## 10-40 seconds or more to complete. For usability, I'm no longer
  ## supporting it. Perhaps there are ways to speed it up, but all the
  ## important count data is still available quickly.
  ## total <-
  ##   do.call(integrate,
  ##           c(substitute(fh$model),
  ##             as.list(coef(fh$nls)),
  ##             SCvals = substitute(fh$data$SCvals),
  ##             lower = lower, upper = upper,
  ##             subdivisions = subdivisions))
  firstPeak <-
    integrate(fh@comps$fA1@func, a1 = coef(fh@nls)["a1"],
              Ma = coef(fh@nls)["Ma"],
              Sa = coef(fh@nls)["Sa"],
              lower = lower, upper = upper,
              subdivisions = 1000)
  if("fB1" %in% names(fh@comps)){
    secondPeak <-
      integrate(fh@comps$fB1@func, b1 = coef(fh@nls)["b1"],
                Mb = coef(fh@nls)["Mb"],
                Sb = coef(fh@nls)["Sb"],
                lower = lower, upper = upper,
                subdivisions = 1000)
  } else {
    secondPeak <- NULL
  }

  fh@counts <- list(firstPeak = firstPeak, secondPeak = secondPeak)

  fh
}  

fhCV <- function(fh){
  CVa <- coef(fh@nls)["Sa"]/coef(fh@nls)["Ma"]
  if("fB1" %in% names(fh@comps)){
    CVb <- coef(fh@nls)["Sb"]/coef(fh@nls)["Mb"]
    CI <- deltaMethod(fh@nls, "Ma/Mb")
  } else {
    CVb <- CI <- NULL
  }
  fh@CV <- list(CVa = CVa, CVb = CVb, CI = CI)
  fh
}

fhRCS <- function(fh){
  #########################################################################
  ## This may not be a useful measure of analysis quality. It is heavily ##
  ## influenced by the highest channels, where the expected value is     ##
  ## close to zero. Consequently, observing 2 or 3 stray events in one   ##
  ## of these channels, where the expected value may be < 0.02, produces ##
  ## a higher value of (obs - exp)^2 / exp than larger absolute          ##
  ## differences in the main region of the histogram.                    ##
  ##                                                                     ##
  ## Rabinovitch 1994:                                                   ##
  ##                                                                     ##
  ## The x2 is affected by a large number of variables, not all related  ##
  ## to goodness of the fit; these include the number of cells acquired  ##
  ## in the histogram and the end points of the analysis region used     ##
  ## within the histogram.                                               ##
  #########################################################################

  obs <- fh@histData$intensity
  exp <- predict(fh@nls)
  zeros <- zapsmall(exp, digits = 5) == 0
  obs <- obs[!zeros]
  exp <- exp[!zeros]
  chi <- sum(((obs - exp)^2) / exp)

  ## n <- length(obs)
  ## m <- length(formals(fh@model)) - length(getSpecialParams(fh))
  ## ## Don't count special parameters (xx, MCvals etc) or intensity in
  ## ## determining the number of parameters fit in the NLS!

  fh@RCS <- chi/summary(fh@nls)$df[2]
  fh
}
