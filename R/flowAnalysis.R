## Functions for analyzing flowHist datasets

#' @importFrom car deltaMethod
NULL

#' @importFrom minpack.lm nlsLM
NULL

#' Complete non-linear regression analysis of flowHist histogram data
#'
#' Completes the NLS analysis, and calculates the modelled events and CVs
#' for the result.
#' 
#' @title fhAnalyze
#' @param fh a flowHist object
#' @return a flowHist object with the analysis (nls, counts, cv, RCS)
#'   slots filled.
#' @seealso \code{\link{flowHist}}
#' @author Tyler Smith
#' @export
fhAnalyze <- function(fh){
  tryVal <- try(fh <- fhNLS(fh))
  if(inherits(tryVal, "try-error")){
    message("\n*** Analysis Failed: ", fh$file, " ***\n")
  } else {
    fh$counts <- fhCount(fh)
    fh$cv <- fhCV(fh)
    fh$RCS <- fhRCS(fh)
  }
  return(fh)
}

##' @rdname fhAnalyze
##' @export
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

  obs <- fh$data$intensity
  exp <- predict(fh$nls)
  zeros <- zapsmall(exp, digits = 5) == 0
  obs <- obs[!zeros]
  exp <- exp[!zeros]
  chi <- sum(((obs - exp)^2) / exp)

  n <- length(obs)
  m <- length(formals(fh$model)) - 2    # don't count xx or intensity as
                                        # parameters

  chi/summary(fh$nls)$df[2]
}
                  
##' @rdname fhAnalyze
##' @export
fhNLS <- function(fh){
  model <- fh$model
  form1 <- paste("intensity ~ model(")
  args <- as.character(names(formals(fh$model)))
  args <- args[!args %in% c("", "SCvals", "xx")]
  args <- paste(args, collapse = ", ")
  form3 <- ", SCvals = SCvals, xx = x)"
  form <- as.formula(paste(form1, args, form3))

  fh$nls <- eval(call("nlsLM", form, start = fh$init, data = fh$data,
                      lower = rep(0, length = length(fh$init)),
                      control = list(ftol = .Machine$double.xmin,
                                     ptol = .Machine$double.xmin)))
  return(fh)
}

##' @rdname fhAnalyze
##' @export
fhCount <- function(fh){
  ## lower was originally an argument to fhCount, but I don't think it will
  ## ever be anything other than 0?
  lower = 0
  ## similarly, upper was an argument, but it should always be the number of bins
  upper = nrow(fh$data)
  ## I think anything >= the number of bins should be fine for subdivisions:
  subdivisions = upper * 2
  ## total <-
  ##   do.call(integrate,
  ##           c(substitute(fh$model),
  ##             as.list(coef(fh$nls)),
  ##             SCvals = substitute(fh$data$SCvals),
  ##             lower = lower, upper = upper,
  ##             subdivisions = subdivisions))
  firstPeak <-
    integrate(fA1, a1 = coef(fh$nls)["a1"],
              Ma = coef(fh$nls)["Ma"],
              Sa = coef(fh$nls)["Sa"],
              lower = lower, upper = upper,
              subdivisions = 1000)
  if("fB1" %in% names(fh$comps)){
    secondPeak <-
      integrate(fB1, b1 = coef(fh$nls)["b1"],
                Mb = coef(fh$nls)["Mb"],
                Sb = coef(fh$nls)["Sb"],
                lower = lower, upper = upper,
                subdivisions = 1000)
  } else {
    secondPeak <- NULL
  }
  
  return(list(firstPeak = firstPeak,
              secondPeak = secondPeak)) 
}  

##' @rdname fhAnalyze
##' @export
fhCV <- function(fh){
  CVa <- coef(fh$nls)["Sa"]/coef(fh$nls)["Ma"]
  if("fB1" %in% names(fh$comps)){
    CVb <- coef(fh$nls)["Sb"]/coef(fh$nls)["Mb"]
    CI <- deltaMethod(fh$nls, "Ma/Mb")
  } else {
    CVb <- CI <- NULL
  }
  return(list(CVa = CVa, CVb = CVb, CI = CI))
}

