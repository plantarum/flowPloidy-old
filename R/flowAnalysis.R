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
#' @param fh a \code{\link{FlowHist}} object
#' @return a \code{\link{FlowHist}} object with the analysis (nls, counts,
#'   cv, RCS) slots filled.
#' @seealso \code{\link{FlowHist}}
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData)
#' fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN")
#' fh1 <- fhAnalyze(fh1)
#' @export
fhAnalyze <- function(fh){
  message("analyzing ", fhFile(fh))
  tryVal <- try(fh <- fhDoNLS(fh), silent = TRUE)
  if(inherits(tryVal, "try-error")){
    message("    *** Model Fit Needs Attention: ", fhFile(fh), " ***")
  } else {
    fh <- fhDoCounts(fh)
    fh <- fhDoCV(fh)
    fh <- fhDoRCS(fh)
  }
  return(fh)
}

#' Fit the NLS model for a \code{\link{FlowHist}} model
#'
#' Constructs the call to \code{\link{nlsLM}} for the
#' \code{\link{FlowHist}} object.
#' 
#' @param fh a \code{\link{FlowHist}} object
#' @return The \code{\link{FlowHist}} object with the \code{NLS} slot
#'   updated to include the results of the analysis
#' @author Tyler Smith
#' @seealso \code{\link{nlsLM}}, \code{\link{fhDoCV}},
#'   \code{\link{fhDoRCS}}, \code{\link{fhDoCounts}}
#' @keywords internal
fhDoNLS <- function(fh){
  model <- fhModel(fh)
  form1 <- paste("intensity ~ model(")
  args <- fhArgs(fh)
  pLims <- fhLimits(fh)
  pLims <- pLims[! names(pLims) %in% fhSpecialParams(fh)]
  lLims <- sapply(pLims, function(x) x[1])
  uLims <- sapply(pLims, function(x) x[2])
  args <- paste(args, collapse = ", ")
  form3 <- paste(", ", getSpecialParamArgs(fh), ")")
  form <- as.formula(paste(form1, args, form3))

  ## Ignore the lowest channels, before we start modeling the debris
  ## component.
  dat <- fhHistData(fh)
  start <- fhStart(dat$intensity)
  dat <- dat[-(1:(start - 1)), ]

  fhNLS(fh) <- nlsLM(formula = form, start = fhInit(fh),
                     data = dat, 
                     lower = lLims, upper = uLims,
                     control = list(ftol = .Machine$double.xmin,
                                    ptol = .Machine$double.xmin,
                                    maxiter = 1024))
  return(fh)
}
  
#' Calculcate \code{\link{FlowHist}} nuclei counts
#'
#' Uses the fitted NLS model for a \code{\link{FlowHist}} object to
#' calculate the cell counts for the G1 peak model components. The actual
#' values are generated by numerical integration of the model components.
#' 
#' @param fh a \code{\link{FlowHist}} object
#' @return The updated \code{\link{FlowHist}} object with the \code{counts}
#'   slot updated.
#' @author Tyler Smith
#' @keywords internal
#' @seealso \code{\link{integrate}}, \code{\link{fhDoCV}}
#'   \code{\link{fhDoNLS}}, \code{\link{fhDoRCS}}.
fhDoCounts <- function(fh){
  ## lower was originally an argument to fhCount, but I don't think it will
  ## ever be anything other than 0?
  lower = 0
  ## similarly, upper was an argument, but it should always be the number
  ## of bins 
  upper = nrow(fhHistData(fh))
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
    integrate(mcFunc(fhComps(fh)$fA1), a1 = coef(fhNLS(fh))["a1"],
              Ma = coef(fhNLS(fh))["Ma"],
              Sa = coef(fhNLS(fh))["Sa"],
              lower = lower, upper = upper,
              subdivisions = 1000)
  if("fB1" %in% names(fhComps(fh))){
    secondPeak <-
      integrate(mcFunc(fhComps(fh)$fB1), b1 = coef(fhNLS(fh))["b1"],
                Mb = coef(fhNLS(fh))["Mb"],
                Sb = coef(fhNLS(fh))["Sb"],
                lower = lower, upper = upper,
                subdivisions = 1000)
  } else {
    secondPeak <- NULL
  }

  if("fC1" %in% names(fhComps(fh))){
    thirdPeak <-
      integrate(mcFunc(fhComps(fh)$fC1), c1 = coef(fhNLS(fh))["c1"],
                Mc = coef(fhNLS(fh))["Mc"],
                Sc = coef(fhNLS(fh))["Sc"],
                lower = lower, upper = upper,
                subdivisions = 1000)
  } else {
    thirdPeak <- NULL
  }

  fhCounts(fh) <- list(firstPeak = firstPeak, secondPeak = secondPeak,
                       thirdPeak = thirdPeak)

  fh
}  

#' Calculate CVs from a \code{\link{FlowHist}} object
#'
#' Extracts the model parameters (G1 peak means and standard deviations)
#' and calculates the CVs. It also calculates the standard errors for
#' the peak ratios.
#'
#' Note that the standard errors here are in fact the SE for the model fit
#' to the particular data set, NOT the SE for the DNA content ratio. In
#' other words, it's a measure of how well the model has fit the data, not
#' a measure of confidence in the actual amount of DNA in the original
#' samples. It's almost always very small, even with very noisy data.
#' 
#' @param fh a \code{\link{FlowHist}} object
#' @return The updated \code{\link{FlowHist}} object.
#' @seealso \code{\link{deltaMethod}}, \code{\link{fhDoCounts}},
#'   \code{\link{fhDoNLS}}, \code{\link{fhDoRCS}}.
#' @author Tyler Smith
#' @aliases PeakRatio
#' @keywords internal
fhDoCV <- function(fh){
  CVa <- coef(fhNLS(fh))["Sa"]/coef(fhNLS(fh))["Ma"]
  if("fB1" %in% names(fhComps(fh))){
    CVb <- coef(fhNLS(fh))["Sb"]/coef(fhNLS(fh))["Mb"]
    AB <- deltaMethod(fhNLS(fh), "Ma/Mb")
  } else {
    CVb <- CI <- AB <- NULL
  }

  if("fC1" %in% names(fhComps(fh))){
    CVc <- coef(fhNLS(fh))["Sc"]/coef(fhNLS(fh))["Mc"]
    AC <- deltaMethod(fhNLS(fh), "Ma/Mc")
    BC <- deltaMethod(fhNLS(fh), "Mb/Mc")
  } else {
    CVc <- AC <- BC <- NULL
  }
  
  fhCV(fh) <- list(CVa = CVa, CVb = CVb, CVc = CVc,
                   AB = AB, AC = AC, BC = BC)
  fh
}

#' Calculate the Residual Chi-Square for a \code{\link{FlowHist}} object
#'
#' Calculate the Residual Chi-Square value for a \code{\link{FlowHist}}
#' model fit. 
#'
#' @section Overview:
#' 
#' The algorithm used to fit the non-linear regression model works by
#' adjusting the model parameters to minimize the Chi-Square value for the
#' resulting fit. The Chi-Square value calculates the departure of observed
#' values from the values predicted by the fitted model:
#'
#' \deqn{\Chi^2 = \sum \frac{(observed(x) -
#'   predicted(x))^2}{observed(x)}}{%
#' Chi^2 = Sum [(observed(x) - predicted(x))^2 / observed(x)]}
#'
#' This would make the Chi-Square value a natural choice for an index to
#' determine the overall goodness-of-fit of the model. However, the
#' Chi-Square value is sensitive to the number of data points in our
#' histogram. We could aggregate the same data into 256, 512 or 1024 bins.
#' All else being equal, the analysis based on 256 bins would give us a
#' lower Chi-Square value than the analyses that use more bins, despite
#' providing essentially identical results.
#'
#' Bagwell (1993) suggested using the Reduced Chi-Square (RCS) value as an
#' superior alternative. It is defined as:
#'
#' \deqn{RCS = \frac{\Chi^2}{n - m}}{RCS = Chi^2/(n - m)}
#'
#' Where n is the number of data points (bins), and m is the number of
#' model parameters. Thus, we correct for the inflation of the Chi-Square
#' value that obtains for higher numbers of bins. At the same time, we
#' introduce a penalty for increasing model complexity, increasing the
#' Chi-Square value proportional to the number of model parameters. This
#' helps us protect against over-fitting the model.
#'
#' @section Guidelines:
#' 
#' As a rule of thumb, RCS values below 0.7 suggest over-fitting, and above
#' 4.0 suggest a poor fit.
#'
#' These are guidelines only, and should not be treated as significance
#' tests. From a statistical perspective, there is limited value to a
#' 'goodness-of-fit' index for a single model. In other contexts we'd
#' compare several competing models to determine which is better. For this
#' application, the RCS is serving as a rough sanity check.
#'
#' Additionally, the absolute value of the RCS is influenced by particular
#' design decisions I made in writing the model-fitting routines.
#' Consequently, other, equally valid approaches may yield slightly
#' different values (Rabinovitch 1994).
#'
#' With this in mind, as long as the values are close to the ideal range
#' 0.7-4.0, we can be reasonably confident that our anlaysis is acceptable.
#' If we get values outside this range, it is a caution that we ought to
#' carefully inspect our model fit, to make sure it appears sensible; the
#' results may still be fine.
#'
#' The most common issue identified by extreme RCS values is poor fitting
#' of the debris component. Occassionally, an otherwise sensible looking
#' model fit will produce extremely high RCS values. Switching from
#' Single-Cut to Multiple-Cut, or vice versa, will often provide a much
#' better fit, with a corresondingly lower RCS value. Visually, the fit may
#' not look much different, and usually the model parameters don't change
#' much either way.
#'
#' @references
#'
#' Bagwell, C.B., 1993. Theoretical aspects of flow cytometry data
#' analysis. pp.41-61 in Clinical flow cytometry: principles and
#' applications. Baltimore: Williams & Wilkins.
#'
#' Rabinovitch, P. S. 1994. DNA content histogram and cell-cycle analysis.
#' Methods in Cell Biology 41:263-296.
#' 
#' @param fh a \code{\link{FlowHist}} object
#' @return The updated \code{\link{FlowHist}} object.
#' @author Tyler Smith
#' @keywords internal
#' @aliases RCS
#' @seealso \code{\link{fhDoCV}}, \code{\link{fhDoNLS}},
#'   \code{\link{fhDoCounts}}, \code{\link{DebrisModels}}
fhDoRCS <- function(fh){
  ## Ignoring the lowest channels, before the debris component starts. We
  ## calculate RCS based on the number of channels fit in the model, not
  ## the full data set, which includes a number of empty/unmodelled
  ## channels at the beginning.
  dat <- fhHistData(fh)
  start <- fhStart(fhHistData(fh)$intensity)
  dat <- dat[-(1:(start - 1)), ]
  obs <- dat$intensity

  exp <- predict(fhNLS(fh))
  chi <- sum(((obs - exp)^2) / exp)

  fhRCS(fh) <- chi/summary(fhNLS(fh))$df[2]
  fh
}
