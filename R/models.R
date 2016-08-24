## Contains the various model components used in the flowPloidy analysis.

setClass(
  Class = "modelComponent",
  representation = representation(
    name = "character",
    desc = "character", 
    color = "character",
    includeTest = "function",
    func = "function",
    initParams = "function"
  )
)

setMethod(
  f = "show",
  signature = "modelComponent",
  def = function(object){
    cat("** flowHist model component: ")
    cat(object@name); cat(" ** \n")
    cat(object@desc); cat(" \n")
    cat("Parameters: ")
    pnames <- names(formals(object@func))
    pnames <- pnames[which(pnames != "xx")]
    cat(paste(pnames, collapse = ", "))
    cat("\n")
 }
)

fhComponents <- list()

## Define new components with the following template:
##
## fhComponents$<name> <-
##   new("modelComponent", name = "<name>", color = "<colour>",
##       desc = "<one-line description>",
##       includeTest = function(fh){
##         test based on the flowHist object. Return TRUE if the component
##         should be included, FALSE otherwise. 
##       },
##       func = function(){
##         a single-line function that returns the value of the component.
##         Can take multiple arguments, usually one of which will be 'xx'
##       },
##       initParams = function(fh){
##         a function that returns a named list of initial parameter
##         estimates, based on the single argument of the flowHist object 
##         list(param1 = param1, ...)
##       }
##       )

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
fhComponents$fA1 <-
  new("modelComponent", name = "fA1", color = "blue",
      desc = "Gaussian curve for G1 peak of sample A",
      includeTest = function(fh) {TRUE},
      func = function(a1, Ma, Sa, xx){
        (a1 / (sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 * Sa^2)))
      },
      initParams = function(fh){
        Ma <- as.numeric(fh@peaks[1, "mean"])
        Sa <- as.numeric(Ma / 20)
        a1 <- as.numeric(fh@peaks[1, "height"] * Sa / 0.45)
        list(Ma = Ma, Sa = Sa, a1 = a1)
      }
      )

fhComponents$fA2 <-
  new("modelComponent", name = "fA2", color = "blue",
      desc = "Gaussian curve for G2 peak of sample A",
      includeTest = function(fh){
        (fh@peaks[1, "mean"] * 2) <= nrow(fh@histData)
      },
      func = function(a2, Ma, Sa, xx){
        (a2 / (sqrt(2 * pi) * Sa * 2) *
         exp(-((xx - Ma * 2)^2)/(2 * (Sa * 2)^2))) 
      },
      initParams = function(fh){
        Ma <- as.numeric(fh@peaks[1, "mean"])
        Sa <- as.numeric(Ma / 20)
        a2 <- as.numeric(fh@histData[Ma * 2, "intensity"] *
                         Sa * 2 / 0.45)
        list(a2 = a2)
      }
      )

fhComponents$fB1 <-
  new("modelComponent", name = "fB1", color = "orange",
      desc = "Gaussian curve for G1 peak of sample B",
      includeTest = function(fh){
        nrow(fh@peaks) > 1
      },
      func = function(b1, Mb, Sb, xx){
        (b1 / (sqrt(2 * pi) * Sb) * exp(-((xx - Mb)^2)/(2 * Sb^2)))
      },
      initParams = function(fh){
        Mb <- as.numeric(fh@peaks[2, "mean"])
        Sb <- as.numeric(Mb / 20)
        b1 <- as.numeric(fh@peaks[2, "height"] * Sb / 0.45)
        list(Mb = Mb, Sb = Sb, b1 = b1)
      }
      )

fhComponents$fB2 <-
  new("modelComponent", name = "fB2", color = "blue",
      desc = "Gaussian curve for G2 peak of sample B",
      includeTest = function(fh){
        if(nrow(fh@peaks) > 1)
          (fh@peaks[2, "mean"] * 2) <= nrow(fh@histData)
        else
          FALSE
      },
      func = function(b2, Mb, Sb, xx){
        (b2 / (sqrt(2 * pi) * Sb * 2) *
         exp(-((xx - Mb * 2)^2)/(2 * (Sb * 2)^2))) 
      },
      initParams = function(fh){
        Mb <- fh@peaks[2, "mean"]
        Sb <- Mb / 20
        b2 <- as.numeric(fh@histData[fh@peaks[2, "mean"] * 2, "intensity"]
          * Sb * 2 / 0.45)
        list(b2 = b2)
      }
      )

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
#' @param SCvals a numeric vector, stored in the \code{flowHist} object
#'   slot `SCvals`. Users shouldn't need this.
#' @return NA
#'
#' @references Bagwell, C. B., Mayo, S. W., Whetstone, S. D., Hitchcox, S.
#'   A., Baker, D. R., Herbert, D. J., Weaver, D. L., Jones, M. A. and
#'   Lovett, E. J. (1991), DNA histogram debris theory and compensation.
#'   Cytometry, 12: 107-118. doi: 10.1002/cyto.990120203
#'
#' @author Tyler Smith

#' @rdname singleCut


fhComponents$SC <-
  new("modelComponent", name = "SC", color = "green",
      desc = "The single-cut debris model.",
      includeTest = function(fh){
        TRUE
      },
      func = function(SCa, SCvals){
        SCa * SCvals
      },
      initParams = function(fh){
        list(SCa = 0.1)
      }
      )

## Broadened rectangles:
## simplified with a fixed sd of 1. Very little change in results with more
## flexible models, so keeping it simple.
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

fhComponents$brA <-
  new("modelComponent", name = "brA", color = "magenta",
      desc = "Broadened rectangle for S-phase of sample A",
      includeTest = function(fh){
        TRUE
      },
      func = function(BRA, Ma, xx){
        ## 2 * 1 is a placeholder for 2 * sd, should we decide it's worth
        ## adding sd as a separate parameter
        BRA * ((erf(((2 * Ma) - xx)/sqrt(2 * 1)) -
                erf((Ma - xx)/sqrt(2 * 1))) / 2)
      },
      initParams = function(fh){
        list(BRA = 10)
      }
      )

fhComponents$brB <-
  new("modelComponent", name = "brB", color = "turquoise",
      desc = "Broadened rectangle for S-phase of sample B",
      includeTest = function(fh){
        nrow(fh@peaks) > 1        
      },
      func = function(BRB, Mb, xx){
        ## 2 * 1 is a placeholder for 2 * sd, should we decide it's worth
        ## adding sd as a separate parameter
        BRB * ((erf(((2 * Mb) - xx)/sqrt(2 * 1)) -
                erf((Mb - xx)/sqrt(2 * 1))) / 2)
      },
      initParams = function(fh){
        list(BRB = 10)
      }
      )

## for testing the influence of sd:
## This isn't used in any other code, retained here for further study if
## needed. 
## brA1 <- function(BRA, Ma, xx, sd){
##   BRA * ((flowPloidy::erf(((2 * Ma) - xx)/sqrt(2 * sd)) -
##           erf((Ma - xx)/sqrt(2 * sd))) / 2)
## }

## The basic broadened trapezoid functions
## Retained here for study, but the complexity doesn't provide much/any
## useful improvement in the model fit.
## broadenedTrapezoid <- function(BTt1, BTt2, BTx1, BTx2, BTs1, BTs2, xx){
##   ((BTt2 - BTt1) / (BTx2 - BTx1) * (xx - BTx2) + BTt2) *
##     ((erf((BTx2 - xx)/sqrt(2 * BTs2)) -
##       erf((BTx1 - xx)/sqrt(2 * BTs1))) / 2)
## }

## ## Translated into model components:
## btA <- function(BTt1A, BTt2A, Ma, xx){
##   ## Simplified to use a fixed sd (5), and bounded to the G1 mean value
##   ## (and by extension the G2 mean value).
##   ((BTt2A - BTt1A) / Ma * (xx - (2 * Ma)) + BTt2A) *
##     ((erf(((2 * Ma) - xx)/sqrt(2 * 5)) -
##       erf((Ma - xx)/sqrt(2 * 5))) / 2)
## }

## btB <- function(BTt1B, BTt2B, Mb, xx){
##   ## Simplified to use a fixed sd (5), and bounded to the G1 mean value
##   ## (and by extension the G2 mean value).
##   ((BTt2B - BTt1B) / Mb * (xx - (2 * Mb)) + BTt2B) *
##     ((erf(((2 * Mb) - xx)/sqrt(2 * 5)) -
##       erf((Mb - xx)/sqrt(2 * 5))) / 2)
## }


addComponents <- function(fh){
  for(i in fhComponents)
    if(i@includeTest(fh))
      fh@comps[[i@name]] <- i
  fh
}

makeModel <- function(fh, env = parent.frame()){
  components <- fh@comps
  names(components) <- NULL
  args <- unlist(lapply(components, FUN = function(x) formals(x@func)))
  args <- args[unique(names(args))]

  bodList <- lapply(components, FUN = function(x) body(x@func))
  bod <- bodList[[1]]
  bodList <- bodList[-1]

  while(length(bodList) > 0){
    bod <- call("+", bod, bodList[[1]])
    bodList <- bodList[-1]
  }

  fh@model <- eval(call("function", as.pairlist(args), bod), env)
  fh@nls <- structure(list(), class = "nls")
  fh
}

getInit <- function(fh){
  fh@init <- list()
  for(i in fh@comps){
    fh@init <- c(fh@init, i@initParams(fh))
  }
  fh
}
