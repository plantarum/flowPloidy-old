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

fhComponents$fA1 <-
  new("modelComponent", name = "fA1", color = "blue",
      desc = "Gaussian curve for G1 peak of sample A",
      includeTest = function(fh) {TRUE},
      func = function(a1, Ma, Sa, xx){
        (a1 / (sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 * Sa^2)))
      },
      initParams = function(fh){
        Ma <- fh@peaks[1, "mean"]
        Sa <- Ma / 20
        a1 <- fh@peaks[1, "height"] * Sa / 0.45
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
        Ma <- fh@peaks[1, "mean"]
        Sa <- Ma / 20
        a1 <- fh@peaks[1, "height"] * Sa / 0.45
        a2 <- fh@histData[fh@peaks[1, "mean"] * 2, "intensity"] *
          Sa * 2 / 0.45
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
        Mb <- fh@peaks[1, "mean"]
        Sb <- Mb / 20
        b1 <- fh@peaks[1, "height"] * Sb / 0.45
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
        Mb <- fh@peaks[1, "mean"]
        Sb <- Mb / 20
        b2 <- fh@histData[fh@peaks[1, "mean"] * 2, "intensity"] *
          Sb * 2 / 0.45
        list(b2 = b2)
      }
      )

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
        list(BRA = 10)
      }
      )

addComponents4 <- function(fh){
  for(i in fhComponents)
    if(i@includeTest(fh))
      fh@comps[[i@name]] <- i
  fh
}

makeModel4 <- function(fh, env = parent.frame()){
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
  fh@nls <- structure(list(), class = "nls.lm")
  fh
}

getInit4 <- function(fh){
  fh@init <- list()
  for(i in fh@comps){
    fh@init <- c(fh@init, i@initParams(fh1S4))
  }
  fh
}
