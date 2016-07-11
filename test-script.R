library(devtools)
load_all()
chan = "FL3.INT.LIN"
files <- list.files(system.file("extdata", package = "flowPloidy"),
                    full.names = TRUE)

i <- 1
#filei <- system.file("extdata", files[i], package = "flowPloidy")
filei <- files[i]
fhi <- flowHist(FILE = filei, CHANNEL = chan)
plot(fhi, init = TRUE)
fhi <- pickInit(fhi)
fhi <- fhAnalyze(fhi)
plot(fhi)
fhi


fhNLS.lm <- function(fh){
  model <- fh$model
  form1 <- paste("intensity ~ model(")
  args <- as.character(names(formals(fh$model)))
  args <- args[!args %in% c("", "intensity", "xx")]
  args <- paste(args, collapse = ", ")
  form3 <- ", intensity = intensity, xx = x)"
  form <- as.formula(paste(form1, args, form3))

  eval(call("nlsLM", form, start = fh$init, data = fh$data,
            lower = rep(0, length = length(fh$init)),
            control = list(ftol = .Machine$double.xmin,
                           ptol = .Machine$double.xmin))) 
}



flowDat <- fhi$dat
flowSC <- singleCutVect(1, flowDat$intensity, flowDat$x)

flowFun <- function (SCa, SCvals, xx, a1, Ma, Sa, a2, b1, Mb, Sb) 
{
  SCvals * SCa
} + {
  (a1/(sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 * Sa^2)))
} + {
  (a2/(sqrt(2 * pi) * Sa * 2) * exp(-((xx - Ma * 2)^2)/(2 * 
                                                        (Sa * 2)^2)))
} + {
  (b1/(sqrt(2 * pi) * Sb) * exp(-((xx - Mb)^2)/(2 * Sb^2)))
}

library(pander)
myReport <- Pandoc$new("Tyler Smith", "flowPloidy Test")
myReport$format <- "html"
res <- list()
for(i in seq_along(files)){
  message("processing ", files[i])
  list[[files[i]]] <- NULL
  filei <- system.file("extdata", files[i], package = "flowPloidy")
  fhi <- flowHist(FILE = filei, CHANNEL = chan)
  try(list[[files[[files[i]]] <- fhAnalyze(fhi))
  myReport$add.paragraph(paste("#", files[i]))
  myReport$add(exportFlowHist(fhi))
  myReport$add(plot(fhi))
}              

histReport <- function(hb, author, title, reportFile = NULL,
                       verbose = TRUE){
  if(is.null(reportFile))
    reportFile <- tempfile(tmpdir = getwd())
  myReport <- Pandoc$new(author, title)
  for(i in seq_along(hb)){
    myReport$add(plot(hb[[i]], init = TRUE))
  }

  myReport$format <- "html"
  myReport$export(reportFile, options = " ")

}
  
histProcess <- function(files, chan, author, title,
                        dataFile = NULL, reportFile = NULL, bins = 256,
                        verbose = TRUE){ 
  if(is.null(reportFile))
    reportFile <- tempfile(tmpdir = getwd())
  myReport <- Pandoc$new(author, title)
  res <- list()

  for(i in seq_along(files)){
    if(verbose) message("processing ", files[i])
    #filei <- system.file("extdata", files[i], package = "flowPloidy")
    res[[files[i]]] <- flowHist(FILE = files[i], CHANNEL = chan)
    tryVal <- try(res[[files[i]]] <- fhAnalyze(res[[files[i]]]))
    if(verbose && inherits(tryVal, "try-error")) message("-- analysis failed")
    myReport$add(plot(res[[files[i]]], init = TRUE))
  }              

  exportFlowHist(res, file = dataFile)
  myReport$format <- "html"
  myReport$export(reportFile, options = " ")

  return(res)
}
  
## 12: 734.LMD requires manually setting the init values
## 8: "337.LMD" requires manually setting the init values
## "240+S.LMD" requires manually setting the init values
## "248+r.LMD" check linearity

file1 <- system.file("extdata", "188-15.LMD", package = "flowPloidy")
fh1 <- flowHist(FILE = file1, CHANNEL = chan)
plot(fh1, init = TRUE)
fh1 <- fhAnalyze(fh1)
plot(fh1)
fh1

file2 <- system.file("extdata", "SM239.LMD", package = "flowPloidy")
fh2 <- flowHist(FILE = file2, CHANNEL = chan)
plot(fh2, init = TRUE)
fh2 <- fhAnalyze(fh2)
plot(fh2)
fh2

file3 <- system.file("extdata", "226.LMD", package = "flowPloidy")
fh3 <- flowHist(FILE = file3, CHANNEL = chan)
plot(fh3, init = TRUE)
fh3 <- fhAnalyze(fh3)
plot(fh3)
fh3


222.LMD

