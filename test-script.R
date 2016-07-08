library(devtools)
load_all()
chan = "FL3.INT.LIN"
files <- list.files(system.file("extdata", package = "flowPloidy"),
                    full.names = TRUE)

i <- 3
filei <- system.file("extdata", files[i], package = "flowPloidy")
fhi <- flowHist(FILE = filei, CHANNEL = chan)
plot(fhi, init = TRUE)
fhi <- pickInit(fhi)
fhi <- fhAnalyze(fhi)
plot(fhi)
fhi

library(pander)
myReport <- Pandoc$new("Tyler Smith", "flowPloidy Test")
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



histProcess <- function(author, title, files, chan, outFile, bins = 256, verbose = TRUE){
  myReport <- Pandoc$new("Tyler Smith", "flowPloidy Test")
  res <- list()

  for(i in seq_along(files)){
    if(verbose) message("processing ", files[i])
    #filei <- system.file("extdata", files[i], package = "flowPloidy")
    res[[files[i]]] <- flowHist(FILE = files[i], CHANNEL = chan)
    tryVal <- try(res[[files[i]]] <- fhAnalyze(res[[files[i]]]))
    if(verbose && inherits(tryVal, "try-error")) message("-- analysis failed")
    myReport$add.paragraph(paste("#", files[i]))
    myReport$add(plot(res[[files[i]]], init = TRUE))
  }              

  myReport$add(exportFlowHist(res, file = outFile))
  myReport$export()
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

