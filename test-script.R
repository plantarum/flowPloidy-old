chan = "FL3.INT.LIN"
files <- list.files(system.file("extdata", package = "flowPloidy"))

i <- 12
filei <- system.file("extdata", files[i], package = "flowPloidy")
fhi <- flowHist(FILE = filei, CHANNEL = chan)
plot(fhi, init = TRUE)
fhi <- pickInit(fhi)
fhi <- fhAnalyze(fhi)
plot(fhi)
fhi


## 734.LMD
## Something wrong with this one -- no fA2 component in the model, but it
## should be??

## "337.LMD" requires manually setting the init values
## "240+2.LMD" requires manually setting the init values
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
