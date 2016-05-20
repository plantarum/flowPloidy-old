chan = "FL3.INT.LIN"

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
