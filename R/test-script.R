library(flowPloidy)

file1 <- system.file("extdata", "188-15.LMD", package = "flowPloidy")
fh1 <- flowHist(FILE = file1, CHANNEL = "FL3.INT.LIN")
