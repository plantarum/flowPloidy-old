% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FlowHist.R
\docType{class}
\name{FlowStandards}
\alias{FlowStandards}
\alias{stdSizes}
\alias{stdSelected}
\alias{stdPeak}
\title{An S4 class to represent internal standard details for
\code{\link{FlowHist}} objects}
\usage{
stdSizes(std)

stdSelected(std)

stdPeak(std)
}
\arguments{
\item{std}{a \code{\link{FlowStandards}} object}
}
\value{
\code{\link{stdSizes}}, \code{\link{stdSelected}} and
  \code{\link{stdPeak}} return the corresponding slot values
}
\description{
The \code{sizes} slot is set in \code{\link{FlowHist}} or
\code{\link{batchFlowHist}}. The other values are updates through
interaction with the \code{\link{browseFlowHist}} GUI.
}
\section{Slots}{

\describe{
\item{\code{sizes}}{numeric, the size (in pg) of the internal size standard. Can
be a vector of multiple values, if the sample is part of a set that
included different standards for different samples.}

\item{\code{selected}}{numeric, the size (in pg) of the internal size standard
actually used for this sample. Must be one of the values in the
\code{sizes} slot.}

\item{\code{peak}}{character, "A" or "B", indicating which of the histogram
peaks is the size standard.}
}}

\examples{
library(flowPloidyData) 
fh1 <- FlowHist(file = flowPloidyFiles[1], channel = "FL3.INT.LIN",
                standards = c(1.96, 5.43))
fhStandards(fh1)  ## display standards included in this object
stdSizes(fhStandards(fh1))  ## list standard sizes
}
