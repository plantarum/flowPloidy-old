% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny.R
\name{browseFlowHist}
\alias{browseFlowHist}
\title{browseFlowHist}
\usage{
browseFlowHist(flowList, debug = FALSE)
}
\arguments{
\item{flowList}{either a \code{\link{FlowHist}} object, or a list of
\code{\link{FlowHist}} objects}

\item{debug}{boolean, turns on debugging messages}
}
\value{
Returns the list of \code{\link{FlowHist}} objects, updated by
  any changes made in the GUI.
}
\description{
Visually assess and correct histogram fits
}
\details{
Visually assess histogram fits, correcting initial values, and selecting
model components.

This function will open a browser tab displaying the first
\code{\link{FlowHist}} object from the argument \code{flowList}. Using
the interface, the user can modify the starting values for the histogram
peaks, select different debris model components, toggle the linearity
option, select which peak to treat as the standard, and, if multiple
standard sizes are available, select which one to apply.

See the "Getting Started" vignette for a tutorial introduction.
}
\examples{
library(flowPloidyData)
batch1 <- batchFlowHist(flowPloidyFiles, channel = "FL3.INT.LIN")
\dontrun{
batch1 <- browseFlowHist(batch1)
}
}
\author{
Tyler Smith
}
