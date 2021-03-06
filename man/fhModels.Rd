% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{fhModels}
\alias{fhModels}
\alias{addComponents}
\alias{flowModels}
\alias{dropComponents}
\alias{setLimits}
\alias{makeModel}
\title{Building Flow Histogram Models}
\usage{
addComponents(fh)

dropComponents(fh, components)

setLimits(fh)

makeModel(fh, env = parent.frame())
}
\arguments{
\item{fh}{a \code{\link{FlowHist}} object}

\item{components}{character, a vector of \code{\link{ModelComponent}}
names.}

\item{env}{an R environment. Don't change this, it's R magic to keep the
appropriate environment in scope when building our model.}
}
\value{
The updated \code{\link{FlowHist}} object.
}
\description{
Functions for assembling non-linear regression models for
\code{\link{FlowHist}} objects.
}
\details{
\code{\link{addComponents}} examines the model components in
\code{fhComponents} and includes the ones that pass their
\code{includeTest}.

\code{\link{dropComponents}} removes a component from the
\code{\link{FlowHist}} model

\code{\link{setLimits}} collates the parameter limits for the model
  components included in a \code{\link{FlowHist}} object. (could be
  called automatically from \code{\link{addComponents}}, as it already
  is from \code{\link{dropComponents}}?)

\code{\link{makeModel}} creates a model out of all the included
components.
}
\author{
Tyler Smith
}
