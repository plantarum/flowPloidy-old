% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{DebrisModels}
\alias{DebrisModels}
\alias{getSingleCutValsBase}
\alias{SingleCut}
\alias{getMultipleCutVals}
\alias{MultipleCut}
\title{Histogram Debris Models}
\usage{
getSingleCutValsBase(intensity, xx, first.channel)

getMultipleCutVals(intensity, first.channel)
}
\arguments{
\item{intensity}{a numeric vector, the histogram intensity in each
channel}

\item{xx}{an integer vector, the ordered channels corresponding to the
values in `intensity'.}

\item{first.channel}{integer, the lowest bin to include in the modelling
process. Determined by the internal function \code{fhStart}.}
}
\value{
\code{getSingleCutVals}, the vectorized function built from
  getSingleCutValsBase, returns the fixed \code{SCvals} for the
  histogram.

\code{getMultipleCutVals}, a vectorized function, returns the
  fixed \code{MCvals} for the histogram.
}
\description{
Implementation of debris models described by Bagwell et al. (1991).
}
\section{Single Cut Model}{


This is the theoretical probability distribution of the size of pieces
formed by a single random cut through an ellipsoid. In other words, we
assume that the debris is composed of nuclei pieces generated by cutting
a subset of the nuclei in a sample into two pieces.

The model is:
\deqn{S(x) = a \sum_{j = x + 1}^{n} \sqrt[3]{j} Y_j P_s(j, x)}

\enumerate{
\item \code{x} the histogram channel that we're estimating the debris
value for.
\item \code{SCa} the amplitude parameter.
\item \code{Y_j} the histogram intensity for channel j.
}

where P_s(j, x) is the probability of a nuclei from channel j falling
into channel x when cut. That is, for j > x, the probability that
fragmenting a nuclei from channel j with a single cut will produce a
fragment of size x. This probability is calculated as:

\deqn{P_s(j, x) = \frac{2}{(\pi j \sqrt{(x/j) (1 - x/j)}}}

This model involves a recursive calculation, since the fitted value
for channel x depends not just on the intensity for channel x, but also
the intensities at all channels > x. I deal with this by pre-calculating
the raw values, which don't actually depend on the only parameter,
\code{SCa}. These raw values are stored in the \code{histData} matrix
(which is a slot in the \code{\link{FlowHist}} object). This must be
accomodated by treating \code{SCvals} as a 'special parameter' in the
\code{\link{ModelComponent}} definition. See that help page for details.
}

\section{Multiple-Cut Model}{


The Multiple-Cut model extends the Single-Cut model by assuming that a
single nuclei may be cut multiple times, thus creating more than two
fragments.

The model is:
\deqn{S(x) = MCa e^{-kx}\sum_{j = x + 1}^{n} Y_j}

\enumerate{
\item \code{x} the histogram channel that we're estimating the debris
value for.
\item \code{k} an exponential fitting parameter
\item \code{MCa} the amplitiude parameter
\item \code{Y_j} the histogram intensity for channel j.
}

This model involves another recursive or "histogram-dependent"
component. Again, the sum is independent of the fitted parameters, so we
can pre-compute that and add it to the \code{histData} slot, in the
column \code{MCvals}. This is treated as a 'special parameter' when the
Multiple-Cut model is applied, so we only need to fit the parameters k
and MCa.
}

\examples{
## This is an internal function, called from setBins()
\dontrun{
  ## ...
  SCvals <- getSingleCutVals(intensity, xx, startBin)
  MCvals <- getMultipleCutVals(intensity, startBin)
  ## ...
  fhHistData(fh) <- data.frame(xx = xx, intensity = intensity,
                           SCvals = SCvals, MCvals = MCvals,
                           DBvals = DBvals, TRvals = TRvals,
                           QDvals = QDvals, gateResid = gateResid)
  ## ...
}
}
\references{
Bagwell, C. B., Mayo, S. W., Whetstone, S. D., Hitchcox, S.
A., Baker, D. R., Herbert, D. J., Weaver, D. L., Jones, M. A. and
Lovett, E. J. (1991), DNA histogram debris theory and compensation.
Cytometry, 12: 107-118. doi: 10.1002/cyto.990120203
}
\author{
Tyler Smith
}
\keyword{internal}
