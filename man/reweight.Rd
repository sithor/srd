% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/srd_functions.R
\name{reweight}
\alias{reweight}
\title{To  plot an srd drawing object}
\usage{
reweight(x, weight = NULL, alpha = 0.5)
}
\arguments{
\item{x}{binary/logical variable}

\item{weight}{existing weights}

\item{alpha}{forced probability for \eqn{x=1} for new weights}
}
\description{
To  plot an srd drawing object
}
\details{
A function to create weights to enable srd() to create diagram with forced
proportion of \eqn{x=1}
 \code{\link{srd}}.
}
\author{
Roger Marshall <rj.marshall@auckland.ac.nz>, The University of Auckland, New Zealand
}
