\name{qualitySample}
\alias{qualitySample}
\title{Sample quality computation in PCA}
\description{
  Sample quality computation in Principal Component Analysis: compute cos2 between samples and each combination of 2 axes
}
\usage{
qualitySample(acp, axes = c(1:3))
}
\arguments{
  \item{acp}{result from PCA or runPCA function}
  \item{axes}{vector of integer, axes number for quality computation, by default 1:3}
}
\details{
  Representation quality for 2 axes is the sum of representation quality for each axe
}
\value{
 A matrix with cos2 between samples and axes
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}


\seealso{\code{\link{runPCA}}, \code{\link{PCA}}}
\keyword{internal}

