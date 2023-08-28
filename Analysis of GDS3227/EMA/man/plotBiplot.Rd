\name{plotBiplot}
\alias{plotBiplot}

\title{Sample and variable representation on a same graph for PCA}
\description{
  Sample and variable representation on a same graph for Principal Component Analysis (PCA)
}
\usage{
plotBiplot(acp, ...)
}

\arguments{
  \item{acp}{result from PCA or do.pca function}
    \item{...}{Arguments to be passed to methods, such as graphical parameters (see 'par').}

}

\value{
 Plot of samples and variables on a same graph
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link{runPCA}},\code{\link{PCA}} }
\examples{
\dontrun{
data(marty)

## PCA on sample
## example dataset
example.subset <- marty[1:100,1:100]
pca <- runPCA(t(example.subset), verbose = FALSE, plotSample = FALSE,
    plotInertia = FALSE)

## Biplot of PCA object
plotBiplot(pca)
}
}

\keyword{plot}
