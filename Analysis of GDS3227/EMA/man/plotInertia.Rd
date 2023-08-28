\name{plotInertia}
\alias{plotInertia}
\title{Barplot of component inertia percentage for PCA}
\description{
  Barplot of component inertia percentage for Principal Component Analysis (PCA)
}
\usage{
plotInertia(acp, ncp = 5, ...)
}

\arguments{
  \item{acp}{result from do.pca or PCA function}
  \item{ncp}{number of components displayed, by default 5}
  \item{...}{Arguments to be passed to methods, such as graphical parameters (see 'par').}
}

\value{
 Barplot of component inertia percentage
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link{runPCA}},\code{\link{PCA}}}
\examples{
data(marty)

## PCA on sample - example set
example.subset <- marty[1:100,]
pca <- runPCA(t(example.subset), verbose = FALSE, plotInertia = FALSE, plotSample = FALSE)

## Inertia plot of PCA object
plotInertia(pca)
}

\keyword{plot}

