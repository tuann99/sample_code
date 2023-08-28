\name{jaccard}
\alias{jaccard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute Jaccard distance on a data matrix}
\description{
  Computes Jaccard distances between the columns of a 0-1-matrix.
}
\usage{
jaccard(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{0-1-matrix.}
}
\details{
  Missing values are allowed, and are excluded from all computations
  involving the cols within which they occur.
}
\value{
  A distance matrix
}
\references{Jaccard, P. (1908), Nouvelles recherches sur la distribution
lorale, Bull. Soc. Vaud. Sci. Nat., 44, 223-270.}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\keyword{cluster}
\keyword{internal}

